"""Time tensorQTL on one exported dataset and save per-pair results.

    python bench/tensorqtl/run_tensorqtl.py <dataset> <device> <threads> <design> <n_pcs> <out_prefix> [--perm]

device:  cpu | gpu        (cpu hides the GPU: CUDA_VISIBLE_DEVICES is set before torch is imported)
threads: torch CPU threads
design:  sample  every sample is a phenotype column (genotypes repeated per sample of an individual);
                 eQTL: map_nominal with condition + PCs (+ dataset covariates) as covariates;
                 response: map_nominal with interaction_df = condition (dummies for a factor), PCs as covariates
         indiv   one column per individual (2-level balanced designs): eQTL on the mean of the two
                 conditions, response on the difference (treated - control), each rank-normalized, PCs of each
--perm   also time map_cis (gene-level permutations, 10,000 by default) on the eQTL phenotype
Writes <out_prefix>.tsv.gz (gene, snp, p_geno, p_interact[, p_interact_<term>]) and <out_prefix>.json.
"""
import sys, os, json, time, resource, glob, shutil
ds, device, threads, design, n_pcs, out = sys.argv[1:7]
perm = '--perm' in sys.argv
threads, n_pcs = int(threads), int(n_pcs)
if device == 'cpu':
    os.environ['CUDA_VISIBLE_DEVICES'] = ''
import numpy as np, pandas as pd, scipy.stats as st, torch
torch.set_num_threads(threads)
from tensorqtl import cis
if device == 'gpu':
    assert torch.cuda.is_available()
    t0 = time.perf_counter(); torch.ones(1, device='cuda').sum().item(); t_cuda_init = time.perf_counter() - t0
else:
    assert not torch.cuda.is_available(); t_cuda_init = 0.0

D = os.path.join('bench/tensorqtl/data', ds)
t0 = time.perf_counter()
pheno = pd.read_csv(f'{D}/pheno.tsv.gz', sep='\t', index_col=0).astype(np.float32)
ppos = pd.read_csv(f'{D}/pheno_pos.tsv', sep='\t', index_col=0)
anno = pd.read_csv(f'{D}/anno.tsv', sep='\t', index_col=0)
inds = open(f'{D}/individuals.txt').read().split()
variants = pd.read_csv(f'{D}/variants.tsv', sep='\t', index_col=0)
G = np.fromfile(f'{D}/geno.i8', dtype=np.int8).reshape(len(variants), len(inds)).astype(np.float32)
geno_ind = pd.DataFrame(G, index=variants.index, columns=inds)
window = int(open(f'{D}/window.txt').read()) if os.path.exists(f'{D}/window.txt') else 100000   # cis, bp
dcov = pd.read_csv(f'{D}/covariates.tsv', sep='\t', index_col=0) if os.path.exists(f'{D}/covariates.tsv') else None

def inv_norm(x):  # rank-based inverse normal per row
    r = st.rankdata(x, axis=1)
    return st.norm.ppf((r - 0.5) / x.shape[1])

def pcs(P, fixed, k):  # top k PCs of genes x samples P after regressing out the fixed columns
    if k == 0: return None
    X = np.c_[np.ones(P.shape[1]), fixed] if fixed is not None else np.ones((P.shape[1], 1))
    R = P.values.T - X @ np.linalg.lstsq(X, P.values.T, rcond=None)[0]
    u, s, _ = np.linalg.svd(R - R.mean(0), full_matrices=False)
    return pd.DataFrame(u[:, :k] * s[:k], index=P.columns, columns=[f'PC{i+1}' for i in range(k)])

cond = anno['condition']
if pd.api.types.is_numeric_dtype(cond) and cond.nunique() > 3:
    inter = pd.DataFrame({'cond': cond.astype(float)}, index=anno.index)          # numeric condition
elif cond.nunique() > 1:
    inter = pd.get_dummies(cond.astype(str), drop_first=True, dtype=float)         # factor: k - 1 dummies
else:
    inter = None                                                                   # eQTL only

if design == 'sample':
    geno = geno_ind[anno['individual'].values]; geno.columns = anno.index
    pc = pcs(pheno, None if inter is None else inter.values, n_pcs)
    cat = lambda xs: pd.concat(xs, axis=1) if xs else None
    cov_e = cat([x for x in (inter, pc, dcov) if x is not None])
    cov_i = cat([x for x in (pc, dcov) if x is not None])
    jobs = [('geno', pheno, cov_e, None)]
    if inter is not None: jobs.append(('interact', pheno, cov_i, inter))
else:
    assert inter is not None and inter.shape[1] == 1 and set(cond.astype(str).unique()) == {'ctl', 'trt'}
    a = anno.reset_index()
    ctl = a[a.condition == 'ctl'].set_index('individual')['sample']; trt = a[a.condition == 'trt'].set_index('individual')['sample']
    assert ctl.index.equals(trt.index)
    pm = pd.DataFrame(inv_norm((pheno[ctl.values].values + pheno[trt.values].values) / 2), index=pheno.index, columns=ctl.index)
    pd_ = pd.DataFrame(inv_norm(pheno[trt.values].values - pheno[ctl.values].values), index=pheno.index, columns=ctl.index)
    geno = geno_ind[list(ctl.index)]
    jobs = [('geno', pm, pcs(pm, None, n_pcs), None), ('interact', pd_, pcs(pd_, None, n_pcs), None)]
t_load = time.perf_counter() - t0

tmp = out + '_tq'; os.makedirs(tmp, exist_ok=True)
timing, res = {}, {}
for name, P, C, I in jobs:
    if device == 'gpu': torch.cuda.reset_peak_memory_stats()
    t0 = time.perf_counter()
    try:
        cis.map_nominal(geno, variants, P, ppos, name, covariates_df=C, interaction_df=I,
                        maf_threshold=0, maf_threshold_interaction=0, window=window, output_dir=tmp,
                        write_top=False, verbose=False)
    except KeyError as e:
        # tensorQTL 1.0.10 bug: with more than one interaction term, map_nominal fails while summarising the
        # top association per gene ('pval_i' does not exist), after every per-pair file has been written
        if not (I is not None and I.shape[1] > 1 and str(e) == "'pval_i'"): raise
        timing['tensorqtl_top_assoc_bug'] = True
    timing[f't_{name}'] = time.perf_counter() - t0
    if device == 'gpu': timing[f'gpu_mem_{name}_mb'] = torch.cuda.max_memory_allocated() / 2**20
    r = pd.concat([pd.read_parquet(f) for f in sorted(glob.glob(f'{tmp}/{name}.cis_qtl_pairs.*.parquet'))])
    if I is None:
        res[name] = r.set_index(['phenotype_id', 'variant_id'])[['pval_nominal', 'slope', 'slope_se']].rename(
            columns={'pval_nominal': 'p_geno' if name == 'geno' else 'p_interact', 'slope': f'b_{name}', 'slope_se': f'se_{name}'})
    else:
        cols = [c for c in r.columns if c.startswith('pval_g-')] or ['pval_gi']   # one term: not renamed
        x = r.set_index(['phenotype_id', 'variant_id'])[cols]
        x.columns = ['p_interact_' + (c[len('pval_g-'):] if c.startswith('pval_g-') else I.columns[0]) for c in cols]
        # one p-value per pair: the single term, or Bonferroni over the k - 1 terms of a factor
        x['p_interact'] = np.minimum(1, x.min(1) * len(cols))
        res[name] = x
if perm:
    name, P, C, _ = jobs[0]
    t0 = time.perf_counter()
    perm_res = cis.map_cis(geno, variants, P, ppos, covariates_df=C, window=window, maf_threshold=0,
                           seed=1, verbose=False)
    timing['t_map_cis'] = time.perf_counter() - t0
    perm_res.to_csv(out + '_map_cis.tsv.gz', sep='\t')
shutil.rmtree(tmp)
df = pd.concat(res.values(), axis=1, join='outer')
df.index.names = ['gene', 'snp']
df.reset_index().to_csv(out + '.tsv.gz', sep='\t', index=False)
info = dict(tool='tensorqtl', dataset=ds, device=device, threads=threads, design=design, n_pcs=n_pcs,
            n_samples=int(jobs[0][1].shape[1]), n_tests=int(len(df)), t_load=t_load, t_cuda_init=t_cuda_init,
            peak_rss_mb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024,
            torch=torch.__version__, **timing)
json.dump(info, open(out + '.json', 'w'), indent=1)
print(json.dumps(info))
