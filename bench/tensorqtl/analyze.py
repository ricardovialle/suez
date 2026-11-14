"""Summarise the suez vs tensorQTL runs and write the HTML report.

    pixi run --manifest-path bench/tensorqtl/pixi.toml -e py python bench/tensorqtl/analyze.py

Reads bench/tensorqtl/out/ (runs.tsv and per-pair results), writes summary tables to
bench/tensorqtl/results/ and the report to doc/suez_vs_tensorqtl.html.
"""
import json, os, subprocess, html, math
import numpy as np, pandas as pd, scipy.stats as st

O, R = 'bench/tensorqtl/out', 'bench/tensorqtl/results'
os.makedirs(R, exist_ok=True)
DESC = {  # dataset: (description, N samples shown from runs)
    'S': '100 individuals x 2 conditions', 'M': '200 x 2', 'L': '500 x 2', 'XL': '1,000 x 2',
    'XXL': '2,000 x 2', 'G': '200 x 2, 4,000 genes', 'M_grm': '200 x 2, GRM kernel',
    'L_grm': '500 x 2, GRM kernel', 'XL_grm': '1,000 x 2, GRM kernel', 'M_fac3': '200 x 3 levels',
    'M_num': '200 x 2, numeric condition, GRM', 'M_cov': '200 x 2, 10 covariates',
    'M_eqtl': '200 x 1, no condition, GRM', 'M_unbal': '200, 1-4 samples each', 'HS': '200 x 2, strong effects, GRM',
    'NULLR': '200 x 2, no cis effects, ICC 0.5'}
SPEED = ['S', 'M', 'L', 'XL', 'XXL', 'G', 'M_grm', 'L_grm', 'XL_grm', 'M_fac3', 'M_num', 'M_cov', 'M_eqtl', 'M_unbal']
STATS = ['M', 'M_fac3', 'M_num', 'M_unbal', 'M_cov', 'M_eqtl', 'HS', 'NULLR']

# ---------------------------------------------------------------------------------------------- runs
runs = pd.read_csv(f'{O}/runs.tsv', sep='\t')
rows = []
for _, r in runs.iterrows():
    j = json.loads(r['json'])
    t = j['t_step2'] if j['tool'] == 'suez' else j.get('t_geno', 0) + j.get('t_interact', 0)
    rows.append(dict(kind=r['kind'], dataset=r['dataset'], method=r['method'], rep=r['rep'], time=t,
                     t_geno=j.get('t_geno'), t_interact=j.get('t_interact'), t_map_cis=j.get('t_map_cis'),
                     t_step1=j.get('t_step1'), n_tests=j['n_tests'], n_samples=j['n_samples'],
                     rss_mb=r['peak_rss_kb'] / 1024,
                     gpu_mb=max(j.get('gpu_mem_geno_mb', 0) or 0, j.get('gpu_mem_interact_mb', 0) or 0) or None))
runs = pd.DataFrame(rows)
sp = runs[runs.kind == 'speed']
speed = sp.groupby(['dataset', 'method']).agg(time=('time', 'median'), tmin=('time', 'min'), tmax=('time', 'max'),
                                              reps=('time', 'size'), n_tests=('n_tests', 'first'),
                                              n_samples=('n_samples', 'first'), rss_mb=('rss_mb', 'median'),
                                              gpu_mb=('gpu_mb', 'max'), t_geno=('t_geno', 'median'),
                                              t_interact=('t_interact', 'median')).reset_index()
speed['tests_per_s'] = speed.n_tests / speed.time
speed.to_csv(f'{R}/speed.tsv', sep='\t', index=False, float_format='%.4g')
perm = runs[runs.kind == 'perm'][['dataset', 'method', 'n_samples', 'n_tests', 't_map_cis', 'rss_mb', 'gpu_mb']]
perm.to_csv(f'{R}/permutations.tsv', sep='\t', index=False, float_format='%.4g')

# ---------------------------------------------------------------------------------------------- statistics
def load(path):
    d = pd.read_csv(path, sep='\t')
    d = d.rename(columns={'cis_snp': 'snp'})
    return d.set_index(['gene', 'snp'])

def lam(p):
    p = p.dropna(); p = p[p > 0]
    return np.median(st.chi2.isf(p, 1)) / st.chi2.ppf(0.5, 1) if len(p) else np.nan

METHODS = [('suez_learned', 'suez, step-1 kernel'), ('suez_dataset', 'suez, repeated-individual kernel'),
           ('tq_sample5', 'tensorQTL, all samples, 5 PCs'), ('tq_sample0', 'tensorQTL, all samples, no PCs'),
           ('tq_indiv5', 'tensorQTL, per individual (mean / difference), 5 PCs')]
cal_rows, pow_rows, qq = [], [], {}
for ds in STATS:
    truth = pd.read_csv(f'bench/tensorqtl/data/{ds}/truth.tsv', sep='\t')
    if ds == 'NULLR': truth = truth.iloc[0:0]
    eg, rg = set(truth.gene[truth.type == 'eQTL']), set(truth.gene[truth.type == 'response'])
    base = load(f'{O}/stats/{ds}_suez_dataset.tsv.gz').index      # common pair universe: pairs suez tested
    for key, label in METHODS:
        f = f'{O}/stats/{ds}_{key}.tsv.gz'.replace('_suez_', '_suez_')
        if not os.path.exists(f): continue
        d = load(f).reindex(base)
        genes = d.index.get_level_values(0)
        null = ~genes.isin(eg | rg)
        rec = dict(dataset=ds, method=label, key=key, n_null_pairs=int(null.sum()))
        for test in ['p_geno', 'p_interact']:
            if test not in d or d[test].isna().all(): continue
            p = d[test]
            pn = p[null]
            rec[f'lambda_{test}'] = lam(pn)
            rec[f'fpr05_{test}'] = (pn < 0.05).mean()
            rec[f'fpr1e3_{test}'] = (pn < 1e-3).mean()
            gp = p.groupby(level=0).agg(lambda x: min(1, x.min() * x.notna().sum()))   # Bonferroni within gene
            gnull = gp[~gp.index.isin(eg | rg)]
            target = eg if test == 'p_geno' else rg
            gt = gp[gp.index.isin(target)]
            thr = np.quantile(gnull.dropna(), 0.05) if len(gnull) else np.nan          # 5% of null genes
            pow_rows.append(dict(dataset=ds, method=label, key=key, test=test, n_causal_genes=len(gt),
                                 n_null_genes=len(gnull), false_gene_rate=(gnull < 0.05).mean(),
                                 power=(gt < 0.05).mean() if len(gt) else np.nan,
                                 power_calibrated=(gt <= thr).mean() if len(gt) else np.nan))
            if ds == 'NULLR': qq[(test, key)] = np.sort(pn.dropna().values)
        cal_rows.append(rec)
cal = pd.DataFrame(cal_rows); pw = pd.DataFrame(pow_rows)
cal.to_csv(f'{R}/calibration.tsv', sep='\t', index=False, float_format='%.4g')
pw.to_csv(f'{R}/power.tsv', sep='\t', index=False, float_format='%.4g')

# concordance on M: suez repeated-individual kernel vs tensorQTL; GPU vs CPU
conc = []
a = load(f'{O}/stats/M_suez_dataset.tsv.gz')
for key in ['suez_learned', 'tq_sample5', 'tq_sample0', 'tq_indiv5']:
    b = load(f'{O}/stats/M_{key}.tsv.gz').reindex(a.index)
    for test in ['p_geno', 'p_interact']:
        x, y = -np.log10(a[test]), -np.log10(b[test]); ok = x.notna() & y.notna()
        conc.append(dict(a='suez_dataset', b=key, test=test, spearman=st.spearmanr(x[ok], y[ok])[0],
                         pearson_log10p=np.corrcoef(x[ok], y[ok])[0, 1]))
g, c = load(f'{O}/stats/M_tq_sample5.tsv.gz'), load(f'{O}/stats/M_tq_sample5_cpu.tsv.gz').reindex(load(f'{O}/stats/M_tq_sample5.tsv.gz').index)
for test in ['p_geno', 'p_interact']:
    dd = (-np.log10(g[test]) + np.log10(c[test])).abs()
    conc.append(dict(a='tq_sample5_gpu', b='tq_sample5_cpu', test=test, spearman=st.spearmanr(g[test], c[test], nan_policy='omit')[0],
                     pearson_log10p=np.nan, max_abs_diff_log10p=dd.max()))
conc = pd.DataFrame(conc); conc.to_csv(f'{R}/concordance.tsv', sep='\t', index=False, float_format='%.5g')
step1 = runs[(runs.kind == 'stats') & (runs.method == 'suez_learned_kernel')][['dataset', 'n_samples', 't_step1', 'time']]

# ---------------------------------------------------------------------------------------------- environment
def sh(c):
    try: return subprocess.run(c, shell=True, capture_output=True, text=True, timeout=60).stdout.strip()
    except Exception: return ''
env = dict(
    cpu=sh("lscpu | grep -m1 'Model name' | sed 's/.*: *//'"), cores=sh('nproc'),
    ram=sh("free -g | awk '/Mem:/{print $2}'") + ' GB',
    gpu=sh("nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader"),
    os=sh("grep PRETTY_NAME /etc/os-release | cut -d'\"' -f2") + ' (WSL2, kernel ' + sh('uname -r') + ')',
    suez=sh("git rev-parse --short HEAD") + ' (' + sh("grep ^Version DESCRIPTION | cut -d' ' -f2") + ')',
    r=sh("pixi run --manifest-path bench/tensorqtl/pixi.toml -e r Rscript -e 'cat(R.version.string, basename(La_library()))' 2>/dev/null"),
    tq=sh("pixi run --manifest-path bench/tensorqtl/pixi.toml -e py python -c \"import importlib.metadata as m, torch; print('tensorQTL', m.version('tensorqtl'), '/ PyTorch', torch.__version__, '/ CUDA', torch.version.cuda)\" 2>/dev/null | tail -1"),
    date=sh('date +%Y-%m-%d'))
open(f'{R}/environment.txt', 'w').write('\n'.join(f'{k}: {v}' for k, v in env.items()) + '\n')

# ---------------------------------------------------------------------------------------------- HTML helpers
E = html.escape
def fmt_t(x):
    if x is None or (isinstance(x, float) and math.isnan(x)): return '-'
    return f'{x:.2f}' if x < 10 else f'{x:.1f}' if x < 100 else f'{x:,.0f}'
def fmt_n(x): return '-' if pd.isna(x) else f'{int(x):,}'
def fmt_p(x, d=3): return '-' if pd.isna(x) else f'{x:.{d}f}'
def fmt_rate(x): return '-' if pd.isna(x) else (f'{x:.1e}' if 0 < x < 1e-3 else f'{x:.4f}' if x < 0.01 else f'{x:.3f}')

SPEED_M = [('suez_1core', 'suez, 1 core'), ('suez_8cores', 'suez, 8 cores'),
           ('tensorqtl_gpu', 'tensorQTL, GPU'), ('tensorqtl_cpu1', 'tensorQTL, CPU 1 thread'),
           ('tensorqtl_cpu8', 'tensorQTL, CPU 8 threads')]
# marks: (css color var, shape, filled)
MARK = {'suez_1core': ('--s1', 'circle', False), 'suez_8cores': ('--s1', 'circle', True),
        'tensorqtl_gpu': ('--s2', 'square', True), 'tensorqtl_cpu1': ('--s3', 'triangle', True)}

def marker(shape, x, y, var, filled, r=5):
    fill = f'var({var})' if filled else 'var(--surface)'
    common = f'fill="{fill}" stroke="var({var})" stroke-width="2"'
    if shape == 'circle': return f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" {common}/>'
    if shape == 'square': return f'<rect x="{x-r+0.5:.1f}" y="{y-r+0.5:.1f}" width="{2*r-1}" height="{2*r-1}" rx="1.5" {common}/>'
    h = r * 1.15
    return f'<path d="M{x:.1f},{y-h:.1f} L{x+h:.1f},{y+h*0.8:.1f} L{x-h:.1f},{y+h*0.8:.1f} Z" {common} stroke-linejoin="round"/>'

def speed_svg():
    W, left, right, top, rowh = 760, 200, 24, 34, 30
    dss = [d for d in SPEED if d in set(speed.dataset)]
    H = top + rowh * len(dss) + 40
    lo, hi = 0.05, 2000
    X = lambda v: left + (math.log10(v) - math.log10(lo)) / (math.log10(hi) - math.log10(lo)) * (W - left - right)
    s = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-labelledby="spt" class="chart"><title id="spt">Compute time per dataset and method, log scale</title>']
    for v in [0.1, 1, 10, 100, 1000]:
        s.append(f'<line x1="{X(v):.1f}" x2="{X(v):.1f}" y1="{top-8}" y2="{H-34}" class="grid"/>'
                 f'<text x="{X(v):.1f}" y="{H-18}" class="tick" text-anchor="middle">{v:g} s</text>')
    for i, ds in enumerate(dss):
        y = top + i * rowh + rowh / 2
        n = speed[speed.dataset == ds].n_tests.max(); N = speed[speed.dataset == ds].n_samples.max()
        s.append(f'<text x="{left-12}" y="{y-2:.1f}" class="lab" text-anchor="end">{E(ds)}</text>'
                 f'<text x="{left-12}" y="{y+11:.1f}" class="sub" text-anchor="end">N={N:,} · {n/1000:,.0f}k tests</text>')
        s.append(f'<line x1="{left}" x2="{W-right}" y1="{y:.1f}" y2="{y:.1f}" class="row"/>')
        for m, lab in SPEED_M:
            if m not in MARK: continue
            q = speed[(speed.dataset == ds) & (speed.method == m)]
            if not len(q): continue
            q = q.iloc[0]; var, shp, fil = MARK[m]
            s.append(f'<g class="mk"><title>{E(ds)} · {E(lab)}: {fmt_t(q.time)} s (median of {q.reps}), '
                     f'{q.tests_per_s:,.0f} tests/s</title>'
                     f'<circle cx="{X(q.time):.1f}" cy="{y:.1f}" r="11" fill="transparent"/>'
                     f'{marker(shp, X(q.time), y, var, fil)}</g>')
    s.append(f'<text x="{(left+W-right)/2:.0f}" y="{H-2}" class="tick" text-anchor="middle">compute time, seconds (log scale; both tests: eQTL + G x condition)</text>')
    s.append('</svg>')
    return ''.join(s)

def legend(items):
    out = ['<div class="legend">']
    for var, shp, fil, lab in items:
        out.append(f'<span><svg width="16" height="16" viewBox="0 0 16 16">{marker(shp, 8, 8.5, var, fil, 5)}</svg>{E(lab)}</span>')
    return ''.join(out) + '</div>'

QQM = [('suez_learned', '--s1', '', 'suez, step-1 kernel'), ('suez_dataset', '--s1', '5 4', 'suez, repeated-individual kernel'),
       ('tq_sample5', '--s2', '', 'tensorQTL, all samples, 5 PCs'), ('tq_indiv5', '--s3', '', 'tensorQTL, per individual, 5 PCs')]
def qq_svg(test, title):
    W, H, l, b, t, r = 360, 330, 46, 40, 26, 12
    mx_e = max(-np.log10((0.5) / len(v)) for (tt, k), v in qq.items() if tt == test)
    mx_o = max(-np.log10(v[0]) for (tt, k), v in qq.items() if tt == test and len(v))
    mx = math.ceil(max(mx_e, mx_o) + 0.2)
    X = lambda v: l + v / mx * (W - l - r); Y = lambda v: H - b - v / mx * (H - b - t)
    s = [f'<svg viewBox="0 0 {W} {H}" role="img" class="chart qq"><title>{E(title)}</title>',
         f'<text x="{l}" y="14" class="lab">{E(title)}</text>']
    for v in range(0, mx + 1):
        s.append(f'<line x1="{X(v):.1f}" x2="{X(v):.1f}" y1="{t}" y2="{H-b}" class="grid"/><line x1="{l}" x2="{W-r}" y1="{Y(v):.1f}" y2="{Y(v):.1f}" class="grid"/>'
                 f'<text x="{X(v):.1f}" y="{H-b+15}" class="tick" text-anchor="middle">{v}</text><text x="{l-6}" y="{Y(v)+4:.1f}" class="tick" text-anchor="end">{v}</text>')
    s.append(f'<line x1="{X(0)}" y1="{Y(0)}" x2="{X(mx)}" y2="{Y(mx)}" class="diag"/>')
    for key, var, dash, lab in QQM:
        v = qq.get((test, key))
        if v is None or not len(v): continue
        n = len(v); e = -np.log10((np.arange(1, n + 1) - 0.5) / n); o = -np.log10(np.maximum(v, 1e-300))
        idx = np.unique(np.r_[np.arange(min(n, 300)), np.geomspace(1, n, 400).astype(int) - 1])
        pts = ' '.join(f'{X(e[i]):.1f},{Y(min(o[i], mx)):.1f}' for i in idx)
        da = f' stroke-dasharray="{dash}"' if dash else ''
        s.append(f'<g class="ln"><title>{E(lab)}: λ = {lam(pd.Series(v)):.2f}; smallest null p = {v[0]:.1e}</title>'
                 f'<polyline points="{pts}" fill="none" stroke="var({var})" stroke-width="2"{da} stroke-linejoin="round"/></g>')
    s.append(f'<text x="{(l+W-r)/2:.0f}" y="{H-6}" class="tick" text-anchor="middle">expected −log10 p</text>'
             f'<text transform="translate(12,{(t+H-b)/2:.0f}) rotate(-90)" class="tick" text-anchor="middle">observed −log10 p</text></svg>')
    return ''.join(s)

def table(df, cols, cls=''):
    h = ['<div class="tw"><table class="' + cls + '"><thead><tr>' + ''.join(f'<th>{E(c[0])}</th>' for c in cols) + '</tr></thead><tbody>']
    for _, r in df.iterrows():
        h.append('<tr>' + ''.join(f'<td>{c[1](r)}</td>' for c in cols) + '</tr>')
    return ''.join(h) + '</tbody></table></div>'

# speed table (wide)
def sget(ds, m, col='time'):
    q = speed[(speed.dataset == ds) & (speed.method == m)]
    return q.iloc[0][col] if len(q) else np.nan
wide = []
for ds in SPEED:
    if ds not in set(speed.dataset): continue
    r = dict(dataset=ds, desc=DESC[ds], N=speed[speed.dataset == ds].n_samples.max(),
             tests=sget(ds, 'suez_1core', 'n_tests'), tq_tests=sget(ds, 'tensorqtl_gpu', 'n_tests'))
    for m, _ in SPEED_M: r[m] = sget(ds, m)
    r['boot'] = sget(ds, 'suez_1core_boot')
    r['rss_suez'] = sget(ds, 'suez_1core', 'rss_mb'); r['rss_gpu'] = sget(ds, 'tensorqtl_gpu', 'rss_mb')
    r['gpu_mb'] = sget(ds, 'tensorqtl_gpu', 'gpu_mb')
    best_tq = np.nanmin([r['tensorqtl_gpu'], r['tensorqtl_cpu1'], r['tensorqtl_cpu8']])
    r['ratio1'] = r['suez_1core'] / best_tq; r['ratio8'] = r['suez_8cores'] / best_tq
    wide.append(r)
wide = pd.DataFrame(wide); wide.to_csv(f'{R}/speed_wide.tsv', sep='\t', index=False, float_format='%.4g')

def ratio_txt(x): return '-' if pd.isna(x) else (f'{x:.1f}x slower' if x >= 1.05 else f'{1/x:.1f}x faster' if x <= 0.95 else 'same')
speed_tbl = table(wide, [
    ('dataset', lambda r: f'<b>{E(r.dataset)}</b><br><span class="mut">{E(r.desc)}</span>'),
    ('N', lambda r: fmt_n(r.N)), ('tests (suez / tensorQTL)', lambda r: f'{fmt_n(r.tests)} / {fmt_n(r.tq_tests)}'),
    ('suez 1 core', lambda r: fmt_t(r.suez_1core)), ('suez 8 cores', lambda r: fmt_t(r.suez_8cores)),
    ('suez 1 core + bootstrap', lambda r: fmt_t(r.boot)),
    ('tensorQTL GPU', lambda r: fmt_t(r.tensorqtl_gpu)), ('tensorQTL CPU 1 thread', lambda r: fmt_t(r.tensorqtl_cpu1)),
    ('tensorQTL CPU 8 threads', lambda r: fmt_t(r.tensorqtl_cpu8)),
    ('suez 8 cores vs fastest tensorQTL', lambda r: ratio_txt(r.ratio8)),
    ('peak RSS MB (suez 1 core / tensorQTL GPU)', lambda r: f'{r.rss_suez:,.0f} / {r.rss_gpu:,.0f}'),
    ('GPU MB', lambda r: '-' if pd.isna(r.gpu_mb) else f'{r.gpu_mb:,.0f}')], 'num')

pm = perm.pivot_table(index='dataset', columns='method', values='t_map_cis').reindex([d for d in SPEED if d in set(perm.dataset)])
perm_tbl = table(pm.reset_index(), [('dataset', lambda r: f'<b>{E(r.dataset)}</b>'),
    ('genes', lambda r: fmt_n(len(pd.read_csv(f"bench/tensorqtl/data/{r.dataset}/pheno_pos.tsv", sep="\t")))),
    ('tensorQTL map_cis, GPU (s)', lambda r: fmt_t(r.get('tensorqtl_gpu_mapcis'))),
    ('tensorQTL map_cis, CPU 1 thread (s)', lambda r: fmt_t(r.get('tensorqtl_cpu1_mapcis')))], 'num')

LAB = dict(METHODS)
cal_tbl = table(cal, [
    ('dataset', lambda r: f'<b>{E(r.dataset)}</b>'), ('method', lambda r: E(r.method)),
    ('λ eQTL', lambda r: fmt_p(r.get('lambda_p_geno'), 2)), ('P(p<0.05) eQTL', lambda r: fmt_rate(r.get('fpr05_p_geno'))),
    ('P(p<1e-3) eQTL', lambda r: fmt_rate(r.get('fpr1e3_p_geno'))),
    ('λ response', lambda r: fmt_p(r.get('lambda_p_interact'), 2)), ('P(p<0.05) response', lambda r: fmt_rate(r.get('fpr05_p_interact'))),
    ('P(p<1e-3) response', lambda r: fmt_rate(r.get('fpr1e3_p_interact')))], 'num cal')
pw2 = pw[pw.dataset != 'NULLR']
def pw_cell(ds, key, test, col):
    q = pw[(pw.dataset == ds) & (pw.key == key) & (pw.test == test)]
    return np.nan if not len(q) else q.iloc[0][col]
prow = []
for ds in [d for d in STATS if d != 'NULLR']:
    for key, label in METHODS:
        if not len(pw[(pw.dataset == ds) & (pw.key == key)]): continue
        prow.append(dict(dataset=ds, method=label, **{f'{t}_{c}': pw_cell(ds, key, t, c) for t in ['p_geno', 'p_interact']
                                                     for c in ['power', 'power_calibrated', 'false_gene_rate']}))
prow = pd.DataFrame(prow)
pow_tbl = table(prow, [
    ('dataset', lambda r: f'<b>{E(r.dataset)}</b>'), ('method', lambda r: E(r.method)),
    ('eQTL genes found', lambda r: fmt_p(r.p_geno_power, 2)), ('at calibrated threshold', lambda r: fmt_p(r.p_geno_power_calibrated, 2)),
    ('null genes called', lambda r: fmt_p(r.p_geno_false_gene_rate, 3)),
    ('response genes found', lambda r: fmt_p(r.get('p_interact_power'), 2)), ('at calibrated threshold', lambda r: fmt_p(r.get('p_interact_power_calibrated'), 2)),
    ('null genes called', lambda r: fmt_p(r.get('p_interact_false_gene_rate'), 3))], 'num pow')
conc_tbl = table(conc, [('comparison', lambda r: f'{E(r.a)} vs {E(r.b)}'), ('test', lambda r: E(r.test)),
                        ('Spearman ρ of p', lambda r: fmt_p(r.spearman, 4)),
                        ('max |Δ −log10 p|', lambda r: '-' if pd.isna(r.get('max_abs_diff_log10p')) else f'{r.max_abs_diff_log10p:.1e}')], 'num')
st1_tbl = table(step1, [('dataset', lambda r: f'<b>{E(r.dataset)}</b>'), ('N', lambda r: fmt_n(r.n_samples)),
                        ('step 1, learn_covariance P = 10 (s)', lambda r: fmt_t(r.t_step1)),
                        ('step 2, 8 cores (s)', lambda r: fmt_t(r.time))], 'num')

# headline numbers
def W(ds, col): q = wide[wide.dataset == ds]; return q.iloc[0][col] if len(q) else np.nan
def C(ds, key, col): q = cal[(cal.dataset == ds) & (cal.key == key)]; return q.iloc[0].get(col, np.nan) if len(q) else np.nan

body = open('bench/tensorqtl/report_template.html').read()
subs = dict(
    ENV_CPU=E(env['cpu']), ENV_CORES=E(env['cores']), ENV_RAM=E(env['ram']), ENV_GPU=E(env['gpu']), ENV_OS=E(env['os']),
    ENV_SUEZ=E(env['suez']), ENV_R=E(env['r']), ENV_TQ=E(env['tq']), ENV_DATE=E(env['date']),
    SPEED_SVG=speed_svg(),
    SPEED_LEGEND=legend([('--s1', 'circle', False, 'suez, 1 core'), ('--s1', 'circle', True, 'suez, 8 cores'),
                         ('--s2', 'square', True, 'tensorQTL, GPU'), ('--s3', 'triangle', True, 'tensorQTL, CPU 1 thread')]),
    SPEED_TABLE=speed_tbl, PERM_TABLE=perm_tbl, CAL_TABLE=cal_tbl, POW_TABLE=pow_tbl, CONC_TABLE=conc_tbl, STEP1_TABLE=st1_tbl,
    QQ_GENO=qq_svg('p_geno', 'eQTL test, null data (NULLR)'), QQ_INT=qq_svg('p_interact', 'Response test, null data (NULLR)'),
    QQ_LEGEND=''.join(f'<span><svg width="26" height="10"><line x1="1" y1="5" x2="25" y2="5" stroke="var({v})" stroke-width="2"'
                      f'{" stroke-dasharray=%s" % chr(34) + d + chr(34) if d else ""}/></svg>{E(l)}</span>' for _, v, d, l in QQM),
)
for k, v in subs.items():
    body = body.replace('{{' + k + '}}', str(v))
# numeric placeholders: {{num:expr}}
import re
def PM(ds, m):
    q = perm[(perm.dataset == ds) & (perm.method == m)]
    return q.t_map_cis.iloc[0] if len(q) else np.nan
def ev(m):
    try: return str(eval(m.group(1), {'W': W, 'C': C, 'PM': PM, 'fmt_t': fmt_t, 'fmt_p': fmt_p, 'np': np, 'wide': wide, 'pw_cell': pw_cell, 'speed': speed}))
    except Exception as e: return f'[{e}]'
body = re.sub(r'\{\{num:(.*?)\}\}', ev, body)
open('doc/suez_vs_tensorqtl.html', 'w').write(body)
print('wrote doc/suez_vs_tensorqtl.html')
