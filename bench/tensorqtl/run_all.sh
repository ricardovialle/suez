#!/usr/bin/env bash
# suez vs tensorQTL: all timing and statistics runs (run from the repository root, machine otherwise idle).
#   bash bench/tensorqtl/run_all.sh            # REPS=3 by default
# Needs: pixi (envs in bench/tensorqtl/pixi.toml), suez built into bench/lib/cur, datasets simulated and
# exported (sim_tq.R, export.R). Per-pair outputs go to bench/tensorqtl/out/, run log to out/runs.tsv.
set -uo pipefail
cd "$(git rev-parse --show-toplevel)"
PIXI="pixi run --manifest-path bench/tensorqtl/pixi.toml"
REPS=${REPS:-3}; O=bench/tensorqtl/out; mkdir -p $O/speed $O/stats
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
RUNS=$O/runs.tsv; [ -f $RUNS ] || printf 'kind\tdataset\tmethod\trep\tpeak_rss_kb\tjson\n' > $RUNS
run() {   # run <kind> <dataset> <method> <rep> <out_prefix> <command...>
  local kind=$1 ds=$2 m=$3 rep=$4 out=$5; shift 5
  [ -f "$out.json" ] && return 0
  /usr/bin/time -f '%M' -o "$out.rss" "$@" > "$out.log" 2>&1 || { echo "FAILED $kind $ds $m $rep (see $out.log)"; return 0; }
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$kind" "$ds" "$m" "$rep" "$(tail -1 "$out.rss")" "$(tr -d '\n' < "$out.json")" >> $RUNS
  echo "[$(date +%H:%M:%S)] $kind $ds $m rep$rep: $(grep -o '"t_[a-z0-9_]*": *[0-9.]*' "$out.json" | tr '\n' ' ')"
}
suez()  { run "$1" "$2" "$3" "$4" "$5" $PIXI -e r  Rscript bench/tensorqtl/run_suez.R "$2" "${@:6}" "$5"; }
tq()    { run "$1" "$2" "$3" "$4" "$5" $PIXI -e py python bench/tensorqtl/run_tensorqtl.py "$2" "${@:6:4}" "$5" "${@:10}"; }

SPEED="S M L XL XXL G M_grm L_grm XL_grm M_fac3 M_num M_cov M_eqtl M_unbal"
LADDER=" S M L XL XXL G "
# ---- speed: interleaved repetitions -------------------------------------------------------------------------
for r in $(seq $REPS); do for ds in $SPEED; do
  p=$O/speed/${ds}
  suez speed $ds suez_1core      $r ${p}_suez1_r$r      1 dataset none
  suez speed $ds suez_8cores     $r ${p}_suez8_r$r      8 dataset none
  tq   speed $ds tensorqtl_gpu   $r ${p}_tqgpu_r$r      gpu 1 sample 5
  tq   speed $ds tensorqtl_cpu1  $r ${p}_tqcpu1_r$r     cpu 1 sample 5
  if [ $r = 1 ]; then
    OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 tq speed $ds tensorqtl_cpu8 $r ${p}_tqcpu8_r$r cpu 8 sample 5
    [[ "$LADDER" == *" $ds "* ]] && suez speed $ds suez_1core_boot $r ${p}_suez1boot_r$r 1 dataset boot
  fi
done; done
# ---- realistic scale (REAL: +-1 Mb windows, 500 x 2, GRM), once: rotate-once memory cap 2 GB (default) vs 50 MB
p=$O/speed/REAL
suez speed REAL suez_1core          1 ${p}_suez1_r1   1 dataset none
suez speed REAL suez_8cores         1 ${p}_suez8_r1   8 dataset none
SUEZ_ROTATE_MEM=5e7 suez speed REAL suez_1core_mem5e7  1 ${p}_suez1m_r1 1 dataset none
SUEZ_ROTATE_MEM=5e7 suez speed REAL suez_8cores_mem5e7 1 ${p}_suez8m_r1 8 dataset none
tq   speed REAL tensorqtl_gpu       1 ${p}_tqgpu_r1   gpu 1 sample 5
tq   speed REAL tensorqtl_cpu1      1 ${p}_tqcpu1_r1  cpu 1 sample 5
# ---- gene-level permutations (tensorQTL map_cis, 10,000 permutations), once -------------------------------
for ds in $LADDER; do
  p=$O/speed/${ds}
  tq perm $ds tensorqtl_gpu_mapcis  1 ${p}_tqgpuperm  gpu 1 sample 5 --perm
  tq perm $ds tensorqtl_cpu1_mapcis 1 ${p}_tqcpu1perm cpu 1 sample 5 --perm
done
# ---- statistics ---------------------------------------------------------------------------------------------
for ds in M M_fac3 M_num M_unbal M_cov M_eqtl HS NULLR; do
  p=$O/stats/${ds}
  suez stats $ds suez_dataset_kernel 1 ${p}_suez_dataset 8 dataset none
  suez stats $ds suez_learned_kernel 1 ${p}_suez_learned 8 learned none
  tq   stats $ds tensorqtl_samples_5pc 1 ${p}_tq_sample5 gpu 1 sample 5
  tq   stats $ds tensorqtl_samples_0pc 1 ${p}_tq_sample0 gpu 1 sample 0
  case $ds in M|M_cov|HS|NULLR) tq stats $ds tensorqtl_individual_5pc 1 ${p}_tq_indiv5 gpu 1 indiv 5 ;; esac
done
tq stats M tensorqtl_samples_5pc_cpu 1 $O/stats/M_tq_sample5_cpu cpu 1 sample 5   # CPU vs GPU agreement
echo done
