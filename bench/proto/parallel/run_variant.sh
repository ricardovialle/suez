#!/bin/bash
# run_variant.sh <impl> <ds> <cores> <rep>: bench.sh with the prototype switches implied by the impl name
# (bench/lib/par_* are symlinks to the one prototype build bench/lib/par_threads)
cd /Users/rvialle/suez
impl=$1
case $impl in
  par_static) export SUEZ_SCHEDULE=static;;
  par_dyn)    export SUEZ_SCHEDULE=dynamic;;
  par_chunk4) export SUEZ_SCHEDULE=chunk SUEZ_CPW=4;;
  par_chunk8) export SUEZ_SCHEDULE=chunk SUEZ_CPW=8;;
  par_thr*)   export SUEZ_SCHEDULE=chunk SUEZ_THREADS=${impl#par_thr};;
esac
bench/bench.sh "$@"
