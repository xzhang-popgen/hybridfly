#!/usr/bin/env bash
set -euo pipefail

SLIM=${SLIM:-slim}
MODEL=slim/fly_hybrid_3R_fullhistory.slim
STATE=results/simulations/parental_13MB.bin
OUT=results/simulations
mkdir -p "$OUT"

if [[ ! -f "$STATE" ]]; then
  echo "run bash/make_parental_state.sh first"
  exit 1
fi

run_sims () {
  tag=$1
  reps=$2
  seed0=$3
  n=$4
  rec=$5
  sa=$6
  sb=$7
  sc=$8
  sd=$9
  se=${10}
  hneg=${11:-0.5}
  hpos=${12:-0.5}

  tmpdir="$OUT/${tag}_tmp"
  mkdir -p "$tmpdir"

  for rep in $(seq 1 $reps); do
    $SLIM -s $((seed0 + rep)) \
      -d REP=$rep -d N=$n -d NGEN=45 \
      -d CROSS_REC_MULT=1.0 -d REC_MULT=$rec \
      -d H_NEG=$hneg -d H_POS=$hpos \
      -d 'STATEFILE="results/simulations/parental_13MB.bin"' \
      -d "SCENARIO=\"$tag\"" \
      -d S_A=$sa -d S_B=$sb -d S_C=$sc -d S_D=$sd -d S_E=$se \
      "$MODEL" | awk '/^rep,scenario,/ {data=1} data' > "$tmpdir/rep${rep}.csv"
  done

  awk '
    /^rep,scenario,/ {
      if (!header++) print;
      data=1;
      next
    }
    data
  ' "$tmpdir"/rep*.csv > "$OUT/${tag}_${reps}.csv"

  rm -rf "$tmpdir"
}

#neutral drift
run_sims neutral 1000 12000000 50 1.0 0 0 0 0 0

#recombination sensitivity, AE20 model
for rec in 0 0.25 0.5 1 2; do
  tag=$(echo $rec | tr '.' '_')
  run_sims AE20_rec_${tag} 500 15000000 50 $rec -0.20 0.145 0 0.078 -0.20
 done

#population size sensitivity; N=50 is AE20_rec_1_500.csv
for n in 10 20 30; do
  run_sims AE20_N${n} 500 $((17000000 + n*10000)) $n 1.0 -0.20 0.145 0 0.078 -0.20
 done
