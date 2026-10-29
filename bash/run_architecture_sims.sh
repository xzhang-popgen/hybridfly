#!/usr/bin/env bash
set -euo pipefail

SLIM=${SLIM:-slim}
MODEL=slim/fly_hybrid_3R_fullhistory.slim
STATE=results/simulations/parental_13MB.bin
OUT=results/simulations/architecture
mkdir -p "$OUT"

if [[ ! -f "$STATE" ]]; then
  echo "run bash/make_parental_state.sh first"
  exit 1
fi

run_sims () {
  tag=$1
  reps=$2
  seed0=$3
  sa=$4
  sb=$5
  sc=$6
  sd=$7
  se=$8
  hneg=${9:-0.5}
  hpos=${10:-0.5}

  tmpdir="$OUT/${tag}_tmp"
  mkdir -p "$tmpdir"

  for rep in $(seq 1 $reps); do
    $SLIM -s $((seed0 + rep)) \
      -d REP=$rep -d N=50 -d NGEN=45 \
      -d CROSS_REC_MULT=1.0 -d REC_MULT=1.0 \
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

#initial selected model
run_sims selected_original 1000 13000000 -0.14 0.07 -0.20 0.04 -0.14

#dominance diagnostic: beneficial alleles dominant, deleterious alleles additive
run_sims selected_hpos1 500 13500000 -0.14 0.07 -0.20 0.04 -0.14 0.5 1.0

#positive effects converted to diploid SLiM coefficients
run_sims selected_corrected 500 14000000 -0.14 0.145 -0.20 0.078 -0.14

#one linked deleterious locus
run_sims oneC 500 14100000 0 0.145 -0.20 0.078 0

#two-deleterious-locus arrangements
run_sims AC 300 14200000 -0.14 0.145 -0.20 0.078 0
run_sims CE 300 14300000 0 0.145 -0.20 0.078 -0.14
run_sims AE 300 14400000 -0.14 0.145 0 0.078 -0.14

#AE linked-load sensitivity
run_sims AE18 300 14500000 -0.18 0.145 0 0.078 -0.18
run_sims AE20 300 14600000 -0.20 0.145 0 0.078 -0.20
run_sims AE22 300 14700000 -0.22 0.145 0 0.078 -0.22
run_sims AE24 300 14800000 -0.24 0.145 0 0.078 -0.24
