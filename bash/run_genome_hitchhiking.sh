#!/usr/bin/env bash
set -euo pipefail

SLIM=${SLIM:-slim}
MODEL=slim/fly_hybrid_genome_hitchhiking.slim
OUT=results/simulations/genome_hitchhiking
NREP=${1:-1000}
mkdir -p "$OUT"

run_sims () {
  n=$1
  nsel=$2
  sdel=$3
  tag=N${n}_nsel${nsel}_s${sdel}
  outfile="$OUT/${tag}_${NREP}.csv"

  for rep in $(seq 1 $NREP); do
    seed=$((2000000 + rep + n*10000 + nsel*1000))
    $SLIM -s $seed -d REP=$rep -d N=$n -d NSEL=$nsel -d SDEL=$sdel \
      -d REC_MULT=1.0 -d "SCENARIO=\"$tag\"" "$MODEL" \
      | awk '/^rep,scenario,/ {data=1} data' > "$OUT/tmp.csv"

    if [[ $rep == 1 ]]; then
      cat "$OUT/tmp.csv" > "$outfile"
    else
      tail -n +2 "$OUT/tmp.csv" >> "$outfile"
    fi
  done
  rm -f "$OUT/tmp.csv"
}

for n in 10 20 30 50; do
  run_sims $n 0 -0.20
  run_sims $n 1 -0.20
  run_sims $n 3 -0.20
  run_sims $n 5 -0.20
 done
