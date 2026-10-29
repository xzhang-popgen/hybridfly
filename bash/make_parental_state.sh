#!/usr/bin/env bash
set -euo pipefail

SLIM=${SLIM:-slim}
mkdir -p results/simulations

$SLIM -d 'STATEFILE="results/simulations/parental_13MB.bin"' \
  slim/make_parental_state_13MB.slim
