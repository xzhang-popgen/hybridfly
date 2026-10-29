#!/usr/bin/env bash
set -euo pipefail

bash bash/make_parental_state.sh
bash bash/run_main_sims.sh
bash bash/run_architecture_sims.sh
