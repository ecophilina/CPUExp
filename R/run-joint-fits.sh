#!/bin/bash
# Fit the joint models in `model_specs` (R/cpue-index.R) as parallel R
# processes, in phases of at most two: each joint fit needs ~9 GB of memory.
#
# Usage (from the project root): bash R/run-joint-fits.sh
# Logs go to output/logs/fit-<model>.log. Models with a cached fit are skipped.

cd "$(dirname "$0")/.." || exit 1
mkdir -p output/logs

fit() {
  Rscript R/fit-joint.R "$1" > "output/logs/fit-$1.log" 2>&1
  echo "[$(date +%H:%M)] $1 exit $?: $(grep -hE 'FIT DONE|all_ok|Error|already exists' "output/logs/fit-$1.log" | tr '\n' ' ')"
}

echo "[$(date +%H:%M)] phase 1: delta-Gamma"
fit constant_xi_expected_dgamma &
fit zero_xi_expected_dgamma &
wait

echo "[$(date +%H:%M)] phase 2: delta-lognormal"
fit constant_xi_expected_dlnorm &
fit zero_xi_expected_dlnorm &
wait

echo "[$(date +%H:%M)] phase 3: delta-lognormal, rw baseline"
fit constant_rwbase_xi_expected_dlnorm

echo "[$(date +%H:%M)] done"
