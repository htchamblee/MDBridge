#!/bin/bash
# Master script of plan apt1 - written by MDATK <VERSION> on <DATE>.
# Runs the main script of every step in plan order.  Scripts marked 'run separately'
# (e.g. production runs for a cluster) are listed but not started.

set -eo pipefail
cd "$(dirname "$0")"

# 01 · Amber · system build
bash "01_amber-build/apt1_01_amber-build_20261125-1200.sh"

# 02 · Amber · helix restraints
bash "02_amber-restraints/apt1_02_amber-restraints_20261125-1200.sh"

# 03 · Amber · MD stages
bash "03_amber-md/apt1_03_amber-md_20261125-1200.sh"
#   run separately: bash 03_amber-md/apt1_03_amber-prod_20261125-1200.sh

# 04 · VMD · visualisation
bash "04_vmd/apt1_04_vmd_20261125-1200.sh"
