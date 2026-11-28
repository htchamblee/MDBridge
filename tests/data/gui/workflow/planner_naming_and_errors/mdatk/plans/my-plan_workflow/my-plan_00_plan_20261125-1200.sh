#!/bin/bash
# Master script of plan my-plan - written by MDATK <VERSION> on <DATE>.
# Runs the main script of every step in plan order.  Scripts marked 'run separately'
# (e.g. production runs for a cluster) are listed but not started.

set -eo pipefail
cd "$(dirname "$0")"

# 01 · Amber · system build
#   (not written: the step has errors)

# 02 · Amber · helix restraints
#   (not written: the step has errors)

# 03 · Amber · MD stages
bash "my-plan_03_amber-md_20261125-1200.sh"
#   run separately: bash my-plan_03_amber-prod_20261125-1200.sh
