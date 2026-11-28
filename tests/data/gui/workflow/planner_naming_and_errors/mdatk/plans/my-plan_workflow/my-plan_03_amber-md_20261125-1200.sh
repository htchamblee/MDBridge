#!/bin/bash
# Amber MD stages
# Project my-plan, step 03 - written by MDATK <VERSION> on <DATE>.
# MDATK only writes this script; review it, then run it with: bash my-plan_03_amber-md_20261125-1200.sh
# Amber 26; engine pmemd.cuda.
# Please cite (AmberTools, pmemd (Amber), CPPTRAJ):
#   - Case D. A., Aktulga H. M., Belfon K., Cerutti D. S., Cisneros G. A., Cruzeiro V. W. D.,
#     Forouzesh N., Giese T. J., Götz A. W., Gohlke H., et al. J. Chem. Inf. Model. 2023, 63,
#     6183-6191. doi:10.1021/acs.jcim.3c01153
#   - Salomon-Ferrer R., Götz A. W., Poole D., Le Grand S., Walker R. C. J. Chem. Theory Comput.
#     2013, 9, 3878-3888. doi:10.1021/ct400314y
#   - Götz A. W., Williamson M. J., Xu D., Poole D., Le Grand S., Walker R. C. J. Chem. Theory
#     Comput. 2012, 8, 1542-1555. doi:10.1021/ct200909j
#   - Roe D. R., Cheatham T. E. J. Chem. Theory Comput. 2013, 9, 3084-3095. doi:10.1021/ct400341p

set -eo pipefail
cd "$(dirname "$0")"

PMEMD="pmemd.cuda"
TOP="my-plan_01_parmed_20261125-1200.parm7"

# a stage counts as finished when its output has the final timing lines
# (pmemd: "Total wall time", sander: "Run   done at")
finished() { [ -s "$1" ] && grep -qE "Total wall time|Run +done at" "$1"; }

# min1: Minimize water/ions, RNA restrained
if finished my-plan_03_min1_20261125-1200.out; then
  echo "min1: already finished - skipping"
else
  $PMEMD -O -i my-plan_03_min1_20261125-1200.in -p "$TOP" -c my-plan_01_tleap_20261125-1200.rst7 -ref my-plan_01_tleap_20261125-1200.rst7 -o my-plan_03_min1_20261125-1200.out -r my-plan_03_min1_20261125-1200.rst7 -inf my-plan_03_min1_20261125-1200.mdinfo
fi
[ -s "my-plan_03_min1_20261125-1200.rst7" ] || { echo "ERROR: pmemd (min1) did not write my-plan_03_min1_20261125-1200.rst7" >&2; exit 1; }

# min2: Unrestrained min, helix pairs held
if finished my-plan_03_min2_20261125-1200.out; then
  echo "min2: already finished - skipping"
else
  $PMEMD -O -i my-plan_03_min2_20261125-1200.in -p "$TOP" -c my-plan_03_min1_20261125-1200.rst7 -ref my-plan_03_min1_20261125-1200.rst7 -o my-plan_03_min2_20261125-1200.out -r my-plan_03_min2_20261125-1200.rst7 -inf my-plan_03_min2_20261125-1200.mdinfo
fi
[ -s "my-plan_03_min2_20261125-1200.rst7" ] || { echo "ERROR: pmemd (min2) did not write my-plan_03_min2_20261125-1200.rst7" >&2; exit 1; }

# heat: Heat 0 -> 300 K, helix pairs held
if finished my-plan_03_heat_20261125-1200.out; then
  echo "heat: already finished - skipping"
else
  $PMEMD -O -i my-plan_03_heat_20261125-1200.in -p "$TOP" -c my-plan_03_min2_20261125-1200.rst7 -ref my-plan_03_min2_20261125-1200.rst7 -o my-plan_03_heat_20261125-1200.out -r my-plan_03_heat_20261125-1200.rst7 -x my-plan_03_heat_20261125-1200.nc -inf my-plan_03_heat_20261125-1200.mdinfo
fi
[ -s "my-plan_03_heat_20261125-1200.rst7" ] || { echo "ERROR: pmemd (heat) did not write my-plan_03_heat_20261125-1200.rst7" >&2; exit 1; }

# eq: NPT eq, light restraint
if finished my-plan_03_eq_20261125-1200.out; then
  echo "eq: already finished - skipping"
else
  $PMEMD -O -i my-plan_03_eq_20261125-1200.in -p "$TOP" -c my-plan_03_heat_20261125-1200.rst7 -ref my-plan_03_heat_20261125-1200.rst7 -o my-plan_03_eq_20261125-1200.out -r my-plan_03_eq_20261125-1200.rst7 -x my-plan_03_eq_20261125-1200.nc -inf my-plan_03_eq_20261125-1200.mdinfo
fi
[ -s "my-plan_03_eq_20261125-1200.rst7" ] || { echo "ERROR: pmemd (eq) did not write my-plan_03_eq_20261125-1200.rst7" >&2; exit 1; }

# eq-free: Unrestrained NPT eq
if finished my-plan_03_eq-free_20261125-1200.out; then
  echo "eq-free: already finished - skipping"
else
  $PMEMD -O -i my-plan_03_eq-free_20261125-1200.in -p "$TOP" -c my-plan_03_eq_20261125-1200.rst7 -ref my-plan_03_eq_20261125-1200.rst7 -o my-plan_03_eq-free_20261125-1200.out -r my-plan_03_eq-free_20261125-1200.rst7 -x my-plan_03_eq-free_20261125-1200.nc -inf my-plan_03_eq-free_20261125-1200.mdinfo
fi
[ -s "my-plan_03_eq-free_20261125-1200.rst7" ] || { echo "ERROR: pmemd (eq-free) did not write my-plan_03_eq-free_20261125-1200.rst7" >&2; exit 1; }

# Strip a dry PDB of eq-free for visual checks
cpptraj << EOF
parm "$TOP"
trajin my-plan_03_eq-free_20261125-1200.rst7
strip :WAT,HOH,Na+,Cl-,K+,MG
trajout my-plan_03_dry_20261125-1200.pdb pdb
run
EOF
[ -s "my-plan_03_dry_20261125-1200.pdb" ] || { echo "ERROR: cpptraj did not write my-plan_03_dry_20261125-1200.pdb" >&2; exit 1; }
