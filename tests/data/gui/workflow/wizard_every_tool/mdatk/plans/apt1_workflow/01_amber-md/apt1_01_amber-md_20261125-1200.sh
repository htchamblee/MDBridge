#!/bin/bash
# Amber MD stages
# Project apt1, step 01 - written by MDATK <VERSION> on <DATE>.
# MDATK only writes this script; review it, then run it with: bash apt1_01_amber-md_20261125-1200.sh
# Amber 26; engine pmemd.MPI.
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

PMEMD="mpirun -np 8 pmemd.MPI"
TOP="/x/sys.parm7"

# a stage counts as finished when its output has the final timing lines
# (pmemd: "Total wall time", sander: "Run   done at")
finished() { [ -s "$1" ] && grep -qE "Total wall time|Run +done at" "$1"; }

# min1: Minimize water/ions, RNA restrained
if finished apt1_01_min1_20261125-1200.out; then
  echo "min1: already finished - skipping"
else
  $PMEMD -O -i apt1_01_min1_20261125-1200.in -p "$TOP" -c /x/sys.rst7 -ref /x/sys.rst7 -o apt1_01_min1_20261125-1200.out -r apt1_01_min1_20261125-1200.rst7 -inf apt1_01_min1_20261125-1200.mdinfo
fi
[ -s "apt1_01_min1_20261125-1200.rst7" ] || { echo "ERROR: pmemd (min1) did not write apt1_01_min1_20261125-1200.rst7" >&2; exit 1; }

# min2: Unrestrained min, helix pairs held
if finished apt1_01_min2_20261125-1200.out; then
  echo "min2: already finished - skipping"
else
  $PMEMD -O -i apt1_01_min2_20261125-1200.in -p "$TOP" -c apt1_01_min1_20261125-1200.rst7 -ref apt1_01_min1_20261125-1200.rst7 -o apt1_01_min2_20261125-1200.out -r apt1_01_min2_20261125-1200.rst7 -inf apt1_01_min2_20261125-1200.mdinfo
fi
[ -s "apt1_01_min2_20261125-1200.rst7" ] || { echo "ERROR: pmemd (min2) did not write apt1_01_min2_20261125-1200.rst7" >&2; exit 1; }

# heat: Heat 0 -> 300 K, helix pairs held
if finished apt1_01_heat_20261125-1200.out; then
  echo "heat: already finished - skipping"
else
  $PMEMD -O -i apt1_01_heat_20261125-1200.in -p "$TOP" -c apt1_01_min2_20261125-1200.rst7 -ref apt1_01_min2_20261125-1200.rst7 -o apt1_01_heat_20261125-1200.out -r apt1_01_heat_20261125-1200.rst7 -x apt1_01_heat_20261125-1200.nc -inf apt1_01_heat_20261125-1200.mdinfo
fi
[ -s "apt1_01_heat_20261125-1200.rst7" ] || { echo "ERROR: pmemd (heat) did not write apt1_01_heat_20261125-1200.rst7" >&2; exit 1; }

# eq: NPT eq, light restraint
if finished apt1_01_eq_20261125-1200.out; then
  echo "eq: already finished - skipping"
else
  $PMEMD -O -i apt1_01_eq_20261125-1200.in -p "$TOP" -c apt1_01_heat_20261125-1200.rst7 -ref apt1_01_heat_20261125-1200.rst7 -o apt1_01_eq_20261125-1200.out -r apt1_01_eq_20261125-1200.rst7 -x apt1_01_eq_20261125-1200.nc -inf apt1_01_eq_20261125-1200.mdinfo
fi
[ -s "apt1_01_eq_20261125-1200.rst7" ] || { echo "ERROR: pmemd (eq) did not write apt1_01_eq_20261125-1200.rst7" >&2; exit 1; }

# eq-free: Unrestrained NPT eq
if finished apt1_01_eq-free_20261125-1200.out; then
  echo "eq-free: already finished - skipping"
else
  $PMEMD -O -i apt1_01_eq-free_20261125-1200.in -p "$TOP" -c apt1_01_eq_20261125-1200.rst7 -ref apt1_01_eq_20261125-1200.rst7 -o apt1_01_eq-free_20261125-1200.out -r apt1_01_eq-free_20261125-1200.rst7 -x apt1_01_eq-free_20261125-1200.nc -inf apt1_01_eq-free_20261125-1200.mdinfo
fi
[ -s "apt1_01_eq-free_20261125-1200.rst7" ] || { echo "ERROR: pmemd (eq-free) did not write apt1_01_eq-free_20261125-1200.rst7" >&2; exit 1; }

# Strip a dry PDB of eq-free for visual checks
cpptraj << EOF
parm "$TOP"
trajin apt1_01_eq-free_20261125-1200.rst7
strip :WAT,HOH,Na+,Cl-,K+,MG
trajout apt1_01_dry_20261125-1200.pdb pdb
run
EOF
[ -s "apt1_01_dry_20261125-1200.pdb" ] || { echo "ERROR: cpptraj did not write apt1_01_dry_20261125-1200.pdb" >&2; exit 1; }
