#!/bin/bash
# Helix-preserving distance restraints (DISANG)
# Project apt1, step 02 - written by MDATK <VERSION> on <DATE>.
# MDATK only writes this script; review it, then run it with: bash apt1_02_amber-restraints_20261125-1200.sh
# Please cite (AmberTools):
#   - Case D. A., Aktulga H. M., Belfon K., Cerutti D. S., Cisneros G. A., Cruzeiro V. W. D.,
#     Forouzesh N., Giese T. J., Götz A. W., Gohlke H., et al. J. Chem. Inf. Model. 2023, 63,
#     6183-6191. doi:10.1021/acs.jcim.3c01153

set -eo pipefail
cd "$(dirname "$0")"

# 6 base pairs from the typed secondary structure
python3 apt1_02_helix_20261125-1200.py ../01_amber-build/apt1_01_parmed_20261125-1200.parm7 apt1_02_helix_20261125-1200.dist
[ -s "apt1_02_helix_20261125-1200.dist" ] || { echo "ERROR: the restraint script did not write apt1_02_helix_20261125-1200.dist" >&2; exit 1; }
