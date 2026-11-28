#!/bin/bash
# VMD
# Project apt1, step 04 - written by MDATK <VERSION> on <DATE>.
# MDATK only writes this script; review it, then run it with: bash apt1_04_vmd_20261125-1200.sh
# Please cite (VMD):
#   - Humphrey W., Dalke A., Schulten K. J. Mol. Graph. 1996, 14, 33-38.
#     doi:10.1016/0263-7855(96)00018-5

set -eo pipefail
cd "$(dirname "$0")"

vmd -e apt1_04_vmd_20261125-1200.tcl
