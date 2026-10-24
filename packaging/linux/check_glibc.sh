#!/usr/bin/env bash
# Report the newest glibc symbol version the bundled programs and libraries need.
#
#   packaging/linux/check_glibc.sh BUNDLE_DIR [MAX_VERSION]
#
# MAX_VERSION defaults to 2.34 (RHEL 9; Ubuntu 22.04 has 2.35). Lists every file that needs a
# newer glibc and exits with 1 if there is one. Needs objdump (binutils).
set -euo pipefail
DIR="${1:?usage: check_glibc.sh BUNDLE_DIR [MAX_VERSION]}"
MAX="${2:-2.34}"

newest="0"
bad=0
while IFS= read -r -d '' f; do
    head -c 4 "$f" 2>/dev/null | grep -q $'\x7fELF' || continue
    v="$( (objdump -T "$f" 2>/dev/null || true) | (grep -o 'GLIBC_[0-9.]*' || true) | sed 's/GLIBC_//' | sort -V | tail -n 1)"
    [ -n "$v" ] || continue
    newest="$(printf '%s\n%s\n' "$newest" "$v" | sort -V | tail -n 1)"
    if [ "$(printf '%s\n%s\n' "$MAX" "$v" | sort -V | tail -n 1)" != "$MAX" ]; then
        echo "needs glibc $v: ${f#"$DIR"/}"
        bad=1
    fi
done < <(find "$DIR" -type f \( -name '*.so*' -o -perm -u+x \) -print0)

echo "newest glibc symbol version needed: $newest (limit $MAX)"
exit "$bad"
