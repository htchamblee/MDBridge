#!/usr/bin/env bash
# Build MDBridge's video encoding libraries on Linux: SVT-AV1 and a minimal LGPL-only FFmpeg.
#
#   packaging/ffmpeg/build_ffmpeg.sh [--prefix DIR] [--work DIR] [--jobs N]
#
# Result (default PREFIX: build/ffmpeg/prefix), found by CMake with -DMDBRIDGE_FFMPEG_ROOT=PREFIX:
#   PREFIX/include, PREFIX/lib      libavcodec, libavformat, libavutil, libswscale and libSvtAv1Enc
#                                   (shared libraries, RUNPATH $ORIGIN, no static libraries)
#   PREFIX/share/mdbridge-ffmpeg/   BUILDINFO.txt (versions, commits, the configure line),
#                                   the licences of both projects
#   PREFIX/../ffmpeg-<version>-source.tar.xz
#                                   the corresponding source of the FFmpeg libraries (LGPL 2.1 section
#                                   6): FFmpeg at the pinned commit, this script and versions.env
#
# Sources are git checkouts of the tags in versions.env, checked against the pinned commit hashes.
# Offline: FFMPEG_SRC / SVTAV1_SRC may name existing checkouts (they must be at the pinned commits
# and clean). Releases are built on Rocky Linux 9 (glibc 2.34, the AppImage's floor) in CI; needs
# gcc, make, cmake 3.16+, nasm 2.14+, pkg-config, git, xz, readelf (binutils).
#
# The FFmpeg build is LGPL 2.1 or later only: --disable-gpl --disable-nonfree --disable-version3, no
# external library except SVT-AV1 (--disable-autodetect), only the components MDBridge uses (the
# libsvtav1 and mpeg4 encoders, the mpeg4 decoder and the av1 / mpeg4 parsers for the self-check,
# the mp4 muxer, the mov demuxer, the file protocol). The script fails if FFmpeg's configure reports
# any other licence.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
PREFIX="$ROOT/build/ffmpeg/prefix"
WORK="$ROOT/build/ffmpeg/work"
JOBS="$(nproc 2>/dev/null || echo 4)"

while [ $# -gt 0 ]; do
    case "$1" in
        --prefix) PREFIX="$(realpath -m "$2")"; shift 2 ;;
        --work) WORK="$(realpath -m "$2")"; shift 2 ;;
        --jobs) JOBS="$2"; shift 2 ;;
        -h|--help) sed -n '2,23p' "$0"; exit 0 ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
done

# pins: KEY=VALUE lines of versions.env (no shell code is evaluated)
while IFS='=' read -r key value; do
    case "$key" in ''|\#*) continue ;; esac
    printf -v "$key" '%s' "$value"
done < "$HERE/versions.env"

for tool in git make cmake nasm pkg-config xz tar readelf; do
    command -v "$tool" > /dev/null || { echo "$tool is needed to build FFmpeg (see $0 --help)" >&2; exit 1; }
done
mkdir -p "$WORK"

# checkout(<name> <git url> <tag> <commit> <existing checkout or empty>) -> source directory
checkout() {
    local name="$1" url="$2" tag="$3" commit="$4" dir="${5:-$WORK/$1}"
    if [ ! -d "$dir/.git" ]; then
        rm -rf "$dir"
        git -c advice.detachedHead=false clone --quiet --depth 1 --branch "$tag" "$url" "$dir" >&2
    fi
    local head
    head="$(git -C "$dir" rev-parse HEAD)"
    if [ "$head" != "$commit" ]; then
        echo "$name: $dir is at $head, expected $commit ($tag)" >&2
        exit 1
    fi
    if [ -n "$(git -C "$dir" status --porcelain --untracked-files=no)" ]; then
        echo "$name: $dir has local changes" >&2
        exit 1
    fi
    echo "$dir"
}

echo "== sources (FFmpeg $FFMPEG_VERSION, SVT-AV1 $SVTAV1_VERSION)"
SVT_SRC="$(checkout SVT-AV1 "$SVTAV1_GIT" "$SVTAV1_TAG" "$SVTAV1_COMMIT" "${SVTAV1_SRC:-}")"
FF_SRC="$(checkout FFmpeg "$FFMPEG_GIT" "$FFMPEG_TAG" "$FFMPEG_COMMIT" "${FFMPEG_SRC:-}")"

rm -rf "$PREFIX"
mkdir -p "$PREFIX"

echo "== SVT-AV1 (shared library, no applications)"
SVT_BUILD="$WORK/svt-av1-build"
rm -rf "$SVT_BUILD"
# portable code with run-time CPU dispatch (no -march=native); libraries in PREFIX/lib, not lib64
cmake -S "$SVT_SRC" -B "$SVT_BUILD" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PREFIX" \
    -DCMAKE_INSTALL_LIBDIR=lib -DBUILD_SHARED_LIBS=ON -DBUILD_APPS=OFF -DBUILD_TESTING=OFF \
    -DNATIVE=OFF -DSVT_AV1_LTO=OFF -DREPRODUCIBLE_BUILDS=ON -DEXCLUDE_HASH=ON \
    -DCMAKE_INSTALL_RPATH='$ORIGIN' > "$WORK/svt-av1-configure.log" 2>&1 ||
    { tail -n 40 "$WORK/svt-av1-configure.log"; exit 1; }
cmake --build "$SVT_BUILD" --parallel "$JOBS" > "$WORK/svt-av1-build.log" 2>&1 ||
    { tail -n 40 "$WORK/svt-av1-build.log"; exit 1; }
cmake --install "$SVT_BUILD" > "$WORK/svt-av1-install.log" 2>&1 ||
    { tail -n 40 "$WORK/svt-av1-install.log"; exit 1; }

echo "== FFmpeg (LGPL only, shared libraries)"
FF_BUILD="$WORK/ffmpeg-build"
rm -rf "$FF_BUILD"
mkdir -p "$FF_BUILD"
CONFIGURE_FLAGS=(
    --prefix="$PREFIX"
    --enable-shared --disable-static --enable-pic
    --disable-gpl --disable-nonfree --disable-version3
    --disable-autodetect --disable-programs --disable-doc --disable-network
    --disable-avdevice --disable-avfilter --disable-swresample
    --disable-everything
    --enable-libsvtav1
    --enable-encoder=libsvtav1,mpeg4
    --enable-decoder=mpeg4
    --enable-parser=av1,mpeg4video
    --enable-muxer=mp4
    --enable-demuxer=mov
    --enable-protocol=file
)
(cd "$FF_BUILD" && PKG_CONFIG_PATH="$PREFIX/lib/pkgconfig" "$FF_SRC/configure" "${CONFIGURE_FLAGS[@]}" \
    > "$WORK/ffmpeg-configure.log" 2>&1) ||
    { tail -n 40 "$WORK/ffmpeg-configure.log"; tail -n 40 "$FF_BUILD/ffbuild/config.log" 2>/dev/null; exit 1; }
# RUNPATH $ORIGIN: libavcodec finds libavutil and libSvtAv1Enc next to itself. Added to the
# generated makefile rather than passed to configure, whose eval mangles '$'. library.mak expands it
# twice ($(eval)), so make turns '\$$$$' into '\$', which the shell passes on as '$'
printf '%s\n' 'LDSOFLAGS += -Wl,-rpath,\$$$$ORIGIN' >> "$FF_BUILD/ffbuild/config.mak"
license="$(sed -n 's/^License: //p' "$WORK/ffmpeg-configure.log" | head -n 1)"
if [ "$license" != "LGPL version 2.1 or later" ]; then
    echo "FFmpeg's configure reports the licence '$license'; MDBridge needs 'LGPL version 2.1 or later'" >&2
    exit 1
fi
make -C "$FF_BUILD" -j "$JOBS" > "$WORK/ffmpeg-build.log" 2>&1 || { tail -n 40 "$WORK/ffmpeg-build.log"; exit 1; }
make -C "$FF_BUILD" install > "$WORK/ffmpeg-install.log" 2>&1 || { tail -n 40 "$WORK/ffmpeg-install.log"; exit 1; }

echo "== check"
# the configuration string avcodec_configuration() returns is a literal in libavcodec
config="$(grep -a -o -- '--prefix=[^[:cntrl:]]*' "$PREFIX"/lib/libavcodec.so.* | head -n 1 | cut -d: -f2-)"
case " $config " in
    *" --enable-gpl "*|*" --enable-nonfree "*|*" --enable-version3 "*)
        echo "libavcodec was configured with a GPL / nonfree / version3 flag: $config" >&2; exit 1 ;;
esac
for lib in libavcodec libavformat libavutil libswscale libSvtAv1Enc; do
    so="$(ls "$PREFIX"/lib/$lib.so.* 2>/dev/null | head -n 1)"
    [ -n "$so" ] || { echo "$lib was not built" >&2; exit 1; }
    readelf -d "$so" | grep -q 'R.*PATH.*\$ORIGIN' || { echo "$so has no RUNPATH \$ORIGIN" >&2; exit 1; }
done
rm -rf "$PREFIX/share/doc" "$PREFIX/share/man" "$PREFIX"/share/ffmpeg/examples

echo "== licences and build information"
INFO="$PREFIX/share/mdbridge-ffmpeg"
mkdir -p "$INFO"
cp "$FF_SRC/COPYING.LGPLv2.1" "$INFO/FFmpeg-COPYING.LGPLv2.1.txt"
cp "$FF_SRC/LICENSE.md" "$INFO/FFmpeg-LICENSE.md"
cp "$SVT_SRC/LICENSE.md" "$INFO/SVT-AV1-LICENSE.md"
cp "$SVT_SRC/LICENSE-BSD2.md" "$INFO/SVT-AV1-LICENSE-BSD2.md"
cp "$SVT_SRC/PATENTS.md" "$INFO/SVT-AV1-PATENTS.md"
{
    echo "MDBridge video libraries, built by packaging/ffmpeg/build_ffmpeg.sh"
    echo "FFmpeg $FFMPEG_VERSION ($FFMPEG_GIT tag $FFMPEG_TAG, commit $FFMPEG_COMMIT)"
    echo "SVT-AV1 $SVTAV1_VERSION ($SVTAV1_GIT tag $SVTAV1_TAG, commit $SVTAV1_COMMIT)"
    echo "FFmpeg licence: $license"
    echo "FFmpeg configuration: $config"
    echo "Built on: $( (. /etc/os-release && echo "$PRETTY_NAME") 2>/dev/null || uname -s), $(gcc --version | head -n 1), glibc $(ldd --version | head -n 1 | awk '{print $NF}')"
} > "$INFO/BUILDINFO.txt"

echo "== corresponding source of the FFmpeg libraries"
ARCHIVE="$(dirname "$PREFIX")/ffmpeg-$FFMPEG_VERSION-source.tar.xz"
STAGE="$WORK/source-stage"
rm -rf "$STAGE" && mkdir -p "$STAGE/ffmpeg-$FFMPEG_VERSION/MDBRIDGE-BUILD"
git -C "$FF_SRC" archive --format=tar --prefix="ffmpeg-$FFMPEG_VERSION/" HEAD | tar -x -C "$STAGE"
cp "$HERE/build_ffmpeg.sh" "$HERE/versions.env" "$INFO/BUILDINFO.txt" "$STAGE/ffmpeg-$FFMPEG_VERSION/MDBRIDGE-BUILD/"
tar -C "$STAGE" --sort=name --owner=0 --group=0 --numeric-owner --mtime='2000-01-01 00:00Z' -cf - \
    "ffmpeg-$FFMPEG_VERSION" | xz -6 -T1 > "$ARCHIVE"
rm -rf "$STAGE"

echo "== done"
cat "$INFO/BUILDINFO.txt" | sed 's/^/   /'
(cd "$PREFIX/lib" && ls -1 ./*.so.*) | sed 's|^\./|   |'
du -sh "$ARCHIVE" | sed 's/^/   source: /'
