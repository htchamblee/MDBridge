#!/usr/bin/env bash
# Build the MDBridge AppImage from a release build.
#
#   packaging/linux/build_appimage.sh [--build DIR] [--out DIR] [--work DIR]
#
# Steps: cmake --install into an AppDir -> linuxdeploy and its Qt plugin bundle Qt (xcb and offscreen
# platform plugins) and the libraries that not every Linux system has -> licences of those
# libraries -> appimagetool. Result: OUT/MDBridge-<version>-x86_64.AppImage (default OUT: dist/).
#
# --build: the configured and built CMake tree (default: build/release; docs/BUILDING.md). Build on
# the oldest distribution the AppImage must run on (glibc is the host's): Rocky Linux 9 for
# releases. Needs curl (to fetch the pinned tools), file, and rpm or dpkg for the licences.
#
# Environment:
#   QMAKE                qmake of the Qt the build used (default: found from the build's CMake cache)
#   LINUXDEPLOY, LINUXDEPLOY_PLUGIN_QT, APPIMAGETOOL
#                        local copies of the tools instead of the pinned downloads below
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
BUILD="$ROOT/build/release"
OUT="$ROOT/dist"
WORK="$ROOT/build/appimage"

# pinned tool releases (sha256 of the release files)
LINUXDEPLOY_URL="https://github.com/linuxdeploy/linuxdeploy/releases/download/1-alpha-20251107-1/linuxdeploy-x86_64.AppImage"
LINUXDEPLOY_SHA256="c20cd71e3a4e3b80c3483cef793cda3f4e990aca14014d23c544ca3ce1270b4d"
PLUGIN_QT_URL="https://github.com/linuxdeploy/linuxdeploy-plugin-qt/releases/download/1-alpha-20250213-1/linuxdeploy-plugin-qt-x86_64.AppImage"
PLUGIN_QT_SHA256="15106be885c1c48a021198e7e1e9a48ce9d02a86dd0a1848f00bdbf3c1c92724"
APPIMAGETOOL_URL="https://github.com/AppImage/appimagetool/releases/download/1.9.1/appimagetool-x86_64.AppImage"
APPIMAGETOOL_SHA256="ed4ce84f0d9caff66f50bcca6ff6f35aae54ce8135408b3fa33abfc3cb384eb0"

while [ $# -gt 0 ]; do
    case "$1" in
        --build) BUILD="$(realpath -m "$2")"; shift 2 ;;
        --out) OUT="$(realpath -m "$2")"; shift 2 ;;
        --work) WORK="$(realpath -m "$2")"; shift 2 ;;
        -h|--help) sed -n '2,19p' "$0"; exit 0 ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
done

[ -f "$BUILD/CMakeCache.txt" ] || { echo "no CMake build in $BUILD (see docs/BUILDING.md)" >&2; exit 1; }
cache() { sed -n "s/^$1:[A-Z]*=//p" "$BUILD/CMakeCache.txt" | head -n 1; }
VERSION="$(cache CMAKE_PROJECT_VERSION)"
ARCH="$(uname -m)"
QT_DIR="$(cd "$(cache Qt6_DIR)/../../.." && pwd)"  # <Qt>/lib/cmake/Qt6 -> <Qt>
export QMAKE="${QMAKE:-$QT_DIR/bin/qmake}"
[ -x "$QMAKE" ] || { echo "qmake not found ($QMAKE); set QMAKE" >&2; exit 1; }
echo "== MDBridge $VERSION ($ARCH), Qt $("$QMAKE" -query QT_VERSION) at $QT_DIR"
mkdir -p "$WORK/tools" "$OUT"

# fetch(<url> <sha256> <file>): download a pinned tool once and check it
fetch() {
    local target="$WORK/tools/$3"
    if [ ! -x "$target" ]; then
        curl -fsSL --retry 3 -o "$target.part" "$1"
        echo "$2  $target.part" | sha256sum -c --quiet - || { echo "checksum mismatch: $1" >&2; exit 1; }
        chmod +x "$target.part" && mv "$target.part" "$target"
    fi
    echo "$target"
}
LINUXDEPLOY="${LINUXDEPLOY:-$(fetch "$LINUXDEPLOY_URL" "$LINUXDEPLOY_SHA256" linuxdeploy-x86_64.AppImage)}"
PLUGIN_QT="${LINUXDEPLOY_PLUGIN_QT:-$(fetch "$PLUGIN_QT_URL" "$PLUGIN_QT_SHA256" linuxdeploy-plugin-qt-x86_64.AppImage)}"
APPIMAGETOOL="${APPIMAGETOOL:-$(fetch "$APPIMAGETOOL_URL" "$APPIMAGETOOL_SHA256" appimagetool-x86_64.AppImage)}"
# linuxdeploy finds its plugins on PATH; containers and CI runners have no FUSE
mkdir -p "$WORK/plugins" && ln -sf "$PLUGIN_QT" "$WORK/plugins/linuxdeploy-plugin-qt"
export PATH="$WORK/plugins:$PATH" APPIMAGE_EXTRACT_AND_RUN=1

echo "== install into the AppDir"
APPDIR="$WORK/AppDir"
rm -rf "$APPDIR"
cmake --install "$BUILD" --prefix "$APPDIR/usr" > "$WORK/install.log"

echo "== linuxdeploy (Qt and libraries)"
export EXTRA_PLATFORM_PLUGINS="libqoffscreen.so"  # headless checks; xcb is always deployed
# the build's Qt, not a distribution's Qt that may be installed too
LD_LIBRARY_PATH="$QT_DIR/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" "$LINUXDEPLOY" --appdir "$APPDIR" \
    --executable "$APPDIR/usr/bin/mdbridge" \
    --desktop-file "$APPDIR/usr/share/applications/mdbridge.desktop" \
    --icon-file "$APPDIR/usr/share/icons/hicolor/256x256/apps/mdbridge.png" \
    --custom-apprun "$HERE/AppRun" \
    --plugin qt > "$WORK/linuxdeploy.log" 2>&1 || { tail -n 40 "$WORK/linuxdeploy.log"; exit 1; }
# the Vulkan loader is the graphics driver's: never bundle one
if ls "$APPDIR"/usr/lib/libvulkan.so* > /dev/null 2>&1; then
    echo "the AppDir contains a Vulkan loader; it must come from the host" >&2
    exit 1
fi

echo "== libraries taken from the host"
# what every Linux desktop has; anything else must be bundled (e.g. libOpenGL.so.0 is often missing)
HOST_ALLOWED='^(ld-linux-x86-64\.so\.2|lib(c|m|dl|pthread|rt)\.so\.[0-9]+|libstdc\+\+\.so\.6|libgcc_s\.so\.1|libz\.so\.1|libX11(-xcb)?\.so\.[0-9]+|libxcb\.so\.1|libfontconfig\.so\.1|libfreetype\.so\.6|libGL\.so\.1|libEGL\.so\.1|libgpg-error\.so\.0)$'
host_libs="$(find "$APPDIR/usr" -type f \( -name '*.so*' -o -perm -u+x \) -exec sh -c \
    'readelf -d "$1" 2>/dev/null | sed -n "s/.*(NEEDED).*\[\(.*\)\]/\1/p"' _ {} \; | sort -u |
    while read -r lib; do [ -e "$APPDIR/usr/lib/$lib" ] || echo "$lib"; done)"
echo "   $(echo $host_libs)"
unexpected="$(echo "$host_libs" | grep -vE "$HOST_ALLOWED" || true)"
if [ -n "$unexpected" ]; then
    echo "the AppImage would need these libraries from the host: $(echo $unexpected)" >&2
    exit 1
fi

echo "== licences of the bundled system libraries"
LICENSES="$APPDIR/usr/share/doc/mdbridge/licenses/system"
mkdir -p "$LICENSES"
for lib in "$APPDIR"/usr/lib/*.so*; do
    name="$(basename "$lib")"
    [ -e "$QT_DIR/lib/$name" ] && continue  # Qt's own (and its ICU): Qt's notice covers them
    source="$(ldconfig -p | awk -v n="$name" '$1 == n {print $NF; exit}')"
    [ -n "$source" ] || source="$(find /usr/lib64 /usr/lib /lib64 /lib -name "$name" 2>/dev/null | head -n 1)"
    package=""
    if command -v rpm > /dev/null && [ -n "$source" ]; then
        package="$(rpm -qf --qf '%{NAME}\n' "$source" 2>/dev/null | head -n 1)" || package=""
        for f in $(rpm -qL "$package" 2>/dev/null; rpm -qd "$package" 2>/dev/null | grep -iE 'copying|license' || true); do
            [ -f "$f" ] && mkdir -p "$LICENSES/$package" && cp "$f" "$LICENSES/$package/"
        done
    elif command -v dpkg > /dev/null && [ -n "$source" ]; then
        package="$(dpkg -S "$(readlink -f "$source")" 2>/dev/null | head -n 1 | cut -d: -f1)" || package=""
        [ -f "/usr/share/doc/$package/copyright" ] && mkdir -p "$LICENSES/$package" &&
            cp "/usr/share/doc/$package/copyright" "$LICENSES/$package/"
    fi
    printf '   %-32s %s\n' "$name" "${package:-(package not found)}"
done

echo "== appimagetool"
target="$OUT/MDBridge-$VERSION-$ARCH.AppImage"
rm -f "$target"
ARCH="$ARCH" "$APPIMAGETOOL" --no-appstream "$APPDIR" "$target" > "$WORK/appimagetool.log" 2>&1 ||
    { tail -n 40 "$WORK/appimagetool.log"; exit 1; }
chmod +x "$target"

echo "== done"
du -sh "$APPDIR" | sed 's/^/   AppDir:   /'
du -sh "$target" | sed 's/^/   AppImage: /'
