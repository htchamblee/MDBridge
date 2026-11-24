# Dependencies fetched at configure time, so a source build installs nothing on the system.
#
# Each one is pinned to a GitHub release tarball and its SHA512 (the hashes vcpkg's ports use;
# zlib: its release asset zlib-1.3.1.tar.gz, SHA-256 9a93b2b7...df23 as published).
# Offline or mirrored builds: point FetchContent at a checkout of the same tag with
#   -DFETCHCONTENT_SOURCE_DIR_<NAME>=<dir>   (NAME: TOMLPLUSPLUS, CLI11, GOOGLETEST, VULKANHEADERS,
#                                             ZLIB)
# or use installed packages with -DMDBRIDGE_USE_SYSTEM_DEPS=ON.
# Qt 6 is not fetched (too large to build per checkout): see docs/BUILDING.md.
include(FetchContent)

# mdbridge_dependency(<name> <url> <sha512> [find_package arguments ...])
function(mdbridge_dependency name url sha512)
    set(find_args)
    if(MDBRIDGE_USE_SYSTEM_DEPS)
        set(find_args FIND_PACKAGE_ARGS ${ARGN})
    endif()
    FetchContent_Declare(${name}
        URL "${url}"
        URL_HASH SHA512=${sha512}
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        SYSTEM  # third-party headers: no warnings from them in MDBridge's strict builds
        ${find_args})
endfunction()

# TOML parser of the project files (MIT)
mdbridge_dependency(tomlplusplus
    https://github.com/marzer/tomlplusplus/archive/v3.4.0.tar.gz
    c227fc8147c9459b29ad24002aaf6ab2c42fac22ea04c1c52b283a0172581ccd4527b33c1931e0ef0d1db6b6a53f9e9882c6d4231c7f3494cf070d0220741aa5
    3.4)

# command line parser of the CLI (BSD-3-Clause)
mdbridge_dependency(CLI11
    https://github.com/CLIUtils/CLI11/archive/v2.7.2.tar.gz
    437cd9d704c1b0de0516bdc8cb115befb2e7170417375a64845a2bc2385fa07948c4af724ec615d27fbc102dbe13a7117f33873d6776fe9b17f293d42fe64dd9
    2.4)

# Vulkan API headers of the 3D view (Apache-2.0 / MIT); the loader is the user's (graphics driver)
mdbridge_dependency(VulkanHeaders
    https://github.com/KhronosGroup/Vulkan-Headers/archive/vulkan-sdk-1.4.363.0.tar.gz
    262e4fabdcea0e13c18aac831b2ad6d99d17c83ae0e740d658792a5f48985befe128e9cd6e69309b21580eda17aecf564359379344300bf2a352a95830f9f944
    CONFIG)

# gzip-compressed density maps (zlib licence). zlib's own CMake project is not used (it builds a
# shared library, installs itself and renames zconf.h in its source tree): its sources are compiled
# into the static mdbridge_zlib, a name that cannot clash with the zlib Qt bundles or finds.
# MDBRIDGE_ZLIB_TARGET names the library to link.
if(MDBRIDGE_USE_SYSTEM_DEPS)
    find_package(ZLIB 1.2.11)
endif()
if(ZLIB_FOUND)
    set(MDBRIDGE_ZLIB_TARGET ZLIB::ZLIB)
else()
    FetchContent_Declare(zlib
        URL "https://github.com/madler/zlib/releases/download/v1.3.1/zlib-1.3.1.tar.gz"
        URL_HASH SHA512=580677aad97093829090d4b605ac81c50327e74a6c2de0b85dd2e8525553f3ddde17556ea46f8f007f89e435493c9a20bc997d1ef1c1c2c23274528e3c46b94f
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        SOURCE_SUBDIR no-cmake-project)  # no CMakeLists.txt there: nothing is added
    FetchContent_MakeAvailable(zlib)
    enable_language(C)
    set(zlib_sources adler32.c compress.c crc32.c deflate.c infback.c inffast.c inflate.c inftrees.c
                     trees.c uncompr.c zutil.c)
    list(TRANSFORM zlib_sources PREPEND "${zlib_SOURCE_DIR}/")
    add_library(mdbridge_zlib STATIC ${zlib_sources})
    target_include_directories(mdbridge_zlib SYSTEM PUBLIC "${zlib_SOURCE_DIR}")
    # z_-prefixed symbols: never interposed on the system zlib that libraries Qt loads use
    target_compile_definitions(mdbridge_zlib PUBLIC Z_PREFIX)
    if(MSVC)
        target_compile_definitions(mdbridge_zlib PRIVATE
            _CRT_SECURE_NO_DEPRECATE _CRT_NONSTDC_NO_DEPRECATE)
    endif()
    set_target_properties(mdbridge_zlib PROPERTIES POSITION_INDEPENDENT_CODE ON)
    set(MDBRIDGE_ZLIB_TARGET mdbridge_zlib)
endif()

# unit tests only, never bundled (BSD-3-Clause)
mdbridge_dependency(googletest
    https://github.com/google/googletest/archive/v1.18.0.tar.gz
    ba0f5769ccf34acf1bc72d1f7e9ffb8202176d02b64f6f3d9047accfc0cf9026ff5a653d24935e2705fff8709566676452616c93ca0ca6277f1e21d79b58a10a
    NAMES GTest)
