# Dependencies fetched at configure time, so a source build installs nothing on the system.
#
# Each one is pinned to a GitHub release tarball and its SHA512 (the hashes vcpkg's ports use).
# Offline or mirrored builds: point FetchContent at a checkout of the same tag with
#   -DFETCHCONTENT_SOURCE_DIR_<NAME>=<dir>   (NAME: TOMLPLUSPLUS, CLI11, GOOGLETEST, VULKANHEADERS)
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

# unit tests only, never bundled (BSD-3-Clause)
mdbridge_dependency(googletest
    https://github.com/google/googletest/archive/v1.18.0.tar.gz
    ba0f5769ccf34acf1bc72d1f7e9ffb8202176d02b64f6f3d9047accfc0cf9026ff5a653d24935e2705fff8709566676452616c93ca0ca6277f1e21d79b58a10a
    NAMES GTest)
