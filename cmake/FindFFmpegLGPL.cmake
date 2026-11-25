# FindFFmpegLGPL: the LGPL-only FFmpeg libraries built by packaging/ffmpeg (video encoding).
#
#   find_package(FFmpegLGPL [REQUIRED])
#
# Looks in MDBRIDGE_FFMPEG_ROOT (CMake or environment variable; the --prefix of
# packaging/ffmpeg/build_ffmpeg.sh or -Prefix of packaging/ffmpeg/windows/build_ffmpeg.ps1), then in
# FFmpegLGPL_ROOT and the usual places. Needs libavcodec, libavformat, libavutil and libswscale.
#
# The libraries must be an LGPL-only build: the configuration string that avcodec_configuration()
# returns (a literal in libavcodec's shared library / DLL) must not contain --enable-gpl or
# --enable-nonfree, and the licence string must be LGPL. A GPL or nonfree build is a fatal error,
# not "not found": it must never be bundled. (The same check runs again at run time,
# mdbridge::render::video::ffmpegBuildInfo, and in a unit test.)
#
# Results:
#   FFmpegLGPL_FOUND, FFmpegLGPL_VERSION (FFmpeg's release, e.g. 9.0.2)
#   FFmpegLGPL_CONFIGURATION    FFmpeg's configure flags, FFmpegLGPL_LICENSE its licence string
#   FFmpegLGPL_HAS_SVTAV1       the libsvtav1 encoder is built in
#   FFmpegLGPL_RUNTIME_FILES    the shared libraries to bundle (Linux: lib*.so.* with their soname
#                               links, FFmpeg's and SVT-AV1's; Windows: the DLLs of the prefix's bin)
#   FFmpegLGPL_BUILDINFO        share/mdbridge-ffmpeg/BUILDINFO.txt of the prefix, if there
#   FFmpegLGPL::FFmpeg          imported target linking the four libraries

set(_ffmpeg_hints)
foreach(_root "${MDBRIDGE_FFMPEG_ROOT}" "$ENV{MDBRIDGE_FFMPEG_ROOT}")
    if(_root)
        list(APPEND _ffmpeg_hints "${_root}")
    endif()
endforeach()

find_path(FFmpegLGPL_INCLUDE_DIR libavcodec/avcodec.h HINTS ${_ffmpeg_hints} PATH_SUFFIXES include)
foreach(_lib avcodec avformat avutil swscale)
    find_library(FFmpegLGPL_${_lib}_LIBRARY NAMES ${_lib} HINTS ${_ffmpeg_hints} PATH_SUFFIXES lib)
    mark_as_advanced(FFmpegLGPL_${_lib}_LIBRARY)
endforeach()
mark_as_advanced(FFmpegLGPL_INCLUDE_DIR)

set(FFmpegLGPL_VERSION "")
if(FFmpegLGPL_INCLUDE_DIR AND EXISTS "${FFmpegLGPL_INCLUDE_DIR}/libavutil/ffversion.h")
    file(STRINGS "${FFmpegLGPL_INCLUDE_DIR}/libavutil/ffversion.h" _ffversion REGEX "#define FFMPEG_VERSION ")
    string(REGEX REPLACE ".*\"n?([^\"]*)\".*" "\\1" FFmpegLGPL_VERSION "${_ffversion}")
endif()

# the shared libraries: what the configuration check reads and what installs bundle
set(FFmpegLGPL_RUNTIME_FILES)
set(_avcodec_binary "")
if(FFmpegLGPL_avcodec_LIBRARY)
    get_filename_component(_libdir "${FFmpegLGPL_avcodec_LIBRARY}" DIRECTORY)
    get_filename_component(_prefix "${_libdir}" DIRECTORY)
    if(WIN32)
        file(GLOB FFmpegLGPL_RUNTIME_FILES "${_prefix}/bin/*.dll")
        file(GLOB _avcodec_binary "${_prefix}/bin/avcodec-*.dll")
    else()
        file(GLOB FFmpegLGPL_RUNTIME_FILES "${_libdir}/libavcodec.so.*" "${_libdir}/libavformat.so.*"
             "${_libdir}/libavutil.so.*" "${_libdir}/libswscale.so.*" "${_libdir}/libSvtAv1Enc.so.*")
        get_filename_component(_avcodec_binary "${FFmpegLGPL_avcodec_LIBRARY}" REALPATH)
    endif()
    list(GET _avcodec_binary 0 _avcodec_binary)
    if(EXISTS "${_prefix}/share/mdbridge-ffmpeg/BUILDINFO.txt")
        set(FFmpegLGPL_BUILDINFO "${_prefix}/share/mdbridge-ffmpeg/BUILDINFO.txt")
    endif()
endif()

set(FFmpegLGPL_CONFIGURATION "")
set(FFmpegLGPL_LICENSE "")
set(FFmpegLGPL_HAS_SVTAV1 FALSE)
if(_avcodec_binary AND EXISTS "${_avcodec_binary}")
    file(STRINGS "${_avcodec_binary}" _configuration REGEX "--prefix=" LENGTH_MINIMUM 20 LIMIT_COUNT 1)
    file(STRINGS "${_avcodec_binary}" _license REGEX "^(libavcodec license: )?(L?GPL|nonfree)" LIMIT_COUNT 4)
    file(STRINGS "${_avcodec_binary}" _svtav1 REGEX "^libsvtav1$" LIMIT_COUNT 1)
    string(STRIP "${_configuration}" FFmpegLGPL_CONFIGURATION)
    foreach(_line IN LISTS _license)
        if(_line MATCHES "license: (.*)$")
            set(FFmpegLGPL_LICENSE "${CMAKE_MATCH_1}")
        endif()
    endforeach()
    if(_svtav1)
        set(FFmpegLGPL_HAS_SVTAV1 TRUE)
    endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFmpegLGPL
    REQUIRED_VARS FFmpegLGPL_INCLUDE_DIR FFmpegLGPL_avcodec_LIBRARY FFmpegLGPL_avformat_LIBRARY
                  FFmpegLGPL_avutil_LIBRARY FFmpegLGPL_swscale_LIBRARY FFmpegLGPL_CONFIGURATION
    VERSION_VAR FFmpegLGPL_VERSION)

if(FFmpegLGPL_FOUND)
    separate_arguments(_flags UNIX_COMMAND "${FFmpegLGPL_CONFIGURATION}")
    foreach(_forbidden --enable-gpl --enable-nonfree)
        if(_forbidden IN_LIST _flags)
            message(FATAL_ERROR "The FFmpeg libraries at ${_libdir} were configured with ${_forbidden}; MDBridge "
                                "bundles only LGPL builds of FFmpeg (packaging/ffmpeg). Configuration: "
                                "${FFmpegLGPL_CONFIGURATION}")
        endif()
    endforeach()
    if(FFmpegLGPL_LICENSE AND NOT FFmpegLGPL_LICENSE MATCHES "^LGPL")
        message(FATAL_ERROR "The FFmpeg libraries at ${_libdir} are licensed '${FFmpegLGPL_LICENSE}'; MDBridge "
                            "bundles only LGPL builds of FFmpeg (packaging/ffmpeg)")
    endif()
    if(NOT TARGET FFmpegLGPL::FFmpeg)
        add_library(FFmpegLGPL::FFmpeg INTERFACE IMPORTED)
        foreach(_lib avcodec avformat avutil swscale)
            add_library(FFmpegLGPL::${_lib} UNKNOWN IMPORTED)
            set_target_properties(FFmpegLGPL::${_lib} PROPERTIES
                IMPORTED_LOCATION "${FFmpegLGPL_${_lib}_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${FFmpegLGPL_INCLUDE_DIR}")
            target_link_libraries(FFmpegLGPL::FFmpeg INTERFACE FFmpegLGPL::${_lib})
        endforeach()
    endif()
endif()
