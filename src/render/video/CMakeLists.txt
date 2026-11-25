# video: MP4 encoding with the bundled LGPL FFmpeg libraries (packaging/ffmpeg). No Qt, so it is
# built and tested without the Qt front ends too.
#
#   -DMDBRIDGE_FFMPEG=AUTO|ON|OFF   AUTO (default): use the libraries when found; ON: required;
#                                   OFF: build without them (video then fails with "video encoding
#                                   is not available in this build")
#   -DMDBRIDGE_FFMPEG_ROOT=<prefix> where packaging/ffmpeg's script installed them
set(MDBRIDGE_FFMPEG AUTO CACHE STRING "Video encoding with packaging/ffmpeg's LGPL FFmpeg libraries: AUTO, ON or OFF")
set_property(CACHE MDBRIDGE_FFMPEG PROPERTY STRINGS AUTO ON OFF)
set(MDBRIDGE_FFMPEG_ROOT "" CACHE PATH "Prefix of the FFmpeg libraries built by packaging/ffmpeg")

if(MDBRIDGE_FFMPEG STREQUAL "ON")
    find_package(FFmpegLGPL REQUIRED)
elseif(MDBRIDGE_FFMPEG STREQUAL "AUTO")
    find_package(FFmpegLGPL)
endif()

add_library(mdbridge_video STATIC
    encoder.cpp
    encoder.hpp
    ffmpeg_support.hpp
    probe.cpp
    probe.hpp)
add_library(mdbridge::video ALIAS mdbridge_video)
# included as "video/encoder.hpp"
target_include_directories(mdbridge_video PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/..")
target_link_libraries(mdbridge_video PUBLIC mdbridge::core)
target_compile_features(mdbridge_video PUBLIC cxx_std_20)
mdbridge_target_warnings(mdbridge_video)

if(FFmpegLGPL_FOUND)
    message(STATUS "Video encoding: FFmpeg ${FFmpegLGPL_VERSION} (${FFmpegLGPL_LICENSE}, libsvtav1: "
                   "${FFmpegLGPL_HAS_SVTAV1}) at ${FFmpegLGPL_INCLUDE_DIR}/..")
    if(NOT FFmpegLGPL_HAS_SVTAV1)
        message(WARNING "The FFmpeg libraries have no libsvtav1 encoder: videos fall back to mpeg4")
    endif()
    target_link_libraries(mdbridge_video PRIVATE FFmpegLGPL::FFmpeg)
    target_compile_definitions(mdbridge_video PRIVATE MDBRIDGE_HAVE_FFMPEG=1)
    set(MDBRIDGE_VIDEO_RUNTIME_FILES "${FFmpegLGPL_RUNTIME_FILES}" CACHE INTERNAL "")
    # bundled with the programs: Linux in lib/ next to Qt, as linuxdeploy lays out the AppImage
    # (usr/lib, not GNUInstallDirs' lib64 on RHEL), Windows next to mdbridge.exe; the build
    # information goes with the licence texts
    if(MDBRIDGE_BUILD_APP)
        if(WIN32)
            install(FILES ${FFmpegLGPL_RUNTIME_FILES} DESTINATION ${CMAKE_INSTALL_BINDIR})
        else()
            install(FILES ${FFmpegLGPL_RUNTIME_FILES} DESTINATION lib)
        endif()
        if(FFmpegLGPL_BUILDINFO)
            # MDBRIDGE_DOC_DIR is set after the subdirectories: same rule as the top level
            if(WIN32)
                set(_doc .)
            else()
                set(_doc ${CMAKE_INSTALL_DATADIR}/doc/mdbridge)
            endif()
            install(FILES "${FFmpegLGPL_BUILDINFO}" DESTINATION ${_doc}/licenses RENAME FFmpeg-BUILDINFO.txt)
        endif()
    endif()
else()
    message(STATUS "Video encoding: not available (no LGPL FFmpeg libraries; docs/BUILDING.md, section Video encoding)")
    target_compile_definitions(mdbridge_video PRIVATE MDBRIDGE_HAVE_FFMPEG=0)
    set(MDBRIDGE_VIDEO_RUNTIME_FILES "" CACHE INTERNAL "")
endif()

# the frame loop (frame_loop.hpp: MDATK's render_video) on the figure's panels and the encoder; no
# Qt (frames are composed through FrameCompositor, whose QPainter implementation is in mdbridge::render)
add_library(mdbridge_video_loop STATIC
    frame_loop.cpp
    frame_loop.hpp)
add_library(mdbridge::video_loop ALIAS mdbridge_video_loop)
find_package(Threads REQUIRED)  # the composing threads (std::thread)
target_link_libraries(mdbridge_video_loop PUBLIC mdbridge::video mdbridge::render_panels mdbridge::core PRIVATE Threads::Threads)
target_compile_features(mdbridge_video_loop PUBLIC cxx_std_20)
mdbridge_target_warnings(mdbridge_video_loop)
