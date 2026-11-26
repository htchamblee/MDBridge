#pragma once
// Paths as the interface shows them (milestone 5 wave 3; docs/DIFFERENCES.md WN-12).
//
// MDBridge stores and compares paths in one generic form with '/' (the Windows path rules of
// docs/DIFFERENCES.md "Project files"), and std::filesystem joins a name to a folder with '\' on
// Windows, so a path built from the two reads "C:/Users/me/out\video.mp4". On Windows every path
// the interface shows (window title, recent list, notices, log lines, dialogs, status bar,
// problems) is shown with the native '\' throughout instead; storage (project files, settings,
// parameter logs, the cache) and the command line keep their forms. On Linux and macOS these
// functions return their input unchanged.
//
//   displayPath(path)     one path known as a path
//   displayPathsIn(text)  a message that may quote paths (the pipeline's log lines, problems,
//                         error texts): each absolute Windows path in it (a drive letter and a
//                         separator, or a UNC / device prefix) gets '\' up to the next white space,
//                         quote or ':' of the text; where a root folder (the home folder, the
//                         temporary folder, and the folders the window names with
//                         setDisplayRoots(): the project's and its output and cache folders) starts
//                         such a path, the root is converted as a whole, spaces in it included.

#include <filesystem>

#include <QString>
#include <QStringList>

namespace mdbridge::app::gui {

enum class PathStyle { Posix, Windows };

#ifdef _WIN32
inline constexpr PathStyle kNativePathStyle = PathStyle::Windows;
#else
inline constexpr PathStyle kNativePathStyle = PathStyle::Posix;
#endif

// The rules above for a given style (tests use the Windows rules on every platform).
QString displayPathText(const QString& path, PathStyle style);
QString displayPathsInText(const QString& text, PathStyle style, const QStringList& roots);

// With the platform's style and the roots below.
QString displayPath(const std::filesystem::path& path);
QString displayPath(const QString& path);
QString displayPathsIn(const QString& text);

// Folders whose names may hold spaces, besides the home and temporary folders (set by the window
// when a project is opened; any thread may read them).
void setDisplayRoots(const QStringList& roots);
QStringList displayRoots();

}  // namespace mdbridge::app::gui
