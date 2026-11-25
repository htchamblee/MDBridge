#pragma once
// Windows path rules for the paths stored in project files, on UTF-8 strings: MDBridge's own rules,
// modelled on Python's ntpath. MDATK 0.0.4.1 has a Windows build too (a beta, packaging/windows:
// MSYS2's Python in a PyInstaller bundle), where Project.resolve_paths uses ntpath as it is:
// expanduser, then the path itself when isabs, else normpath(join(folder, path)), with '\' in the
// result and absolute paths left unnormalised. CI cannot run it; MDBridge's rules point at the same
// file as Python 3.13's ntpath except where docs/DIFFERENCES.md ("Windows paths") lists otherwise
// (a drive-relative path on another drive, "~user"; Python 3.12 and older keep "\x" without a
// drive). Pure string functions: the same results on every host (the project code uses them on
// Windows only; Linux and macOS keep posixpath, as MDATK).
//
// The rules:
// - '/' and '\' both separate.
// - A drive is an ASCII letter and ':' ("C:"), kept as given (no case change). A UNC drive is
//   "\\server\share" or "//server/share" (any mix of separators; also "\\?\UNC\server\share").
// - Absolute (isabs): a path with a root ("C:\x", "C:/x", "\x", "/x") or a UNC drive
//   ("\\server\share\x"). A root without a drive ("\x") takes the drive of the base directory.
// - Drive-relative ("C:x"): against the base directory when that is on the same drive (drive
//   letters compared without case), otherwise against that drive's root ("C:/x").
// - "~" and "~\x" (or "~/x"): the profile folder, %USERPROFILE%, else %HOMEDRIVE%%HOMEPATH%, as
//   Python's ntpath; unchanged when neither is set. "~user" stays as it is.
// - Results are normalised and in generic form: '/' separators, no "." and no ".." components
//   (a ".." never goes above the drive or UNC root), no doubled or trailing separator except the
//   root itself ("C:/").
// - Device paths ("\\?\C:\x", "\\.\pipe\x", either separator) are absolute and left exactly as
//   given: Windows does not normalise them, so neither does MDBridge.
// - Comparisons of stored text stay case-sensitive, as in MDATK 0.0.4.1 (it never calls normcase);
//   only the "same drive?" test of a drive-relative path ignores the case of the letter, as
//   ntpath.join does.

#include <optional>
#include <string>
#include <string_view>

namespace mdbridge::core::ntpath {

// '/' or '\'
constexpr bool isSep(char c) { return c == '/' || c == '\\'; }

// The parts of a path, as views into it (Python 3.12's ntpath.splitroot):
// "C:\a\b" -> ("C:", "\", "a\b"); "C:a" -> ("C:", "", "a"); "\a" -> ("", "\", "a");
// "//srv/share/a" -> ("//srv/share", "/", "a"); "a/b" -> ("", "", "a/b").
struct Parts {
    std::string_view drive;
    std::string_view root;
    std::string_view tail;
};
Parts splitroot(std::string_view path);

// "\\?\..." or "\\.\..." (either separator): a Win32 device path, never normalised.
bool isDevicePath(std::string_view path);

// A root, or a UNC drive (see the rules above). "C:x" and "x" are relative.
bool isabs(std::string_view path);

// '\' replaced by '/'.
std::string generic(std::string_view path);

// ntpath.join(base, path), separators kept as given: a path with a root replaces the base's root
// and tail (keeping the base's drive when it has none); a path on another drive replaces the base;
// otherwise the path is appended with a separator when needed.
std::string join(std::string_view base, std::string_view path);

// Normalised, generic form (see the rules); "." for an empty result without drive or root.
// Device paths unchanged.
std::string normpath(std::string_view path);

// The absolute form of a stored path against `base` (the project's folder, or the current
// directory): normpath(join(base, path)), a drive-relative result put at its drive's root.
// "~" is not expanded here (see resolve).
std::string abspath(std::string_view path, std::string_view base);

// os.path.dirname: everything before the last component, trailing separators removed (but not
// the root): "C:/a/b" -> "C:/a", "C:/a" -> "C:/", "//srv/share/a" -> "//srv/share/".
std::string dirname(std::string_view path);

// The profile folder from the variables ntpath.expanduser reads: USERPROFILE when set, else
// join(HOMEDRIVE or "", HOMEPATH) when HOMEPATH is set, else none.
std::optional<std::string> profileHome(const std::optional<std::string>& userProfile,
                                       const std::optional<std::string>& homeDrive,
                                       const std::optional<std::string>& homePath);
// profileHome of this process's environment.
std::optional<std::string> environmentHome();

// "~" or "~" followed by a separator: `home` in place of "~" (path unchanged when home is none);
// anything else, "~user" included, unchanged.
std::string expanduser(std::string_view path, const std::optional<std::string>& home);
// With the home of the environment.
std::string expanduser(std::string_view path);

// A stored project path made absolute against the project's folder: "" stays "", "~" expanded
// with `home`, then abspath(path, baseDir).
std::string resolve(std::string_view path, std::string_view baseDir, const std::optional<std::string>& home);
// With the home of the environment.
std::string resolve(std::string_view path, std::string_view baseDir);

}  // namespace mdbridge::core::ntpath
