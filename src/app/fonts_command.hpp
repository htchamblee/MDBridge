#pragma once
// `mdbridge fonts`: the bundled figure fonts, their metrics and a text image drawn without a display
// (the check that figures render the same on a bare server as on the build machine).

#include <string>

namespace mdbridge::app {

// Prints the report (--json: as JSON) and with `renderPath` writes the text image as PNG. Exit status
// 0, or 1 when the fonts cannot be loaded or the image not written.
int fontsCommand(const char* program, bool json, const std::string& renderPath);

}  // namespace mdbridge::app
