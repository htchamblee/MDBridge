#pragma once
// `dict pubchem`'s transport in the `mdbridge` program: one HTTP(S) GET through Qt Network
// (QNetworkAccessManager in a local event loop; the system's proxy settings, redirects followed),
// only when the command runs (owner decision D7). The request logic and the messages are the
// commands library's (commands/pubchem.hpp).

#include "pubchem.hpp"

namespace mdbridge::app {

// A NetworkBackend whose get() runs a QCoreApplication when none exists yet.
NetworkBackend qtNetworkBackend();

}  // namespace mdbridge::app
