#pragma once
// Reader warnings: what a reader defaulted, skipped or guessed while reading a file (a placeholder
// cell, unknown elements, malformed CONECT records, a time step the file does not have ...).
// MDATK loads with every Python warning ignored, so its users never learn about these; MDBridge
// collects them and the front ends print each once (CLI: "WARNING: <message>" on stderr, GUI: the
// log). docs/DIFFERENCES.md, "Reader warnings".
//
// For reader authors (any format):
//
//     io::readerWarning(fileText + ": 3 CONECT records are malformed and were skipped (first: line 12)");
//
// The message starts with the file as given (detail::pathText) and names the line or record where
// one is known, and says what was done. Issue one warning per kind and file, with a count and the
// first place, not one per atom or frame. readerWarning() may be called from any thread; it goes to
// the innermost ReaderWarningScope open on the calling thread. io::loadSystem() and
// io::readTopologyFile() open one while they open the files, so warnings issued then end up in
// core::System::warnings / TopologyFile::warnings, which the front ends print (the `info` command,
// analysis::Pipeline::system(), `export-pdb`). Outside any scope (a frame decoded later, during an
// analysis) the warning goes to core::logWarning() at once. Readers that decode frames lazily
// should remember what they already reported, so that a warning is not repeated for every frame.

#include <string>
#include <vector>

namespace mdbridge::io {

// Reports a reader warning (see above). Exact repeats within one scope are dropped.
void readerWarning(std::string message);

// Collects the reader warnings issued on this thread while it exists. Scopes nest: a warning goes
// to the innermost one only. Destroy scopes in reverse order of creation (use them as locals).
class ReaderWarningScope {
public:
    ReaderWarningScope();
    ~ReaderWarningScope();
    ReaderWarningScope(const ReaderWarningScope&) = delete;
    ReaderWarningScope& operator=(const ReaderWarningScope&) = delete;

    // The warnings collected so far, in the order they were issued; the scope is emptied.
    std::vector<std::string> take();

private:
    friend void readerWarning(std::string message);
    std::vector<std::string> m_warnings;
    ReaderWarningScope* m_outer;
};

}  // namespace mdbridge::io
