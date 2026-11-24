// TEMPORARY STUB (replaced by the project model)
#include "mdbridge/core/project.hpp"

#include <stdexcept>

namespace mdbridge::core {

Value toPlain(const Project&) { throw std::logic_error("toPlain: not implemented yet"); }

Project projectFromPlain(const Value&, std::string_view, std::vector<std::string>*) {
    throw std::logic_error("projectFromPlain: not implemented yet");
}

std::optional<std::string_view> modelComment(std::string_view) { return std::nullopt; }

}  // namespace mdbridge::core
