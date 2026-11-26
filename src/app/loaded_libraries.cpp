#include "loaded_libraries.hpp"

#if defined(__linux__)
#include <link.h>
#endif

namespace mdbridge::app {

std::vector<std::string> loadedLibraries() {
    std::vector<std::string> libraries;
#if defined(__linux__)
    dl_iterate_phdr(
        [](dl_phdr_info* info, size_t, void* data) {
            if (info->dlpi_name && info->dlpi_name[0] == '/') {  // the program itself and the vDSO have no path
                static_cast<std::vector<std::string>*>(data)->emplace_back(info->dlpi_name);
            }
            return 0;
        },
        &libraries);
#endif
    return libraries;
}

}  // namespace mdbridge::app
