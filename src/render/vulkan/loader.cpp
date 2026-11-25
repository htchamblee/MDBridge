#include "mdbridge/render/vulkan/loader.hpp"

#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif

namespace mdbridge::render::vulkan {

namespace {

#ifdef _WIN32
std::wstring wide(const std::string& utf8) {
    const int n = MultiByteToWideChar(CP_UTF8, 0, utf8.data(), static_cast<int>(utf8.size()), nullptr, 0);
    std::wstring out(static_cast<std::size_t>(n), L'\0');
    MultiByteToWideChar(CP_UTF8, 0, utf8.data(), static_cast<int>(utf8.size()), out.data(), n);
    return out;
}

std::string lastError() {
    const DWORD code = GetLastError();
    char* text = nullptr;
    const DWORD length = FormatMessageA(
        FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS, nullptr, code, 0,
        reinterpret_cast<LPSTR>(&text), 0, nullptr);
    std::string message = length ? std::string(text, length) : "error " + std::to_string(code);
    LocalFree(text);
    while (!message.empty() && (message.back() == '\n' || message.back() == '\r' || message.back() == '.')) {
        message.pop_back();
    }
    return message;
}
#endif

}  // namespace

const char* defaultLoaderName() {
#ifdef _WIN32
    return "vulkan-1.dll";
#else
    return "libvulkan.so.1";
#endif
}

std::unique_ptr<Loader> Loader::open(const std::string& path, std::string* problem) {
    const std::string name = path.empty() ? defaultLoaderName() : path;
#ifdef _WIN32
    HMODULE handle = LoadLibraryW(wide(name).c_str());
    if (!handle) {
        *problem = "no Vulkan loader (" + name + ": " + lastError() + ")";
        return nullptr;
    }
    auto entry = reinterpret_cast<PFN_vkGetInstanceProcAddr>(
        reinterpret_cast<void (*)()>(GetProcAddress(handle, "vkGetInstanceProcAddr")));
    if (!entry) {
        *problem = "no Vulkan loader (" + name + " has no vkGetInstanceProcAddr)";
        FreeLibrary(handle);
        return nullptr;
    }
    return std::unique_ptr<Loader>(new Loader(handle, entry));
#else
    void* handle = dlopen(name.c_str(), RTLD_NOW | RTLD_LOCAL);
    if (!handle) {
        const char* error = dlerror();
        *problem = "no Vulkan loader (" + std::string(error ? error : name + ": cannot be opened") + ")";
        return nullptr;
    }
    auto entry = reinterpret_cast<PFN_vkGetInstanceProcAddr>(dlsym(handle, "vkGetInstanceProcAddr"));
    if (!entry) {
        *problem = "no Vulkan loader (" + name + " has no vkGetInstanceProcAddr)";
        dlclose(handle);
        return nullptr;
    }
    return std::unique_ptr<Loader>(new Loader(handle, entry));
#endif
}

Loader::~Loader() {
#ifdef _WIN32
    FreeLibrary(static_cast<HMODULE>(m_handle));
#else
    dlclose(m_handle);
#endif
}

std::uint32_t Loader::instanceVersion() const {
    auto enumerate = reinterpret_cast<PFN_vkEnumerateInstanceVersion>(
        m_getInstanceProcAddr(nullptr, "vkEnumerateInstanceVersion"));
    std::uint32_t version = VK_API_VERSION_1_0;  // a 1.0 loader has no vkEnumerateInstanceVersion
    if (enumerate && enumerate(&version) != VK_SUCCESS) {
        version = VK_API_VERSION_1_0;
    }
    return version;
}

std::string resultText(VkResult result) {
    switch (result) {
#define MDBRIDGE_VK_RESULT(name) \
    case name:                   \
        return #name;
        MDBRIDGE_VK_RESULT(VK_SUCCESS)
        MDBRIDGE_VK_RESULT(VK_NOT_READY)
        MDBRIDGE_VK_RESULT(VK_TIMEOUT)
        MDBRIDGE_VK_RESULT(VK_INCOMPLETE)
        MDBRIDGE_VK_RESULT(VK_ERROR_OUT_OF_HOST_MEMORY)
        MDBRIDGE_VK_RESULT(VK_ERROR_OUT_OF_DEVICE_MEMORY)
        MDBRIDGE_VK_RESULT(VK_ERROR_INITIALIZATION_FAILED)
        MDBRIDGE_VK_RESULT(VK_ERROR_DEVICE_LOST)
        MDBRIDGE_VK_RESULT(VK_ERROR_MEMORY_MAP_FAILED)
        MDBRIDGE_VK_RESULT(VK_ERROR_LAYER_NOT_PRESENT)
        MDBRIDGE_VK_RESULT(VK_ERROR_EXTENSION_NOT_PRESENT)
        MDBRIDGE_VK_RESULT(VK_ERROR_FEATURE_NOT_PRESENT)
        MDBRIDGE_VK_RESULT(VK_ERROR_INCOMPATIBLE_DRIVER)
        MDBRIDGE_VK_RESULT(VK_ERROR_TOO_MANY_OBJECTS)
        MDBRIDGE_VK_RESULT(VK_ERROR_FORMAT_NOT_SUPPORTED)
        MDBRIDGE_VK_RESULT(VK_ERROR_FRAGMENTED_POOL)
        MDBRIDGE_VK_RESULT(VK_ERROR_UNKNOWN)
        MDBRIDGE_VK_RESULT(VK_ERROR_OUT_OF_POOL_MEMORY)
#undef MDBRIDGE_VK_RESULT
    default:
        return "VkResult " + std::to_string(static_cast<int>(result));
    }
}

}  // namespace mdbridge::render::vulkan
