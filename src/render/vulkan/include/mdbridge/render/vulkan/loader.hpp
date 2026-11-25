#pragma once
// The system's Vulkan loader (libvulkan.so.1 / vulkan-1.dll), opened at run time. MDBridge never
// links or bundles a loader: it comes with the graphics driver, and machines without one get the
// software renderer.

#include <cstdint>
#include <memory>
#include <string>

#ifndef VK_NO_PROTOTYPES
#define VK_NO_PROTOTYPES
#endif
#include <vulkan/vulkan.h>

namespace mdbridge::render::vulkan {

// "vulkan-1.dll" on Windows, "libvulkan.so.1" elsewhere: the name the driver installs.
const char* defaultLoaderName();

class Loader {
public:
    // Opens `path` (a file name searched as the system searches libraries, or a full path; empty:
    // defaultLoaderName()). nullptr with `problem` set ("no Vulkan loader (<system message>)") if
    // it is missing or has no vkGetInstanceProcAddr.
    static std::unique_ptr<Loader> open(const std::string& path, std::string* problem);
    ~Loader();
    Loader(const Loader&) = delete;
    Loader& operator=(const Loader&) = delete;

    PFN_vkGetInstanceProcAddr getInstanceProcAddr() const { return m_getInstanceProcAddr; }

    // The highest instance API version the loader supports (1.0 for a loader without
    // vkEnumerateInstanceVersion).
    std::uint32_t instanceVersion() const;

private:
    Loader(void* handle, PFN_vkGetInstanceProcAddr getInstanceProcAddr)
        : m_handle(handle), m_getInstanceProcAddr(getInstanceProcAddr) {}

    void* m_handle;
    PFN_vkGetInstanceProcAddr m_getInstanceProcAddr;
};

// "VK_ERROR_OUT_OF_DEVICE_MEMORY" for the VkResult values a renderer meets, else "VkResult <n>".
std::string resultText(VkResult result);

}  // namespace mdbridge::render::vulkan
