#pragma once
// A Vulkan 1.3 device and queue without a window: what the offscreen renderer (batch `render`,
// tests, the video loop) draws with. No Qt; the GUI's window keeps QVulkanWindow.
//
//   std::string problem;
//   auto context = vulkan::Context::create({}, &problem);   // nullptr + problem: use the software renderer
//   auto target = vulkan::OffscreenTarget::create(*context, 1920, 1080, {}, &problem);
//   vulkan::Image image;
//   target->render({0, 0, 0, 1}, [&](VkCommandBuffer cmd) { ... }, &image, &problem);
//
// The device is the one probeVulkan() would choose (discrete > integrated > virtual > CPU, Vulkan
// 1.3 required) among those with a graphics queue. Dynamic rendering and synchronization2 (core in
// 1.3) are enabled; there are no render pass objects.

#include <cstdint>
#include <filesystem>
#include <functional>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <string_view>

#include "mdbridge/render/vulkan/functions.hpp"
#include "mdbridge/render/vulkan_probe.hpp"

namespace mdbridge::render::vulkan {

struct ContextOptions {
    std::string applicationName = "MDBridge";
    // the loader to open (empty: the system's, defaultLoaderName())
    std::string loaderPath;
    // only devices whose name contains this text (empty: any), e.g. "llvmpipe" for lavapipe
    std::string deviceName;
    // Khronos validation layer when installed (silently without it); also switched on by the
    // environment variable MDBRIDGE_VULKAN_VALIDATION=1. Messages go to stderr.
    bool validation = false;
    // where the pipeline cache is kept between runs (empty: in memory only). A file written for
    // another device or driver, or a damaged one, is ignored; savePipelineCache() writes it.
    std::filesystem::path pipelineCachePath;
};

class Context {
public:
    // nullptr with `problem` set when there is no loader, no driver, no Vulkan 1.3 device with a
    // graphics queue, or device creation fails; the problem reads like probeVulkan()'s.
    static std::unique_ptr<Context> create(const ContextOptions& options, std::string* problem);
    ~Context();
    Context(const Context&) = delete;
    Context& operator=(const Context&) = delete;

    // context.vk.vkCmdDraw(...): the device functions; vki: the instance functions
    const DeviceFunctions& vk = m_vk;
    const InstanceFunctions& vki = m_vki;

    VkInstance instance() const { return m_instance; }
    VkPhysicalDevice physicalDevice() const { return m_physicalDevice; }
    VkDevice device() const { return m_device; }
    VkQueue queue() const { return m_queue; }
    std::uint32_t queueFamily() const { return m_queueFamily; }
    const VkPhysicalDeviceProperties& properties() const { return m_properties; }
    const VulkanDevice& deviceInfo() const { return m_info; }  // name, type, API version
    bool validationEnabled() const { return m_debugMessenger != VK_NULL_HANDLE; }

    // index of a memory type allowed by `typeBits` with all `required` flags, preferring one that
    // also has `preferred`; nullopt if none
    std::optional<std::uint32_t> memoryType(std::uint32_t typeBits, VkMemoryPropertyFlags required,
                                            VkMemoryPropertyFlags preferred = 0) const;
    // whether optimal-tiling images of `format` support `features`
    bool supportsFormat(VkFormat format, VkFormatFeatureFlags features) const;

    // Records commands with `record` into a one-time command buffer, submits it to the queue and
    // waits for it. Thread-safe (submissions are serialised). VK_SUCCESS or the first failure.
    VkResult submitAndWait(const std::function<void(VkCommandBuffer)>& record);

    // the pipeline cache every pipeline of this context is created with
    VkPipelineCache pipelineCache() const { return m_pipelineCache; }
    // whether ContextOptions::pipelineCachePath held a cache this device accepted
    bool pipelineCacheLoaded() const { return m_pipelineCacheLoaded; }
    // writes the pipeline cache to ContextOptions::pipelineCachePath (temporary file, then
    // rename); false with `problem` set on failure. Nothing to do (true) without a path.
    bool savePipelineCache(std::string* problem) const;

private:
    Context() = default;
    std::string init(const ContextOptions& options);  // "" or the problem

    DeviceFunctions m_vk;
    InstanceFunctions m_vki;
    std::unique_ptr<Loader> m_loader;
    VkInstance m_instance = VK_NULL_HANDLE;
    VkDebugUtilsMessengerEXT m_debugMessenger = VK_NULL_HANDLE;
    PFN_vkDestroyDebugUtilsMessengerEXT m_destroyDebugMessenger = nullptr;
    VkPhysicalDevice m_physicalDevice = VK_NULL_HANDLE;
    VkPhysicalDeviceProperties m_properties{};
    VkPhysicalDeviceMemoryProperties m_memory{};
    VulkanDevice m_info;
    VkDevice m_device = VK_NULL_HANDLE;
    VkQueue m_queue = VK_NULL_HANDLE;
    std::uint32_t m_queueFamily = 0;
    VkCommandPool m_commandPool = VK_NULL_HANDLE;
    VkFence m_fence = VK_NULL_HANDLE;
    std::mutex m_submitMutex;
    VkPipelineCache m_pipelineCache = VK_NULL_HANDLE;
    bool m_pipelineCacheLoaded = false;
    std::filesystem::path m_pipelineCachePath;
};

// The pipeline cache file: MDBridge's header (magic "MDBVKPC1", payload size, SHA-256 of the
// payload) followed by the driver's cache data. Drivers check their own header but not always
// their payload, and a truncated cache can crash a driver; so the file is checked before the
// driver sees it.
std::string wrapPipelineCache(std::string_view driverData);
// The driver data of a cache file if it is intact and was written for this device and driver
// (Vulkan's header: vendor id, device id, pipelineCacheUUID); else nullopt with `why` set.
std::optional<std::string> unwrapPipelineCache(std::string_view file, const VkPhysicalDeviceProperties& device,
                                               std::string* why);

}  // namespace mdbridge::render::vulkan
