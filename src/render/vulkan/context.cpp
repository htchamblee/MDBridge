#include "mdbridge/render/vulkan/context.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iterator>
#include <system_error>
#include <vector>

#include "mdbridge/core/hash.hpp"

namespace mdbridge::render::vulkan {

namespace {

const std::uint32_t kMinimum = VK_MAKE_API_VERSION(0, kVulkanMajor, kVulkanMinor, 0);
const char* const kValidationLayer = "VK_LAYER_KHRONOS_validation";

bool environmentFlag(const char* name) {
    const char* value = std::getenv(name);  // NOLINT(concurrency-mt-unsafe): read before any thread starts
    return value && *value && std::strcmp(value, "0") != 0;
}

VKAPI_ATTR VkBool32 VKAPI_CALL debugMessage(VkDebugUtilsMessageSeverityFlagBitsEXT severity,
                                            VkDebugUtilsMessageTypeFlagsEXT /*types*/,
                                            const VkDebugUtilsMessengerCallbackDataEXT* data, void* /*user*/) {
    const char* level = (severity & VK_DEBUG_UTILS_MESSAGE_SEVERITY_ERROR_BIT_EXT) ? "error" : "warning";
    std::fprintf(stderr, "MDBridge: Vulkan validation %s: %s\n", level, data && data->pMessage ? data->pMessage : "");
    return VK_FALSE;
}

bool hasInstanceLayer(const GlobalFunctions& g, const char* name) {
    std::uint32_t count = 0;
    if (!g.vkEnumerateInstanceLayerProperties || g.vkEnumerateInstanceLayerProperties(&count, nullptr) != VK_SUCCESS) {
        return false;
    }
    std::vector<VkLayerProperties> layers(count);
    g.vkEnumerateInstanceLayerProperties(&count, layers.data());
    for (const VkLayerProperties& layer : layers) {
        if (std::strcmp(layer.layerName, name) == 0) {
            return true;
        }
    }
    return false;
}

bool hasInstanceExtension(const GlobalFunctions& g, const char* layer, const char* name) {
    std::uint32_t count = 0;
    if (!g.vkEnumerateInstanceExtensionProperties ||
        g.vkEnumerateInstanceExtensionProperties(layer, &count, nullptr) != VK_SUCCESS) {
        return false;
    }
    std::vector<VkExtensionProperties> extensions(count);
    g.vkEnumerateInstanceExtensionProperties(layer, &count, extensions.data());
    for (const VkExtensionProperties& extension : extensions) {
        if (std::strcmp(extension.extensionName, name) == 0) {
            return true;
        }
    }
    return false;
}

template <typename T>
void putLittle(std::string& out, T value) {
    for (std::size_t i = 0; i < sizeof(T); ++i) {
        out.push_back(static_cast<char>((static_cast<std::uint64_t>(value) >> (8 * i)) & 0xff));
    }
}

template <typename T>
T getLittle(std::string_view in, std::size_t at) {
    std::uint64_t value = 0;
    for (std::size_t i = 0; i < sizeof(T); ++i) {
        value |= static_cast<std::uint64_t>(static_cast<unsigned char>(in[at + i])) << (8 * i);
    }
    return static_cast<T>(value);
}

constexpr std::string_view kCacheMagic = "MDBVKPC1";
constexpr std::size_t kCacheHeader = 8 + 8 + 64;  // magic, payload size, SHA-256 (hex) of the payload

}  // namespace

// ------------------------------------------------------------------------------- pipeline cache

std::string wrapPipelineCache(std::string_view driverData) {
    std::string out(kCacheMagic);
    putLittle<std::uint64_t>(out, driverData.size());
    out += core::sha256Hex(driverData);
    out += driverData;
    return out;
}

std::optional<std::string> unwrapPipelineCache(std::string_view file, const VkPhysicalDeviceProperties& device,
                                               std::string* why) {
    if (file.size() < kCacheHeader || file.substr(0, kCacheMagic.size()) != kCacheMagic) {
        *why = "not an MDBridge pipeline cache";
        return std::nullopt;
    }
    const auto size = getLittle<std::uint64_t>(file, 8);
    const std::string_view payload = file.substr(kCacheHeader);
    if (payload.size() != size) {
        *why = "truncated (" + std::to_string(payload.size()) + " of " + std::to_string(size) + " bytes)";
        return std::nullopt;
    }
    if (core::sha256Hex(payload) != file.substr(16, 64)) {
        *why = "damaged (checksum mismatch)";
        return std::nullopt;
    }
    // Vulkan's pipeline cache header, version one: header size, version, vendor id, device id, UUID
    constexpr std::size_t kVulkanHeader = 16 + VK_UUID_SIZE;
    if (payload.size() < kVulkanHeader || getLittle<std::uint32_t>(payload, 0) < kVulkanHeader ||
        getLittle<std::uint32_t>(payload, 4) != VK_PIPELINE_CACHE_HEADER_VERSION_ONE) {
        *why = "no Vulkan pipeline cache header";
        return std::nullopt;
    }
    if (getLittle<std::uint32_t>(payload, 8) != device.vendorID ||
        getLittle<std::uint32_t>(payload, 12) != device.deviceID ||
        std::memcmp(payload.data() + 16, device.pipelineCacheUUID, VK_UUID_SIZE) != 0) {
        *why = "written for another device or driver";
        return std::nullopt;
    }
    return std::string(payload);
}

// ----------------------------------------------------------------------------------- context

std::unique_ptr<Context> Context::create(const ContextOptions& options, std::string* problem) {
    std::unique_ptr<Context> context(new Context());
    *problem = context->init(options);
    if (!problem->empty()) {
        return nullptr;  // the destructor releases what init() created
    }
    return context;
}

std::string Context::init(const ContextOptions& options) {
    std::string problem;
    m_loader = Loader::open(options.loaderPath, &problem);
    if (!m_loader) {
        return problem;
    }
    const PFN_vkGetInstanceProcAddr gipa = m_loader->getInstanceProcAddr();
    GlobalFunctions g;
#define MDBRIDGE_LOAD(fn) g.fn = reinterpret_cast<PFN_##fn>(gipa(nullptr, #fn));
    MDBRIDGE_VULKAN_GLOBAL_FUNCTIONS(MDBRIDGE_LOAD)
#undef MDBRIDGE_LOAD
    if (!g.vkCreateInstance) {
        return "the Vulkan loader has no vkCreateInstance";
    }
    const std::uint32_t loaderVersion = m_loader->instanceVersion();
    if (loaderVersion < kMinimum) {
        return "the Vulkan loader supports only Vulkan " + vulkanVersionText(loaderVersion) +
               " (update the graphics driver)";
    }

    // instance, with the validation layer and its messenger when asked for and installed
    std::vector<const char*> layers;
    std::vector<const char*> extensions;
    const bool wantValidation = options.validation || environmentFlag("MDBRIDGE_VULKAN_VALIDATION");
    if (wantValidation && hasInstanceLayer(g, kValidationLayer)) {
        layers.push_back(kValidationLayer);
        if (hasInstanceExtension(g, nullptr, VK_EXT_DEBUG_UTILS_EXTENSION_NAME) ||
            hasInstanceExtension(g, kValidationLayer, VK_EXT_DEBUG_UTILS_EXTENSION_NAME)) {
            extensions.push_back(VK_EXT_DEBUG_UTILS_EXTENSION_NAME);
        }
    }
    VkApplicationInfo app{};
    app.sType = VK_STRUCTURE_TYPE_APPLICATION_INFO;
    app.pApplicationName = options.applicationName.c_str();
    app.pEngineName = "MDBridge";
    app.apiVersion = kMinimum;
    VkInstanceCreateInfo instanceInfo{};
    instanceInfo.sType = VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO;
    instanceInfo.pApplicationInfo = &app;
    instanceInfo.enabledLayerCount = static_cast<std::uint32_t>(layers.size());
    instanceInfo.ppEnabledLayerNames = layers.data();
    instanceInfo.enabledExtensionCount = static_cast<std::uint32_t>(extensions.size());
    instanceInfo.ppEnabledExtensionNames = extensions.data();
    if (const VkResult r = g.vkCreateInstance(&instanceInfo, nullptr, &m_instance); r != VK_SUCCESS) {
        m_instance = VK_NULL_HANDLE;
        return "the Vulkan loader could not create an instance (" + resultText(r) + "; is a Vulkan driver installed?)";
    }
#define MDBRIDGE_LOAD(fn) m_vki.fn = reinterpret_cast<PFN_##fn>(gipa(m_instance, #fn));
    MDBRIDGE_VULKAN_INSTANCE_FUNCTIONS(MDBRIDGE_LOAD)
#undef MDBRIDGE_LOAD
#define MDBRIDGE_CHECK(fn) \
    if (!m_vki.fn) {       \
        return "the Vulkan loader has no " #fn; \
    }
    MDBRIDGE_VULKAN_INSTANCE_FUNCTIONS(MDBRIDGE_CHECK)
#undef MDBRIDGE_CHECK
    if (!extensions.empty()) {
        auto createMessenger = reinterpret_cast<PFN_vkCreateDebugUtilsMessengerEXT>(
            gipa(m_instance, "vkCreateDebugUtilsMessengerEXT"));
        m_destroyDebugMessenger = reinterpret_cast<PFN_vkDestroyDebugUtilsMessengerEXT>(
            gipa(m_instance, "vkDestroyDebugUtilsMessengerEXT"));
        VkDebugUtilsMessengerCreateInfoEXT messenger{};
        messenger.sType = VK_STRUCTURE_TYPE_DEBUG_UTILS_MESSENGER_CREATE_INFO_EXT;
        messenger.messageSeverity =
            VK_DEBUG_UTILS_MESSAGE_SEVERITY_WARNING_BIT_EXT | VK_DEBUG_UTILS_MESSAGE_SEVERITY_ERROR_BIT_EXT;
        messenger.messageType = VK_DEBUG_UTILS_MESSAGE_TYPE_GENERAL_BIT_EXT |
                                VK_DEBUG_UTILS_MESSAGE_TYPE_VALIDATION_BIT_EXT |
                                VK_DEBUG_UTILS_MESSAGE_TYPE_PERFORMANCE_BIT_EXT;
        messenger.pfnUserCallback = debugMessage;
        if (!createMessenger || !m_destroyDebugMessenger ||
            createMessenger(m_instance, &messenger, nullptr, &m_debugMessenger) != VK_SUCCESS) {
            m_debugMessenger = VK_NULL_HANDLE;
        }
    }

    // the device probeVulkan() would choose, among those with a graphics queue (and the name asked for)
    std::uint32_t count = 0;
    m_vki.vkEnumeratePhysicalDevices(m_instance, &count, nullptr);
    std::vector<VkPhysicalDevice> handles(count);
    if (count) {
        m_vki.vkEnumeratePhysicalDevices(m_instance, &count, handles.data());
    }
    handles.resize(count);
    std::vector<VulkanDevice> candidates;
    std::vector<VkPhysicalDevice> candidateHandles;
    std::vector<std::uint32_t> candidateQueues;
    std::vector<VulkanDevice> all;
    for (VkPhysicalDevice handle : handles) {
        VkPhysicalDeviceProperties p{};
        m_vki.vkGetPhysicalDeviceProperties(handle, &p);
        const VulkanDevice info = vulkanDevice(p.deviceName, p.deviceType, p.apiVersion, p.vendorID, p.deviceID);
        all.push_back(info);
        if (!options.deviceName.empty() && info.name.find(options.deviceName) == std::string::npos) {
            continue;
        }
        std::uint32_t families = 0;
        m_vki.vkGetPhysicalDeviceQueueFamilyProperties(handle, &families, nullptr);
        std::vector<VkQueueFamilyProperties> queueFamilies(families);
        m_vki.vkGetPhysicalDeviceQueueFamilyProperties(handle, &families, queueFamilies.data());
        for (std::uint32_t q = 0; q < families; ++q) {
            if (queueFamilies[q].queueFlags & VK_QUEUE_GRAPHICS_BIT) {
                candidates.push_back(info);
                candidateHandles.push_back(handle);
                candidateQueues.push_back(q);
                break;
            }
        }
    }
    const int chosen = chooseDevice(candidates);
    if (chosen < 0) {
        if (all.empty()) {
            return "no Vulkan device (is a Vulkan driver installed?)";
        }
        std::string names;
        for (const VulkanDevice& d : all) {
            names += (names.empty() ? "" : ", ") + d.name + " (Vulkan " + vulkanVersionText(d.apiVersion) + ")";
        }
        if (!options.deviceName.empty()) {
            return "no Vulkan " + std::to_string(kVulkanMajor) + "." + std::to_string(kVulkanMinor) +
                   " device with a graphics queue named like '" + options.deviceName + "' (devices: " + names + ")";
        }
        return "no device supports Vulkan " + std::to_string(kVulkanMajor) + "." + std::to_string(kVulkanMinor) +
               " with a graphics queue (devices: " + names + ")";
    }
    const auto index = static_cast<std::size_t>(chosen);
    m_physicalDevice = candidateHandles[index];
    m_queueFamily = candidateQueues[index];
    m_info = candidates[index];
    m_vki.vkGetPhysicalDeviceProperties(m_physicalDevice, &m_properties);
    m_vki.vkGetPhysicalDeviceMemoryProperties(m_physicalDevice, &m_memory);

    // device: one graphics queue, dynamic rendering and synchronization2 (both core in 1.3)
    VkPhysicalDeviceVulkan13Features supported13{};
    supported13.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_3_FEATURES;
    VkPhysicalDeviceFeatures2 supported{};
    supported.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_FEATURES_2;
    supported.pNext = &supported13;
    m_vki.vkGetPhysicalDeviceFeatures2(m_physicalDevice, &supported);
    if (!supported13.dynamicRendering || !supported13.synchronization2) {
        return m_info.name + " lacks Vulkan 1.3's dynamic rendering or synchronization2";
    }
    VkPhysicalDeviceVulkan13Features enabled13{};
    enabled13.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_3_FEATURES;
    enabled13.dynamicRendering = VK_TRUE;
    enabled13.synchronization2 = VK_TRUE;
    VkPhysicalDeviceFeatures2 enabled{};
    enabled.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_FEATURES_2;
    enabled.pNext = &enabled13;
    enabled.features.fillModeNonSolid = supported.features.fillModeNonSolid;
    enabled.features.wideLines = supported.features.wideLines;
    enabled.features.sampleRateShading = supported.features.sampleRateShading;
    const float priority = 1.0f;
    VkDeviceQueueCreateInfo queueInfo{};
    queueInfo.sType = VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO;
    queueInfo.queueFamilyIndex = m_queueFamily;
    queueInfo.queueCount = 1;
    queueInfo.pQueuePriorities = &priority;
    VkDeviceCreateInfo deviceInfo{};
    deviceInfo.sType = VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO;
    deviceInfo.pNext = &enabled;
    deviceInfo.queueCreateInfoCount = 1;
    deviceInfo.pQueueCreateInfos = &queueInfo;
    if (const VkResult r = m_vki.vkCreateDevice(m_physicalDevice, &deviceInfo, nullptr, &m_device); r != VK_SUCCESS) {
        m_device = VK_NULL_HANDLE;
        return "could not create a Vulkan device on " + m_info.name + " (" + resultText(r) + ")";
    }
#define MDBRIDGE_LOAD(fn) m_vk.fn = reinterpret_cast<PFN_##fn>(m_vki.vkGetDeviceProcAddr(m_device, #fn));
    MDBRIDGE_VULKAN_DEVICE_FUNCTIONS(MDBRIDGE_LOAD)
#undef MDBRIDGE_LOAD
#define MDBRIDGE_CHECK(fn) \
    if (!m_vk.fn) {        \
        return "the Vulkan driver of " + m_info.name + " has no " #fn; \
    }
    MDBRIDGE_VULKAN_DEVICE_FUNCTIONS(MDBRIDGE_CHECK)
#undef MDBRIDGE_CHECK
    m_vk.vkGetDeviceQueue(m_device, m_queueFamily, 0, &m_queue);

    VkCommandPoolCreateInfo poolInfo{};
    poolInfo.sType = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO;
    poolInfo.flags = VK_COMMAND_POOL_CREATE_TRANSIENT_BIT;
    poolInfo.queueFamilyIndex = m_queueFamily;
    if (const VkResult r = m_vk.vkCreateCommandPool(m_device, &poolInfo, nullptr, &m_commandPool); r != VK_SUCCESS) {
        m_commandPool = VK_NULL_HANDLE;
        return "could not create a Vulkan command pool (" + resultText(r) + ")";
    }
    VkFenceCreateInfo fenceInfo{};
    fenceInfo.sType = VK_STRUCTURE_TYPE_FENCE_CREATE_INFO;
    if (const VkResult r = m_vk.vkCreateFence(m_device, &fenceInfo, nullptr, &m_fence); r != VK_SUCCESS) {
        m_fence = VK_NULL_HANDLE;
        return "could not create a Vulkan fence (" + resultText(r) + ")";
    }

    // pipeline cache: the file's data when it is intact and for this device, else empty
    m_pipelineCachePath = options.pipelineCachePath;
    std::string initial;
    if (!m_pipelineCachePath.empty()) {
        std::ifstream in(m_pipelineCachePath, std::ios::binary);
        if (in) {
            const std::string file((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
            std::string why;
            if (std::optional<std::string> data = unwrapPipelineCache(file, m_properties, &why)) {
                initial = std::move(*data);
            }
        }
    }
    VkPipelineCacheCreateInfo cacheInfo{};
    cacheInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_CACHE_CREATE_INFO;
    cacheInfo.initialDataSize = initial.size();
    cacheInfo.pInitialData = initial.empty() ? nullptr : initial.data();
    VkResult cache = m_vk.vkCreatePipelineCache(m_device, &cacheInfo, nullptr, &m_pipelineCache);
    if (cache != VK_SUCCESS && !initial.empty()) {  // the driver refused the data: start empty
        cacheInfo.initialDataSize = 0;
        cacheInfo.pInitialData = nullptr;
        initial.clear();
        cache = m_vk.vkCreatePipelineCache(m_device, &cacheInfo, nullptr, &m_pipelineCache);
    }
    if (cache != VK_SUCCESS) {
        m_pipelineCache = VK_NULL_HANDLE;
        return "could not create a Vulkan pipeline cache (" + resultText(cache) + ")";
    }
    m_pipelineCacheLoaded = !initial.empty();
    return {};
}

Context::~Context() {
    if (m_device) {
        m_vk.vkDeviceWaitIdle(m_device);
        if (m_pipelineCache) {
            m_vk.vkDestroyPipelineCache(m_device, m_pipelineCache, nullptr);
        }
        if (m_fence) {
            m_vk.vkDestroyFence(m_device, m_fence, nullptr);
        }
        if (m_commandPool) {
            m_vk.vkDestroyCommandPool(m_device, m_commandPool, nullptr);
        }
        if (m_vk.vkDestroyDevice) {
            m_vk.vkDestroyDevice(m_device, nullptr);
        }
    }
    if (m_instance) {
        if (m_debugMessenger && m_destroyDebugMessenger) {
            m_destroyDebugMessenger(m_instance, m_debugMessenger, nullptr);
        }
        if (m_vki.vkDestroyInstance) {
            m_vki.vkDestroyInstance(m_instance, nullptr);
        }
    }
}

std::optional<std::uint32_t> Context::memoryType(std::uint32_t typeBits, VkMemoryPropertyFlags required,
                                                 VkMemoryPropertyFlags preferred) const {
    std::optional<std::uint32_t> fallback;
    for (std::uint32_t i = 0; i < m_memory.memoryTypeCount; ++i) {
        const VkMemoryPropertyFlags flags = m_memory.memoryTypes[i].propertyFlags;
        if (!(typeBits & (1u << i)) || (flags & required) != required) {
            continue;
        }
        if ((flags & preferred) == preferred) {
            return i;
        }
        if (!fallback) {
            fallback = i;
        }
    }
    return fallback;
}

bool Context::supportsFormat(VkFormat format, VkFormatFeatureFlags features) const {
    VkFormatProperties properties{};
    m_vki.vkGetPhysicalDeviceFormatProperties(m_physicalDevice, format, &properties);
    return (properties.optimalTilingFeatures & features) == features;
}

VkResult Context::submitAndWait(const std::function<void(VkCommandBuffer)>& record) {
    const std::lock_guard lock(m_submitMutex);
    VkCommandBufferAllocateInfo allocate{};
    allocate.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO;
    allocate.commandPool = m_commandPool;
    allocate.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
    allocate.commandBufferCount = 1;
    VkCommandBuffer cmd = VK_NULL_HANDLE;
    VkResult r = m_vk.vkAllocateCommandBuffers(m_device, &allocate, &cmd);
    if (r != VK_SUCCESS) {
        return r;
    }
    VkCommandBufferBeginInfo begin{};
    begin.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
    begin.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
    r = m_vk.vkBeginCommandBuffer(cmd, &begin);
    if (r == VK_SUCCESS) {
        record(cmd);
        r = m_vk.vkEndCommandBuffer(cmd);
    }
    if (r == VK_SUCCESS) {
        VkCommandBufferSubmitInfo command{};
        command.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_SUBMIT_INFO;
        command.commandBuffer = cmd;
        VkSubmitInfo2 submit{};
        submit.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO_2;
        submit.commandBufferInfoCount = 1;
        submit.pCommandBufferInfos = &command;
        r = m_vk.vkQueueSubmit2(m_queue, 1, &submit, m_fence);
        if (r == VK_SUCCESS) {
            r = m_vk.vkWaitForFences(m_device, 1, &m_fence, VK_TRUE, UINT64_MAX);
            m_vk.vkResetFences(m_device, 1, &m_fence);
        }
    }
    m_vk.vkFreeCommandBuffers(m_device, m_commandPool, 1, &cmd);
    m_vk.vkResetCommandPool(m_device, m_commandPool, 0);
    return r;
}

bool Context::savePipelineCache(std::string* problem) const {
    if (m_pipelineCachePath.empty()) {
        return true;
    }
    std::size_t size = 0;
    VkResult r = m_vk.vkGetPipelineCacheData(m_device, m_pipelineCache, &size, nullptr);
    std::string data(size, '\0');
    if (r == VK_SUCCESS) {
        r = m_vk.vkGetPipelineCacheData(m_device, m_pipelineCache, &size, data.data());
        data.resize(size);
    }
    if (r != VK_SUCCESS) {
        *problem = "could not read the Vulkan pipeline cache (" + resultText(r) + ")";
        return false;
    }
    std::error_code error;
    if (m_pipelineCachePath.has_parent_path()) {
        std::filesystem::create_directories(m_pipelineCachePath.parent_path(), error);
    }
    std::filesystem::path temporary = m_pipelineCachePath;
    temporary += ".tmp";
    {
        std::ofstream out(temporary, std::ios::binary | std::ios::trunc);
        const std::string file = wrapPipelineCache(data);
        out.write(file.data(), static_cast<std::streamsize>(file.size()));
        if (!out.flush()) {
            *problem = "could not write " + temporary.string();
            return false;
        }
    }
    std::filesystem::rename(temporary, m_pipelineCachePath, error);
    if (error) {
        std::filesystem::remove(temporary, error);
        *problem = "could not write " + m_pipelineCachePath.string();
        return false;
    }
    return true;
}

}  // namespace mdbridge::render::vulkan
