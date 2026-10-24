#include "mdbridge/render/vulkan_probe.hpp"

#include <algorithm>

#include <QLibrary>

#ifndef VK_NO_PROTOTYPES
#define VK_NO_PROTOTYPES
#endif
#include <vulkan/vulkan.h>

namespace mdbridge::render {

namespace {

const std::uint32_t kMinimum = VK_MAKE_API_VERSION(0, kVulkanMajor, kVulkanMinor, 0);

std::pair<const char*, int> deviceType(VkPhysicalDeviceType type) {
    switch (type) {
    case VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU:
        return {"discrete GPU", 4};
    case VK_PHYSICAL_DEVICE_TYPE_INTEGRATED_GPU:
        return {"integrated GPU", 3};
    case VK_PHYSICAL_DEVICE_TYPE_VIRTUAL_GPU:
        return {"virtual GPU", 2};
    case VK_PHYSICAL_DEVICE_TYPE_CPU:
        return {"CPU", 1};
    default:
        return {"other", 0};
    }
}

std::string required() { return std::to_string(kVulkanMajor) + "." + std::to_string(kVulkanMinor); }

}  // namespace

bool VulkanDevice::meetsMinimum() const { return apiVersion >= kMinimum; }

VulkanDevice vulkanDevice(const char* name, int physicalDeviceType, std::uint32_t apiVersion, std::uint32_t vendorId,
                          std::uint32_t deviceId) {
    const auto [type, rank] = deviceType(static_cast<VkPhysicalDeviceType>(physicalDeviceType));
    return {name, type, apiVersion, vendorId, deviceId, rank};
}

std::string vulkanVersionText(std::uint32_t version) {
    return std::to_string(VK_API_VERSION_MAJOR(version)) + "." + std::to_string(VK_API_VERSION_MINOR(version)) + "." +
           std::to_string(VK_API_VERSION_PATCH(version));
}

int chooseDevice(const std::vector<VulkanDevice>& devices) {
    int best = -1;
    for (int i = 0; i < static_cast<int>(devices.size()); ++i) {
        if (devices[i].meetsMinimum() && (best < 0 || devices[i].rank > devices[best].rank)) {
            best = i;
        }
    }
    return best;
}

VulkanReport probeVulkan() {
    VulkanReport report;
#ifdef Q_OS_WIN
    QLibrary loader(QStringLiteral("vulkan-1"));
#else
    QLibrary loader(QStringLiteral("vulkan"), 1);
#endif
    auto getInstanceProcAddr = reinterpret_cast<PFN_vkGetInstanceProcAddr>(loader.resolve("vkGetInstanceProcAddr"));
    if (!getInstanceProcAddr) {
        report.problem = "no Vulkan loader (" + loader.errorString().toStdString() + ")";
        return report;
    }
    report.loaderFound = true;

    auto enumerateInstanceVersion = reinterpret_cast<PFN_vkEnumerateInstanceVersion>(
        getInstanceProcAddr(nullptr, "vkEnumerateInstanceVersion"));
    report.loaderVersion = VK_API_VERSION_1_0;  // a 1.0 loader has no vkEnumerateInstanceVersion
    if (enumerateInstanceVersion) {
        enumerateInstanceVersion(&report.loaderVersion);
    }

    auto createInstance = reinterpret_cast<PFN_vkCreateInstance>(getInstanceProcAddr(nullptr, "vkCreateInstance"));
    VkApplicationInfo app{};
    app.sType = VK_STRUCTURE_TYPE_APPLICATION_INFO;
    app.pApplicationName = "MDBridge (Vulkan probe)";
    app.apiVersion = report.loaderVersion < kMinimum ? report.loaderVersion : kMinimum;
    VkInstanceCreateInfo info{};
    info.sType = VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO;
    info.pApplicationInfo = &app;
    VkInstance instance = VK_NULL_HANDLE;
    const VkResult created = createInstance ? createInstance(&info, nullptr, &instance) : VK_ERROR_INITIALIZATION_FAILED;
    if (created != VK_SUCCESS) {
        report.problem = "the Vulkan loader could not create an instance (VkResult " + std::to_string(created) +
                         "; is a Vulkan driver installed?)";
        return report;
    }

    auto destroyInstance = reinterpret_cast<PFN_vkDestroyInstance>(getInstanceProcAddr(instance, "vkDestroyInstance"));
    auto enumerateDevices = reinterpret_cast<PFN_vkEnumeratePhysicalDevices>(
        getInstanceProcAddr(instance, "vkEnumeratePhysicalDevices"));
    auto deviceProperties = reinterpret_cast<PFN_vkGetPhysicalDeviceProperties>(
        getInstanceProcAddr(instance, "vkGetPhysicalDeviceProperties"));
    std::uint32_t count = 0;
    if (enumerateDevices && deviceProperties && enumerateDevices(instance, &count, nullptr) == VK_SUCCESS) {
        std::vector<VkPhysicalDevice> handles(count);
        enumerateDevices(instance, &count, handles.data());
        for (VkPhysicalDevice handle : handles) {
            VkPhysicalDeviceProperties properties{};
            deviceProperties(handle, &properties);
            report.devices.push_back(vulkanDevice(properties.deviceName, properties.deviceType, properties.apiVersion,
                                                  properties.vendorID, properties.deviceID));
        }
    }
    if (destroyInstance) {
        destroyInstance(instance, nullptr);
    }

    report.chosen = chooseDevice(report.devices);
    if (report.devices.empty()) {
        report.problem = "no Vulkan device (is a Vulkan driver installed?)";
    } else if (report.chosen < 0) {
        const auto newest = std::max_element(report.devices.begin(), report.devices.end(),
                                             [](const auto& a, const auto& b) { return a.apiVersion < b.apiVersion; });
        report.problem = "no device supports Vulkan " + required() + " (newest: " + newest->name + ", Vulkan " +
                         vulkanVersionText(newest->apiVersion) + ")";
    } else if (report.loaderVersion < kMinimum) {
        report.chosen = -1;
        report.problem = "the Vulkan loader supports only Vulkan " + vulkanVersionText(report.loaderVersion) +
                         " (update the graphics driver)";
    }
    return report;
}

}  // namespace mdbridge::render
