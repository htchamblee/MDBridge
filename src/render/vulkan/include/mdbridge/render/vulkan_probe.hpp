#pragma once
// What Vulkan the machine offers, found without opening a window (works over ssh, in CI).
// The 3D view needs Vulkan 1.3; machines without it get the software view.

#include <cstdint>
#include <string>
#include <vector>

namespace mdbridge::render {

inline constexpr int kVulkanMajor = 1;
inline constexpr int kVulkanMinor = 3;

struct VulkanDevice {
    std::string name;
    std::string type;  // "discrete GPU", "integrated GPU", "virtual GPU", "CPU" or "other"
    std::uint32_t apiVersion = 0;
    std::uint32_t vendorId = 0;
    std::uint32_t deviceId = 0;
    int rank = 0;  // preference: discrete 4, integrated 3, virtual 2, CPU 1, other 0

    bool meetsMinimum() const;
};

struct VulkanReport {
    bool loaderFound = false;
    std::uint32_t loaderVersion = 0;  // highest API version the loader supports
    std::vector<VulkanDevice> devices;
    int chosen = -1;      // index of the device the 3D view uses, -1: none qualifies
    std::string problem;  // why the 3D view cannot use Vulkan ("" if a device qualifies)

    bool usable() const { return chosen >= 0; }
};

// Loads the system's Vulkan loader (the graphics driver's; never bundled), creates an instance
// and lists the devices.
VulkanReport probeVulkan();

// A device from its VkPhysicalDeviceProperties fields (`physicalDeviceType`: a VkPhysicalDeviceType).
VulkanDevice vulkanDevice(const char* name, int physicalDeviceType, std::uint32_t apiVersion, std::uint32_t vendorId,
                          std::uint32_t deviceId);

// Index of the preferred device that supports Vulkan 1.3, or -1.
int chooseDevice(const std::vector<VulkanDevice>& devices);

// "1.3.275" for a packed Vulkan version number.
std::string vulkanVersionText(std::uint32_t version);

}  // namespace mdbridge::render
