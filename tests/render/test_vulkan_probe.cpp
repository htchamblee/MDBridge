#include <gtest/gtest.h>

#include <vulkan/vulkan_core.h>

#include "mdbridge/render/view.hpp"
#include "mdbridge/render/vulkan_probe.hpp"

using namespace mdbridge::render;

namespace {
const std::uint32_t k12 = VK_MAKE_API_VERSION(0, 1, 2, 198);
const std::uint32_t k13 = VK_MAKE_API_VERSION(0, 1, 3, 275);
const std::uint32_t k14 = VK_MAKE_API_VERSION(0, 1, 4, 303);
}  // namespace

TEST(VulkanProbe, VersionText) {
    EXPECT_EQ(vulkanVersionText(k13), "1.3.275");
    EXPECT_EQ(vulkanVersionText(VK_API_VERSION_1_0), "1.0.0");
}

TEST(VulkanProbe, DevicesNeedVulkan13) {
    EXPECT_FALSE(vulkanDevice("old GPU", VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU, k12, 0, 0).meetsMinimum());
    EXPECT_TRUE(vulkanDevice("GPU", VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU, k13, 0, 0).meetsMinimum());
    EXPECT_TRUE(vulkanDevice("GPU", VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU, k14, 0, 0).meetsMinimum());
    EXPECT_EQ(vulkanDevice("llvmpipe", VK_PHYSICAL_DEVICE_TYPE_CPU, k13, 0, 0).type, "CPU");
}

TEST(VulkanProbe, PrefersHardwareDevicesWithVulkan13) {
    const std::vector<VulkanDevice> devices = {
        vulkanDevice("llvmpipe", VK_PHYSICAL_DEVICE_TYPE_CPU, k13, 0, 0),
        vulkanDevice("old discrete", VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU, k12, 0, 0),
        vulkanDevice("integrated", VK_PHYSICAL_DEVICE_TYPE_INTEGRATED_GPU, k13, 0, 0),
    };
    EXPECT_EQ(chooseDevice(devices), 2);
    EXPECT_EQ(chooseDevice({devices[0], devices[1]}), 0);  // a CPU implementation beats no Vulkan 1.3
    EXPECT_EQ(chooseDevice({devices[1]}), -1);
    EXPECT_EQ(chooseDevice({}), -1);
}

TEST(VulkanProbe, ReportIsConsistent) {
    const VulkanReport report = probeVulkan();  // whatever this machine has
    EXPECT_EQ(report.usable(), report.problem.empty()) << report.problem;
    if (report.usable()) {
        EXPECT_TRUE(report.devices.at(static_cast<std::size_t>(report.chosen)).meetsMinimum());
    }
}

TEST(Renderer, Names) {
    for (Renderer r : {Renderer::Auto, Renderer::Vulkan, Renderer::Software}) {
        EXPECT_EQ(rendererFromName(rendererName(r)), r);
    }
    EXPECT_FALSE(rendererFromName("vtk").has_value());
}
