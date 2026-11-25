// Vulkan without a window (src/render/vulkan): loader, device and queue, offscreen target and
// read-back, pipelines and the pipeline cache. The GPU tests run on whatever Vulkan 1.3 driver the
// machine has (lavapipe in CI) and are skipped without one (vulkan_harness.hpp).
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <thread>

#include "mdbridge/render/vulkan/context.hpp"
#include "mdbridge/render/vulkan/loader.hpp"
#include "mdbridge/render/vulkan/offscreen.hpp"
#include "mdbridge/render/vulkan/pipeline.hpp"
#include "shaders/spirv.hpp"
#include "vulkan_harness.hpp"

using namespace mdbridge::render;
using namespace mdbridge::render::vulkan;

namespace {

VkPhysicalDeviceProperties deviceProperties(std::uint32_t vendor, std::uint32_t device, unsigned char uuidByte) {
    VkPhysicalDeviceProperties p{};
    p.vendorID = vendor;
    p.deviceID = device;
    for (unsigned char& b : p.pipelineCacheUUID) {
        b = uuidByte;
    }
    return p;
}

// a driver's cache data: Vulkan's header (version one) and some payload
std::string driverCache(const VkPhysicalDeviceProperties& p, std::string_view payload) {
    std::string out;
    const auto put = [&](std::uint32_t v) {
        for (int i = 0; i < 4; ++i) {
            out.push_back(static_cast<char>((v >> (8 * i)) & 0xff));
        }
    };
    put(16 + VK_UUID_SIZE);
    put(VK_PIPELINE_CACHE_HEADER_VERSION_ONE);
    put(p.vendorID);
    put(p.deviceID);
    out.append(reinterpret_cast<const char*>(p.pipelineCacheUUID), VK_UUID_SIZE);
    out += payload;
    return out;
}

std::uint8_t unorm8(float v) { return static_cast<std::uint8_t>(std::lround(v * 255.0f)); }

std::filesystem::path tempPath(const char* name) {
    return std::filesystem::temp_directory_path() / (std::string("mdbridge-test-") + name);
}

}  // namespace

// ------------------------------------------------------------------------- without a GPU

TEST(VulkanLoader, MissingLoaderIsAProblemNotACrash) {
    std::string problem;
    EXPECT_EQ(Loader::open("mdbridge-no-such-vulkan-loader-library", &problem), nullptr);
    EXPECT_EQ(problem.rfind("no Vulkan loader (", 0), 0u) << problem;
    EXPECT_NE(problem.find("mdbridge-no-such-vulkan-loader-library"), std::string::npos) << problem;

    ContextOptions options;
    options.loaderPath = "mdbridge-no-such-vulkan-loader-library";
    problem.clear();
    EXPECT_EQ(Context::create(options, &problem), nullptr);
    EXPECT_EQ(problem.rfind("no Vulkan loader (", 0), 0u) << problem;
}

TEST(VulkanLoader, ResultText) {
    EXPECT_EQ(resultText(VK_SUCCESS), "VK_SUCCESS");
    EXPECT_EQ(resultText(VK_ERROR_OUT_OF_DEVICE_MEMORY), "VK_ERROR_OUT_OF_DEVICE_MEMORY");
    EXPECT_EQ(resultText(static_cast<VkResult>(-123456)), "VkResult -123456");
}

TEST(VulkanPipelineCacheFile, RoundTrip) {
    const VkPhysicalDeviceProperties device = deviceProperties(0x10005, 7, 0xab);
    const std::string data = driverCache(device, "driver blob");
    const std::string file = wrapPipelineCache(data);
    EXPECT_EQ(file.substr(0, 8), "MDBVKPC1");
    std::string why;
    const std::optional<std::string> back = unwrapPipelineCache(file, device, &why);
    ASSERT_TRUE(back) << why;
    EXPECT_EQ(*back, data);
}

TEST(VulkanPipelineCacheFile, RejectsForeignDamagedAndTruncatedFiles) {
    const VkPhysicalDeviceProperties device = deviceProperties(0x10de, 0x2204, 1);
    const std::string file = wrapPipelineCache(driverCache(device, std::string(100, 'x')));
    std::string why;
    EXPECT_FALSE(unwrapPipelineCache("", device, &why));
    EXPECT_EQ(why, "not an MDBridge pipeline cache");
    EXPECT_FALSE(unwrapPipelineCache(driverCache(device, "raw driver data"), device, &why));  // no MDBridge header
    EXPECT_EQ(why, "not an MDBridge pipeline cache");
    EXPECT_FALSE(unwrapPipelineCache(file.substr(0, file.size() - 1), device, &why));
    EXPECT_EQ(why.rfind("truncated", 0), 0u) << why;
    std::string damaged = file;
    damaged.back() = 'y';
    EXPECT_FALSE(unwrapPipelineCache(damaged, device, &why));
    EXPECT_EQ(why, "damaged (checksum mismatch)");
    EXPECT_FALSE(unwrapPipelineCache(file, deviceProperties(0x10de, 0x2204, 2), &why));  // driver update: new UUID
    EXPECT_EQ(why, "written for another device or driver");
    EXPECT_FALSE(unwrapPipelineCache(file, deviceProperties(0x10de, 0x2205, 1), &why));
    EXPECT_FALSE(unwrapPipelineCache(wrapPipelineCache("short"), device, &why));
    EXPECT_EQ(why, "no Vulkan pipeline cache header");
}

// ------------------------------------------------------------------------------ on a GPU

TEST(VulkanContext, DeviceWithVulkan13) {
    MDBRIDGE_VULKAN_CONTEXT_OR_SKIP(context);
    EXPECT_FALSE(context->deviceInfo().name.empty());
    EXPECT_TRUE(context->deviceInfo().meetsMinimum());
    EXPECT_NE(context->device(), VK_NULL_HANDLE);
    EXPECT_NE(context->queue(), VK_NULL_HANDLE);
    EXPECT_NE(context->pipelineCache(), VK_NULL_HANDLE);
    std::printf("Vulkan device: %s (%s, Vulkan %s)%s\n", context->deviceInfo().name.c_str(),
                context->deviceInfo().type.c_str(), vulkanVersionText(context->deviceInfo().apiVersion).c_str(),
                context->validationEnabled() ? ", validation layer on" : "");
}

TEST(VulkanContext, UnknownDeviceNameIsAProblem) {
    MDBRIDGE_VULKAN_CONTEXT_OR_SKIP(context);
    ContextOptions options = vulkan_test::testOptions();
    options.deviceName = "no device is called like this";
    std::string problem;
    EXPECT_EQ(Context::create(options, &problem), nullptr);
    EXPECT_NE(problem.find("named like 'no device is called like this'"), std::string::npos) << problem;
    EXPECT_NE(problem.find(context->deviceInfo().name), std::string::npos) << problem;  // lists the devices
}

TEST(VulkanOffscreen, ClearColourReadBack) {
    MDBRIDGE_VULKAN_CONTEXT_OR_SKIP(context);
    std::string problem;
    for (const auto& [w, h] : {std::pair{64u, 48u}, std::pair{37u, 23u}, std::pair{1u, 1u}}) {
        auto target = OffscreenTarget::create(*context, w, h, {}, &problem);
        ASSERT_TRUE(target) << problem;
        // exact in 8 bits: k / 255
        const std::array<float, 4> clear = {51.0f / 255, 102.0f / 255, 204.0f / 255, 1.0f};
        Image image;
        ASSERT_TRUE(target->render(clear, {}, &image, &problem)) << problem;
        ASSERT_EQ(image.width, w);
        ASSERT_EQ(image.height, h);
        ASSERT_EQ(image.rgba.size(), std::size_t{w} * h * 4);
        std::size_t wrong = 0;
        for (std::size_t i = 0; i < image.rgba.size(); i += 4) {
            const bool differs = image.rgba[i] != 51 || image.rgba[i + 1] != 102 || image.rgba[i + 2] != 204 ||
                                 image.rgba[i + 3] != 255;
            wrong += differs ? 1u : 0u;
        }
        EXPECT_EQ(wrong, 0u) << w << " x " << h << ": first pixel " << int(image.rgba[0]) << " " << int(image.rgba[1])
                             << " " << int(image.rgba[2]) << " " << int(image.rgba[3]);
    }
}

TEST(VulkanOffscreen, RendersAgainIntoTheSameTarget) {
    MDBRIDGE_VULKAN_CONTEXT_OR_SKIP(context);
    std::string problem;
    auto target = OffscreenTarget::create(*context, 16, 8, {}, &problem);
    ASSERT_TRUE(target) << problem;
    Image image;
    for (const float red : {0.0f, 1.0f, 128.0f / 255}) {
        ASSERT_TRUE(target->render({red, 0.0f, 0.0f, 0.0f}, {}, &image, &problem)) << problem;
        EXPECT_EQ(image.pixel(15, 7)[0], unorm8(red));
        EXPECT_EQ(image.pixel(0, 0)[3], 0);  // transparent clear kept as is
    }
}

TEST(VulkanOffscreen, TargetWithoutDepth) {
    MDBRIDGE_VULKAN_CONTEXT_OR_SKIP(context);
    std::string problem;
    OffscreenOptions options;
    options.depth = false;
    auto target = OffscreenTarget::create(*context, 8, 8, options, &problem);
    ASSERT_TRUE(target) << problem;
    EXPECT_EQ(target->depthFormat(), VK_FORMAT_UNDEFINED);
    Image image;
    ASSERT_TRUE(target->render({0.0f, 1.0f, 0.0f, 1.0f}, {}, &image, &problem)) << problem;
    EXPECT_EQ(image.pixel(3, 3)[1], 255);
}

TEST(VulkanOffscreen, RefusesImpossibleSizes) {
    MDBRIDGE_VULKAN_CONTEXT_OR_SKIP(context);
    std::string problem;
    EXPECT_EQ(OffscreenTarget::create(*context, 0, 10, {}, &problem), nullptr);
    EXPECT_NE(problem.find("0 x 10 px"), std::string::npos) << problem;
    const std::uint32_t tooLarge = context->properties().limits.maxImageDimension2D + 1;
    EXPECT_EQ(OffscreenTarget::create(*context, tooLarge, 10, {}, &problem), nullptr);
    EXPECT_NE(problem.find("must be between 1 and"), std::string::npos) << problem;
}

// a full-screen pass: red grows to the right, green downwards; the read-back rows go top to bottom
TEST(VulkanOffscreen, DrawnPatternHasTheRightOrientation) {
    MDBRIDGE_VULKAN_CONTEXT_OR_SKIP(context);
    std::string problem;
    const std::uint32_t w = 40;
    const std::uint32_t h = 30;
    auto target = OffscreenTarget::create(*context, w, h, {}, &problem);
    ASSERT_TRUE(target) << problem;
    GraphicsPipelineDesc desc;
    desc.stages = {{VK_SHADER_STAGE_VERTEX_BIT, spirv::kFullscreenVertSpirv},
                   {VK_SHADER_STAGE_FRAGMENT_BIT, spirv::kUvGradientFragSpirv}};
    desc.colorFormat = target->colorFormat();
    desc.depthFormat = target->depthFormat();
    desc.pushConstantBytes = 2 * sizeof(float);
    desc.pushConstantStages = VK_SHADER_STAGE_FRAGMENT_BIT;
    auto pipeline = Pipeline::create(*context, desc, &problem);
    ASSERT_TRUE(pipeline) << problem;
    const float pattern[2] = {153.0f / 255, 1.0f};  // blue, alpha
    Image image;
    ASSERT_TRUE(target->render(
        {0.0f, 0.0f, 0.0f, 0.0f},
        [&](VkCommandBuffer cmd) {
            pipeline->bind(cmd);
            pipeline->push(cmd, pattern, sizeof pattern);
            context->vk.vkCmdDraw(cmd, 3, 1, 0, 0);
        },
        &image, &problem))
        << problem;
    int worst = 0;
    for (std::uint32_t y = 0; y < h; ++y) {
        for (std::uint32_t x = 0; x < w; ++x) {
            const std::uint8_t* p = image.pixel(x, y);
            const int red = unorm8((static_cast<float>(x) + 0.5f) / static_cast<float>(w));
            const int green = unorm8((static_cast<float>(y) + 0.5f) / static_cast<float>(h));
            worst = std::max({worst, std::abs(p[0] - red), std::abs(p[1] - green)});
            ASSERT_EQ(p[2], 153) << x << ", " << y;
            ASSERT_EQ(p[3], 255) << x << ", " << y;
        }
    }
    EXPECT_LE(worst, 1);  // interpolation and float -> unorm rounding
    EXPECT_LT(image.pixel(0, 0)[0], 10);  // top-left: dark
    EXPECT_GT(image.pixel(w - 1, h - 1)[1], 245);  // bottom-right: full green
}

// the milestone-1 placeholder shaders (vertex-index triangle, push constants) offscreen, with depth
TEST(VulkanOffscreen, PlaceholderTriangle) {
    MDBRIDGE_VULKAN_CONTEXT_OR_SKIP(context);
    std::string problem;
    auto target = OffscreenTarget::create(*context, 64, 64, {}, &problem);
    ASSERT_TRUE(target) << problem;
    GraphicsPipelineDesc desc;
    desc.stages = {{VK_SHADER_STAGE_VERTEX_BIT, spirv::kPlaceholderVertSpirv},
                   {VK_SHADER_STAGE_FRAGMENT_BIT, spirv::kPlaceholderFragSpirv}};
    desc.colorFormat = target->colorFormat();
    desc.depthFormat = target->depthFormat();
    desc.depthTest = true;
    desc.depthWrite = true;
    desc.pushConstantBytes = 2 * sizeof(float);
    desc.pushConstantStages = VK_SHADER_STAGE_VERTEX_BIT;
    auto pipeline = Pipeline::create(*context, desc, &problem);
    ASSERT_TRUE(pipeline) << problem;
    const float scene[2] = {0.0f, 1.0f};  // angle, aspect
    Image image;
    ASSERT_TRUE(target->render(
        {1.0f, 1.0f, 1.0f, 1.0f},
        [&](VkCommandBuffer cmd) {
            pipeline->bind(cmd);
            pipeline->push(cmd, scene, sizeof scene);
            context->vk.vkCmdDraw(cmd, 3, 1, 0, 0);
        },
        &image, &problem))
        << problem;
    const std::uint8_t* centre = image.pixel(32, 32);
    EXPECT_FALSE(centre[0] == 255 && centre[1] == 255 && centre[2] == 255) << "the triangle covers the centre";
    const std::uint8_t* corner = image.pixel(0, 0);
    EXPECT_TRUE(corner[0] == 255 && corner[1] == 255 && corner[2] == 255) << "the corner keeps the clear colour";
    // the first corner (blue, #2a78d6) points up: near the top centre the blue channel dominates
    const std::uint8_t* top = image.pixel(32, 15);
    EXPECT_GT(top[2], top[0]) << int(top[0]) << " " << int(top[1]) << " " << int(top[2]);
}

TEST(VulkanPipelineCache, SavedAndReloaded) {
    MDBRIDGE_VULKAN_CONTEXT_OR_SKIP(context);
    const std::filesystem::path path = tempPath("pipeline-cache.bin");
    std::error_code ignored;
    std::filesystem::remove(path, ignored);
    ContextOptions options = vulkan_test::testOptions();
    options.pipelineCachePath = path;
    std::string problem;
    {
        auto first = Context::create(options, &problem);
        ASSERT_TRUE(first) << problem;
        EXPECT_FALSE(first->pipelineCacheLoaded());  // no file yet
        GraphicsPipelineDesc desc;
        desc.stages = {{VK_SHADER_STAGE_VERTEX_BIT, spirv::kFullscreenVertSpirv},
                       {VK_SHADER_STAGE_FRAGMENT_BIT, spirv::kUvGradientFragSpirv}};
        desc.pushConstantBytes = 2 * sizeof(float);
        desc.pushConstantStages = VK_SHADER_STAGE_FRAGMENT_BIT;
        ASSERT_TRUE(Pipeline::create(*first, desc, &problem)) << problem;
        ASSERT_TRUE(first->savePipelineCache(&problem)) << problem;
    }
    ASSERT_TRUE(std::filesystem::exists(path));
    EXPECT_FALSE(std::filesystem::exists(path.string() + ".tmp"));
    {
        auto second = Context::create(options, &problem);
        ASSERT_TRUE(second) << problem;
        EXPECT_TRUE(second->pipelineCacheLoaded());
    }
    {  // a damaged file is ignored, not handed to the driver
        std::fstream file(path, std::ios::in | std::ios::out | std::ios::binary);
        file.seekp(-1, std::ios::end);
        file.put('\x5a');
    }
    {
        auto third = Context::create(options, &problem);
        ASSERT_TRUE(third) << problem;
        EXPECT_FALSE(third->pipelineCacheLoaded());
    }
    std::filesystem::remove(path, ignored);
}

TEST(VulkanContext, SubmitsFromSeveralThreads) {
    MDBRIDGE_VULKAN_CONTEXT_OR_SKIP(context);
    std::string problem;
    auto a = OffscreenTarget::create(*context, 32, 32, {}, &problem);
    auto b = OffscreenTarget::create(*context, 32, 32, {}, &problem);
    ASSERT_TRUE(a && b) << problem;
    bool okA = true;
    bool okB = true;
    std::thread other([&] {
        Image image;
        std::string p;
        for (int i = 0; i < 20; ++i) {
            okB = okB && b->render({0.0f, 0.0f, 1.0f, 1.0f}, {}, &image, &p) && image.pixel(5, 5)[2] == 255;
        }
    });
    Image image;
    for (int i = 0; i < 20; ++i) {
        okA = okA && a->render({1.0f, 0.0f, 0.0f, 1.0f}, {}, &image, &problem) && image.pixel(5, 5)[0] == 255;
    }
    other.join();
    EXPECT_TRUE(okA);
    EXPECT_TRUE(okB);
}
