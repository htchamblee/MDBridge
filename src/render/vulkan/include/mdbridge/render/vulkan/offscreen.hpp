#pragma once
// An offscreen render target (colour + optional depth image, no surface or swapchain) and the
// read-back of what was drawn into it: works without a display, over ssh and in CI (lavapipe).

#include <array>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "mdbridge/render/vulkan/context.hpp"

namespace mdbridge::render::vulkan {

// A read-back image: RGBA, 8 bits per channel, rows top to bottom, no padding (the layout of a
// QImage::Format_RGBA8888 and of the software renderer's frames).
struct Image {
    std::uint32_t width = 0;
    std::uint32_t height = 0;
    std::vector<std::uint8_t> rgba;  // width * height * 4 bytes

    const std::uint8_t* pixel(std::uint32_t x, std::uint32_t y) const {
        return rgba.data() + (static_cast<std::size_t>(y) * width + x) * 4;
    }
};

struct OffscreenOptions {
    VkFormat colorFormat = VK_FORMAT_R8G8B8A8_UNORM;  // read back as is: must be 4 bytes RGBA
    bool depth = true;  // a depth attachment (the first of D32_SFLOAT, X8_D24, D16 the device supports)
};

class OffscreenTarget {
public:
    // nullptr with `problem` set if the size is 0 or above the device's limit, the format is not
    // renderable, or memory runs out
    static std::unique_ptr<OffscreenTarget> create(Context& context, std::uint32_t width, std::uint32_t height,
                                                   const OffscreenOptions& options, std::string* problem);
    ~OffscreenTarget();
    OffscreenTarget(const OffscreenTarget&) = delete;
    OffscreenTarget& operator=(const OffscreenTarget&) = delete;

    // Clears the target to `clear` (RGBA, 0-1; depth to 1), lets `draw` record its commands inside
    // the rendering (viewport and scissor already set to the whole target; may be empty), then
    // copies the colour image into `out`. Waits for the GPU. False with `problem` set on failure.
    bool render(const std::array<float, 4>& clear, const std::function<void(VkCommandBuffer)>& draw, Image* out,
                std::string* problem);

    std::uint32_t width() const { return m_width; }
    std::uint32_t height() const { return m_height; }
    VkFormat colorFormat() const { return m_colorFormat; }
    VkFormat depthFormat() const { return m_depthFormat; }  // VK_FORMAT_UNDEFINED without depth
    Context& context() const { return m_context; }

private:
    explicit OffscreenTarget(Context& context) : m_context(context) {}
    std::string init(std::uint32_t width, std::uint32_t height, const OffscreenOptions& options);
    std::string allocateImage(VkFormat format, VkImageUsageFlags usage, VkImageAspectFlags aspect, VkImage* image,
                              VkDeviceMemory* memory, VkImageView* view);

    Context& m_context;
    std::uint32_t m_width = 0;
    std::uint32_t m_height = 0;
    VkFormat m_colorFormat = VK_FORMAT_UNDEFINED;
    VkFormat m_depthFormat = VK_FORMAT_UNDEFINED;
    VkImage m_color = VK_NULL_HANDLE;
    VkDeviceMemory m_colorMemory = VK_NULL_HANDLE;
    VkImageView m_colorView = VK_NULL_HANDLE;
    VkImage m_depth = VK_NULL_HANDLE;
    VkDeviceMemory m_depthMemory = VK_NULL_HANDLE;
    VkImageView m_depthView = VK_NULL_HANDLE;
    VkBuffer m_readback = VK_NULL_HANDLE;
    VkDeviceMemory m_readbackMemory = VK_NULL_HANDLE;
    bool m_readbackCoherent = false;
    void* m_mapped = nullptr;
};

}  // namespace mdbridge::render::vulkan
