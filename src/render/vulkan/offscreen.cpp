#include "mdbridge/render/vulkan/offscreen.hpp"

#include <cstring>

namespace mdbridge::render::vulkan {

namespace {

VkImageMemoryBarrier2 imageBarrier(VkImage image, VkImageAspectFlags aspect, VkImageLayout from, VkImageLayout to,
                                   VkPipelineStageFlags2 srcStage, VkAccessFlags2 srcAccess,
                                   VkPipelineStageFlags2 dstStage, VkAccessFlags2 dstAccess) {
    VkImageMemoryBarrier2 b{};
    b.sType = VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER_2;
    b.srcStageMask = srcStage;
    b.srcAccessMask = srcAccess;
    b.dstStageMask = dstStage;
    b.dstAccessMask = dstAccess;
    b.oldLayout = from;
    b.newLayout = to;
    b.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    b.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    b.image = image;
    b.subresourceRange = {aspect, 0, 1, 0, 1};
    return b;
}

}  // namespace

std::unique_ptr<OffscreenTarget> OffscreenTarget::create(Context& context, std::uint32_t width, std::uint32_t height,
                                                         const OffscreenOptions& options, std::string* problem) {
    std::unique_ptr<OffscreenTarget> target(new OffscreenTarget(context));
    *problem = target->init(width, height, options);
    if (!problem->empty()) {
        return nullptr;
    }
    return target;
}

std::string OffscreenTarget::allocateImage(VkFormat format, VkImageUsageFlags usage, VkImageAspectFlags aspect,
                                           VkImage* image, VkDeviceMemory* memory, VkImageView* view) {
    const DeviceFunctions& vk = m_context.vk;
    const VkDevice device = m_context.device();
    VkImageCreateInfo info{};
    info.sType = VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO;
    info.imageType = VK_IMAGE_TYPE_2D;
    info.format = format;
    info.extent = {m_width, m_height, 1};
    info.mipLevels = 1;
    info.arrayLayers = 1;
    info.samples = VK_SAMPLE_COUNT_1_BIT;
    info.tiling = VK_IMAGE_TILING_OPTIMAL;
    info.usage = usage;
    info.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
    info.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
    if (const VkResult r = vk.vkCreateImage(device, &info, nullptr, image); r != VK_SUCCESS) {
        *image = VK_NULL_HANDLE;
        return "could not create a " + std::to_string(m_width) + " x " + std::to_string(m_height) +
               " Vulkan image (" + resultText(r) + ")";
    }
    VkMemoryRequirements needs{};
    vk.vkGetImageMemoryRequirements(device, *image, &needs);
    const std::optional<std::uint32_t> type = m_context.memoryType(needs.memoryTypeBits, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
    const std::optional<std::uint32_t> any = type ? type : m_context.memoryType(needs.memoryTypeBits, 0);
    if (!any) {
        return "no Vulkan memory type for the render target";
    }
    VkMemoryAllocateInfo allocate{};
    allocate.sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO;
    allocate.allocationSize = needs.size;
    allocate.memoryTypeIndex = *any;
    if (const VkResult r = vk.vkAllocateMemory(device, &allocate, nullptr, memory); r != VK_SUCCESS) {
        *memory = VK_NULL_HANDLE;
        return "could not allocate " + std::to_string(needs.size) + " bytes for the render target (" + resultText(r) +
               ")";
    }
    vk.vkBindImageMemory(device, *image, *memory, 0);
    VkImageViewCreateInfo viewInfo{};
    viewInfo.sType = VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO;
    viewInfo.image = *image;
    viewInfo.viewType = VK_IMAGE_VIEW_TYPE_2D;
    viewInfo.format = format;
    viewInfo.subresourceRange = {aspect, 0, 1, 0, 1};
    if (const VkResult r = vk.vkCreateImageView(device, &viewInfo, nullptr, view); r != VK_SUCCESS) {
        *view = VK_NULL_HANDLE;
        return "could not create a Vulkan image view (" + resultText(r) + ")";
    }
    return {};
}

std::string OffscreenTarget::init(std::uint32_t width, std::uint32_t height, const OffscreenOptions& options) {
    const DeviceFunctions& vk = m_context.vk;
    const VkDevice device = m_context.device();
    const std::uint32_t limit = m_context.properties().limits.maxImageDimension2D;
    if (width == 0 || height == 0 || width > limit || height > limit) {
        return "render target " + std::to_string(width) + " x " + std::to_string(height) +
               " px: each side must be between 1 and " + std::to_string(limit) + " px on " + m_context.deviceInfo().name;
    }
    m_width = width;
    m_height = height;
    m_colorFormat = options.colorFormat;
    if (!m_context.supportsFormat(m_colorFormat, VK_FORMAT_FEATURE_COLOR_ATTACHMENT_BIT | VK_FORMAT_FEATURE_TRANSFER_SRC_BIT)) {
        return "the Vulkan device cannot render into format " + std::to_string(static_cast<int>(m_colorFormat));
    }
    if (std::string p = allocateImage(m_colorFormat,
                                      VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT | VK_IMAGE_USAGE_TRANSFER_SRC_BIT |
                                          VK_IMAGE_USAGE_SAMPLED_BIT,
                                      VK_IMAGE_ASPECT_COLOR_BIT, &m_color, &m_colorMemory, &m_colorView);
        !p.empty()) {
        return p;
    }
    if (options.depth) {
        for (VkFormat format : {VK_FORMAT_D32_SFLOAT, VK_FORMAT_X8_D24_UNORM_PACK32, VK_FORMAT_D16_UNORM}) {
            if (m_context.supportsFormat(format, VK_FORMAT_FEATURE_DEPTH_STENCIL_ATTACHMENT_BIT)) {
                m_depthFormat = format;
                break;
            }
        }
        if (m_depthFormat == VK_FORMAT_UNDEFINED) {
            return "the Vulkan device has no depth format";
        }
        if (std::string p = allocateImage(m_depthFormat, VK_IMAGE_USAGE_DEPTH_STENCIL_ATTACHMENT_BIT,
                                          VK_IMAGE_ASPECT_DEPTH_BIT, &m_depth, &m_depthMemory, &m_depthView);
            !p.empty()) {
            return p;
        }
    }

    // host-visible buffer the colour image is copied into, mapped for the target's lifetime
    const VkDeviceSize bytes = VkDeviceSize{m_width} * m_height * 4;
    VkBufferCreateInfo bufferInfo{};
    bufferInfo.sType = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO;
    bufferInfo.size = bytes;
    bufferInfo.usage = VK_BUFFER_USAGE_TRANSFER_DST_BIT;
    bufferInfo.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
    if (const VkResult r = vk.vkCreateBuffer(device, &bufferInfo, nullptr, &m_readback); r != VK_SUCCESS) {
        m_readback = VK_NULL_HANDLE;
        return "could not create the Vulkan read-back buffer (" + resultText(r) + ")";
    }
    VkMemoryRequirements needs{};
    vk.vkGetBufferMemoryRequirements(device, m_readback, &needs);
    const std::optional<std::uint32_t> type =
        m_context.memoryType(needs.memoryTypeBits, VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT,
                             VK_MEMORY_PROPERTY_HOST_CACHED_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT);
    if (!type) {
        return "no host-visible Vulkan memory for the read-back buffer";
    }
    VkPhysicalDeviceMemoryProperties memory{};
    m_context.vki.vkGetPhysicalDeviceMemoryProperties(m_context.physicalDevice(), &memory);
    m_readbackCoherent = (memory.memoryTypes[*type].propertyFlags & VK_MEMORY_PROPERTY_HOST_COHERENT_BIT) != 0;
    VkMemoryAllocateInfo allocate{};
    allocate.sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO;
    allocate.allocationSize = needs.size;
    allocate.memoryTypeIndex = *type;
    if (const VkResult r = vk.vkAllocateMemory(device, &allocate, nullptr, &m_readbackMemory); r != VK_SUCCESS) {
        m_readbackMemory = VK_NULL_HANDLE;
        return "could not allocate the Vulkan read-back buffer (" + resultText(r) + ")";
    }
    vk.vkBindBufferMemory(device, m_readback, m_readbackMemory, 0);
    if (const VkResult r = vk.vkMapMemory(device, m_readbackMemory, 0, VK_WHOLE_SIZE, 0, &m_mapped); r != VK_SUCCESS) {
        m_mapped = nullptr;
        return "could not map the Vulkan read-back buffer (" + resultText(r) + ")";
    }
    return {};
}

OffscreenTarget::~OffscreenTarget() {
    const DeviceFunctions& vk = m_context.vk;
    const VkDevice device = m_context.device();
    if (m_mapped) {
        vk.vkUnmapMemory(device, m_readbackMemory);
    }
    if (m_readback) {
        vk.vkDestroyBuffer(device, m_readback, nullptr);
    }
    if (m_readbackMemory) {
        vk.vkFreeMemory(device, m_readbackMemory, nullptr);
    }
    for (VkImageView view : {m_colorView, m_depthView}) {
        if (view) {
            vk.vkDestroyImageView(device, view, nullptr);
        }
    }
    for (VkImage image : {m_color, m_depth}) {
        if (image) {
            vk.vkDestroyImage(device, image, nullptr);
        }
    }
    for (VkDeviceMemory memory : {m_colorMemory, m_depthMemory}) {
        if (memory) {
            vk.vkFreeMemory(device, memory, nullptr);
        }
    }
}

bool OffscreenTarget::render(const std::array<float, 4>& clear, const std::function<void(VkCommandBuffer)>& draw,
                             Image* out, std::string* problem) {
    const DeviceFunctions& vk = m_context.vk;
    const VkResult r = m_context.submitAndWait([&](VkCommandBuffer cmd) {
        // the previous read-back (copy) must finish before the image is drawn into again
        VkImageMemoryBarrier2 before[2] = {
            imageBarrier(m_color, VK_IMAGE_ASPECT_COLOR_BIT, VK_IMAGE_LAYOUT_UNDEFINED,
                         VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL, VK_PIPELINE_STAGE_2_COPY_BIT, 0,
                         VK_PIPELINE_STAGE_2_COLOR_ATTACHMENT_OUTPUT_BIT,
                         VK_ACCESS_2_COLOR_ATTACHMENT_WRITE_BIT | VK_ACCESS_2_COLOR_ATTACHMENT_READ_BIT),
            imageBarrier(m_depth, VK_IMAGE_ASPECT_DEPTH_BIT, VK_IMAGE_LAYOUT_UNDEFINED,
                         VK_IMAGE_LAYOUT_DEPTH_ATTACHMENT_OPTIMAL,
                         VK_PIPELINE_STAGE_2_EARLY_FRAGMENT_TESTS_BIT | VK_PIPELINE_STAGE_2_LATE_FRAGMENT_TESTS_BIT,
                         VK_ACCESS_2_DEPTH_STENCIL_ATTACHMENT_WRITE_BIT,
                         VK_PIPELINE_STAGE_2_EARLY_FRAGMENT_TESTS_BIT | VK_PIPELINE_STAGE_2_LATE_FRAGMENT_TESTS_BIT,
                         VK_ACCESS_2_DEPTH_STENCIL_ATTACHMENT_READ_BIT |
                             VK_ACCESS_2_DEPTH_STENCIL_ATTACHMENT_WRITE_BIT),
        };
        VkDependencyInfo dependency{};
        dependency.sType = VK_STRUCTURE_TYPE_DEPENDENCY_INFO;
        dependency.imageMemoryBarrierCount = m_depth ? 2 : 1;
        dependency.pImageMemoryBarriers = before;
        vk.vkCmdPipelineBarrier2(cmd, &dependency);

        VkRenderingAttachmentInfo color{};
        color.sType = VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO;
        color.imageView = m_colorView;
        color.imageLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;
        color.loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
        color.storeOp = VK_ATTACHMENT_STORE_OP_STORE;
        std::memcpy(color.clearValue.color.float32, clear.data(), sizeof(float) * 4);
        VkRenderingAttachmentInfo depth{};
        depth.sType = VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO;
        depth.imageView = m_depthView;
        depth.imageLayout = VK_IMAGE_LAYOUT_DEPTH_ATTACHMENT_OPTIMAL;
        depth.loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
        depth.storeOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
        depth.clearValue.depthStencil = {1.0f, 0};
        VkRenderingInfo rendering{};
        rendering.sType = VK_STRUCTURE_TYPE_RENDERING_INFO;
        rendering.renderArea = {{0, 0}, {m_width, m_height}};
        rendering.layerCount = 1;
        rendering.colorAttachmentCount = 1;
        rendering.pColorAttachments = &color;
        rendering.pDepthAttachment = m_depth ? &depth : nullptr;
        vk.vkCmdBeginRendering(cmd, &rendering);
        const VkViewport viewport{0.0f, 0.0f, static_cast<float>(m_width), static_cast<float>(m_height), 0.0f, 1.0f};
        const VkRect2D scissor{{0, 0}, {m_width, m_height}};
        vk.vkCmdSetViewport(cmd, 0, 1, &viewport);
        vk.vkCmdSetScissor(cmd, 0, 1, &scissor);
        if (draw) {
            draw(cmd);
        }
        vk.vkCmdEndRendering(cmd);

        VkImageMemoryBarrier2 toCopy = imageBarrier(
            m_color, VK_IMAGE_ASPECT_COLOR_BIT, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL,
            VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL, VK_PIPELINE_STAGE_2_COLOR_ATTACHMENT_OUTPUT_BIT,
            VK_ACCESS_2_COLOR_ATTACHMENT_WRITE_BIT, VK_PIPELINE_STAGE_2_COPY_BIT, VK_ACCESS_2_TRANSFER_READ_BIT);
        dependency.imageMemoryBarrierCount = 1;
        dependency.pImageMemoryBarriers = &toCopy;
        vk.vkCmdPipelineBarrier2(cmd, &dependency);
        VkBufferImageCopy region{};
        region.imageSubresource = {VK_IMAGE_ASPECT_COLOR_BIT, 0, 0, 1};
        region.imageExtent = {m_width, m_height, 1};
        vk.vkCmdCopyImageToBuffer(cmd, m_color, VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL, m_readback, 1, &region);

        VkBufferMemoryBarrier2 toHost{};
        toHost.sType = VK_STRUCTURE_TYPE_BUFFER_MEMORY_BARRIER_2;
        toHost.srcStageMask = VK_PIPELINE_STAGE_2_COPY_BIT;
        toHost.srcAccessMask = VK_ACCESS_2_TRANSFER_WRITE_BIT;
        toHost.dstStageMask = VK_PIPELINE_STAGE_2_HOST_BIT;
        toHost.dstAccessMask = VK_ACCESS_2_HOST_READ_BIT;
        toHost.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
        toHost.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
        toHost.buffer = m_readback;
        toHost.size = VK_WHOLE_SIZE;
        VkDependencyInfo hostDependency{};
        hostDependency.sType = VK_STRUCTURE_TYPE_DEPENDENCY_INFO;
        hostDependency.bufferMemoryBarrierCount = 1;
        hostDependency.pBufferMemoryBarriers = &toHost;
        vk.vkCmdPipelineBarrier2(cmd, &hostDependency);
    });
    if (r != VK_SUCCESS) {
        *problem = "Vulkan rendering failed (" + resultText(r) + ")";
        return false;
    }
    if (!m_readbackCoherent) {
        VkMappedMemoryRange range{};
        range.sType = VK_STRUCTURE_TYPE_MAPPED_MEMORY_RANGE;
        range.memory = m_readbackMemory;
        range.size = VK_WHOLE_SIZE;
        vk.vkInvalidateMappedMemoryRanges(m_context.device(), 1, &range);
    }
    out->width = m_width;
    out->height = m_height;
    const std::size_t bytes = static_cast<std::size_t>(m_width) * m_height * 4;
    out->rgba.resize(bytes);
    std::memcpy(out->rgba.data(), m_mapped, bytes);
    return true;
}

}  // namespace mdbridge::render::vulkan
