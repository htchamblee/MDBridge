#pragma once
// Graphics pipelines for dynamic rendering, created through the context's pipeline cache.
// Viewport and scissor are dynamic state (OffscreenTarget::render sets them to the whole target).

#include <cstdint>
#include <memory>
#include <span>
#include <string>
#include <vector>

#include "mdbridge/render/vulkan/context.hpp"

namespace mdbridge::render::vulkan {

struct ShaderStage {
    VkShaderStageFlagBits stage = VK_SHADER_STAGE_VERTEX_BIT;
    std::span<const std::uint32_t> spirv;  // e.g. spirv::kFullscreenVertSpirv (src/render/shaders/spirv.hpp)
};

struct GraphicsPipelineDesc {
    std::vector<ShaderStage> stages;
    VkFormat colorFormat = VK_FORMAT_R8G8B8A8_UNORM;
    VkFormat depthFormat = VK_FORMAT_UNDEFINED;  // the target's depth format, or none
    VkSampleCountFlagBits samples = VK_SAMPLE_COUNT_1_BIT;
    VkPrimitiveTopology topology = VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST;
    VkCullModeFlags cullMode = VK_CULL_MODE_NONE;
    VkFrontFace frontFace = VK_FRONT_FACE_COUNTER_CLOCKWISE;
    bool depthTest = false;
    bool depthWrite = false;
    VkCompareOp depthCompare = VK_COMPARE_OP_LESS_OR_EQUAL;
    bool blend = false;  // premultiplied "over": src + (1 - src.a) * dst
    std::uint32_t pushConstantBytes = 0;
    VkShaderStageFlags pushConstantStages = VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT;
    std::vector<VkVertexInputBindingDescription> vertexBindings;
    std::vector<VkVertexInputAttributeDescription> vertexAttributes;
    std::vector<VkDescriptorSetLayout> setLayouts;  // not owned
    // a render pass object (subpass 0) instead of dynamic rendering, e.g. a QVulkanWindow's
    // default render pass; colorFormat / depthFormat are then ignored
    VkRenderPass renderPass = VK_NULL_HANDLE;
};

class Pipeline {
public:
    // nullptr with `problem` set if a shader module or the pipeline cannot be created
    static std::unique_ptr<Pipeline> create(Context& context, const GraphicsPipelineDesc& desc, std::string* problem);
    ~Pipeline();
    Pipeline(const Pipeline&) = delete;
    Pipeline& operator=(const Pipeline&) = delete;

    VkPipeline handle() const { return m_pipeline; }
    VkPipelineLayout layout() const { return m_layout; }

    // vkCmdBindPipeline for graphics
    void bind(VkCommandBuffer cmd) const;
    // vkCmdPushConstants of `bytes` at offset 0 with the pipeline's push-constant stages
    void push(VkCommandBuffer cmd, const void* data, std::uint32_t bytes) const;

private:
    explicit Pipeline(Context& context) : m_context(context) {}

    Context& m_context;
    VkPipelineLayout m_layout = VK_NULL_HANDLE;
    VkPipeline m_pipeline = VK_NULL_HANDLE;
    VkShaderStageFlags m_pushStages = 0;
};

struct ComputePipelineDesc {
    std::span<const std::uint32_t> spirv;  // a compute shader
    std::uint32_t pushConstantBytes = 0;
    std::vector<VkDescriptorSetLayout> setLayouts;  // not owned
};

// A compute pipeline through the context's pipeline cache.
class ComputePipeline {
public:
    // nullptr with `problem` set if the shader module or the pipeline cannot be created
    static std::unique_ptr<ComputePipeline> create(Context& context, const ComputePipelineDesc& desc,
                                                   std::string* problem);
    ~ComputePipeline();
    ComputePipeline(const ComputePipeline&) = delete;
    ComputePipeline& operator=(const ComputePipeline&) = delete;

    VkPipeline handle() const { return m_pipeline; }
    VkPipelineLayout layout() const { return m_layout; }
    // vkCmdBindPipeline for compute
    void bind(VkCommandBuffer cmd) const;
    // vkCmdPushConstants of `bytes` at offset 0 (compute stage)
    void push(VkCommandBuffer cmd, const void* data, std::uint32_t bytes) const;

private:
    explicit ComputePipeline(Context& context) : m_context(context) {}

    Context& m_context;
    VkPipelineLayout m_layout = VK_NULL_HANDLE;
    VkPipeline m_pipeline = VK_NULL_HANDLE;
};

}  // namespace mdbridge::render::vulkan
