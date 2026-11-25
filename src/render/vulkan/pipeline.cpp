#include "mdbridge/render/vulkan/pipeline.hpp"

namespace mdbridge::render::vulkan {

std::unique_ptr<Pipeline> Pipeline::create(Context& context, const GraphicsPipelineDesc& desc, std::string* problem) {
    const DeviceFunctions& vk = context.vk;
    const VkDevice device = context.device();
    std::unique_ptr<Pipeline> pipeline(new Pipeline(context));
    pipeline->m_pushStages = desc.pushConstantStages;

    VkPushConstantRange range{};
    range.stageFlags = desc.pushConstantStages;
    range.size = desc.pushConstantBytes;
    VkPipelineLayoutCreateInfo layoutInfo{};
    layoutInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
    layoutInfo.setLayoutCount = static_cast<std::uint32_t>(desc.setLayouts.size());
    layoutInfo.pSetLayouts = desc.setLayouts.data();
    layoutInfo.pushConstantRangeCount = desc.pushConstantBytes ? 1 : 0;
    layoutInfo.pPushConstantRanges = &range;
    if (const VkResult r = vk.vkCreatePipelineLayout(device, &layoutInfo, nullptr, &pipeline->m_layout);
        r != VK_SUCCESS) {
        pipeline->m_layout = VK_NULL_HANDLE;
        *problem = "could not create a Vulkan pipeline layout (" + resultText(r) + ")";
        return nullptr;
    }

    std::vector<VkShaderModule> modules;
    std::vector<VkPipelineShaderStageCreateInfo> stages;
    const auto release = [&] {
        for (VkShaderModule module : modules) {
            vk.vkDestroyShaderModule(device, module, nullptr);
        }
    };
    for (const ShaderStage& stage : desc.stages) {
        VkShaderModuleCreateInfo moduleInfo{};
        moduleInfo.sType = VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO;
        moduleInfo.codeSize = stage.spirv.size_bytes();
        moduleInfo.pCode = stage.spirv.data();
        VkShaderModule module = VK_NULL_HANDLE;
        if (const VkResult r = vk.vkCreateShaderModule(device, &moduleInfo, nullptr, &module); r != VK_SUCCESS) {
            release();
            *problem = "could not create a Vulkan shader module (" + resultText(r) + ")";
            return nullptr;
        }
        modules.push_back(module);
        VkPipelineShaderStageCreateInfo info{};
        info.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
        info.stage = stage.stage;
        info.module = module;
        info.pName = "main";
        stages.push_back(info);
    }

    VkPipelineVertexInputStateCreateInfo vertexInput{};
    vertexInput.sType = VK_STRUCTURE_TYPE_PIPELINE_VERTEX_INPUT_STATE_CREATE_INFO;
    vertexInput.vertexBindingDescriptionCount = static_cast<std::uint32_t>(desc.vertexBindings.size());
    vertexInput.pVertexBindingDescriptions = desc.vertexBindings.data();
    vertexInput.vertexAttributeDescriptionCount = static_cast<std::uint32_t>(desc.vertexAttributes.size());
    vertexInput.pVertexAttributeDescriptions = desc.vertexAttributes.data();
    VkPipelineInputAssemblyStateCreateInfo assembly{};
    assembly.sType = VK_STRUCTURE_TYPE_PIPELINE_INPUT_ASSEMBLY_STATE_CREATE_INFO;
    assembly.topology = desc.topology;
    VkPipelineViewportStateCreateInfo viewport{};
    viewport.sType = VK_STRUCTURE_TYPE_PIPELINE_VIEWPORT_STATE_CREATE_INFO;
    viewport.viewportCount = 1;
    viewport.scissorCount = 1;
    VkPipelineRasterizationStateCreateInfo raster{};
    raster.sType = VK_STRUCTURE_TYPE_PIPELINE_RASTERIZATION_STATE_CREATE_INFO;
    raster.polygonMode = VK_POLYGON_MODE_FILL;
    raster.cullMode = desc.cullMode;
    raster.frontFace = desc.frontFace;
    raster.lineWidth = 1.0f;
    VkPipelineMultisampleStateCreateInfo multisample{};
    multisample.sType = VK_STRUCTURE_TYPE_PIPELINE_MULTISAMPLE_STATE_CREATE_INFO;
    multisample.rasterizationSamples = desc.samples;
    VkPipelineDepthStencilStateCreateInfo depth{};
    depth.sType = VK_STRUCTURE_TYPE_PIPELINE_DEPTH_STENCIL_STATE_CREATE_INFO;
    depth.depthTestEnable = desc.depthTest ? VK_TRUE : VK_FALSE;
    depth.depthWriteEnable = desc.depthWrite ? VK_TRUE : VK_FALSE;
    depth.depthCompareOp = desc.depthCompare;
    VkPipelineColorBlendAttachmentState attachment{};
    attachment.colorWriteMask =
        VK_COLOR_COMPONENT_R_BIT | VK_COLOR_COMPONENT_G_BIT | VK_COLOR_COMPONENT_B_BIT | VK_COLOR_COMPONENT_A_BIT;
    if (desc.blend) {
        attachment.blendEnable = VK_TRUE;
        attachment.srcColorBlendFactor = VK_BLEND_FACTOR_ONE;
        attachment.dstColorBlendFactor = VK_BLEND_FACTOR_ONE_MINUS_SRC_ALPHA;
        attachment.colorBlendOp = VK_BLEND_OP_ADD;
        attachment.srcAlphaBlendFactor = VK_BLEND_FACTOR_ONE;
        attachment.dstAlphaBlendFactor = VK_BLEND_FACTOR_ONE_MINUS_SRC_ALPHA;
        attachment.alphaBlendOp = VK_BLEND_OP_ADD;
    }
    VkPipelineColorBlendStateCreateInfo blend{};
    blend.sType = VK_STRUCTURE_TYPE_PIPELINE_COLOR_BLEND_STATE_CREATE_INFO;
    blend.attachmentCount = 1;
    blend.pAttachments = &attachment;
    const VkDynamicState dynamicStates[] = {VK_DYNAMIC_STATE_VIEWPORT, VK_DYNAMIC_STATE_SCISSOR};
    VkPipelineDynamicStateCreateInfo dynamic{};
    dynamic.sType = VK_STRUCTURE_TYPE_PIPELINE_DYNAMIC_STATE_CREATE_INFO;
    dynamic.dynamicStateCount = 2;
    dynamic.pDynamicStates = dynamicStates;
    VkPipelineRenderingCreateInfo rendering{};
    rendering.sType = VK_STRUCTURE_TYPE_PIPELINE_RENDERING_CREATE_INFO;
    rendering.colorAttachmentCount = 1;
    rendering.pColorAttachmentFormats = &desc.colorFormat;
    rendering.depthAttachmentFormat = desc.depthFormat;

    VkGraphicsPipelineCreateInfo info{};
    info.sType = VK_STRUCTURE_TYPE_GRAPHICS_PIPELINE_CREATE_INFO;
    info.pNext = desc.renderPass ? nullptr : &rendering;
    info.renderPass = desc.renderPass;
    info.stageCount = static_cast<std::uint32_t>(stages.size());
    info.pStages = stages.data();
    info.pVertexInputState = &vertexInput;
    info.pInputAssemblyState = &assembly;
    info.pViewportState = &viewport;
    info.pRasterizationState = &raster;
    info.pMultisampleState = &multisample;
    info.pDepthStencilState = &depth;
    info.pColorBlendState = &blend;
    info.pDynamicState = &dynamic;
    info.layout = pipeline->m_layout;
    const VkResult r =
        vk.vkCreateGraphicsPipelines(device, context.pipelineCache(), 1, &info, nullptr, &pipeline->m_pipeline);
    release();
    if (r != VK_SUCCESS) {
        pipeline->m_pipeline = VK_NULL_HANDLE;
        *problem = "could not create a Vulkan graphics pipeline (" + resultText(r) + ")";
        return nullptr;
    }
    return pipeline;
}

Pipeline::~Pipeline() {
    if (m_pipeline) {
        m_context.vk.vkDestroyPipeline(m_context.device(), m_pipeline, nullptr);
    }
    if (m_layout) {
        m_context.vk.vkDestroyPipelineLayout(m_context.device(), m_layout, nullptr);
    }
}

void Pipeline::bind(VkCommandBuffer cmd) const {
    m_context.vk.vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_pipeline);
}

void Pipeline::push(VkCommandBuffer cmd, const void* data, std::uint32_t bytes) const {
    m_context.vk.vkCmdPushConstants(cmd, m_layout, m_pushStages, 0, bytes, data);
}

std::unique_ptr<ComputePipeline> ComputePipeline::create(Context& context, const ComputePipelineDesc& desc,
                                                         std::string* problem) {
    const DeviceFunctions& vk = context.vk;
    const VkDevice device = context.device();
    std::unique_ptr<ComputePipeline> pipeline(new ComputePipeline(context));
    VkPushConstantRange range{};
    range.stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
    range.size = desc.pushConstantBytes;
    VkPipelineLayoutCreateInfo layoutInfo{};
    layoutInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
    layoutInfo.setLayoutCount = static_cast<std::uint32_t>(desc.setLayouts.size());
    layoutInfo.pSetLayouts = desc.setLayouts.data();
    layoutInfo.pushConstantRangeCount = desc.pushConstantBytes ? 1 : 0;
    layoutInfo.pPushConstantRanges = &range;
    if (const VkResult r = vk.vkCreatePipelineLayout(device, &layoutInfo, nullptr, &pipeline->m_layout);
        r != VK_SUCCESS) {
        pipeline->m_layout = VK_NULL_HANDLE;
        *problem = "could not create a Vulkan pipeline layout (" + resultText(r) + ")";
        return nullptr;
    }
    VkShaderModuleCreateInfo moduleInfo{};
    moduleInfo.sType = VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO;
    moduleInfo.codeSize = desc.spirv.size_bytes();
    moduleInfo.pCode = desc.spirv.data();
    VkShaderModule module = VK_NULL_HANDLE;
    if (const VkResult r = vk.vkCreateShaderModule(device, &moduleInfo, nullptr, &module); r != VK_SUCCESS) {
        *problem = "could not create a Vulkan shader module (" + resultText(r) + ")";
        return nullptr;
    }
    VkComputePipelineCreateInfo info{};
    info.sType = VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO;
    info.stage.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
    info.stage.stage = VK_SHADER_STAGE_COMPUTE_BIT;
    info.stage.module = module;
    info.stage.pName = "main";
    info.layout = pipeline->m_layout;
    const VkResult r =
        vk.vkCreateComputePipelines(device, context.pipelineCache(), 1, &info, nullptr, &pipeline->m_pipeline);
    vk.vkDestroyShaderModule(device, module, nullptr);
    if (r != VK_SUCCESS) {
        pipeline->m_pipeline = VK_NULL_HANDLE;
        *problem = "could not create a Vulkan compute pipeline (" + resultText(r) + ")";
        return nullptr;
    }
    return pipeline;
}

ComputePipeline::~ComputePipeline() {
    if (m_pipeline) {
        m_context.vk.vkDestroyPipeline(m_context.device(), m_pipeline, nullptr);
    }
    if (m_layout) {
        m_context.vk.vkDestroyPipelineLayout(m_context.device(), m_layout, nullptr);
    }
}

void ComputePipeline::bind(VkCommandBuffer cmd) const {
    m_context.vk.vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, m_pipeline);
}

void ComputePipeline::push(VkCommandBuffer cmd, const void* data, std::uint32_t bytes) const {
    m_context.vk.vkCmdPushConstants(cmd, m_layout, VK_SHADER_STAGE_COMPUTE_BIT, 0, bytes, data);
}

}  // namespace mdbridge::render::vulkan
