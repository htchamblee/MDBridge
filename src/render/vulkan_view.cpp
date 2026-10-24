#include "vulkan_view.hpp"

#include <QElapsedTimer>
#include <QGuiApplication>
#include <QTimer>
#include <QVulkanDeviceFunctions>
#include <QVulkanFunctions>

#include "mdbridge/render/vulkan_probe.hpp"
#include "scene.hpp"
#include "shaders/spirv.hpp"

namespace mdbridge::render {

namespace {

// Time from the first expose to the check whether the window started (swapchain creation can
// fail late, e.g. over a remote X display).
constexpr int kStartCheckMs = 2000;

struct Scene {  // push constants of shaders/placeholder.vert
    float angle;
    float aspect;
};

class PlaceholderRenderer : public QVulkanWindowRenderer {
public:
    explicit PlaceholderRenderer(VulkanWindow* window) : m_window(window) {}

    void initResources() override;
    void releaseResources() override;
    void startNextFrame() override;

private:
    VkShaderModule shaderModule(const std::uint32_t* code, std::size_t bytes);

    VulkanWindow* m_window;
    QVulkanDeviceFunctions* m_f = nullptr;
    VkPipelineLayout m_layout = VK_NULL_HANDLE;
    VkPipeline m_pipeline = VK_NULL_HANDLE;
    QElapsedTimer m_clock;
};

VkShaderModule PlaceholderRenderer::shaderModule(const std::uint32_t* code, std::size_t bytes) {
    const VkShaderModuleCreateInfo info{
        .sType = VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO,
        .codeSize = bytes,
        .pCode = code,
    };
    VkShaderModule module = VK_NULL_HANDLE;
    if (m_f->vkCreateShaderModule(m_window->device(), &info, nullptr, &module) != VK_SUCCESS) {
        qWarning("MDBridge: could not create a Vulkan shader module");
    }
    return module;
}

void PlaceholderRenderer::initResources() {
    const VkDevice device = m_window->device();
    m_f = m_window->vulkanInstance()->deviceFunctions(device);
    m_clock.start();

    const VkPushConstantRange range{VK_SHADER_STAGE_VERTEX_BIT, 0, sizeof(Scene)};
    const VkPipelineLayoutCreateInfo layoutInfo{
        .sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO,
        .pushConstantRangeCount = 1,
        .pPushConstantRanges = &range,
    };
    if (m_f->vkCreatePipelineLayout(device, &layoutInfo, nullptr, &m_layout) != VK_SUCCESS) {
        qWarning("MDBridge: could not create the Vulkan pipeline layout");
        return;
    }

    const VkShaderModule vertex = shaderModule(spirv::kPlaceholderVertSpirv, sizeof spirv::kPlaceholderVertSpirv);
    const VkShaderModule fragment = shaderModule(spirv::kPlaceholderFragSpirv, sizeof spirv::kPlaceholderFragSpirv);
    const VkPipelineShaderStageCreateInfo stages[] = {
        {.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO,
         .stage = VK_SHADER_STAGE_VERTEX_BIT,
         .module = vertex,
         .pName = "main"},
        {.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO,
         .stage = VK_SHADER_STAGE_FRAGMENT_BIT,
         .module = fragment,
         .pName = "main"},
    };
    const VkPipelineVertexInputStateCreateInfo vertexInput{
        .sType = VK_STRUCTURE_TYPE_PIPELINE_VERTEX_INPUT_STATE_CREATE_INFO};  // corners are built into the shader
    const VkPipelineInputAssemblyStateCreateInfo assembly{
        .sType = VK_STRUCTURE_TYPE_PIPELINE_INPUT_ASSEMBLY_STATE_CREATE_INFO,
        .topology = VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST,
    };
    const VkPipelineViewportStateCreateInfo viewport{
        .sType = VK_STRUCTURE_TYPE_PIPELINE_VIEWPORT_STATE_CREATE_INFO,
        .viewportCount = 1,
        .scissorCount = 1,
    };
    const VkPipelineRasterizationStateCreateInfo rasterization{
        .sType = VK_STRUCTURE_TYPE_PIPELINE_RASTERIZATION_STATE_CREATE_INFO,
        .polygonMode = VK_POLYGON_MODE_FILL,
        .cullMode = VK_CULL_MODE_NONE,
        .frontFace = VK_FRONT_FACE_COUNTER_CLOCKWISE,
        .lineWidth = 1.0f,
    };
    const VkPipelineMultisampleStateCreateInfo multisample{
        .sType = VK_STRUCTURE_TYPE_PIPELINE_MULTISAMPLE_STATE_CREATE_INFO,
        .rasterizationSamples = m_window->sampleCountFlagBits(),
    };
    const VkPipelineDepthStencilStateCreateInfo depthStencil{
        .sType = VK_STRUCTURE_TYPE_PIPELINE_DEPTH_STENCIL_STATE_CREATE_INFO,
        .depthTestEnable = VK_FALSE,
    };
    const VkPipelineColorBlendAttachmentState blendAttachment{
        .colorWriteMask = VK_COLOR_COMPONENT_R_BIT | VK_COLOR_COMPONENT_G_BIT | VK_COLOR_COMPONENT_B_BIT |
                          VK_COLOR_COMPONENT_A_BIT,
    };
    const VkPipelineColorBlendStateCreateInfo blend{
        .sType = VK_STRUCTURE_TYPE_PIPELINE_COLOR_BLEND_STATE_CREATE_INFO,
        .attachmentCount = 1,
        .pAttachments = &blendAttachment,
    };
    const VkDynamicState dynamicStates[] = {VK_DYNAMIC_STATE_VIEWPORT, VK_DYNAMIC_STATE_SCISSOR};
    const VkPipelineDynamicStateCreateInfo dynamic{
        .sType = VK_STRUCTURE_TYPE_PIPELINE_DYNAMIC_STATE_CREATE_INFO,
        .dynamicStateCount = 2,
        .pDynamicStates = dynamicStates,
    };
    const VkGraphicsPipelineCreateInfo info{
        .sType = VK_STRUCTURE_TYPE_GRAPHICS_PIPELINE_CREATE_INFO,
        .stageCount = 2,
        .pStages = stages,
        .pVertexInputState = &vertexInput,
        .pInputAssemblyState = &assembly,
        .pViewportState = &viewport,
        .pRasterizationState = &rasterization,
        .pMultisampleState = &multisample,
        .pDepthStencilState = &depthStencil,
        .pColorBlendState = &blend,
        .pDynamicState = &dynamic,
        .layout = m_layout,
        .renderPass = m_window->defaultRenderPass(),
    };
    if (vertex && fragment &&
        m_f->vkCreateGraphicsPipelines(device, VK_NULL_HANDLE, 1, &info, nullptr, &m_pipeline) != VK_SUCCESS) {
        qWarning("MDBridge: could not create the Vulkan graphics pipeline");
        m_pipeline = VK_NULL_HANDLE;
    }
    m_f->vkDestroyShaderModule(device, vertex, nullptr);
    m_f->vkDestroyShaderModule(device, fragment, nullptr);
}

void PlaceholderRenderer::releaseResources() {
    const VkDevice device = m_window->device();
    m_f->vkDestroyPipeline(device, m_pipeline, nullptr);
    m_f->vkDestroyPipelineLayout(device, m_layout, nullptr);
    m_pipeline = VK_NULL_HANDLE;
    m_layout = VK_NULL_HANDLE;
}

void PlaceholderRenderer::startNextFrame() {
    const QSize size = m_window->swapChainImageSize();
    const auto width = static_cast<std::uint32_t>(size.width());
    const auto height = static_cast<std::uint32_t>(size.height());
    VkClearValue clear[3]{};
    clear[0].color = {{scene::kBackground[0], scene::kBackground[1], scene::kBackground[2], 1.0f}};
    clear[1].depthStencil = {1.0f, 0};
    clear[2].color = clear[0].color;  // multisample resolve target
    const VkRenderPassBeginInfo begin{
        .sType = VK_STRUCTURE_TYPE_RENDER_PASS_BEGIN_INFO,
        .renderPass = m_window->defaultRenderPass(),
        .framebuffer = m_window->currentFramebuffer(),
        .renderArea = {{0, 0}, {width, height}},
        .clearValueCount = m_window->sampleCountFlagBits() > VK_SAMPLE_COUNT_1_BIT ? 3u : 2u,
        .pClearValues = clear,
    };
    const VkCommandBuffer commands = m_window->currentCommandBuffer();
    m_f->vkCmdBeginRenderPass(commands, &begin, VK_SUBPASS_CONTENTS_INLINE);
    if (m_pipeline) {
        const VkViewport viewport{0.0f, 0.0f, static_cast<float>(width), static_cast<float>(height), 0.0f, 1.0f};
        const VkRect2D scissor{{0, 0}, {width, height}};
        m_f->vkCmdSetViewport(commands, 0, 1, &viewport);
        m_f->vkCmdSetScissor(commands, 0, 1, &scissor);
        m_f->vkCmdBindPipeline(commands, VK_PIPELINE_BIND_POINT_GRAPHICS, m_pipeline);
        const Scene push{static_cast<float>(scene::angleAt(m_clock.elapsed())),
                         height ? static_cast<float>(width) / static_cast<float>(height) : 1.0f};
        m_f->vkCmdPushConstants(commands, m_layout, VK_SHADER_STAGE_VERTEX_BIT, 0, static_cast<std::uint32_t>(sizeof push),
                               &push);
        m_f->vkCmdDraw(commands, 3, 1, 0, 0);
    }
    m_f->vkCmdEndRenderPass(commands);
    m_window->frameReady();
    m_window->frameDone();
    m_window->requestUpdate();  // animate: the next frame follows the display (FIFO present)
}

}  // namespace

VulkanWindow::VulkanWindow(std::unique_ptr<QVulkanInstance> instance) : m_instance(std::move(instance)) {
    setVulkanInstance(m_instance.get());
}

// Releases the swapchain and device (QVulkanWindow does on SurfaceAboutToBeDestroyed) while the
// instance still exists.
VulkanWindow::~VulkanWindow() { destroy(); }

QVulkanWindowRenderer* VulkanWindow::createRenderer() { return new PlaceholderRenderer(this); }

void VulkanWindow::exposeEvent(QExposeEvent* event) {
    QVulkanWindow::exposeEvent(event);
    if (isExposed() && !m_checkScheduled) {
        m_checkScheduled = true;
        QTimer::singleShot(kStartCheckMs, this, [this] {
            if (!isValid()) {
                emit failed(tr("the Vulkan window could not start (Qt's messages in the program output say why)"));
            }
        });
    }
}

VulkanView::VulkanView(VulkanWindow* window, const QString& description)
    : View(window, Renderer::Vulkan, description) {
    connect(window, &VulkanWindow::frameRendered, this, [this] { countFrame(); });
    connect(window, &VulkanWindow::failed, this, &View::failed);
}

QImage VulkanView::grabFrame() {
    auto* w = static_cast<VulkanWindow*>(window());
    if (!w || !w->isValid() || !w->supportsGrab()) {
        return {};
    }
    return w->grab();
}

std::unique_ptr<View> createVulkanView(QString* reason) {
    const VulkanReport report = probeVulkan();
    if (!report.usable()) {
        *reason = QString::fromStdString(report.problem);
        return nullptr;
    }
    auto instance = std::make_unique<QVulkanInstance>();
    instance->setApiVersion(QVersionNumber(kVulkanMajor, kVulkanMinor));
    if (qEnvironmentVariableIntValue("MDBRIDGE_VULKAN_VALIDATION")) {
        instance->setLayers({"VK_LAYER_KHRONOS_validation"});
    }
    if (!instance->create()) {
        *reason = instance->errorCode() == VK_SUCCESS  // the platform plugin has no Vulkan support
                      ? QObject::tr("the display (Qt platform \"%1\") does not support Vulkan").arg(QGuiApplication::platformName())
                      : QObject::tr("Qt could not create a Vulkan instance (VkResult %1)").arg(instance->errorCode());
        return nullptr;
    }
    auto window = std::make_unique<VulkanWindow>(std::move(instance));

    // Qt lists the devices in the same order as the probe; choose among them again in case they differ
    const QList<VkPhysicalDeviceProperties> properties = window->availablePhysicalDevices();
    std::vector<VulkanDevice> devices;
    for (const VkPhysicalDeviceProperties& p : properties) {
        devices.push_back(vulkanDevice(p.deviceName, p.deviceType, p.apiVersion, p.vendorID, p.deviceID));
    }
    const int chosen = chooseDevice(devices);
    if (chosen < 0) {
        *reason = QObject::tr("no device supports Vulkan %1.%2").arg(kVulkanMajor).arg(kVulkanMinor);
        return nullptr;
    }
    window->setPhysicalDeviceIndex(chosen);
    // colours are sRGB-encoded like everything QPainter draws, so a UNORM swapchain shows them as
    // the software view does (an _SRGB format would encode them a second time)
    window->setPreferredColorFormats({VK_FORMAT_B8G8R8A8_UNORM, VK_FORMAT_R8G8B8A8_UNORM});
    const VulkanDevice& device = devices[static_cast<std::size_t>(chosen)];
    const QString description = QObject::tr("Vulkan %1 · %2%3")
                                    .arg(QString::fromStdString(vulkanVersionText(device.apiVersion)),
                                         QString::fromStdString(device.name),
                                         device.type == "CPU" ? QObject::tr(" (CPU)") : QString());
    return std::make_unique<VulkanView>(window.release(), description);
}

}  // namespace mdbridge::render
