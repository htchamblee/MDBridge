#pragma once
// The Vulkan 1.3 view: the molecule scene drawn by the GPU scene renderer
// (render/vulkan/scene_renderer.hpp) on Qt's device, straight into the swapchain image (no
// read-back); milestone 1's placeholder without a molecule.

#include <atomic>
#include <memory>
#include <mutex>
#include <utility>

#include <QVulkanInstance>
#include <QVulkanWindow>

#include "mdbridge/render/view.hpp"

namespace mdbridge::render {

// A QVulkanWindow that owns its Vulkan instance, so it can be handed to a window container.
class VulkanWindow : public QVulkanWindow {
    Q_OBJECT

public:
    explicit VulkanWindow(std::unique_ptr<QVulkanInstance> instance);
    ~VulkanWindow() override;

    QVulkanWindowRenderer* createRenderer() override;
    void frameDone() { emit frameRendered(); }

    // Shows `frame` through `camera` from the next frame on (nullptr: the placeholder).
    // Thread-safe.
    void setMoleculeFrame(std::shared_ptr<const MoleculeFrame> frame, const ViewCamera& camera = {});
    std::shared_ptr<const MoleculeFrame> moleculeFrame() const;
    std::pair<std::shared_ptr<const MoleculeFrame>, ViewCamera> moleculeAndCamera() const;
    // Whether the last frame drew the molecule set then (set by the renderer).
    void setMoleculeDrawn(const MoleculeFrame* frame) { m_drawn = frame; }
    bool moleculeDrawn() const;

signals:
    void frameRendered();
    void failed(const QString& reason);

protected:
    void exposeEvent(QExposeEvent* event) override;

private:
    std::unique_ptr<QVulkanInstance> m_instance;
    bool m_checkScheduled = false;
    mutable std::mutex m_frameMutex;
    std::shared_ptr<const MoleculeFrame> m_frame;
    ViewCamera m_camera;
    std::atomic<const MoleculeFrame*> m_drawn{nullptr};
};

class VulkanView : public View {
    Q_OBJECT

public:
    VulkanView(VulkanWindow* window, const QString& description);
    QImage grabFrame() override;
    bool moleculeShown() const override;

protected:
    bool showMolecule(std::shared_ptr<const MoleculeFrame> frame, const ViewCamera& camera) override;
};

// The Vulkan view on the preferred device with Vulkan 1.3, or nullptr with the reason.
std::unique_ptr<View> createVulkanView(QString* reason);

}  // namespace mdbridge::render
