#pragma once
// The Vulkan 1.3 view.

#include <memory>

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

signals:
    void frameRendered();
    void failed(const QString& reason);

protected:
    void exposeEvent(QExposeEvent* event) override;

private:
    std::unique_ptr<QVulkanInstance> m_instance;
    bool m_checkScheduled = false;
};

class VulkanView : public View {
    Q_OBJECT

public:
    VulkanView(VulkanWindow* window, const QString& description);
    QImage grabFrame() override;
};

// The Vulkan view on the preferred device with Vulkan 1.3, or nullptr with the reason.
std::unique_ptr<View> createVulkanView(QString* reason);

}  // namespace mdbridge::render
