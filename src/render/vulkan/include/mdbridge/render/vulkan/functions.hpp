#pragma once
// The Vulkan functions MDBridge calls, loaded from the loader at run time (no link-time Vulkan).
// Add a function to the list of its level; the tables are filled by Context::create().

#include "mdbridge/render/vulkan/loader.hpp"

// functions of the loader itself (vkGetInstanceProcAddr(nullptr, ...))
#define MDBRIDGE_VULKAN_GLOBAL_FUNCTIONS(X) \
    X(vkCreateInstance)                     \
    X(vkEnumerateInstanceLayerProperties)   \
    X(vkEnumerateInstanceExtensionProperties)

// instance functions (vkGetInstanceProcAddr(instance, ...))
#define MDBRIDGE_VULKAN_INSTANCE_FUNCTIONS(X)  \
    X(vkDestroyInstance)                       \
    X(vkEnumeratePhysicalDevices)              \
    X(vkGetPhysicalDeviceProperties)           \
    X(vkGetPhysicalDeviceFeatures2)            \
    X(vkGetPhysicalDeviceQueueFamilyProperties) \
    X(vkGetPhysicalDeviceMemoryProperties)     \
    X(vkGetPhysicalDeviceFormatProperties)     \
    X(vkCreateDevice)                          \
    X(vkGetDeviceProcAddr)

// device functions (vkGetDeviceProcAddr(device, ...)); Vulkan 1.3 core only
#define MDBRIDGE_VULKAN_DEVICE_FUNCTIONS(X) \
    X(vkDestroyDevice)                      \
    X(vkGetDeviceQueue)                     \
    X(vkDeviceWaitIdle)                     \
    X(vkQueueSubmit2)                       \
    X(vkQueueWaitIdle)                      \
    X(vkCreateCommandPool)                  \
    X(vkDestroyCommandPool)                 \
    X(vkResetCommandPool)                   \
    X(vkAllocateCommandBuffers)             \
    X(vkFreeCommandBuffers)                 \
    X(vkBeginCommandBuffer)                 \
    X(vkEndCommandBuffer)                   \
    X(vkCreateFence)                        \
    X(vkDestroyFence)                       \
    X(vkResetFences)                        \
    X(vkWaitForFences)                      \
    X(vkCreateImage)                        \
    X(vkDestroyImage)                       \
    X(vkCreateImageView)                    \
    X(vkDestroyImageView)                   \
    X(vkCreateBuffer)                       \
    X(vkDestroyBuffer)                      \
    X(vkCreateSampler)                      \
    X(vkDestroySampler)                     \
    X(vkGetImageMemoryRequirements)         \
    X(vkGetBufferMemoryRequirements)        \
    X(vkAllocateMemory)                     \
    X(vkFreeMemory)                         \
    X(vkBindImageMemory)                    \
    X(vkBindBufferMemory)                   \
    X(vkMapMemory)                          \
    X(vkUnmapMemory)                        \
    X(vkFlushMappedMemoryRanges)            \
    X(vkInvalidateMappedMemoryRanges)       \
    X(vkCreateShaderModule)                 \
    X(vkDestroyShaderModule)                \
    X(vkCreatePipelineCache)                \
    X(vkDestroyPipelineCache)               \
    X(vkGetPipelineCacheData)               \
    X(vkCreatePipelineLayout)               \
    X(vkDestroyPipelineLayout)              \
    X(vkCreateGraphicsPipelines)            \
    X(vkCreateComputePipelines)             \
    X(vkDestroyPipeline)                    \
    X(vkCreateDescriptorSetLayout)          \
    X(vkDestroyDescriptorSetLayout)         \
    X(vkCreateDescriptorPool)               \
    X(vkDestroyDescriptorPool)              \
    X(vkResetDescriptorPool)                \
    X(vkAllocateDescriptorSets)             \
    X(vkUpdateDescriptorSets)               \
    X(vkCmdPipelineBarrier2)                \
    X(vkCmdBeginRendering)                  \
    X(vkCmdEndRendering)                    \
    X(vkCmdBindPipeline)                    \
    X(vkCmdBindDescriptorSets)              \
    X(vkCmdBindVertexBuffers)               \
    X(vkCmdBindIndexBuffer)                 \
    X(vkCmdSetViewport)                     \
    X(vkCmdSetScissor)                      \
    X(vkCmdPushConstants)                   \
    X(vkCmdDraw)                            \
    X(vkCmdDrawIndexed)                     \
    X(vkCmdDispatch)                        \
    X(vkCmdCopyBuffer)                      \
    X(vkCmdCopyBufferToImage)               \
    X(vkCmdCopyImageToBuffer)               \
    X(vkCmdBlitImage)

namespace mdbridge::render::vulkan {

#define MDBRIDGE_VULKAN_MEMBER(name) PFN_##name name = nullptr;

struct GlobalFunctions {
    MDBRIDGE_VULKAN_GLOBAL_FUNCTIONS(MDBRIDGE_VULKAN_MEMBER)
};

struct InstanceFunctions {
    MDBRIDGE_VULKAN_INSTANCE_FUNCTIONS(MDBRIDGE_VULKAN_MEMBER)
};

struct DeviceFunctions {
    MDBRIDGE_VULKAN_DEVICE_FUNCTIONS(MDBRIDGE_VULKAN_MEMBER)
};

#undef MDBRIDGE_VULKAN_MEMBER

}  // namespace mdbridge::render::vulkan
