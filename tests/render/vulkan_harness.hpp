#pragma once
// Test harness for the Vulkan renderer: one shared offscreen context per test program.
//
//   TEST(MyVulkanTest, Draws) {
//       MDBRIDGE_VULKAN_CONTEXT_OR_SKIP(context);  // vulkan::Context* context
//       ...
//   }
//
// Without a usable Vulkan 1.3 driver the tests are skipped with the reason (a machine without a
// GPU or lavapipe), unless MDBRIDGE_REQUIRE_VULKAN=1 is set: CI sets it where lavapipe is
// installed, so a broken driver setup fails instead of passing silently. MDBRIDGE_VULKAN_DEVICE
// picks a device by name (e.g. "llvmpipe"); MDBRIDGE_VULKAN_VALIDATION=1 adds the Khronos
// validation layer when installed.

#include <gtest/gtest.h>

#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>

#include "mdbridge/render/vulkan/context.hpp"

namespace mdbridge::render::vulkan::vulkan_test {

inline bool requireVulkan() {
    const char* value = std::getenv("MDBRIDGE_REQUIRE_VULKAN");
    return value && *value && std::strcmp(value, "0") != 0;
}

inline ContextOptions testOptions() {
    ContextOptions options;
    options.applicationName = "MDBridge tests";
    if (const char* device = std::getenv("MDBRIDGE_VULKAN_DEVICE")) {
        options.deviceName = device;
    }
    return options;
}

struct SharedContext {
    std::unique_ptr<Context> context;
    std::string problem;
};

// created on first use, destroyed when the tests end (ContextEnvironment)
inline SharedContext& sharedContext() {
    static SharedContext shared = [] {
        SharedContext s;
        s.context = Context::create(testOptions(), &s.problem);
        return s;
    }();
    return shared;
}

// Destroys the shared context after the last test, inside main(): from a static destructor the
// driver or validation layer may already be torn down (lavapipe with the validation layer hangs).
class ContextEnvironment : public ::testing::Environment {
public:
    void TearDown() override { sharedContext().context.reset(); }
};

inline const bool kContextEnvironmentRegistered =
    (::testing::AddGlobalTestEnvironment(new ContextEnvironment), true);  // gtest owns it

}  // namespace mdbridge::render::vulkan::vulkan_test

// Declares `Context* name` = the shared context, or ends the test: skipped (failed with
// MDBRIDGE_REQUIRE_VULKAN=1) when Vulkan is unusable.
#define MDBRIDGE_VULKAN_CONTEXT_OR_SKIP(name)                                                              \
    auto& name##_shared_ = ::mdbridge::render::vulkan::vulkan_test::sharedContext();                          \
    if (!name##_shared_.context) {                                                                         \
        if (::mdbridge::render::vulkan::vulkan_test::requireVulkan()) {                                       \
            FAIL() << "Vulkan is required (MDBRIDGE_REQUIRE_VULKAN) but unusable: " << name##_shared_.problem; \
        }                                                                                                  \
        GTEST_SKIP() << "no usable Vulkan: " << name##_shared_.problem;                                   \
    }                                                                                                      \
    [[maybe_unused]] ::mdbridge::render::vulkan::Context* const name = name##_shared_.context.get()
