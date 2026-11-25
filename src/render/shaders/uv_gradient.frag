#version 450
// Test pattern for the offscreen target (tests/render/test_vulkan_offscreen.cpp): red grows to the
// right, green downwards, blue and alpha from push constants, so a read-back image shows its
// orientation and row order.

layout(push_constant) uniform Pattern {
    float blue;
    float alpha;
} pattern;

layout(location = 0) in vec2 uv;
layout(location = 0) out vec4 fragment;

void main() {
    fragment = vec4(uv.x, uv.y, pattern.blue, pattern.alpha);
}
