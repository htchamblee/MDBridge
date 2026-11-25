#version 450
#extension GL_GOOGLE_include_directive : require
// Draws the molecule scene renderer's final image into a window (with fullscreen.vert): the
// straight-alpha pixels over the background colour, each image pixel 1 x 1 or, for "fast"
// quality, 2 x 2 window pixels. Push constants (molecule_common.glsl's block): W, H = image size,
// base = 1 for 2 x 2 pixels, fogR / fogG / fogB = background.

#include "molecule_common.glsl"

layout(std430, set = 0, binding = 7) readonly buffer Output { uint outImage[]; };

layout(location = 0) in vec2 uv;
layout(location = 0) out vec4 colour;

void main() {
    ivec2 p = ivec2(gl_FragCoord.xy) >> int(pc.base);
    vec3 bg = vec3(pc.fogR, pc.fogG, pc.fogB);
    if (p.x >= int(pc.W) || p.y >= int(pc.H)) {
        colour = vec4(bg, 1.0);
        return;
    }
    uint v = outImage[p.y * int(pc.W) + p.x];
    vec4 c = vec4(float(v & 255u), float((v >> 8) & 255u), float((v >> 16) & 255u), float(v >> 24)) / 255.0;
    colour = vec4(c.rgb * c.a + bg * (1.0 - c.a), 1.0);
}
