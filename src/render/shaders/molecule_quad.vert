#version 450
#extension GL_GOOGLE_include_directive : require
// The molecule scene renderer's geometry passes: one screen-aligned quad per primitive (instance)
// covering the pixel box its record starts with (x0, y0, x1, y1, inclusive); the fragment shaders
// decide coverage per pixel centre as the software renderer does. No vertex buffer: 6 vertices
// per instance.

#include "molecule_common.glsl"

layout(location = 0) flat out uint record;

void main() {
    uint r = pc.base + uint(gl_InstanceIndex) * pc.stride;
    int x0 = primI(r), y0 = primI(r + 1u), x1 = primI(r + 2u), y1 = primI(r + 3u);
    int v = gl_VertexIndex;  // two triangles: (0, 1, 2) and (3, 4, 5)
    bool right = v == 1 || v == 3 || v == 4;
    bool bottom = v == 2 || v == 4 || v == 5;
    vec2 p = vec2(float(right ? x1 + 1 : x0), float(bottom ? y1 + 1 : y0));
    gl_Position = vec4(p / vec2(float(pc.W), float(pc.H)) * 2.0 - 1.0, 0.0, 1.0);
    record = r;
}
