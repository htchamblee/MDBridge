#version 450
// One triangle that covers the whole render target, for full-screen passes (offscreen tests,
// later the outline, glow and anti-aliasing passes). No vertex buffer: draw 3 vertices.
// uv: (0, 0) at the top-left corner of the target, (1, 1) at the bottom-right.

layout(location = 0) out vec2 uv;

void main() {
    vec2 corner = vec2(float((gl_VertexIndex << 1) & 2), float(gl_VertexIndex & 2));  // (0,0) (2,0) (0,2)
    uv = corner;
    gl_Position = vec4(corner * 2.0 - 1.0, 0.0, 1.0);
}
