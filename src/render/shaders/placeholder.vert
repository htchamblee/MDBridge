#version 450
// Placeholder scene of the 3D view (milestone 1): a turning triangle in the residue-group colours.
// The vertices are built in, so the pipeline needs no vertex buffer.

layout(push_constant) uniform Scene {
    float angle;   // rotation in radians
    float aspect;  // width / height of the view
} scene;

layout(location = 0) out vec3 colour;

const vec2 kCorners[3] = vec2[](vec2(0.0, -0.62), vec2(0.54, 0.31), vec2(-0.54, 0.31));
const vec3 kColours[3] = vec3[](vec3(0.165, 0.471, 0.839),   // #2a78d6, MDATK's first group
                                vec3(0.922, 0.408, 0.204),   // #eb6834, its second group
                                vec3(0.298, 0.686, 0.314));  // #4caf50

void main() {
    float c = cos(scene.angle);
    float s = sin(scene.angle);
    vec2 p = mat2(c, s, -s, c) * kCorners[gl_VertexIndex];
    p.x /= max(scene.aspect, 1e-3);  // undistorted in any window shape
    gl_Position = vec4(p, 0.0, 1.0);
    colour = kColours[gl_VertexIndex];
}
