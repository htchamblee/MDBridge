#version 450
#extension GL_GOOGLE_include_directive : require
// Mesh triangle (software.cpp meshFrag / resolvePixel): pixel centres inside the slightly widened
// triangle and its box, plane-linear depth kept within the triangle's depth range, the normal
// interpolated, the colour of the residue with the largest barycentric share (its nearest vertex).
// Record: box(4), half-planes(3 kinds, 3 b, 3 c), a0, b0, c0, a1, b1, c1, za, zb, zc, zlo, zhi,
// v0, v1, v2, all relative to the box's corner (x0, y0). Vertex records (pc.vbase): normal(3),
// owner, rgb(3), unused.

#include "molecule_common.glsl"

layout(location = 0) flat in uint record;
layout(location = 0) out uvec4 gbuffer;

vec3 vertexNormal(uint v) {
    uint at = pc.vbase + v * 8u;
    return vec3(primF(at), primF(at + 1u), primF(at + 2u));
}
int vertexOwner(uint v) { return primI(pc.vbase + v * 8u + 3u); }
vec3 vertexRgb(uint v) {
    uint at = pc.vbase + v * 8u + 4u;
    return vec3(primF(at), primF(at + 1u), primF(at + 2u));
}

void main() {
    uint r = record;
    int x = int(gl_FragCoord.x) - primI(r), y = int(gl_FragCoord.y) - primI(r + 1u);
    precise float fx = float(x) + 0.5;
    precise float fy = float(y) + 0.5;
    if (!inRegion(r + 4u, 3u, float(x), fy)) {
        discard;
    }
    precise float l0 = primF(r + 13u) * fx + primF(r + 14u) * fy + primF(r + 15u);
    precise float l1 = primF(r + 16u) * fx + primF(r + 17u) * fy + primF(r + 18u);
    precise float zp = primF(r + 19u) * fx + primF(r + 20u) * fy + primF(r + 21u);
    float z = min(max(zp, primF(r + 22u)), primF(r + 23u));
    float a = clamp(l0, 0.0, 1.0), b = clamp(l1, 0.0, 1.0);
    float s = max(a + b, 1.0);
    a = a / s;
    b = b / s;
    precise float c = 1.0 - a - b;
    uint v0 = prim[r + 24u], v1 = prim[r + 25u], v2 = prim[r + 26u];
    vec3 n = vertexNormal(v0) * a + vertexNormal(v1) * b + vertexNormal(v2) * c;
    int o0 = vertexOwner(v0), o1 = vertexOwner(v1), o2 = vertexOwner(v2);
    float w0 = a + b * float(o1 == o0) + c * float(o2 == o0);
    float w1 = b + a * float(o0 == o1) + c * float(o2 == o1);
    float w2 = c + a * float(o0 == o2) + b * float(o1 == o2);
    int best = w0 >= w1 ? (w0 >= w2 ? o0 : o2) : (w1 >= w2 ? o1 : o2);
    float s0 = o0 == best ? a : -1.0, s1 = o1 == best ? b : -1.0, s2 = o2 == best ? c : -1.0;
    uint vert = s0 >= s1 ? (s0 >= s2 ? v0 : v2) : (s1 >= s2 ? v1 : v2);
    vec3 rgb = shade(n, vertexRgb(vert));
    gbuffer = packFragment(rgb, pc.opacity, z, vertexOwner(vert));
    gl_FragDepth = depthKey(z);
}
