#version 450
#extension GL_GOOGLE_include_directive : require
// Sphere impostor (software.cpp sphereFrag): the sphere centred a quarter pixel off the pixel grid
// (its sub-pixel phase), covering the pixel centres inside its screen disc; normal and depth of the
// front surface. Record: box(4), cx, cy, sx, sy, rp, rp^2, Z, rgb(3), owner.

#include "molecule_common.glsl"

layout(location = 0) flat in uint record;
layout(location = 0) out uvec4 gbuffer;

void main() {
    uint r = record;
    int x = int(gl_FragCoord.x), y = int(gl_FragCoord.y);
    precise float dx = (float(x - primI(r + 4u)) + 0.5) - primF(r + 6u);
    precise float dy = (float(y - primI(r + 5u)) + 0.5) - primF(r + 7u);
    precise float d2 = dx * dx + dy * dy;
    float rp2 = primF(r + 9u);
    if (!(d2 < rp2)) {
        discard;
    }
    float rp = primF(r + 8u);
    precise float dz = sqrt(rp2 - d2);
    vec3 n = vec3(dx / rp, -dy / rp, dz / rp);
    precise float z = primF(r + 10u) + dz / pc.s;
    vec3 rgb = shade(n, vec3(primF(r + 11u), primF(r + 12u), primF(r + 13u)));
    gbuffer = packFragment(rgb, pc.opacity, z, primI(r + 14u));
    gl_FragDepth = depthKey(z);
}
