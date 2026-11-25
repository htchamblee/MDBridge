#version 450
#extension GL_GOOGLE_include_directive : require
// Stick impostor (software.cpp stickFrag): the ray through the pixel centre against the infinite
// cylinder, kept between the stick's ends (open ends: the balls cover them), inside the screen
// capsule bound. Record: box(4), half-planes(4 kinds, 4 b, 4 c), Ax, Ay, Az, ux, uy, uz, a, R, L,
// rgb(3), owner (pixel units, y up).

#include "molecule_common.glsl"

layout(location = 0) flat in uint record;
layout(location = 0) out uvec4 gbuffer;

void main() {
    uint r = record;
    int x = int(gl_FragCoord.x), y = int(gl_FragCoord.y);
    if (!inRegion(r + 4u, 4u, float(x), float(y) + 0.5)) {
        discard;
    }
    float Ax = primF(r + 16u), Ay = primF(r + 17u), Az = primF(r + 18u);
    float ux = primF(r + 19u), uy = primF(r + 20u), uz = primF(r + 21u);
    float a = primF(r + 22u), R = primF(r + 23u), L = primF(r + 24u);
    precise float p0x = (float(x) + 0.5) - Ax;
    precise float p0y = (pc.hHalf - float(y)) - Ay;
    precise float m = p0x * ux + p0y * uy;
    precise float pu = m - Az * uz;
    precise float b = -uz * m - Az * a;
    precise float c = p0x * p0x + p0y * p0y + Az * Az - pu * pu - R * R;
    precise float disc = b * b - a * c;
    precise float tt = (-b + sqrt(max(disc, 0.0))) / a;
    precise float sp = pu + tt * uz;
    if (!(disc >= 0.0 && sp >= 0.0 && sp <= L)) {
        discard;
    }
    vec3 n = vec3(p0x - sp * ux, p0y - sp * uy, -Az + tt - sp * uz);
    precise float z = tt / pc.s;
    vec3 rgb = shade(n, vec3(primF(r + 25u), primF(r + 26u), primF(r + 27u)));
    gbuffer = packFragment(rgb, pc.opacity, z, primI(r + 28u));
    gl_FragDepth = depthKey(z);
}
