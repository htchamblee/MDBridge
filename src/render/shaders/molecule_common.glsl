// Shared by the molecule scene renderer's shaders (src/render/vulkan/scene_renderer.cpp): the
// primitive records, the push constants (one block for every pipeline) and the software
// renderer's arithmetic (src/render/molecule/software.cpp) that the GPU repeats per pixel in
// float32. The geometry passes see only the records; molecule_compute.glsl adds the rest of the
// descriptor set for the compute passes.

// primitive records written by the CPU (layouts: scene_renderer.cpp, "records")
layout(std430, set = 0, binding = 0) readonly buffer Prims { uint prim[]; };

layout(push_constant, std430) uniform Push {
    uint W, H;          // image size (pixels)
    uint base, vbase;   // first record of the draw (words), first vertex record (meshes, words)
    float zmin, zspan;  // depth key range (float32 of the drawn geometry's depth range)
    float opacity;      // of the drawn item
    float hHalf;        // float32(H - 0.5) (sticks)
    float s;            // pixels per Angstrom, float32 (sticks)
    float fogLo, fogSpan, fogStrength;
    float fogR, fogG, fogB;
    float glowR, glowG, glowB;
    uint nGlow;         // glow values
    int pad;            // halo padding around the molecule's box (pixels)
    int radius;         // of the Gaussian weights
    uint flags;         // FLAG_*
    float norm;         // 1 / peak of the blurred glow disc
    float outline;      // ink outline strength (0: none)
    float zFar;         // float32(zmin - 1e3): depth of uncovered pixels for the outline test
    uint cw, ch;        // coarse grid size (row stride and rows)
    uint stride;        // words per record of the drawn kind (geometry passes)
} pc;

const uint FLAG_FOG = 1u;
const uint FLAG_AA = 2u;

// state[] entries
const int STATE_OP_BOX = 0;     // opaque pixels: min x, min y, max x, max y
const int STATE_HALO_BOX = 4;   // opaque + front translucent pixels
const int STATE_GMAX_OP = 8;    // largest glow of an opaque pixel (float bits; >= 0)
const int STATE_GMAX_HALO = 9;  // ... of an opaque or front translucent pixel
const int STATE_GLASS_GLOW = 10;  // 1 if a front translucent pixel has glow != 0

// ---- the software renderer's arithmetic ---------------------------------------------------------

const float kAmbient = 0.36;
const float kDiffuse = 0.64;
const float kSpecular = 0.22;
const float kInk = 0.12;
const float kGlowTint = 0.25;

// LIGHT = [-0.35, 0.55, 0.76] / norm in float32 (software.cpp: kLight)
vec3 lightDir() { return vec3(uintBitsToFloat(0xbeb2f7aeu), uintBitsToFloat(0x3f0c9e08u), uintBitsToFloat(0x3f424ebdu)); }

float primF(uint i) { return uintBitsToFloat(prim[i]); }
int primI(uint i) { return int(prim[i]); }

// Pixel (x, y) is in the convex region of K half-planes stored at `at` (K kinds: +1, -1 or 0 for a
// horizontal constraint, then K values b, then K values c; software.cpp halfPlanes / spanRow): the
// span of row centre yc runs from ceil(lo - 0.5) to floor(hi - 0.5), so an integer x is in it when
// x >= bound - 0.5 for every left bound and x <= bound - 0.5 for every right bound.
bool inRegion(uint at, uint K, float x, float yc) {
    for (uint k = 0u; k < K; ++k) {
        int kind = primI(at + k);
        precise float v = primF(at + K + k) * yc + primF(at + 2u * K + k);
        if (kind == 0) {
            if (v < 0.0) {
                return false;
            }
            continue;
        }
        precise float edge = -v - 0.5;
        if (kind > 0 ? !(x >= edge) : !(x <= edge)) {
            return false;
        }
    }
    return true;
}

// Ambient + Lambert + Phong highlight (shininess 24) from a light near the camera (MDATK's shade)
vec3 shade(vec3 n, vec3 rgb) {
    const vec3 L = lightDir();
    precise float len = sqrt(n.x * n.x + n.y * n.y + n.z * n.z);
    precise float inv = 1.0 / max(len, 1e-6);
    precise float ndlRaw = (n.x * L.x + n.y * L.y + n.z * L.z) * inv;
    float ndl = clamp(ndlRaw, 0.0, 1.0);
    precise float rv = clamp(2.0 * ndlRaw * n.z * inv - L.z, 0.0, 1.0) * (ndlRaw > 0.0 ? 1.0 : 0.0);
    precise float x2 = rv * rv;
    precise float x4 = x2 * x2;
    precise float x8 = x4 * x4;
    precise float diffuse = kAmbient + kDiffuse * ndl;
    precise float spec = kSpecular * x8 * x8 * x8;
    return clamp(rgb * diffuse + spec, 0.0, 1.0);
}

// astype(uint8) of a value clipped to 0..255 (NaN: 0)
uint toByte(float v) { return v >= 0.0 ? uint(min(v, 255.0)) : 0u; }
uint unitToByte(float v) { return toByte(clamp(v * 255.0 + 0.5, 0.0, 255.0)); }

uint packRgba(uvec4 c) { return c.r | (c.g << 8) | (c.b << 16) | (c.a << 24); }
uvec4 unpackRgba(uint v) { return uvec4(v & 255u, (v >> 8) & 255u, (v >> 16) & 255u, v >> 24); }

uvec4 packFragment(vec3 rgb, float alpha, float depth, int owner) {
    return uvec4(packUnorm2x16(rgb.rg), packUnorm2x16(vec2(rgb.b, alpha)), floatBitsToUint(depth), uint(owner + 2));
}
// depth test value of a fragment (the software renderer's depth key before quantisation)
float depthKey(float z) { return clamp((z - pc.zmin) / pc.zspan, 0.0, 1.0); }

