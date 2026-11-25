// The molecule scene renderer's compute passes: every resource of the descriptor set and the
// compositing arithmetic of the software renderer (src/render/molecule/software.cpp composite).

#include "molecule_common.glsl"

// glow per residue row (0..1)
layout(std430, set = 0, binding = 1) readonly buffer Glow { float glow[]; };
// the halo's Gaussian weights (radius + 1 + radius values)
layout(std430, set = 0, binding = 2) readonly buffer Weights { float weight[]; };
// per-render state found on the GPU: bounding boxes and glow maxima (STATE_* below)
layout(std430, set = 0, binding = 3) coherent buffer State { int state[]; };
// the halo's coarse grid (two ping-pong buffers, row stride push.cw)
layout(std430, set = 0, binding = 4) buffer CoarseA { float coarseA[]; };
layout(std430, set = 0, binding = 5) buffer CoarseB { float coarseB[]; };
// the composited image before the edge anti-aliasing and the translucent layer (RGBA8 per uint)
layout(std430, set = 0, binding = 6) buffer Composite { uint comp[]; };
// the final image (RGBA8 per uint, row 0 at the top): read back or drawn into a window
layout(std430, set = 0, binding = 7) buffer Output { uint outImage[]; };
// the front-most fragment per pixel of the opaque and of the translucent items:
// x = rgb.rg (unorm16), y = (rgb.b, opacity) (unorm16), z = view depth (float bits), w = owner + 2
// (0: no fragment, 1: no owner)
layout(set = 0, binding = 8, rgba32ui) uniform readonly uimage2D gOpaque;
layout(set = 0, binding = 9, rgba32ui) uniform readonly uimage2D gGlass;


float glowOf(int owner) { return owner >= 0 && uint(owner) < pc.nGlow ? glow[owner] : 0.0; }

// glowing residues lighten toward the glow colour (screen blend)
vec3 tint(vec3 rgb, float gi) {
    precise float t = kGlowTint * gi;
    return vec3(1.0) - (vec3(1.0) - rgb) * (vec3(1.0) - t * vec3(pc.glowR, pc.glowG, pc.glowB));
}

// depth cue toward the panel colour over the panel's depth range
vec3 fogged(vec3 rgb, float depth) {
    precise float f = pc.fogStrength * (1.0 - clamp((depth - pc.fogLo) / pc.fogSpan, 0.0, 1.0));
    return rgb * (1.0 - f) + vec3(pc.fogR, pc.fogG, pc.fogB) * f;
}

// a G-buffer texel
struct Fragment {
    bool covered;
    vec3 rgb;
    float alpha;
    float depth;
    int owner;
};
Fragment fragmentOf(uvec4 g) {
    Fragment f;
    f.covered = g.w != 0u;
    vec2 rg = unpackUnorm2x16(g.x);
    vec2 ba = unpackUnorm2x16(g.y);
    f.rgb = vec3(rg, ba.x);
    f.alpha = ba.y;
    f.depth = uintBitsToFloat(g.z);
    f.owner = int(g.w) - 2;
    return f;
}
// the translucent fragment at a pixel is drawn when it lies in front of the opaque one
bool frontGlass(Fragment op, Fragment gl) { return gl.covered && (!op.covered || gl.depth > op.depth + 1e-3); }

// ---- the halo's region (software.cpp glowHalo): the molecule's box + 2 px, padded, aligned to the
// coarse cells
struct Region {
    bool any;
    int gx0, gy0, gx1, gy1;  // pixels (end exclusive)
    int ww, hh;              // coarse cells
};
Region haloRegion() {
    Region r;
    int minx = state[STATE_HALO_BOX + 0], miny = state[STATE_HALO_BOX + 1];
    int maxx = state[STATE_HALO_BOX + 2], maxy = state[STATE_HALO_BOX + 3];
    r.any = maxx >= 0;
    int W = int(pc.W), H = int(pc.H);
    int y0 = max(miny - 2, 0), y1 = min(maxy + 3, H);
    int x0 = max(minx - 2, 0), x1 = min(maxx + 3, W);
    r.gy0 = (max(y0 - pc.pad, 0) / 4) * 4;
    r.gx0 = (max(x0 - pc.pad, 0) / 4) * 4;
    r.gy1 = min(y1 + pc.pad, H);
    r.gx1 = min(x1 + pc.pad, W);
    r.hh = (r.gy1 - r.gy0 + 3) / 4;
    r.ww = (r.gx1 - r.gx0 + 3) / 4;
    return r;
}
bool haloOn() { return uintBitsToFloat(uint(state[STATE_GMAX_HALO])) > 0.01; }
