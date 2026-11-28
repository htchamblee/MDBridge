"""Base-pair-frame helical/step parameters (Olson et al. 2001 standard reference frame,
3DNA-style definitions written from the published definitions, not from 3DNA code).

For each duplex in each crystal copy / NMR model:
  base frame  = least-squares fit of the standard base ring atoms onto the observed ring atoms
  bp frame    = middle frame of base I and (base II rotated 180 deg about x)
  step k->k+1 = local helical twist/rise (screw of bp frame k onto bp frame k+1),
                inclination, tip, x-displacement; plus 3DNA-style step twist/rise/roll (CEHS)
Also computes the rev18 method (Kabsch of C1'+base atoms by name, screw) for comparison.

usage: python3 -I -B bpframes.py OUT.json file.cif [file.cif ...]
"""
import sys, json, shlex
import numpy as np

# Olson et al. 2001 (J Mol Biol 313:229) standard base coordinates (base frame), ring atoms + C1'
STD = {
 "A": {"C1'": (-2.479, 5.346, 0), "N9": (-1.291, 4.498, 0), "C8": (0.024, 4.897, 0), "N7": (0.877, 3.902, 0),
       "C5": (0.071, 2.771, 0), "C6": (0.369, 1.398, 0), "N1": (-0.668, 0.532, 0), "C2": (-1.912, 1.023, 0),
       "N3": (-2.320, 2.290, 0), "C4": (-1.267, 3.124, 0)},
 "G": {"C1'": (-2.477, 5.399, 0), "N9": (-1.289, 4.551, 0), "C8": (0.023, 4.962, 0), "N7": (0.870, 3.969, 0),
       "C5": (0.071, 2.833, 0), "C6": (0.424, 1.460, 0), "N1": (-0.700, 0.641, 0), "C2": (-1.999, 1.087, 0),
       "N3": (-2.342, 2.364, 0.001), "C4": (-1.265, 3.177, 0)},
 "C": {"C1'": (-2.477, 5.402, 0), "N1": (-1.285, 4.542, 0), "C2": (-1.472, 3.158, 0), "N3": (-0.391, 2.344, 0),
       "C4": (0.837, 2.868, 0), "C5": (1.056, 4.275, 0), "C6": (-0.023, 5.068, 0)},
 "T": {"C1'": (-2.481, 5.354, 0), "N1": (-1.284, 4.500, 0), "C2": (-1.462, 3.135, 0), "N3": (-0.298, 2.407, 0),
       "C4": (0.994, 2.897, 0), "C5": (1.106, 4.338, 0), "C6": (-0.024, 5.057, 0)},
 "U": {"C1'": (-2.481, 5.354, 0), "N1": (-1.284, 4.500, 0), "C2": (-1.462, 3.131, 0), "N3": (-0.302, 2.397, 0),
       "C4": (0.989, 2.884, 0), "C5": (1.089, 4.311, 0), "C6": (-0.024, 5.053, 0)},
}
RING = {b: [n for n in STD[b] if n != "C1'"] for b in STD}
COMP = {"A": "A", "G": "G", "C": "C", "U": "U", "DA": "A", "DG": "G", "DC": "C", "DT": "T",
        "6HA": "A", "6HG": "G", "6HC": "C", "6HT": "T"}
# rev18 method atom lists (C1' + base atoms; T adds C7)
BASE18 = {"A": ("N9","C8","N7","C5","C6","N6","N1","C2","N3","C4"),
          "G": ("N9","C8","N7","C5","C6","O6","N1","C2","N2","N3","C4"),
          "C": ("N1","C2","O2","N3","C4","N4","C5","C6"),
          "U": ("N1","C2","O2","N3","C4","O4","C5","C6"),
          "T": ("N1","C2","O2","N3","C4","O4","C5","C6")}

def read_cif(path):
    cols, rows, inloop = [], [], False
    for line in open(path):
        if line.startswith("_atom_site."):
            cols.append(line.split(".", 1)[1].strip()); inloop = True; continue
        if inloop and (line.startswith("ATOM") or line.startswith("HETATM")):
            rows.append(shlex.split(line)); continue
        if inloop and rows and not line.startswith(("ATOM", "HETATM")):
            break
    ix = {c: i for i, c in enumerate(cols)}
    models = {}
    for r in rows:
        if r[ix["label_alt_id"]] not in (".", "A"): continue
        comp = r[ix["label_comp_id"]]
        if comp not in COMP: continue
        m = int(r[ix["pdbx_PDB_model_num"]]); ch = r[ix["auth_asym_id"]]; seq = int(r[ix["auth_seq_id"]])
        nm = r[ix["label_atom_id"]]
        x = np.array([float(r[ix["Cartn_x"]]), float(r[ix["Cartn_y"]]), float(r[ix["Cartn_z"]])])
        res = models.setdefault(m, {}).setdefault((ch, seq), {"comp": comp, "b": COMP[comp], "at": {}})
        res["at"][nm] = x
    return models

def fit(P, Q):
    """R,t minimising |R P + t - Q|"""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc

def rot(axis, ang):
    a = axis / np.linalg.norm(axis); c, s = np.cos(ang), np.sin(ang)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + s * K + (1 - c) * K @ K

def base_frame(res):
    b = res["b"]; names = [n for n in RING[b] if n in res["at"]]
    if len(names) < len(RING[b]) - 1: return None
    R, t = fit(np.array([STD[b][n] for n in names]), np.array([res["at"][n] for n in names]))
    return R, t  # columns of R = x,y,z axes; t = origin

def mid_frame(F1, F2):
    (R1, o1), (R2, o2) = F1, F2
    z1, z2 = R1[:, 2], R2[:, 2]
    G = np.arccos(np.clip(z1 @ z2, -1, 1))
    if G < 1e-8: R1p, R2p, rt = R1, R2, R1[:, 1]
    else:
        rt = np.cross(z1, z2); rt /= np.linalg.norm(rt)
        R1p, R2p = rot(rt, G / 2) @ R1, rot(rt, -G / 2) @ R2
    x = R1p[:, 0] + R2p[:, 0]; x /= np.linalg.norm(x)
    y = R1p[:, 1] + R2p[:, 1]; y /= np.linalg.norm(y)
    z = np.cross(x, y); z /= np.linalg.norm(z); y = np.cross(z, x)
    return np.c_[x, y, z], (o1 + o2) / 2, (R1p, R2p, G, rt)

def signed(a, b, n):
    return np.degrees(np.arctan2(np.cross(a, b) @ n, a @ b))

def bp_frame(f1, f2):
    R2, o2 = f2
    R2f = R2 @ np.diag([1, -1, -1])
    Rm, om, _ = mid_frame(f1, (R2f, o2))
    return Rm, om

def step_params(F1, F2):
    (R1, o1), (R2, o2) = F1, F2
    out = {}
    # CEHS step parameters
    Rm, om, (R1p, R2p, G, rt) = mid_frame(F1, F2)
    zm = Rm[:, 2]
    out["step_twist"] = signed(R1p[:, 1], R2p[:, 1], zm)
    d = o2 - o1
    out["step_rise"] = d @ zm; out["shift"] = d @ Rm[:, 0]; out["slide"] = d @ Rm[:, 1]
    phi = np.radians(signed(rt, Rm[:, 1], zm)) if G > 1e-8 else 0.0
    out["roll"] = np.degrees(G) * np.cos(phi); out["tilt"] = np.degrees(G) * np.sin(phi)
    # local helical parameters: screw that maps frame 1 onto frame 2
    R = R2 @ R1.T
    th = np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1))
    w = np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]]); h = w / np.linalg.norm(w)
    if h @ (R1[:, 2] + R2[:, 2]) < 0: h, th = -h, -th
    out["h_twist"] = np.degrees(th)
    out["h_rise"] = h @ (o2 - o1)
    t = o2 - R @ o1
    p = np.linalg.pinv(np.eye(3) - R) @ (t - (t @ h) * h)
    inc, tip, xd, yd = [], [], [], []
    for Ri, oi in ((R1, o1), (R2, o2)):
        z = Ri[:, 2]; ph = np.arccos(np.clip(z @ h, -1, 1)); u = np.cross(z, h)
        if np.linalg.norm(u) > 1e-9:
            u /= np.linalg.norm(u); Rh = rot(u, ph) @ Ri
        else: Rh = Ri; u = Ri[:, 0]
        # inclination: rotation about x; tip: about y (sign chosen so that A-form inclination > 0)
        inc.append(-np.degrees(ph) * (u @ Ri[:, 0])); tip.append(np.degrees(ph) * (u @ Ri[:, 1]))
        dd = oi - p; dd -= (dd @ h) * h
        xd.append(dd @ Rh[:, 0]); yd.append(dd @ Rh[:, 1])
    out["inclination"] = float(np.mean(inc)); out["tip"] = float(np.mean(tip))
    out["x_disp"] = float(np.mean(xd)); out["y_disp"] = float(np.mean(yd))
    return {k: float(v) for k, v in out.items()}

def rev18_step(pairA, pairB):
    """rev18/ana.py method: Kabsch of C1'+base atoms matched by name, strand by strand; screw."""
    P, Q = [], []
    for ra, rb in zip(pairA, pairB):
        na = ("C1'",) + BASE18[ra["b"]]; nb = set(("C1'",) + BASE18[rb["b"]])
        for n in na:
            if n in nb and n in ra["at"] and n in rb["at"]: P.append(ra["at"][n]); Q.append(rb["at"][n])
    R, t = fit(np.array(P), np.array(Q))
    ang = np.degrees(np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1)))
    w = np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]]); ax = w / np.linalg.norm(w)
    rms = np.sqrt((((np.array(P) @ R.T + t) - np.array(Q)) ** 2).sum(1).mean())
    return float(ang), float(ax @ t), float(rms)

def pair_up(model):
    keys = sorted(model); fr = {}
    for k in keys:
        f = base_frame(model[k])
        if f is not None: fr[k] = f
    cand = []
    for i in fr:
        for j in fr:
            if i[0] >= j[0]: continue
            (Ri, oi), (Rj, oj) = fr[i], fr[j]
            dz = Ri[:, 2] @ Rj[:, 2]; do = np.linalg.norm(oi - oj)
            c1 = np.linalg.norm(model[i]["at"].get("C1'", oi + 99) - model[j]["at"].get("C1'", oj))
            if dz < -0.6 and do < 3.0 and 8.5 < c1 < 12.0: cand.append((do, i, j))
    cand.sort(); used = set(); pairs = []
    for do, i, j in cand:
        if i in used or j in used: continue
        used |= {i, j}; pairs.append((i, j))
    # group into duplexes by chain pair; order by strand-1 residue number
    dup = {}
    for i, j in pairs: dup.setdefault((i[0], j[0]), []).append((i, j))
    out = []
    for cp, ps in sorted(dup.items()):
        ps.sort(key=lambda p: p[0][1])
        if len(ps) < 4: continue
        out.append((cp, ps, fr))
    return out

def analyse(path):
    models = read_cif(path); res = []
    for m, model in sorted(models.items()):
        for cp, ps, fr in pair_up(model):
            bps = [bp_frame(fr[i], fr[j]) for i, j in ps]
            steps = []
            for k in range(len(ps) - 1):
                (i1, j1), (i2, j2) = ps[k], ps[k + 1]
                contiguous = (i2[1] - i1[1] == 1) and (j1[1] - j2[1] == 1)
                sp = step_params(bps[k], bps[k + 1])
                sp["contiguous"] = contiguous
                sp["seq"] = model[i1]["b"] + model[i2]["b"] + "/" + model[j2]["b"] + model[j1]["b"]
                sp["res"] = [i1[1], i2[1], j2[1], j1[1]]
                tw, ri, rms = rev18_step((model[i1], model[j1]), (model[i2], model[j2]))
                sp["rev18_twist"], sp["rev18_rise"], sp["rev18_rms"] = tw, ri, rms
                steps.append(sp)
            seq1 = "".join(model[i]["b"] for i, _ in ps); seq2 = "".join(model[j]["b"] for _, j in ps)
            res.append({"model": m, "chains": cp, "seq1": seq1, "seq2_rev": seq2,
                        "res1": [i[1] for i, _ in ps], "res2": [j[1] for _, j in ps], "steps": steps})
    return res

if __name__ == "__main__":
    out = {}
    for p in sys.argv[2:]:
        key = p.rsplit("/", 1)[-1].split(".")[0]
        out[key] = analyse(p)
        for d in out[key][:12]:
            print(key, "model", d["model"], "chains", d["chains"], d["seq1"], d["seq2_rev"], len(d["steps"]), "steps",
                  "noncontig", sum(not s["contiguous"] for s in d["steps"]))
        print(key, "total duplex instances", len(out[key]))
    json.dump(out, open(sys.argv[1], "w"), indent=0)
