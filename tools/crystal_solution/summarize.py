"""Summarise params.json: per structure mean +- SD over interior / all / end steps; crystal - NMR differences.
usage: python3 -I -B summarize.py params.json"""
import sys, json
import numpy as np
D = json.load(open(sys.argv[1]))
PAIRS = [("hybrid r(GAAGAGAAGC)·d(GCTTCTCTTC) [same molecule]", ["404D"], ["1RRD"]),
         ("HNA·RNA (different sequences: 10-mer xtal vs 8-mer NMR)", ["2BJ6"], ["1EJZ"]),
         ("B-DNA Dickerson-Drew CGCGAATTCGCG [same molecule]", ["355D", "1BNA"], ["1NAJ", "1DUF"]),
         ("A-RNA r(GGCGUGCC) tandem G·U [same molecule]", ["2G3S"], ["1EKD"]),
         ("(ref) DNA·DNA d(GAAGAGAAGC)·d(GCTTCTCTTC) NMR", [], ["1AXP"])]
KEYS = ["h_twist", "h_rise", "inclination", "x_disp", "step_twist", "step_rise", "roll"]

def arr(entry, key):
    inst = D[entry]
    n = min(len(d["steps"]) for d in inst)
    return np.array([[s[key] for s in d["steps"][:n]] for d in inst])  # instances x steps

def stats(entry, key):
    A = arr(entry, key); n = A.shape[1]
    inner = slice(2, n - 2); ends = np.r_[0:2, n - 2:n]
    r = {"n_inst": A.shape[0], "n_steps": n,
         "all": A.mean(1).mean(), "int": A[:, inner].mean(1).mean(), "end": A[:, ends].mean(1).mean(),
         # step-to-step spread within one copy/model, averaged over copies/models
         "sd_all": A.std(1, ddof=1).mean(), "sd_int": A[:, inner].std(1, ddof=1).mean(),
         # spread of the copy/model means (copy-to-copy or model-to-model)
         "inst_sd_int": A[:, inner].mean(1).std(ddof=1) if A.shape[0] > 1 else float("nan"),
         # per-step SD across copies/models, averaged over steps (precision / copy variation)
         "perstep_sd": A.std(0, ddof=1).mean() if A.shape[0] > 1 else float("nan"),
         "profile": A.mean(0), "profile_sd": A.std(0, ddof=1) if A.shape[0] > 1 else np.zeros(n)}
    return r

def pooled(entries, key):
    ss = [stats(e, key) for e in entries]
    return {k: np.mean([s[k] for s in ss]) for k in ("all", "int", "end", "sd_all", "sd_int")}

f = lambda v: "%6.2f" % v
print("Per structure (twist deg, rise A, inclination deg, x-disp A). all = all steps; int = interior steps "
      "(2 terminal steps at each end excluded); end = the 4 terminal steps; sdS = step-to-step SD within one "
      "copy/model averaged over copies/models; sdC = SD across copies/models of the interior mean; "
      "sdP = per-step SD across copies/models averaged over steps")
for title, X, N in PAIRS:
    print("\n##", title)
    for e in X + N:
        print(" ", e, "instances", stats(e, "h_twist")["n_inst"], "steps", stats(e, "h_twist")["n_steps"])
        for k in KEYS:
            s = stats(e, k)
            print("    %-11s all %s  int %s  end %s  sdS(all) %s sdS(int) %s  sdC %s  sdP %s" % (
                k, f(s["all"]), f(s["int"]), f(s["end"]), f(s["sd_all"]), f(s["sd_int"]), f(s["inst_sd_int"]), f(s["perstep_sd"])))
    if X and N:
        print("  crystal - NMR (pooled over listed entries):")
        for k in KEYS[:4] + ["step_twist", "step_rise"]:
            a, b = pooled(X, k), pooled(N, k)
            print("    %-11s d(all) %s  d(int) %s  d(end) %s  d(sdS all) %s  d(sdS int) %s   [end-int: xtal %s nmr %s]" % (
                k, f(a["all"] - b["all"]), f(a["int"] - b["int"]), f(a["end"] - b["end"]),
                f(a["sd_all"] - b["sd_all"]), f(a["sd_int"] - b["sd_int"]), f(a["end"] - a["int"]), f(b["end"] - b["int"])))
    # per-step profile of helical twist and rise
    for k in ("h_twist", "h_rise"):
        print("  profile", k)
        for e in X + N:
            s = stats(e, k)
            print("    %-5s" % e, " ".join("%5.1f" % v if k == "h_twist" else "%5.2f" % v for v in s["profile"]),
                  "   +-", " ".join("%4.1f" % v if k == "h_twist" else "%4.2f" % v for v in s["profile_sd"]))
    e0 = (X + N)[0]
    print("  step labels", [d for d in [s["seq"] for s in D[e0][0]["steps"]]])

# rev18 method cross-check (only meaningful where both strands keep purine/pyrimidine class across the step)
print("\n## rev18 Kabsch-screw method vs bp-frame helical twist/rise, homopurine steps of 404D/1RRD/1AXP (steps 1-8)")
for e in ("404D", "1RRD", "1AXP"):
    a = arr(e, "rev18_twist")[:, :8]; b = arr(e, "h_twist")[:, :8]
    c = arr(e, "rev18_rise")[:, :8]; d = arr(e, "h_rise")[:, :8]
    print("  %s twist rev18 %.2f frame %.2f (max |diff| %.2f); rise rev18 %.2f frame %.2f (max |diff| %.2f)" % (
        e, a.mean(), b.mean(), np.abs(a - b).max(), c.mean(), d.mean(), np.abs(c - d).max()))
