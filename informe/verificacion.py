import sys, json
import numpy as np
sys.path.insert(0, ".")
import scipy.io as sio
from scipy.optimize import curve_fit
import cole_engine as ce

out = {}
# --- A. formula checks
for n, p in ce.TABLE_1_TISSUES.items():
    R0, Ri, fc, a = p["R0"], p["Rinf"], p["fc"], p["alpha"]
    Zlo = ce.cole_impedance(np.array([1e-3]), R0, Ri, fc, a)[0]
    Zhi = ce.cole_impedance(np.array([1e12]), R0, Ri, fc, a)[0]
    Zfc = ce.cole_impedance(np.array([fc]), R0, Ri, fc, a)[0]
    # closed form at fc: (jwt)^a = e^{j a pi/2}
    ref = Ri + (R0 - Ri) / (1 + np.exp(1j * a * np.pi / 2))
    # peak of -Im should be at fc
    f = np.logspace(1, 9, 200001)
    Z = ce.cole_impedance(f, R0, Ri, fc, a)
    fpk = f[np.argmax(-Z.imag)]
    # analytic peak -Im = dR/2 * tan(a*pi/4) / ... : -Im(fc)=dR*sin(a pi/2)/(2(1+cos(a pi/2)))
    dR = R0 - Ri
    xm = dR * np.sin(a * np.pi / 2) / (2 * (1 + np.cos(a * np.pi / 2)))
    out[n] = dict(Zlo=[round(Zlo.real, 3), round(Zlo.imag, 5)], Zhi=[round(Zhi.real, 3), round(Zhi.imag, 5)],
                  err_fc=abs(Zfc - ref), fpk_kHz=fpk / 1e3, fc_kHz=fc / 1e3, peak_im=float((-Z.imag).max()),
                  peak_analytic=float(xm), max_im=float(Z.imag.max()),
                  phase_min=float(np.angle(Z, deg=True).min()))
print("=== A: formula/limits")
for k, v in out.items(): print(k, v)

# --- B. dataset consistency
recs = ce.load_dataset("Data/registros_EIS.mat")
print("\n=== B: dataset  n parsed =", len(recs))
mat = sio.loadmat("Data/registros_EIS.mat")
key = [k for k in mat if not k.startswith("__")][0]
raw = mat[key][0]
print("raw n =", len(raw), "fields", raw.dtype.names)
ids = [str(np.squeeze(r["id"])) for r in raw]
print("ids first/last", ids[:3], ids[-3:])
# index vs id relation and origin per id range
org = {}
for i, r in enumerate(raw):
    o = str(np.squeeze(r["origin"]))
    org.setdefault(o, []).append(i)
for o, ix in org.items(): print("origin", o, "n", len(ix), "idx", min(ix), max(ix))
for i in (1311, 1312, 1350, 1400, 1450):
    r = raw[i]
    print("idx", i, "id", np.squeeze(r["id"]), "origin", np.squeeze(r["origin"]), "sample", np.squeeze(r["sample"]),
          "dev", np.squeeze(r["device"]), "eis", r["eis"].shape)
dropped = [i for i in range(len(raw)) if i not in {r["index"] for r in recs}]
print("dropped by loader:", len(dropped), dropped[:10])
# columns 4,5 consistent with 2,3 ?
bad_mod = bad_ph = 0; signs = {"neg": 0, "pos": 0, "mixed": 0}
for r in recs:
    e = r["eis"]
    if e.shape[1] >= 5:
        if not np.allclose(e[:, 3], np.hypot(e[:, 1], e[:, 2]), rtol=1e-2): bad_mod += 1
        ph = np.degrees(np.arctan2(e[:, 2], e[:, 1]))
        if not (np.allclose(e[:, 4], ph, atol=1.0) or np.allclose(e[:, 4], np.radians(ph), atol=0.02)): bad_ph += 1
    im = e[:, 2]
    signs["neg" if np.all(im <= 0) else "pos" if np.all(im >= 0) else "mixed"] += 1
print("records with col4!=modulus", bad_mod, "col5!=phase", bad_ph, "imag sign", signs)
e = recs[1311]["eis"] if len(recs) > 1311 else None
print("example eis head (idx 1311):\n", raw[1311]["eis"][:3])

# --- C. fit check: app fit vs independent fit (complex, scipy curve_fit, different start/weights)
def cole_c(f, R0, Ri, fc, a):
    return ce.cole_impedance(f, R0, Ri, fc, a)

def indep_fit(f, re, im):
    sgn = -1 if np.mean(im) > 0 else 1
    z = re + 1j * im * sgn
    best = None
    for a0 in (0.3, 0.6, 0.9):
        for fc0 in (f.min() * 3, np.sqrt(f.min() * f.max()), f.max() / 3):
            def fun(f_, R0, Ri, lfc, a):
                zz = cole_c(f_, R0, Ri, 10 ** lfc, a)
                return np.concatenate([zz.real, zz.imag])
            try:
                p, _ = curve_fit(fun, f, np.concatenate([z.real, z.imag]),
                                 p0=[z.real.max(), z.real.min(), np.log10(fc0), a0],
                                 bounds=([0, 0, 0, 0.01], [1e5, 1e5, 9, 1.0]), max_nfev=20000)
                r = np.sum((fun(f, *p) - np.concatenate([z.real, z.imag])) ** 2)
                if best is None or r < best[0]: best = (r, p)
            except Exception:
                pass
    r, p = best
    return dict(R0=p[0], Rinf=p[1], fc=10 ** p[2], alpha=p[3], sse=r)

print("\n=== C: fits (app vs independent) for Dell'Osa (origin 1)")
rows = []
for r in recs:
    if r["origin"] != "1": continue
    e = r["eis"]; f, re, im = e[:, 0], e[:, 1], e[:, 2]
    try:
        a = ce.fit_cole_model(f, re, im)
        b = indep_fit(f, re, im)
    except Exception as ex:
        rows.append((r["index"], r["id"], "ERR", str(ex))); continue
    rows.append((r["index"], r["id"], r["points"], r["sample"], r["device"], a["params"]["R0"], a["params"]["Rinf"],
                 a["params"]["fc"], a["params"]["alpha"], a["metrics"]["r2_total"], a["success"],
                 b["R0"], b["Rinf"], b["fc"], b["alpha"]))
errs = [x for x in rows if x[2] == "ERR"]
ok = [x for x in rows if x[2] != "ERR"]
print("n origin1:", len(rows), "errors:", len(errs), errs[:3])
A = np.array([[x[5], x[6], x[7], x[8]] for x in ok]); B = np.array([[x[11], x[12], x[13], x[14]] for x in ok])
rel = np.abs(A - B) / np.abs(B)
print("rel diff median per param (R0,Rinf,fc,alpha):", np.median(rel, axis=0))
print("rel diff p95:", np.percentile(rel, 95, axis=0))
r2 = np.array([x[9] for x in ok]); print("R2 min/median/max", r2.min(), np.median(r2), r2.max())
print("success False:", sum(1 for x in ok if not x[10]))
print("alpha at bound (>=0.99 or <=0.06):", int(np.sum((A[:, 3] >= 0.99) | (A[:, 3] <= 0.06))))
print("fc outside meas range:", int(sum(1 for x in ok if not (recs[x[0]]["f_min"] <= x[7] <= recs[x[0]]["f_max"]))))
print("worst 5 by R2:", [(x[0], x[1], round(x[9], 4)) for x in sorted(ok, key=lambda t: t[9])[:5]])
print("largest alpha discrepancies:", [(x[0], round(x[8], 3), round(x[14], 3)) for x in sorted(ok, key=lambda t: -abs(t[8] - t[14]))[:5]])
print("\nPoint-4 samples:")
for i in (1311, 1312, 1350, 1400, 1450):
    x = [t for t in ok if t[0] == i]
    print(i, x[0][1:] if x else "not origin1/not parsed")

# --- D. R2 definition check on one record
i0 = ok[0][0]; e = recs[i0]["eis"]; res = ce.fit_cole_model(e[:, 0], e[:, 1], e[:, 2])
print("\nR2 reported vs recomputed:", res["metrics"]["r2_total"], res["metrics"]["r2_re"], res["metrics"]["r2_im"], "sign convention", res["sign_im_convention"])
# test: synthetic recovery
f = np.logspace(2, 7, 60)
for n, p in ce.TABLE_1_TISSUES.items():
    z = ce.cole_impedance(f, p["R0"], p["Rinf"], p["fc"], p["alpha"])
    rng = np.random.default_rng(0)
    z = z + rng.normal(0, 0.2, z.shape) + 1j * rng.normal(0, 0.2, z.shape)
    fr = ce.fit_cole_model(f, z.real, z.imag)["params"]
    print("synthetic", n, {k: round(fr[k], 3) for k in ("R0", "Rinf", "fc_kHz", "alpha")},
          "true", p["R0"], p["Rinf"], p["fc"] / 1e3, p["alpha"])
