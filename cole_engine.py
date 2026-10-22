"""
Cole Model Engine & Utilities for Bioimpedance Spectroscopy (EIS)
Handles Cole equation evaluations, Nyquist, Bode, and Non-Linear Least Squares Fitting.
"""

import numpy as np
import scipy.io as sio
from scipy.optimize import least_squares
from typing import Dict, Tuple, List, Optional, Any

# Table 1 Biological Tissues defined in the assignment
TABLE_1_TISSUES = {
    "Músculo Long.": {
        "R0": 385.0,
        "Rinf": 70.0,
        "fc": 37.0 * 1e3,  # 37 kHz in Hz
        "alpha": 0.27,
        "color": "#2563EB",  # Deep royal sapphire
    },
    "Músculo Trans.": {
        "R0": 730.0,
        "Rinf": 75.0,
        "fc": 27.0 * 1e3,  # 27 kHz in Hz
        "alpha": 0.24,
        "color": "#0D9488",  # Rich teal
    },
    "Hígado": {
        "R0": 540.0,
        "Rinf": 100.0,
        "fc": 72.0 * 1e3,  # 72 kHz in Hz
        "alpha": 0.45,
        "color": "#D97706",  # Warm amber/ochre
    },
    "Pulmón": {
        "R0": 440.0,
        "Rinf": 70.0,
        "fc": 100.0 * 1e3,  # 100 kHz in Hz
        "alpha": 0.59,
        "color": "#7C3AED",  # Royal amethyst
    },
    "Bazo": {
        "R0": 550.0,
        "Rinf": 60.0,
        "fc": 373.0 * 1e3,  # 373 kHz in Hz
        "alpha": 0.36,
        "color": "#E11D48",  # Vibrant crimson rose
    },
}

ORIGIN_LABELS = {
    "1": "Dell'Osa et al. (2019/2025)",
    "2": "Gonzalez Correa (2018)",
    "3": "De Lorenzo et al. (1998)",
    "4": "Ayllon et al. (2016)",
    "5": "Bertemes Filho (Ex-vivo)",
}

SAMPLE_LABELS = {
    "0": "Sin datos",
    "1": "Tejido Humano",
    "2": "Vegetal",
    "3": "Carne",
    "4": "Animal / Ex-vivo",
    "5": "Fantoma Agar-Agar",
}

DEVICE_LABELS = {
    "0": "Desconocido",
    "1": "Solartron",
    "2": "HF2IS Zurich Instruments",
    "3": "Xitron 4001B",
    "4": "SFB7 Impedimed",
}


def cole_impedance(f: np.ndarray, R0: float, Rinf: float, fc: float, alpha: float) -> np.ndarray:
    """
    Computes complex impedance Z(omega) using the Cole Model:
    Z(omega) = Rinf + (R0 - Rinf) / (1 + (j * omega * tau)^alpha)
    where omega = 2*pi*f, tau = 1/(2*pi*fc).
    """
    f = np.asarray(f, dtype=float)
    omega = 2.0 * np.pi * f
    tau = 1.0 / (2.0 * np.pi * fc)
    # j * omega * tau = omega * tau * e^(j*pi/2)
    # (j*omega*tau)^alpha:
    term = (1j * omega * tau) ** alpha
    Z = Rinf + (R0 - Rinf) / (1.0 + term)
    return Z


def compute_cole_components(f: np.ndarray, R0: float, Rinf: float, fc: float, alpha: float) -> Dict[str, np.ndarray]:
    """
    Calculates Real, Imaginary, -Im(Z), Modulus |Z|, and Phase (degrees).
    """
    Z = cole_impedance(f, R0, Rinf, fc, alpha)
    re = np.real(Z)
    im = np.imag(Z)
    modulus = np.abs(Z)
    phase_deg = np.angle(Z, deg=True)
    return {
        "f": f,
        "Z": Z,
        "re": re,
        "im": im,
        "minus_im": -im,
        "modulus": modulus,
        "phase_deg": phase_deg,
    }


def load_dataset(mat_path: str = "Data/registros_EIS.mat") -> List[Dict[str, Any]]:
    """
    Loads and structures the records from registros_EIS.mat into a clean Python list of dicts.
    """
    mat = sio.loadmat(mat_path)
    key = [k for k in mat if not k.startswith("__")][0]
    raw_records = mat[key][0]
    
    parsed = []
    for i, r in enumerate(raw_records):
        try:
            rec_id = str(r["id"].squeeze()).replace("[", "").replace("]", "").replace("'", "").strip()
            dev = str(r["device"].squeeze()).replace("[", "").replace("]", "").replace("'", "").strip()
            sample = str(r["sample"].squeeze()).replace("[", "").replace("]", "").replace("'", "").strip()
            origin = str(r["origin"].squeeze()).replace("[", "").replace("]", "").replace("'", "").strip()
            eis = r["eis"]
            
            if eis.ndim == 2 and eis.shape[1] >= 3 and len(eis) >= 4:
                parsed.append({
                    "index": i,
                    "id": rec_id,
                    "device": dev,
                    "sample": sample,
                    "origin": origin,
                    "device_label": DEVICE_LABELS.get(dev, f"Dev {dev}"),
                    "sample_label": SAMPLE_LABELS.get(sample, f"Tipo {sample}"),
                    "origin_label": ORIGIN_LABELS.get(origin, f"Origen {origin}"),
                    "points": len(eis),
                    "f_min": float(np.min(eis[:, 0])),
                    "f_max": float(np.max(eis[:, 0])),
                    "eis": eis,
                })
        except Exception:
            continue
            
    return parsed


def fit_cole_model(
    f_meas: np.ndarray,
    re_meas: np.ndarray,
    im_meas: np.ndarray,
    p0: Optional[List[float]] = None,
    weights: str = "modulus"
) -> Dict[str, Any]:
    """
    Fits experimental EIS data to the Cole equation using Non-Linear Least Squares (Trust Region Reflective).
    Supports automatic sign alignment for imaginary impedance conventions.
    """
    f = np.asarray(f_meas, dtype=float)
    re = np.asarray(re_meas, dtype=float)
    im = np.asarray(im_meas, dtype=float)

    # Clean non-finite values if any
    mask = np.isfinite(f) & np.isfinite(re) & np.isfinite(im) & (f > 0)
    f, re, im = f[mask], re[mask], im[mask]

    if len(f) < 4:
        raise ValueError("Se necesitan al menos 4 puntos válidos para ajustar el modelo de Cole.")

    # Sort by frequency ascending
    sort_idx = np.argsort(f)
    f, re, im = f[sort_idx], re[sort_idx], im[sort_idx]

    # In Cole theoretical model: Imaginary part is NEGATIVE.
    # If the dataset provided positive imaginary part, we invert it for fitting.
    sign_im = -1.0 if np.mean(im) > 0 else 1.0
    im_aligned = im * sign_im  # Now im_aligned will have negative values (capacitive)

    # Initial guesses
    re_max = float(np.percentile(re, 95))
    re_min = float(np.percentile(re, 5))
    r0_init = max(re_max, 1.0)
    rinf_init = max(min(re_min, r0_init * 0.8), 0.1)

    # Peak imaginary frequency estimate
    idx_peak = np.argmin(im_aligned)  # most negative
    fc_init = float(f[idx_peak]) if 0 < idx_peak < len(f) - 1 else float(np.median(f))
    alpha_init = 0.65

    if p0 is None:
        p0 = [r0_init, rinf_init, fc_init, alpha_init]

    # Bounds: R0 > 0, Rinf > 0, fc > 0, 0 < alpha <= 1.0
    bounds_lower = [0.01, 0.01, float(np.min(f) * 0.01), 0.05]
    bounds_upper = [float(np.max(re) * 10.0 + 1000.0), float(np.max(re) * 5.0 + 500.0), float(np.max(f) * 100.0), 1.0]

    # Ensure p0 is strictly inside bounds
    p0 = [
        np.clip(p0[0], bounds_lower[0] * 1.1, bounds_upper[0] * 0.9),
        np.clip(p0[1], bounds_lower[1] * 1.1, bounds_upper[1] * 0.9),
        np.clip(p0[2], bounds_lower[2] * 1.1, bounds_upper[2] * 0.9),
        np.clip(p0[3], 0.1, 0.95),
    ]

    def residual_fn(params):
        R0, Rinf, fc, alpha = params
        w = 2.0 * np.pi * f
        tau = 1.0 / (2.0 * np.pi * fc)
        denom = 1.0 + (1j * w * tau) ** alpha
        Z_theo = Rinf + (R0 - Rinf) / denom
        
        re_pred = np.real(Z_theo)
        im_pred = np.imag(Z_theo)

        if weights == "modulus":
            w_norm = np.sqrt(re**2 + im_aligned**2) + 1e-4
            res_re = (re_pred - re) / w_norm
            res_im = (im_pred - im_aligned) / w_norm
        else:
            res_re = (re_pred - re) / (np.abs(re) + 1.0)
            res_im = (im_pred - im_aligned) / (np.abs(im_aligned) + 1.0)

        return np.concatenate([res_re, res_im])

    res = least_squares(residual_fn, p0, bounds=(bounds_lower, bounds_upper), method="trf", max_nfev=5000)

    R0_fit, Rinf_fit, fc_fit, alpha_fit = res.x

    # Compute fitted impedance over original f and dense f
    Z_fit_data = cole_impedance(f, R0_fit, Rinf_fit, fc_fit, alpha_fit)
    re_fit = np.real(Z_fit_data)
    im_fit = np.imag(Z_fit_data) * sign_im  # Restore original dataset convention

    # Metrics
    ss_tot_re = np.sum((re - np.mean(re)) ** 2)
    ss_res_re = np.sum((re - re_fit) ** 2)
    r2_re = 1.0 - (ss_res_re / ss_tot_re) if ss_tot_re > 1e-9 else 1.0

    ss_tot_im = np.sum((im - np.mean(im)) ** 2)
    ss_res_im = np.sum((im - im_fit) ** 2)
    r2_im = 1.0 - (ss_res_im / ss_tot_im) if ss_tot_im > 1e-9 else 1.0

    rmse_re = float(np.sqrt(np.mean((re - re_fit) ** 2)))
    rmse_im = float(np.sqrt(np.mean((im - im_fit) ** 2)))
    
    # Combined R2
    y_true = np.concatenate([re, im])
    y_pred = np.concatenate([re_fit, im_fit])
    r2_total = 1.0 - np.sum((y_true - y_pred) ** 2) / np.sum((y_true - np.mean(y_true)) ** 2)

    return {
        "success": bool(res.success),
        "params": {
            "R0": float(R0_fit),
            "Rinf": float(Rinf_fit),
            "fc": float(fc_fit),
            "fc_kHz": float(fc_fit / 1000.0),
            "alpha": float(alpha_fit),
            "tau_us": float((1.0 / (2.0 * np.pi * fc_fit)) * 1e6),
            "delta_R": float(R0_fit - Rinf_fit),
        },
        "metrics": {
            "r2_total": float(r2_total),
            "r2_re": float(r2_re),
            "r2_im": float(r2_im),
            "rmse_re": rmse_re,
            "rmse_im": rmse_im,
        },
        "f_meas": f,
        "re_meas": re,
        "im_meas": im,
        "minus_im_meas": -im if sign_im == 1.0 else im,
        "sign_im_convention": "standard_negative" if sign_im == 1.0 else "inverted_positive",
    }
