#!/usr/bin/env python3
"""Edge-wise desmoplastic conductances on a four-node anatomical graph.

Toy identifiability laboratory for Thesis #21. Research only.
Not a medical device, not a dose, not a patient atlas.

Generator G1 (primary). Each directed edge carries a shedding rate λ and a
desmoplastic conductance κ. Node burdens see only the series flux
    φ = λ κ / (λ + κ).
A barrier reading, when scheduled, is the stalled fraction
    s = λ / (λ + κ),
recorded once per edge. Soil (r, K) is known unless a check frees r.

Generator G2 (check). The same edges are hold-up compartments,
    ċ = λ x_src − κ c,  influx to the target = κ c.
Node-only schedules do not see c. Barrier-aware schedules do.

Seed 20260921 is used only for the one noisy draw. Every other path is
deterministic. Regenerating this file rewrites sim/results.json and
sim/figures/.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parent
FIG = ROOT / "figures"
SEED = 20260921
CHI2_95 = 3.841458820694124  # chi-square, 1 df, 95%

# --- declared toy -----------------------------------------------------------
# Nodes: primary P, lung filter F, liver V, bone B.
# Edges: P→F, F→V (transport-limited), F→B.
EDGES = ((0, 1), (1, 2), (1, 3))
NODE_NAMES = ("P", "F", "V", "B")
EDGE_NAMES = ("P→F", "F→V", "F→B")
R = np.array([0.15, 0.07, 0.05, 0.04])
K = np.array([8.0, 5.0, 4.0, 3.0])
X0 = np.array([1.0, 0.0, 0.0, 0.0])
LAM_TRUE = np.array([0.09, 0.06, 0.04])
KAP_TRUE = np.array([0.25, 0.025, 0.04])
T_OBS = np.linspace(0.0, 36.0, 37)
SIGMA_Y = 0.05
SIGMA_S = 0.02
SIGMA_C = 0.02
# Finite-difference step on log-parameters. Structural rank: singular values
# above 1e-6 of the leading one (above the finite-difference floor, below
# any direction the trajectories actually use). Practical rank: principal
# relative standard errors strictly below 1/2.
FD_STEP = 1e-4
STRUCT_RATIO = 1e-6
PRACTICAL_RSE = 0.5
LOG_LO = np.log(1e-4)
LOG_HI = np.log(20.0)

# Profile grid for λ and κ of the transport-limited edge (index 1), as
# multiples of the truth, log-spaced. The lower end crosses φ, where no
# positive partner exists.
PROFILE_MULT = np.geomspace(0.22, 6.0, 17)


def series_phi(lam: np.ndarray, kap: np.ndarray) -> np.ndarray:
    return lam * kap / (lam + kap)


def stalled_fraction(lam: np.ndarray, kap: np.ndarray) -> np.ndarray:
    return lam / (lam + kap)


def partner_kappa(phi: np.ndarray, lam: np.ndarray) -> np.ndarray:
    """κ that keeps φ when λ > φ. Caller must guarantee the inequality."""
    return phi * lam / (lam - phi)


def _rhs_series(t, x, phi, r, k):
    dx = r * x * (1.0 - x / k)
    for e, (src, tgt) in enumerate(EDGES):
        flux = phi[e] * x[src]
        dx[src] -= flux
        dx[tgt] += flux
    return dx


def _rhs_holdup(t, z, lam, kap, r, k):
    x = z[:4]
    c = z[4:]
    dx = r * x * (1.0 - x / k)
    dc = np.empty(3)
    for e, (src, tgt) in enumerate(EDGES):
        shed = lam[e] * x[src]
        cross = kap[e] * c[e]
        dx[src] -= shed
        dx[tgt] += cross
        dc[e] = shed - cross
    return np.concatenate([dx, dc])


def simulate_series(lam, kap, r=R, k=K, t_obs=T_OBS) -> np.ndarray:
    phi = series_phi(np.asarray(lam, float), np.asarray(kap, float))
    sol = solve_ivp(
        _rhs_series,
        (float(t_obs[0]), float(t_obs[-1])),
        X0,
        t_eval=t_obs,
        args=(phi, np.asarray(r, float), np.asarray(k, float)),
        method="RK45",
        rtol=1e-8,
        atol=1e-9,
        dense_output=False,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y  # (4, n_t)


def simulate_phi(phi, r=R, k=K, t_obs=T_OBS) -> np.ndarray:
    """Single-rate edge model: the node field is given φ directly."""
    # Encode φ as a series pair with κ → the same φ by using λ=φ, κ huge
    # is unnecessary: pass φ straight through a dummy pair identity.
    sol = solve_ivp(
        _rhs_series,
        (float(t_obs[0]), float(t_obs[-1])),
        X0,
        t_eval=t_obs,
        args=(np.asarray(phi, float), np.asarray(r, float), np.asarray(k, float)),
        method="RK45",
        rtol=1e-8,
        atol=1e-9,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y


def simulate_holdup(lam, kap, r=R, k=K, t_obs=T_OBS) -> tuple[np.ndarray, np.ndarray]:
    z0 = np.zeros(7)
    z0[:4] = X0
    sol = solve_ivp(
        _rhs_holdup,
        (float(t_obs[0]), float(t_obs[-1])),
        z0,
        t_eval=t_obs,
        args=(np.asarray(lam, float), np.asarray(kap, float), np.asarray(r, float), np.asarray(k, float)),
        method="RK45",
        rtol=1e-8,
        atol=1e-9,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y[:4], sol.y[4:]


def output_series(nodes: np.ndarray, lam, kap, schedule: str) -> np.ndarray:
    """Unweighted observation vector for a series-model schedule."""
    lam = np.asarray(lam, float)
    kap = np.asarray(kap, float)
    lump = nodes.sum(axis=0)
    parts: list[np.ndarray] = []
    if schedule in ("L", "LB", "LP"):
        parts.append(lump)
    elif schedule in ("S", "SB", "SP"):
        parts.append(nodes.ravel())
    else:
        raise KeyError(schedule)
    s = stalled_fraction(lam, kap)
    if schedule == "LB":
        parts.append(s)
    elif schedule == "LP":
        parts.append(s[1:2])  # transport-limited edge only
    elif schedule == "SB":
        parts.append(s)
    elif schedule == "SP":
        parts.append(s[1:2])
    return np.concatenate(parts)


def sigma_series(schedule: str) -> np.ndarray:
    n_t = T_OBS.size
    if schedule in ("L", "LB", "LP"):
        sig = [np.full(n_t, SIGMA_Y)]
    else:
        sig = [np.full(4 * n_t, SIGMA_Y)]
    if schedule in ("LB", "SB"):
        sig.append(np.full(3, SIGMA_S))
    elif schedule in ("LP", "SP"):
        sig.append(np.full(1, SIGMA_S))
    return np.concatenate(sig)


def output_holdup(nodes, comps, schedule: str) -> np.ndarray:
    lump = nodes.sum(axis=0)
    if schedule == "L":
        return lump.copy()
    if schedule == "S":
        return nodes.ravel().copy()
    if schedule == "LB":
        return np.concatenate([lump, comps.ravel()])
    if schedule == "SB":
        return np.concatenate([nodes.ravel(), comps.ravel()])
    raise KeyError(schedule)


def sigma_holdup(schedule: str) -> np.ndarray:
    n_t = T_OBS.size
    if schedule == "L":
        return np.full(n_t, SIGMA_Y)
    if schedule == "S":
        return np.full(4 * n_t, SIGMA_Y)
    if schedule == "LB":
        return np.concatenate([np.full(n_t, SIGMA_Y), np.full(3 * n_t, SIGMA_C)])
    if schedule == "SB":
        return np.concatenate([np.full(4 * n_t, SIGMA_Y), np.full(3 * n_t, SIGMA_C)])
    raise KeyError(schedule)


def weighted_sensitivity(predict, theta: np.ndarray, sigma: np.ndarray) -> np.ndarray:
    """Columns are ∂y/∂logθ, divided by σ. Central differences."""
    theta = np.asarray(theta, float)
    base_unused = None
    cols = []
    step = FD_STEP
    for j in range(theta.size):
        up = theta.copy()
        dn = theta.copy()
        up[j] *= np.exp(step)
        dn[j] *= np.exp(-step)
        y_up = predict(up)
        y_dn = predict(dn)
        cols.append((y_up - y_dn) / (2.0 * step))
    jac = np.column_stack(cols) / sigma[:, None]
    if base_unused is not None:
        pass
    return jac


def fisher_report(jac: np.ndarray) -> dict:
    singular = np.linalg.svd(jac, compute_uv=False)
    singular = np.maximum(singular, 0.0)
    leading = singular[0] if singular.size else 0.0
    ratio = singular / leading if leading > 0 else singular
    structural = int(np.sum(ratio > STRUCT_RATIO))
    # Principal relative SE is 1/s for each retained component. Components
    # on the structural floor are reported as null, not as huge finite errors.
    principal = []
    for s, r in zip(singular, ratio):
        if r > STRUCT_RATIO and s > 0:
            principal.append(float(1.0 / s))
        else:
            principal.append(None)
    practical = int(sum(p is not None and p < PRACTICAL_RSE for p in principal))
    # Coordinate-wise Cramér–Rao sketch from the Moore–Penrose inverse,
    # infinite where the column is in the numerical kernel.
    gram = jac.T @ jac
    evals, evecs = np.linalg.eigh(gram)
    order = np.argsort(evals)[::-1]
    evals = evals[order]
    evecs = evecs[:, order]
    keep = evals > (leading**2) * (STRUCT_RATIO**2)
    cov_diag = np.full(jac.shape[1], np.inf)
    if np.any(keep):
        inv = (evecs[:, keep] * (1.0 / evals[keep])) @ evecs[:, keep].T
        cov_diag = np.clip(np.diag(inv), 0.0, None)
        # A coordinate with a large projection on the kernel is not a
        # finite standard error. Flag it when half or more of its energy
        # sits on discarded components.
        energy = evecs**2
        kernel_share = energy[:, ~keep].sum(axis=1) if np.any(~keep) else np.zeros(jac.shape[1])
        cov_diag = np.where(kernel_share >= 0.5, np.inf, cov_diag)
    cr = np.sqrt(cov_diag)
    # Right singular vector of the smallest component, for the null check.
    u, svals, vt = np.linalg.svd(jac, full_matrices=False)
    null_vec = vt[-1].tolist() if vt.size else []
    return {
        "n_obs": int(jac.shape[0]),
        "n_param": int(jac.shape[1]),
        "singular_values": singular.tolist(),
        "singular_ratios": ratio.tolist(),
        "structural_rank": structural,
        "practical_rank": practical,
        "principal_rse": principal,
        "cramer_rao_cv": [None if not np.isfinite(v) else float(v) for v in cr],
        "smallest_right_vector": null_vec,
        "condition_structural": (
            float(singular[0] / singular[structural - 1]) if structural else None
        ),
    }


def theoretical_null_vector(lam, kap) -> np.ndarray:
    """One null vector per edge in log-coordinates: (λ, −κ) on that pair, stacked.

    Returned as three columns. Node outputs are invariant along each.
    """
    n = lam.size
    cols = []
    for e in range(n):
        v = np.zeros(2 * n)
        v[e] = lam[e]
        v[n + e] = -kap[e]
        v /= np.linalg.norm(v)
        cols.append(v)
    return np.column_stack(cols)


def null_alignment(report: dict, lam, kap) -> dict:
    """Cosines between the three trailing right singular vectors and the
    theoretical per-edge kernel. Meaningful when structural rank is n_edges.
    """
    # Recompute vectors from stored singular path: we need V. Caller passes
    # the vector of the smallest only; full alignment is computed outside
    # if jac is available. This helper fills a placeholder; see analyze().
    return report


def predict_series_theta(theta, schedule, r=R):
    lam, kap = theta[:3], theta[3:]
    if np.any(lam <= 0) or np.any(kap <= 0):
        raise RuntimeError("nonpositive parameter")
    nodes = simulate_series(lam, kap, r=r)
    return output_series(nodes, lam, kap, schedule)


def analyze_series(schedule: str, r=R, theta=None) -> dict:
    theta = LAM_TRUE if theta is None else theta
    # default theta is the 6-vector
    if theta is None or len(np.asarray(theta)) != 6:
        theta = np.concatenate([LAM_TRUE, KAP_TRUE])
    sigma = sigma_series(schedule)

    def predict(th):
        return predict_series_theta(th, schedule, r=r)

    jac = weighted_sensitivity(predict, np.asarray(theta, float), sigma)
    report = fisher_report(jac)
    report["schedule"] = schedule
    # Kernel alignment against the three theoretical directions.
    _, _, vt = np.linalg.svd(jac, full_matrices=False)
    theory = theoretical_null_vector(LAM_TRUE, KAP_TRUE)
    # Trailing 3 right singular vectors (rows of vt).
    trailing = vt[-3:]
    cosines = np.abs(trailing @ theory)
    # Best matching cosine per theoretical direction.
    report["kernel_cosine_max"] = cosines.max(axis=0).tolist()
    report["jacobian_rank_gap"] = {
        "largest": report["singular_values"][0],
        "fourth": report["singular_values"][3] if len(report["singular_values"]) > 3 else None,
        "smallest": report["singular_values"][-1],
    }
    return report


def analyze_single_rate(schedule: str) -> dict:
    """Fisher on log φ for the single-rate node model. Schedules L and S only."""
    phi = series_phi(LAM_TRUE, KAP_TRUE)
    sigma = sigma_series(schedule)

    def predict(ph):
        nodes = simulate_phi(ph)
        # Reuse the node part of output_series by passing a dummy pair whose
        # stalled fraction is never read for L and S.
        return output_series(nodes, LAM_TRUE, KAP_TRUE, schedule)

    if schedule not in ("L", "S"):
        raise KeyError(schedule)
    jac = weighted_sensitivity(predict, phi, sigma)
    report = fisher_report(jac)
    report["schedule"] = schedule
    report["phi_truth"] = phi.tolist()
    return report


def analyze_holdup(schedule: str) -> dict:
    theta = np.concatenate([LAM_TRUE, KAP_TRUE])
    sigma = sigma_holdup(schedule)

    def predict(th):
        nodes, comps = simulate_holdup(th[:3], th[3:])
        return output_holdup(nodes, comps, schedule)

    jac = weighted_sensitivity(predict, theta, sigma)
    report = fisher_report(jac)
    report["schedule"] = schedule
    return report


def analyze_soil_free(schedule: str) -> dict:
    """θ = (λ×3, κ×3, r×4). Soil rates freed; carrying capacities stay known."""
    theta = np.concatenate([LAM_TRUE, KAP_TRUE, R])
    sigma = sigma_series(schedule)

    def predict(th):
        return predict_series_theta(th[:6], schedule, r=th[6:])

    jac = weighted_sensitivity(predict, theta, sigma)
    report = fisher_report(jac)
    report["schedule"] = schedule
    report["names"] = [f"λ{e+1}" for e in range(3)] + [f"κ{e+1}" for e in range(3)] + [
        f"r{n}" for n in NODE_NAMES
    ]
    return report


def hyperbola_twin() -> dict:
    """Move every edge along φ = constant. λ' = 1.8 λ, κ' the partner."""
    phi = series_phi(LAM_TRUE, KAP_TRUE)
    lam2 = 1.8 * LAM_TRUE
    if np.any(lam2 <= phi):
        raise RuntimeError("twin left the φ < λ branch")
    kap2 = partner_kappa(phi, lam2)
    nodes = simulate_series(LAM_TRUE, KAP_TRUE)
    nodes2 = simulate_series(lam2, kap2)
    s = stalled_fraction(LAM_TRUE, KAP_TRUE)
    s2 = stalled_fraction(lam2, kap2)
    return {
        "lambda_twin": lam2.tolist(),
        "kappa_twin": kap2.tolist(),
        "phi": phi.tolist(),
        "stalled_truth": s.tolist(),
        "stalled_twin": s2.tolist(),
        "node_rmse": float(np.sqrt(np.mean((nodes - nodes2) ** 2))),
        "node_max_abs": float(np.max(np.abs(nodes - nodes2))),
        "lump_rmse": float(np.sqrt(np.mean((nodes.sum(0) - nodes2.sum(0)) ** 2))),
        "stalled_rmse": float(np.sqrt(np.mean((s - s2) ** 2))),
    }


def level_set_edge2() -> dict:
    """Cost along the edge-2 hyperbola, other edges held at the truth."""
    phi = series_phi(LAM_TRUE, KAP_TRUE)
    phi2 = float(phi[1])
    # Stay strictly above φ so κ is finite, and include one point below φ.
    multiples = np.array([0.5, 0.8, 1.0, 1.15, 1.4, 1.8, 2.5, 3.5, 5.0])
    # multiples are of λ_true, not of φ. Record λ itself.
    rows = []
    nodes_ref = simulate_series(LAM_TRUE, KAP_TRUE)
    y_ref = {
        "L": output_series(nodes_ref, LAM_TRUE, KAP_TRUE, "L"),
        "S": output_series(nodes_ref, LAM_TRUE, KAP_TRUE, "S"),
        "SB": output_series(nodes_ref, LAM_TRUE, KAP_TRUE, "SB"),
    }
    sig = {sch: sigma_series(sch) for sch in y_ref}
    for m in multiples:
        lam = LAM_TRUE.copy()
        kap = KAP_TRUE.copy()
        lam[1] = float(m) * LAM_TRUE[1]
        feasible = lam[1] > phi2
        if feasible:
            kap[1] = float(partner_kappa(np.array([phi2]), np.array([lam[1]]))[0])
        else:
            # No series partner. Park κ at the truth so the point is simply
            # off the level set, and say so.
            kap[1] = KAP_TRUE[1]
        nodes = simulate_series(lam, kap)
        entry = {
            "lambda2": float(lam[1]),
            "kappa2": float(kap[1]),
            "on_level_set": bool(feasible and abs(m - 0.5) > 0),
            "feasible_partner": bool(feasible),
        }
        # 0.5 and 0.8 of λ_true=0.06 are 0.03 and 0.048. φ2 ≈ 0.0176, so
        # both are still above φ. Use an explicit infeasible marker below.
        for sch in y_ref:
            y = output_series(nodes, lam, kap, sch)
            wr = (y - y_ref[sch]) / sig[sch]
            entry[f"weighted_rss_{sch}"] = float(np.dot(wr, wr))
        rows.append(entry)
    # Explicit infeasible point: λ2 = 0.5 φ2.
    lam = LAM_TRUE.copy()
    kap = KAP_TRUE.copy()
    lam[1] = 0.5 * phi2
    nodes = simulate_series(lam, kap)
    entry = {
        "lambda2": float(lam[1]),
        "kappa2": float(kap[1]),
        "on_level_set": False,
        "feasible_partner": False,
    }
    for sch in y_ref:
        y = output_series(nodes, lam, kap, sch)
        wr = (y - y_ref[sch]) / sig[sch]
        entry[f"weighted_rss_{sch}"] = float(np.dot(wr, wr))
    rows.append(entry)
    return {"phi2": phi2, "points": rows}


def _fit_free(schedule: str, theta0: np.ndarray, free_mask: np.ndarray, r=R) -> dict:
    """Least squares on the free log-coordinates. Fixed entries stay at theta0."""
    theta0 = np.asarray(theta0, float)
    sigma = sigma_series(schedule)
    y_obs = predict_series_theta(np.concatenate([LAM_TRUE, KAP_TRUE]), schedule, r=R)
    idx = np.where(free_mask)[0]

    def residual(u_free):
        theta = theta0.copy()
        theta[idx] = np.exp(u_free)
        try:
            y = predict_series_theta(theta, schedule, r=r)
        except Exception:
            return np.full(y_obs.shape, 1e3)
        return (y - y_obs) / sigma

    theta0 = np.clip(theta0, np.exp(LOG_LO), np.exp(LOG_HI))
    u0 = np.log(theta0[idx])
    sol = least_squares(
        residual,
        u0,
        bounds=(np.full(idx.size, LOG_LO), np.full(idx.size, LOG_HI)),
        method="trf",
        ftol=1e-12,
        xtol=1e-12,
        gtol=1e-12,
        max_nfev=120,
    )
    theta_hat = theta0.copy()
    theta_hat[idx] = np.exp(sol.x)
    return {
        "success": bool(sol.success),
        "nfev": int(sol.nfev),
        "weighted_rss": float(np.dot(sol.fun, sol.fun)),
        "theta_hat": theta_hat.tolist(),
        "lambda_hat": theta_hat[:3].tolist(),
        "kappa_hat": theta_hat[3:].tolist(),
    }


def profiles() -> dict:
    """Profile λ2 and κ2 under L, S, and SB. Other coordinates re-fit."""
    truth = np.concatenate([LAM_TRUE, KAP_TRUE])
    out = {}
    targets = {"lambda2": 1, "kappa2": 4}
    for name, index in targets.items():
        out[name] = {}
        grid = PROFILE_MULT * truth[index]
        for schedule in ("L", "S", "SB"):
            rows = []
            free = np.ones(6, dtype=bool)
            free[index] = False
            for value in grid:
                theta0 = truth.copy()
                theta0[index] = float(value)
                # Warm start: if the profiled rate is λ2 and a partner exists,
                # put κ2 on the hyperbola so node-only fits start at cost 0.
                phi = series_phi(LAM_TRUE, KAP_TRUE)
                if index == 1 and value > phi[1]:
                    theta0[4] = float(partner_kappa(np.array([phi[1]]), np.array([value]))[0])
                if index == 4 and value > phi[1]:
                    # κ2 profile: partner λ2 = φ κ / (κ − φ) when κ > φ.
                    if value > phi[1]:
                        theta0[1] = float(phi[1] * value / (value - phi[1]))
                fit = _fit_free(schedule, theta0, free)
                rows.append(
                    {
                        "value": float(value),
                        "multiple": float(value / truth[index]),
                        "weighted_rss": fit["weighted_rss"],
                        "inside_95": bool(fit["weighted_rss"] <= CHI2_95),
                        "lambda_hat": fit["lambda_hat"],
                        "kappa_hat": fit["kappa_hat"],
                    }
                )
            out[name][schedule] = rows
    return out


def refits(twin: dict) -> dict:
    truth = np.concatenate([LAM_TRUE, KAP_TRUE])
    twin_theta = np.concatenate([np.array(twin["lambda_twin"]), np.array(twin["kappa_twin"])])
    off = truth * np.array([2.0, 0.40, 1.50, 0.50, 3.0, 0.70])
    free = np.ones(6, dtype=bool)
    report = {}
    for label, start in (("hyperbola_twin", twin_theta), ("off_level", off)):
        report[label] = {"start": start.tolist()}
        for schedule in ("L", "S", "SB"):
            fit = _fit_free(schedule, start, free)
            hat = np.array(fit["theta_hat"])
            rel = np.abs(np.log(hat / truth))
            report[label][schedule] = {
                **fit,
                "max_abs_log_error": float(rel.max()),
                "mean_abs_log_error": float(rel.mean()),
            }
    return report


def noisy_multistart() -> dict:
    rng = np.random.default_rng(SEED)
    truth = np.concatenate([LAM_TRUE, KAP_TRUE])
    phi = series_phi(LAM_TRUE, KAP_TRUE)
    lam_t = 1.8 * LAM_TRUE
    kap_t = partner_kappa(phi, lam_t)
    starts = np.ones((8, 6))
    starts[1] = (1.8, 0.5, 2.0, 0.4, 2.5, 0.6)
    starts[2, :3] = lam_t / LAM_TRUE
    starts[2, 3:] = kap_t / KAP_TRUE
    starts[3] = (0.5, 2.0, 0.5, 2.0, 0.4, 1.5)
    starts[4] = (3.0, 3.0, 3.0, 0.2, 0.2, 0.2)
    starts[5] = (0.3, 0.3, 0.3, 3.0, 3.0, 3.0)
    starts[6] = (2.0, 0.8, 1.2, 1.5, 0.3, 2.0)
    starts[7] = (0.7, 1.4, 0.6, 0.8, 1.6, 0.9)
    out = {"seed": SEED, "n_starts": int(starts.shape[0]), "schedules": {}}
    for schedule in ("L", "SB"):
        sigma = sigma_series(schedule)
        y_clean = predict_series_theta(truth, schedule)
        noise = rng.normal(0.0, 1.0, size=y_clean.shape)
        y_obs = y_clean + sigma * noise
        rows = []
        for m in starts:
            theta0 = truth * m

            def residual(u, theta0=theta0):
                theta = np.exp(u)
                try:
                    y = predict_series_theta(theta, schedule)
                except Exception:
                    return np.full(y_obs.shape, 1e3)
                return (y - y_obs) / sigma

            sol = least_squares(
                residual,
                np.log(theta0),
                bounds=(np.full(6, LOG_LO), np.full(6, LOG_HI)),
                method="trf",
                ftol=1e-12,
                xtol=1e-12,
                gtol=1e-12,
                max_nfev=100,
            )
            hat = np.exp(sol.x)
            rows.append(
                {
                    "start_multiple": m.tolist(),
                    "weighted_rss": float(np.dot(sol.fun, sol.fun)),
                    "lambda_hat": hat[:3].tolist(),
                    "kappa_hat": hat[3:].tolist(),
                    "max_abs_log_error": float(np.max(np.abs(np.log(hat / truth)))),
                }
            )
        rss = np.array([row["weighted_rss"] for row in rows])
        # Retain fits within χ²_6,0.95 of the best. 12.59.
        window = 12.591587243743977
        retain = rss <= rss.min() + window
        lam = np.array([row["lambda_hat"] for row, keep in zip(rows, retain) if keep])
        kap = np.array([row["kappa_hat"] for row, keep in zip(rows, retain) if keep])
        out["schedules"][schedule] = {
            "noise_rss_of_truth": float(np.dot(noise, noise)),
            "best_weighted_rss": float(rss.min()),
            "n_retained": int(retain.sum()),
            "lambda_min": lam.min(axis=0).tolist(),
            "lambda_max": lam.max(axis=0).tolist(),
            "kappa_min": kap.min(axis=0).tolist(),
            "kappa_max": kap.max(axis=0).tolist(),
            "fits": rows,
        }
    return out


def _fmt(x, digits=4) -> str:
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "null"
    return f"{x:.{digits}g}"


def make_figures(nodes, twin, level, profiles_data, series_reports, single_reports) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "font.size": 10,
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "figure.dpi": 160,
            "savefig.dpi": 160,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    colors = {
        "P": "#1b4f72",
        "F": "#117a65",
        "V": "#b03a2e",
        "B": "#7d6608",
    }

    # Figure 1. Node trajectories. The hyperbola twin coincides.
    lam2 = np.array(twin["lambda_twin"])
    kap2 = np.array(twin["kappa_twin"])
    nodes2 = simulate_series(lam2, kap2)
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    for i, name in enumerate(NODE_NAMES):
        ax.plot(T_OBS, nodes[i], color=colors[name], lw=1.8, label=name)
        ax.plot(T_OBS, nodes2[i], color=colors[name], lw=1.0, ls="--", alpha=0.9)
    ax.plot(T_OBS, nodes.sum(0), color="#1c2833", lw=1.6, label="lump")
    ax.plot(T_OBS, nodes2.sum(0), color="#1c2833", lw=1.0, ls="--")
    ax.set_xlabel("Toy time")
    ax.set_ylabel("Burden (toy units)")
    ax.set_title("Node burdens and the lumped sum")
    ax.legend(frameon=False, ncol=5)
    fig.tight_layout()
    fig.savefig(FIG / "fig_trajectories.png")
    plt.close(fig)

    # Figure 2. Level-set costs.
    pts = [p for p in level["points"] if p["feasible_partner"]]
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    lam = [p["lambda2"] for p in pts]
    for sch, style in (("L", "-o"), ("S", "-s"), ("SB", "-D")):
        ax.plot(
            lam,
            [max(p[f"weighted_rss_{sch}"], 1e-18) for p in pts],
            style,
            lw=1.5,
            ms=4.5,
            label=sch,
        )
    ax.axhline(CHI2_95, color="#7b241c", lw=0.8, ls=":", label="χ² 95%, 1 df")
    ax.axvline(LAM_TRUE[1], color="#666", lw=0.7, ls="--")
    ax.set_yscale("log")
    ax.set_xlabel("λ of F→V, κ chosen to hold φ")
    ax.set_ylabel("Weighted residual sum of squares")
    ax.set_title("Edge F→V along the series level set")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "fig_level_set.png")
    plt.close(fig)

    # Figure 3. Spectra.
    order = ["L", "S", "LP", "LB", "SP", "SB"]
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    for sch in order:
        ratios = np.array(series_reports[sch]["singular_ratios"], float)
        ax.plot(np.arange(1, 7), np.maximum(ratios, 1e-16), marker="o", ms=4, lw=1.4, label=sch)
    ax.axhline(STRUCT_RATIO, color="#7b241c", lw=0.8, ls=":", label="structural floor")
    ax.set_yscale("log")
    ax.set_xticks(range(1, 7))
    ax.set_xlabel("Component")
    ax.set_ylabel("Singular value / largest")
    ax.set_title("Weighted sensitivity spectra, series model")
    ax.legend(frameon=False, ncol=3)
    fig.tight_layout()
    fig.savefig(FIG / "fig_spectra.png")
    plt.close(fig)

    # Figure 4. Profiles of λ2 and κ2.
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 4.0), sharey=True)
    for ax, key, truth_v, xlabel in (
        (axes[0], "lambda2", LAM_TRUE[1], "λ of F→V"),
        (axes[1], "kappa2", KAP_TRUE[1], "κ of F→V"),
    ):
        for sch, marker in (("L", "o"), ("S", "s"), ("SB", "D")):
            rows = profiles_data[key][sch]
            ax.plot(
                [row["value"] for row in rows],
                [max(row["weighted_rss"], 1e-18) for row in rows],
                marker=marker,
                ms=3.5,
                lw=1.3,
                label=sch,
            )
        ax.axhline(CHI2_95, color="#7b241c", lw=0.8, ls=":")
        ax.axvline(truth_v, color="#666", lw=0.7, ls="--")
        ax.set_yscale("log")
        ax.set_xlabel(xlabel)
        ax.set_title("Profile, " + key)
    axes[0].set_ylabel("Profiled weighted RSS")
    axes[1].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "fig_profiles.png")
    plt.close(fig)

    # Figure 5. Finite principal relative SE only. Structural nulls are omitted
    # rather than drawn as a capped finite error.
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    labels = []
    heights = []
    bar_colors = []
    palette = {"L": "#1b4f72", "S": "#117a65", "SB": "#b03a2e"}
    for sch in ("L", "S"):
        rse = single_reports[sch]["principal_rse"]
        for i, val in enumerate(rse, start=1):
            if val is None:
                continue
            labels.append(f"φ {sch}:{i}")
            heights.append(val)
            bar_colors.append(palette[sch])
    for sch in ("L", "S", "SB"):
        rse = series_reports[sch]["principal_rse"]
        for i, val in enumerate(rse, start=1):
            if val is None:
                continue
            labels.append(f"joint {sch}:{i}")
            heights.append(val)
            bar_colors.append(palette.get(sch, "#7d6608"))
    xpos = np.arange(len(heights))
    ax.bar(xpos, heights, color=bar_colors, width=0.8)
    ax.axhline(PRACTICAL_RSE, color="#7b241c", lw=0.8, ls="--")
    ax.set_xticks(xpos)
    ax.set_xticklabels(labels, rotation=90, fontsize=7)
    ax.set_ylabel("Principal relative SE")
    ax.set_title("Finite principal errors (structural nulls omitted)")
    fig.tight_layout()
    fig.savefig(FIG / "fig_principal_se.png")
    plt.close(fig)


def _jsonable(obj):
    if isinstance(obj, dict):
        return {k: _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, float) and not np.isfinite(obj):
        return None
    return obj


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    phi = series_phi(LAM_TRUE, KAP_TRUE)
    stalled = stalled_fraction(LAM_TRUE, KAP_TRUE)
    nodes = simulate_series(LAM_TRUE, KAP_TRUE)
    print("phi", phi)
    print("stalled", stalled)
    print("final nodes", nodes[:, -1], "lump", nodes[:, -1].sum())

    schedules = ("L", "S", "LP", "LB", "SP", "SB")
    series_reports = {}
    for sch in schedules:
        print("series", sch)
        series_reports[sch] = analyze_series(sch)
        rep = series_reports[sch]
        print(
            " ",
            rep["structural_rank"],
            "/",
            rep["n_param"],
            "practical",
            rep["practical_rank"],
            "rse",
            rep["principal_rse"],
            "cos",
            rep["kernel_cosine_max"],
        )

    print("single-rate")
    single_reports = {}
    for sch in ("L", "S"):
        single_reports[sch] = analyze_single_rate(sch)
        rep = single_reports[sch]
        print(" ", sch, rep["structural_rank"], rep["practical_rank"], rep["principal_rse"])

    print("hold-up")
    hold_reports = {}
    for sch in ("L", "S", "LB", "SB"):
        hold_reports[sch] = analyze_holdup(sch)
        rep = hold_reports[sch]
        print(" ", sch, rep["structural_rank"], rep["practical_rank"], rep["principal_rse"])

    print("soil-free")
    soil_reports = {}
    for sch in ("S", "SB"):
        soil_reports[sch] = analyze_soil_free(sch)
        rep = soil_reports[sch]
        print(" ", sch, rep["structural_rank"], "/", rep["n_param"], "practical", rep["practical_rank"])

    print("hyperbola")
    twin = hyperbola_twin()
    print(" node rmse", twin["node_rmse"], "stalled rmse", twin["stalled_rmse"])

    print("level set")
    level = level_set_edge2()

    print("profiles")
    prof = profiles()
    for key in prof:
        for sch in prof[key]:
            inside = sum(row["inside_95"] for row in prof[key][sch])
            rss = [row["weighted_rss"] for row in prof[key][sch]]
            print(f"  {key} {sch} inside {inside}/{len(rss)} rss {_fmt(min(rss))}..{_fmt(max(rss))}")

    print("refits")
    fits = refits(twin)
    for label in fits:
        for sch in ("L", "S", "SB"):
            row = fits[label][sch]
            print(f"  {label} {sch} rss {row['weighted_rss']:.3e} maxlog {row['max_abs_log_error']:.3e}")

    print("multistart")
    multi = noisy_multistart()
    for sch, row in multi["schedules"].items():
        print(
            " ",
            sch,
            "retained",
            row["n_retained"],
            "best",
            row["best_weighted_rss"],
            "λ span",
            row["lambda_min"],
            row["lambda_max"],
        )

    make_figures(nodes, twin, level, prof, series_reports, single_reports)

    # Trajectory snapshot for the manuscript.
    snapshot_t = [0, 8, 16, 24, 36]
    snap = {}
    for t in snapshot_t:
        j = int(np.where(np.isclose(T_OBS, t))[0][0])
        snap[str(t)] = {
            "nodes": nodes[:, j].tolist(),
            "lump": float(nodes[:, j].sum()),
        }

    payload = {
        "seed": SEED,
        "sigma_y": SIGMA_Y,
        "sigma_s": SIGMA_S,
        "sigma_c": SIGMA_C,
        "struct_ratio": STRUCT_RATIO,
        "practical_rse_cutoff": PRACTICAL_RSE,
        "chi2_95_1df": CHI2_95,
        "times": T_OBS.tolist(),
        "r": R.tolist(),
        "K": K.tolist(),
        "x0": X0.tolist(),
        "lambda_true": LAM_TRUE.tolist(),
        "kappa_true": KAP_TRUE.tolist(),
        "phi_true": phi.tolist(),
        "stalled_true": stalled.tolist(),
        "trajectory_snapshot": snap,
        "final_nodes": nodes[:, -1].tolist(),
        "final_lump": float(nodes[:, -1].sum()),
        "series": series_reports,
        "single_rate": single_reports,
        "holdup": hold_reports,
        "soil_free": soil_reports,
        "hyperbola_twin": twin,
        "level_set_edge2": level,
        "profiles": prof,
        "refits": fits,
        "noisy_multistart": multi,
    }
    text = json.dumps(_jsonable(payload), indent=2)
    (ROOT / "results.json").write_text(text + "\n", encoding="utf-8")
    print("wrote", ROOT / "results.json")


if __name__ == "__main__":
    main()
