"""Numerical routines for the frozen CRC16 within-section analysis."""
import numpy as np
from scipy.special import expit
from scipy.stats import chi2, norm, t


def bh(values, planned=None):
    """BH on finite values; optional P=1 padding preserves the planned family.

    Padded values are never returned as observed P/q for unavailable tests.
    """
    p = np.asarray(values, float)
    finite = np.isfinite(p)
    available = p[finite]
    m = int(planned) if planned is not None else len(available)
    if m < len(available):
        raise ValueError("Planned family smaller than finite test count")
    out = np.full(p.shape, np.nan)
    if not len(available):
        return out
    padded = np.r_[available, np.ones(m - len(available))]
    order = np.argsort(padded, kind="stable")
    adjusted = np.minimum.accumulate((padded[order] * m / np.arange(1, m + 1))[::-1])[::-1]
    restored = np.empty(m)
    restored[order] = np.minimum(adjusted, 1)
    out[finite] = restored[:len(available)]
    return out


def logistic(y, design, niche_index=1, tolerance=1e-10, max_iterations=100):
    """Unpenalized Newton fit, observed-information two-sided Wald test.

    No event/non-event in either niche comparison group means that the binary
    niche indicator separates the outcome. Such estimates remain unavailable.
    """
    y = np.asarray(y, float)
    x = np.asarray(design, float)
    inside = x[:, niche_index].astype(bool)
    for mask in (inside, ~inside):
        selected = y[mask]
        if not len(selected) or selected.min() == selected.max():
            return np.nan, np.nan, np.nan, "Non-estimable: separation by niche indicator"
    beta = np.zeros(x.shape[1])
    for iteration in range(max_iterations):
        mu = expit(x @ beta)
        information = x.T @ (x * (mu * (1 - mu))[:, None])
        try:
            update = np.linalg.solve(information, x.T @ (y - mu))
        except np.linalg.LinAlgError:
            return np.nan, np.nan, np.nan, "Non-estimable: singular information"
        beta += update
        if not np.isfinite(beta).all():
            return np.nan, np.nan, np.nan, "Non-estimable: nonfinite fit"
        if np.max(np.abs(update)) < tolerance:
            break
    else:
        return np.nan, np.nan, np.nan, "Non-estimable: convergence not achieved"
    mu = expit(x @ beta)
    information = x.T @ (x * (mu * (1 - mu))[:, None])
    variance = np.linalg.inv(information)[niche_index, niche_index]
    b = float(beta[niche_index])
    probability = float(2 * norm.sf(abs(b / np.sqrt(variance))))
    return b, float(np.exp(b)), probability, "Estimable"


def ols(y, design, niche_index=1):
    x, y = np.asarray(design, float), np.asarray(y, float)
    beta = np.linalg.lstsq(x, y, rcond=None)[0]
    residual = y - x @ beta
    degrees_freedom = len(y) - np.linalg.matrix_rank(x)
    variance = residual @ residual / degrees_freedom
    covariance = np.linalg.inv(x.T @ x) * variance
    b = beta[niche_index]
    p = 2 * t.sf(abs(b / np.sqrt(covariance[niche_index, niche_index])), degrees_freedom)
    return float(b), float(p)


def cmh(inside, positive, blocks):
    """CMH without continuity correction/pseudocount, including boundary ORs."""
    rows, differences, variances, numerator, denominator = [], [], [], 0., 0.
    for block in np.unique(blocks):
        mask = blocks == block
        a = int(np.sum(inside & positive & mask))
        b = int(np.sum(inside & ~positive & mask))
        c = int(np.sum(~inside & positive & mask))
        d = int(np.sum(~inside & ~positive & mask))
        n = a + b + c + d
        expected = (a + b) * (a + c) / n if n else 0.
        variance = (a + b) * (c + d) * (a + c) * (b + d) / (n * n * (n - 1)) if n > 1 else 0.
        differences.append(a - expected)
        variances.append(variance)
        numerator += a * d / n if n else 0.
        denominator += b * c / n if n else 0.
        rows.append({"Block": int(block), "a_niche_positive": a, "b_niche_negative": b,
                     "c_rest_positive": c, "d_rest_negative": d, "n": n,
                     "expected_a": expected, "conditional_variance": variance})
    total_variance = sum(variances)
    statistic = sum(differences) ** 2 / total_variance if total_variance > 0 else np.nan
    probability = float(chi2.sf(statistic, 1)) if total_variance > 0 else np.nan
    if denominator > 0:
        odds = numerator / denominator
        status = "finite positive" if odds > 0 else "boundary zero"
    elif numerator > 0:
        odds, status = np.inf, "boundary infinite"
    else:
        odds, status = np.nan, "undefined"
    return odds, probability, statistic, status, rows
