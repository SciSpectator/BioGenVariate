import math
import numpy as np
import pandas as pd
import torch
from scipy.stats import gaussian_kde, norm
from scipy.signal import find_peaks

# Set the device for torch operations, defaulting to CPU if CUDA is not available.
current_device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

def to_numeric_array(vals):
    """Converts input values to a numpy float64 array for numerical operations."""
    return np.array(vals, dtype=np.float64)

def is_bimodal(values, min_prominence_ratio=0.05):
    """
    Checks if a distribution is bimodal using kernel density estimation (KDE) and peak finding.
    Returns True if at least two prominent peaks are found.
    """
    arr = to_numeric_array(values)
    arr = arr[~np.isnan(arr)]
    if arr.size < 2 or arr.min() == arr.max():
        return False
    try:
        kde = gaussian_kde(arr)
        grid = np.linspace(arr.min(), arr.max(), 200)
        pdf = kde(grid)
        peaks, _ = find_peaks(pdf, prominence=min_prominence_ratio * pdf.max())
        return len(peaks) >= 2
    except Exception:
        return False

# ------------------------------------------------------------------------------
# Functions for fitting distributions and calculating log-likelihoods
# ------------------------------------------------------------------------------

def fit_normal(x_t):
    """Fits a normal distribution to the data tensor."""
    mu = x_t.mean()
    sigma = x_t.std(ddof=1)
    if sigma.item() <= 1e-9:
        sigma = torch.tensor(1e-9, device=current_device)
    return mu, sigma

def loglike_normal(x_t, mu, sigma):
    """Calculates the log-likelihood of data given a normal distribution."""
    n = x_t.numel()
    if sigma.item() <= 1e-9:
        return torch.tensor(float('-inf'), device=current_device)
    return -0.5 * n * torch.log(torch.tensor(2 * math.pi, device=current_device)) \
           - n * torch.log(sigma) \
           - 0.5 * (((x_t - mu) / sigma) ** 2).sum()

def fit_cauchy(x_t):
    """Fits a Cauchy distribution to the data tensor."""
    x0 = x_t.median()
    q1 = x_t.quantile(0.25)
    q3 = x_t.quantile(0.75)
    gamma_ = (q3 - q1) / 2
    if gamma_.item() <= 1e-9:
        gamma_ = torch.tensor(1e-9, device=current_device)
    return x0, gamma_

def loglike_cauchy(x_t, x0, gamma_):
    """Calculates the log-likelihood of data given a Cauchy distribution."""
    n = x_t.numel()
    if gamma_.item() <= 1e-9:
        return torch.tensor(float('-inf'), device=current_device)
    z = (x_t - x0) / gamma_
    return -n * torch.log(gamma_) \
           - n * torch.log(torch.tensor(math.pi, device=current_device)) \
           - torch.log(1 + z * z).sum()

def fit_lognormal(x_t):
    """Fits a log-normal distribution to the data tensor."""
    xp = x_t[x_t > 0]
    if xp.numel() < 2:
        return torch.tensor(0., device=current_device), torch.tensor(1e-9, device=current_device)
    ln_x = torch.log(xp)
    mu = ln_x.mean()
    sigma = ln_x.std(ddof=1)
    if sigma.item() <= 1e-9:
        sigma = torch.tensor(1e-9, device=current_device)
    return mu, sigma

def loglike_lognormal(x_t, mu, sigma):
    """Calculates the log-likelihood of data given a log-normal distribution."""
    xp = x_t[x_t > 0]
    if xp.numel() == 0 or sigma.item() <= 1e-9:
        return torch.tensor(float('-inf'), device=current_device)
    n = xp.numel()
    ln_x = torch.log(xp)
    return -0.5 * n * torch.log(torch.tensor(2 * math.pi, device=current_device)) \
           - n * torch.log(sigma) \
           - ln_x.sum() \
           - 0.5 * (((ln_x - mu) / sigma) ** 2).sum()

def fit_gamma(x_t):
    """Fits a Gamma distribution to the data tensor."""
    xp = x_t[x_t > 0]
    if xp.numel() < 2:
        return torch.tensor(1., device=current_device), torch.tensor(1., device=current_device)
    mean_x = xp.mean()
    var_x = xp.var(ddof=1)
    if mean_x.item() <= 1e-9 or var_x.item() <= 1e-9:
        return torch.tensor(1., device=current_device), torch.tensor(1., device=current_device)
    k = mean_x * mean_x / var_x
    theta = var_x / mean_x
    if k.item() <= 1e-9:
        k = torch.tensor(1e-9, device=current_device)
    if theta.item() <= 1e-9:
        theta = torch.tensor(1e-9, device=current_device)
    return k, theta

def loglike_gamma(x_t, k, theta):
    """Calculates the log-likelihood of data given a Gamma distribution."""
    xp = x_t[x_t > 0]
    if xp.numel() == 0 or k.item() <= 1e-9 or theta.item() <= 1e-9:
        return torch.tensor(float('-inf'), device=current_device)
    n = xp.numel()
    return (k - 1) * torch.log(xp).sum() \
           - (xp / theta).sum() \
           - n * (torch.lgamma(k) + k * torch.log(theta))

def fit_pareto(x_t):
    """Fits a Pareto distribution to the data tensor."""
    xp = x_t[x_t > 0]
    if xp.numel() < 2:
        return torch.tensor(1e-9, device=current_device), torch.tensor(1., device=current_device)
    xm = xp.min()
    if xm.item() <= 1e-9:
        return torch.tensor(1e-9, device=current_device), torch.tensor(1., device=current_device)
    denom = torch.log(xp / xm).sum()
    alpha = xp.numel() / denom if denom.item() > 1e-9 else torch.tensor(1., device=current_device)
    if alpha.item() <= 1e-9:
        alpha = torch.tensor(1., device=current_device)
    return xm, alpha

def loglike_pareto(x_t, xm, alpha):
    """Calculates the log-likelihood of data given a Pareto distribution."""
    if xm.item() <= 1e-9 or alpha.item() <= 1e-9:
        return torch.tensor(float('-inf'), device=current_device)
    xp = x_t[x_t >= xm]
    if xp.numel() == 0:
        return torch.tensor(float('-inf'), device=current_device)
    n = xp.numel()
    return n * torch.log(alpha) \
           + n * alpha * torch.log(xm) \
           - (alpha + 1) * torch.log(xp).sum()

def analyze_gene_distribution(series: pd.Series) -> str:
    """
    Analyzes a gene expression distribution to find the best-fitting statistical model.
    Checks for bimodality first, then compares several common distributions.
    """
    vals = series.dropna().astype(np.float32).values
    if vals.size == 0:
        return "No Data"
    if vals.size > 1 and (vals == vals[0]).all():
        return "Constant"
    if is_bimodal(vals):
        return "Bimodal"
    
    x_t = torch.tensor(vals, device=current_device)
    if x_t.numel() < 2:
        return "Insufficient Data"
    
    scores = {}
    distributions_to_test = [
        ("Normal", fit_normal, loglike_normal),
        ("Cauchy", fit_cauchy, loglike_cauchy),
        ("Lognormal", fit_lognormal, loglike_lognormal),
        ("Gamma", fit_gamma, loglike_gamma),
        ("Pareto", fit_pareto, loglike_pareto),
    ]

    for name, fit_fn, ll_fn in distributions_to_test:
        try:
            params = fit_fn(x_t)
            scores[name] = float(ll_fn(x_t, *params))
        except Exception:
            scores[name] = float('-inf')
    
    valid_scores = {k: v for k, v in scores.items() if v != float('-inf') and not np.isnan(v)}
    
    return max(valid_scores, key=valid_scores.get) if valid_scores else "Undetermined"
