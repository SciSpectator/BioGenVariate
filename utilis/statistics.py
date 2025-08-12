# utils/statistics.py
import math
import numpy as np
import pandas as pd
import torch
from scipy.stats import gaussian_kde, norm
from scipy.signal import find_peaks

current_device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

def to_numeric_array(vals):
    return np.array(vals, dtype=np.float64)

def is_bimodal(values, min_prominence_ratio=0.05):
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
    except:
        return False

def fit_normal(x_t):
    # ... (code for fit_normal)
    pass

def loglike_normal(x_t, mu, sigma):
    # ... (code for loglike_normal)
    pass

def fit_cauchy(x_t):
    # ... (code for fit_cauchy)
    pass

def loglike_cauchy(x_t, x0, gamma_):
    # ... (code for loglike_cauchy)
    pass

def fit_lognormal(x_t):
    # ... (code for fit_lognormal)
    pass

def loglike_lognormal(x_t, mu, sigma):
    # ... (code for loglike_lognormal)
    pass

def fit_gamma(x_t):
    # ... (code for fit_gamma)
    pass

def loglike_gamma(x_t, k, theta):
    # ... (code for loglike_gamma)
    pass

def fit_pareto(x_t):
    # ... (code for fit_pareto)
    pass

def loglike_pareto(x_t, xm, alpha):
    # ... (code for loglike_pareto)
    pass

def analyze_gene_distribution(series: pd.Series) -> str:
    # ... (code for analyze_gene_distribution)
    pass
