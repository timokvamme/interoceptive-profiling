"""
================================================================================
SIMPSON'S PARADOX: BAYES FACTOR ANALYSIS
================================================================================

Analyzes the IAS-IATS correlation paradox using Bayes factors.
Uses k-means cluster assignments from 06_cluster_analysis.R

Outputs:
  - 1_s_6_simpsons_paradox_bayes.png
================================================================================
"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.integrate import quad
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os

# Set working directory
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)
plots_dir = "C:/code/projects/intero_mod/plots"
supplementary_plots_dir = os.path.join(plots_dir, "supplementary_plots")
analysis_output_dir = "C:/code/projects/intero_mod/analysis_output"
os.makedirs(supplementary_plots_dir, exist_ok=True)
os.makedirs(analysis_output_dir, exist_ok=True)

# Output redirection class
import sys
class TeeOutput:
    def __init__(self, filename):
        self.file = open(filename, 'w')
        self.stdout = sys.stdout
        sys.stdout = self
    def write(self, text):
        self.file.write(text)
        self.stdout.write(text)
    def flush(self):
        self.file.flush()
        self.stdout.flush()
    def close(self):
        sys.stdout = self.stdout
        self.file.close()

tee = TeeOutput(os.path.join(analysis_output_dir, '06_output.txt'))

def jeffreys_bf_correlation(r, n):
    """
    Calculate BF01 (support for null) using the method from JASP/BayesFactor package.
    Based on Ly, Verhagen, & Wagenmakers (2016).
    """
    def log_likelihood(rho, r, n):
        """Log likelihood of data given rho"""
        if abs(rho) >= 1:
            return -np.inf
        return ((n - 1) / 2) * np.log(1 - rho**2) - ((n - 1) / 2) * np.log(1 - 2*rho*r + rho**2)

    def prior(rho):
        """Uniform prior on rho"""
        return 1/2  # uniform on [-1, 1]

    # Marginal likelihood under H1
    def integrand(rho):
        if abs(rho) >= 0.9999:
            return 0
        ll = log_likelihood(rho, r, n)
        return np.exp(ll) * prior(rho)

    try:
        ml_h1, _ = quad(integrand, -0.9999, 0.9999, limit=100)
    except:
        ml_h1 = 1e-10

    # Likelihood under H0 (rho = 0)
    ml_h0 = np.exp(log_likelihood(0, r, n))

    # BF01 = P(data|H0) / P(data|H1)
    bf01 = ml_h0 / ml_h1 if ml_h1 > 0 else np.nan
    bf10 = 1 / bf01 if bf01 > 0 else np.nan

    return bf01, bf10

# Load k-means cluster data (from 04_cluster_analysis.R)
data = pd.read_csv('supplementary_data/data_with_k_means_clusters.csv')

# Get profile names from cluster_ordered column
# Extract the short profile name (Efficient, Hypervigilant, Uncertain)
def extract_profile(label):
    if 'Efficient' in str(label):
        return 'Efficient'
    elif 'Hypervigilant' in str(label):
        return 'Hypervigilant'
    elif 'Uncertain' in str(label):
        return 'Uncertain'
    return str(label)

data['profile'] = data['cluster_ordered'].apply(extract_profile)

# Colors matching 06_cluster_analysis.R
colors = {
    'Efficient': '#4DAF4A',
    'Hypervigilant': '#FF7F00',
    'Uncertain': '#E41A1C'
}

# Calculate correlations and Bayes factors
overall_r, overall_p = stats.pearsonr(data['ias'], data['iats'])
overall_bf01, overall_bf10 = jeffreys_bf_correlation(overall_r, len(data))

print("=" * 70)
print("SIMPSON'S PARADOX: BAYES FACTOR ANALYSIS")
print("=" * 70)
print(f"\nData source: data_with_clusters.csv (k-means)")
print(f"N = {len(data)}")
print(f"\nProfiles: {data['profile'].value_counts().to_dict()}")

print(f"\n=== OVERALL CORRELATION ===")
print(f"r = {overall_r:.3f}, p = {overall_p:.4f}, n = {len(data)}")
print(f"BF01 = {overall_bf01:.2f} (support for null)")

results = {}
print(f"\n=== WITHIN-PROFILE CORRELATIONS ===")
for profile in ['Efficient', 'Hypervigilant', 'Uncertain']:
    subset = data[data['profile'] == profile]
    if len(subset) < 3:
        print(f"{profile}: Not enough data (n={len(subset)})")
        continue
    r, p = stats.pearsonr(subset['ias'], subset['iats'])
    n = len(subset)
    bf01, bf10 = jeffreys_bf_correlation(r, n)

    results[profile] = {'r': r, 'p': p, 'n': n, 'bf01': bf01, 'bf10': bf10}

    sig = '*' if p < 0.05 else 'ns'
    if bf01 > 10:
        strength = "Strong evidence for null"
    elif bf01 > 3:
        strength = "Moderate evidence for null"
    elif bf01 > 1:
        strength = "Anecdotal evidence for null"
    else:
        strength = "Evidence favors alternative"
    print(f"{profile}: r = {r:.3f} ({sig}), n = {n}, BF01 = {bf01:.1f} ({strength})")

# Create visualization with boxes on right side
fig, ax = plt.subplots(figsize=(14, 8))

# Plot each profile
for profile in ['Efficient', 'Hypervigilant', 'Uncertain']:
    if profile not in results:
        continue
    subset = data[data['profile'] == profile]
    ax.scatter(subset['ias'], subset['iats'],
               c=colors[profile], alpha=0.6, s=50,
               label=profile, edgecolors='white', linewidth=0.5)

# Add regression lines for each cluster
for profile in ['Efficient', 'Hypervigilant', 'Uncertain']:
    if profile not in results:
        continue
    subset = data[data['profile'] == profile]
    z = np.polyfit(subset['ias'], subset['iats'], 1)
    p_line = np.poly1d(z)
    x_range = np.linspace(subset['ias'].min(), subset['ias'].max(), 100)
    ax.plot(x_range, p_line(x_range), c=colors[profile],
            linewidth=2.5, linestyle='--', alpha=0.8)

# Add overall regression line
z_overall = np.polyfit(data['ias'], data['iats'], 1)
p_overall = np.poly1d(z_overall)
x_range_all = np.linspace(data['ias'].min(), data['ias'].max(), 100)
ax.plot(x_range_all, p_overall(x_range_all), 'k-', linewidth=3.5, alpha=0.9)

# Labels (no title)
ax.set_xlabel('Interoceptive Accuracy (IAS)', fontsize=14, fontweight='bold')
ax.set_ylabel('Interoceptive Attention (IATS)', fontsize=14, fontweight='bold')

# Frequentist stats box (right side, top)
textstr = "FREQUENTIST\n"
textstr += f"Overall: r = {overall_r:.3f}, p = {overall_p:.3f}*\n\n"
textstr += "Within-Profile:\n"
for profile in ['Efficient', 'Hypervigilant', 'Uncertain']:
    if profile not in results:
        continue
    r = results[profile]
    sig = '*' if r['p'] < 0.05 else 'ns'
    textstr += f"  {profile}: r = {r['r']:.3f} ({sig})\n"

# Bayesian stats box (right side, middle)
bf_text = "BAYESIAN\n"
bf_text += "BF01 = support for null (r = 0)\n\n"
bf_text += f"Overall: BF01 = {overall_bf01:.2f}\n"
if overall_bf01 < 1:
    bf_text += f"(Evidence FOR correlation)\n\n"
elif overall_bf01 < 3:
    bf_text += f"(Anecdotal)\n\n"
else:
    bf_text += f"(Inconclusive)\n\n"
bf_text += "Within-Profile:\n"
for profile in ['Efficient', 'Hypervigilant', 'Uncertain']:
    if profile not in results:
        continue
    bf = results[profile]['bf01']
    if bf > 10:
        strength = "Strong"
    elif bf > 3:
        strength = "Moderate"
    else:
        strength = "Anecdotal"
    bf_text += f"  {profile}: BF01 = {bf:.1f}\n"
    bf_text += f"    ({strength} evidence for null)\n"

# Legend with k-means profile names (right side, bottom)
legend_elements = [
    Line2D([0], [0], color='k', linewidth=3.5, label=f'Overall (r = {overall_r:.3f}*)')
]
for profile in ['Efficient', 'Hypervigilant', 'Uncertain']:
    if profile not in results:
        continue
    r = results[profile]
    legend_elements.append(
        Line2D([0], [0], marker='o', color='w', markerfacecolor=colors[profile],
               markersize=10, linestyle='--', markeredgecolor='white',
               label=f'{profile} (r = {r["r"]:.3f}, BF01 = {r["bf01"]:.1f})')
    )

plt.tight_layout()
plt.subplots_adjust(right=0.65)  # Make room for boxes on right

# Place all boxes on the right side, stacked vertically with proper spacing
props = dict(boxstyle='round', facecolor='white', alpha=0.95, edgecolor='gray')
props2 = dict(boxstyle='round', facecolor='lightyellow', alpha=0.95, edgecolor='orange')

# Frequentist box (top right)
fig.text(0.67, 0.98, textstr, fontsize=8, verticalalignment='top',
         fontfamily='monospace', bbox=props)

# Bayesian box (middle right) - positioned lower to avoid overlap
fig.text(0.67, 0.52, bf_text, fontsize=8, verticalalignment='top',
         fontfamily='monospace', bbox=props2)

# Legend (bottom right) - positioned at bottom
ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.02, 0.18),
          fontsize=8, frameon=True)

plt.savefig(os.path.join(supplementary_plots_dir, "1_s_6_simpsons_paradox_bayes.png"),
            dpi=150, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()

print("\nFigure saved: plots/supplementary_plots/1_s_6_simpsons_paradox_bayes.png")

# Close output file
tee.close()
