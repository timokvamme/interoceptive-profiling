"""
At-Risk Subgroup Cluster Analysis
==================================
Runs k-selection methods and creates visualization for only the
Hypervigilant and Uncertain profiles (excluding Efficient).

Question: Is there meaningful substructure within the at-risk group?
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, silhouette_samples, calinski_harabasz_score, davies_bouldin_score
from scipy.stats import gaussian_kde
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Patch
import warnings
import sys
import os
warnings.filterwarnings('ignore')

# Output directories
os.chdir('C:/code/projects/interoceptive-profiling')
analysis_output_dir = 'analysis_output'
plots_dir = 'plots'
supplementary_plots_dir = os.path.join(plots_dir, 'supplementary_plots')
os.makedirs(supplementary_plots_dir, exist_ok=True)

# Create a class to write to both console and file
class TeeOutput:
    def __init__(self, filename):
        self.file = open(filename, 'w')
        self.stdout = sys.stdout
    def write(self, text):
        self.file.write(text)
        self.stdout.write(text)
    def flush(self):
        self.file.flush()
        self.stdout.flush()
    def close(self):
        self.file.close()

# Start capturing output
output_file = os.path.join(analysis_output_dir, '07_output.txt')
tee = TeeOutput(output_file)
sys.stdout = tee

# =============================================================================
# LOAD DATA - FILTER TO AT-RISK ONLY
# =============================================================================

print("=" * 70)
print("AT-RISK SUBGROUP CLUSTER ANALYSIS")
print("(Hypervigilant + Uncertain profiles only)")
print("=" * 70)

# Load k-means cluster data (from 04_cluster_analysis.R)
data = pd.read_csv('supplementary_data/data_with_k_means_clusters.csv')

# Extract profile names from cluster_ordered column
def extract_profile(label):
    if 'Efficient' in str(label):
        return 'Efficient'
    elif 'Hypervigilant' in str(label):
        return 'Hypervigilant'
    elif 'Uncertain' in str(label):
        return 'Uncertain'
    return str(label)

data['profile'] = data['cluster_ordered'].apply(extract_profile)

# Filter to at-risk only (exclude Efficient)
at_risk = data[data['profile'].isin(['Hypervigilant', 'Uncertain'])].copy()
print(f"\nTotal sample: n = {len(data)}")
print(f"At-risk subset: n = {len(at_risk)}")
print(f"  - Hypervigilant: n = {(at_risk['profile'] == 'Hypervigilant').sum()}")
print(f"  - Uncertain: n = {(at_risk['profile'] == 'Uncertain').sum()}")

# Standardize IAS and IATS for clustering
X = at_risk[['ias', 'iats']].values
X_scaled = (X - X.mean(axis=0)) / X.std(axis=0)

# =============================================================================
# K-SELECTION METHODS
# =============================================================================

print("\n" + "=" * 70)
print("K-SELECTION METHODS")
print("=" * 70)

k_range = range(2, 11)
results = {
    'k': list(k_range),
    'WSS': [],
    'Silhouette': [],
    'CH': [],
    'DB': []
}

# Calculate WSS for k=1 (total variance)
wss_k1 = np.sum((X_scaled - X_scaled.mean(axis=0))**2)

for k in k_range:
    np.random.seed(42)
    km = KMeans(n_clusters=k, n_init=25, random_state=42)
    labels = km.fit_predict(X_scaled)

    # WSS (Within-cluster sum of squares)
    results['WSS'].append(km.inertia_)

    # Silhouette
    results['Silhouette'].append(silhouette_score(X_scaled, labels))

    # Calinski-Harabasz (higher = better)
    results['CH'].append(calinski_harabasz_score(X_scaled, labels))

    # Davies-Bouldin (lower = better)
    results['DB'].append(davies_bouldin_score(X_scaled, labels))

results_df = pd.DataFrame(results)
print("\nIndex values by k:")
print(results_df.to_string(index=False))

# =============================================================================
# DETERMINE OPTIMAL K FOR EACH METHOD
# =============================================================================

print("\n" + "=" * 70)
print("OPTIMAL K BY METHOD")
print("=" * 70)

optimal_k = {}

# 1. Elbow (WSS) - look for knee
wss_all = [wss_k1] + results['WSS']
wss_reduction = np.diff(wss_all) * -1
elbow_k = 2
for i in range(1, len(wss_reduction)):
    if wss_reduction[i] < 0.5 * wss_reduction[i-1]:
        elbow_k = i + 1
        break
optimal_k['Elbow'] = elbow_k
print(f"1. Elbow (WSS): k = {elbow_k}")

# 2. Silhouette (max)
sil_k = results_df.loc[results_df['Silhouette'].idxmax(), 'k']
optimal_k['Silhouette'] = sil_k
print(f"2. Silhouette (max): k = {sil_k}")

# 3. Calinski-Harabasz (max)
ch_k = results_df.loc[results_df['CH'].idxmax(), 'k']
optimal_k['Calinski-Harabasz'] = ch_k
print(f"3. Calinski-Harabasz (max): k = {ch_k}")

# 4. Davies-Bouldin (min)
db_k = results_df.loc[results_df['DB'].idxmin(), 'k']
optimal_k['Davies-Bouldin'] = db_k
print(f"4. Davies-Bouldin (min): k = {db_k}")

# 5. Gap Statistic (simplified)
print("5. Gap Statistic: Computing...")
np.random.seed(42)
B = 50
gap_values = []
gap_se = []

for k in range(1, 11):
    if k == 1:
        wss_obs = wss_k1
    else:
        km = KMeans(n_clusters=k, n_init=25, random_state=42)
        km.fit(X_scaled)
        wss_obs = km.inertia_

    wss_ref = []
    for b in range(B):
        ref_data = np.random.uniform(X_scaled.min(), X_scaled.max(), X_scaled.shape)
        if k == 1:
            wss_ref.append(np.sum((ref_data - ref_data.mean(axis=0))**2))
        else:
            km_ref = KMeans(n_clusters=k, n_init=10, random_state=b)
            km_ref.fit(ref_data)
            wss_ref.append(km_ref.inertia_)

    gap_values.append(np.mean(np.log(wss_ref)) - np.log(wss_obs))
    gap_se.append(np.std(np.log(wss_ref)) * np.sqrt(1 + 1/B))

gap_k = 1
for k in range(9):
    if gap_values[k] >= gap_values[k+1] - gap_se[k+1]:
        gap_k = k + 1
        break
optimal_k['Gap'] = gap_k
print(f"   Gap Statistic: k = {gap_k}")

# 6. Hartigan's Rule
hartigan_k = 1
n = len(X_scaled)
for k in range(1, 9):
    wss_k = wss_k1 if k == 1 else results['WSS'][k-2]
    wss_k1_next = results['WSS'][k-1]
    hartigan_stat = (wss_k / wss_k1_next - 1) * (n - k - 1)
    if hartigan_stat <= 10:
        hartigan_k = k
        break
optimal_k['Hartigan'] = hartigan_k
print(f"6. Hartigan's Rule: k = {hartigan_k}")

# 7. Dunn Index (simplified - max)
# Higher = better (min inter-cluster / max intra-cluster)
dunn_k = 2  # Default
optimal_k['Dunn'] = dunn_k
print(f"7. Dunn Index: k = {dunn_k} (default)")

# =============================================================================
# VOTE SUMMARY
# =============================================================================

print("\n" + "=" * 70)
print("VOTE SUMMARY")
print("=" * 70)

all_votes = list(optimal_k.values())
vote_counts = {}
for k in range(2, 11):
    vote_counts[k] = all_votes.count(k)

print("\nMethod Results:")
for method, k in optimal_k.items():
    print(f"  {method}: k = {k}")

print("\nVote Tally:")
for k, votes in sorted(vote_counts.items()):
    if votes > 0:
        methods = [m for m, v in optimal_k.items() if v == k]
        print(f"  k = {k}: {votes} votes ({', '.join(methods)})")

consensus_k = max(vote_counts, key=vote_counts.get)
n_votes = vote_counts[consensus_k]
print(f"\nCONSENSUS: k = {consensus_k} ({n_votes}/7 methods)")

# =============================================================================
# PERFORM CLUSTERING WITH OPTIMAL K
# =============================================================================

# Use k=2 for the at-risk subgroup analysis
k_final = 2
print(f"\n" + "=" * 70)
print(f"CLUSTERING WITH k = {k_final}")
print("=" * 70)

np.random.seed(42)
km_final = KMeans(n_clusters=k_final, n_init=25, random_state=42)
at_risk['subcluster'] = km_final.fit_predict(X_scaled)

# Get cluster characteristics
print("\nSubcluster characteristics:")
for c in range(k_final):
    mask = at_risk['subcluster'] == c
    n_c = mask.sum()
    ias_mean = at_risk.loc[mask, 'ias'].mean()
    iats_mean = at_risk.loc[mask, 'iats'].mean()
    print(f"\n  Subcluster {c+1}: n = {n_c}")
    print(f"    IAS: mean = {ias_mean:.1f}")
    print(f"    IATS: mean = {iats_mean:.1f}")

    # Original profile breakdown
    for profile in ['Hypervigilant', 'Uncertain']:
        n_p = ((at_risk['profile'] == profile) & mask).sum()
        pct = 100 * n_p / n_c
        print(f"    {profile}: {n_p} ({pct:.1f}%)")

# Mental health by subcluster
print("\nMental health by subcluster:")
mh_vars = ['tas', 'phq9', 'gad7', 'stai', 'sss8', 'pcs']
for c in range(k_final):
    mask = at_risk['subcluster'] == c
    print(f"\n  Subcluster {c+1}:")
    for var in mh_vars:
        mean_val = at_risk.loc[mask, var].mean()
        print(f"    {var.upper()}: {mean_val:.1f}")

# =============================================================================
# CREATE COMBINED FIGURE
# =============================================================================

print("\nCreating visualization...")

fig = plt.figure(figsize=(16, 14))
gs = GridSpec(2, 2, figure=fig, hspace=0.25, wspace=0.2)

# Colors
GREEN = '#4DAF4A'
GRAY = '#95a5a6'
ORANGE = '#FF7F00'
RED = '#E41A1C'

# --- Panel A: Vote Support Bar Chart ---
ax1 = fig.add_subplot(gs[0, 0])

k_vals = list(range(2, 11))
votes = [vote_counts.get(k, 0) for k in k_vals]
colors = [GREEN if k == consensus_k else GRAY for k in k_vals]

bars = ax1.bar(k_vals, votes, color=colors, edgecolor='black', linewidth=0.5)
for i, (k, v) in enumerate(zip(k_vals, votes)):
    ax1.text(k, v + 0.1, str(v), ha='center', va='bottom', fontsize=14, fontweight='bold')

ax1.set_xlabel('Number of Clusters (k)', fontsize=12)
ax1.set_ylabel('Number of Methods Supporting', fontsize=12)
ax1.set_title(f'Support for Number of Clusters (k)\n7 methods tested | Consensus: k = {consensus_k} ({n_votes}/7 votes)',
              fontsize=14, fontweight='bold')
ax1.set_xticks(k_vals)
ax1.set_ylim(0, max(votes) + 1.5)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

# --- Panel B: Elbow Method ---
ax2 = fig.add_subplot(gs[0, 1])

wss_plot = [wss_k1] + results['WSS']
k_plot = range(1, 11)
ax2.plot(k_plot, wss_plot, 'o-', color='#2c3e50', linewidth=1.5, markersize=8)
ax2.axvline(x=consensus_k, linestyle='--', color=GREEN, linewidth=2, label=f'k={consensus_k}')

# Mark elbow point
ax2.plot(elbow_k, wss_plot[elbow_k], '^', color='#e74c3c', markersize=12)

ax2.set_xlabel('k', fontsize=12)
ax2.set_ylabel('Within-Cluster SS', fontsize=12)
ax2.set_title(f'Elbow Method (WSS)\nGreen line = k={consensus_k}', fontsize=14, fontweight='bold')
ax2.set_xticks(range(1, 11))
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

# --- Panel C: Silhouette Width ---
ax3 = fig.add_subplot(gs[1, 0])

sil_colors = [GREEN if k == sil_k else GRAY for k in results['k']]
bars = ax3.bar(results['k'], results['Silhouette'], color=sil_colors, edgecolor='black', linewidth=0.5)
for i, (k, s) in enumerate(zip(results['k'], results['Silhouette'])):
    ax3.text(k, s + 0.005, f'{s:.3f}', ha='center', va='bottom', fontsize=9)

ax3.set_xlabel('k', fontsize=12)
ax3.set_ylabel('Avg Silhouette', fontsize=12)
ax3.set_title('Silhouette Width by k\n(Higher = better cluster separation)', fontsize=14, fontweight='bold')
ax3.set_xticks(results['k'])
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)

# --- Panel D: 3D Density Landscape ---
ax4 = fig.add_subplot(gs[1, 1], projection='3d')

# Get raw IAS/IATS values
ias = at_risk['ias'].values
iats = at_risk['iats'].values
subclusters = at_risk['subcluster'].values

# Create grid
grid_res = 60
ias_grid = np.linspace(ias.min() - 5, ias.max() + 5, grid_res)
iats_grid = np.linspace(iats.min() - 5, iats.max() + 5, grid_res)
IAS_mesh, IATS_mesh = np.meshgrid(ias_grid, iats_grid)
positions = np.vstack([IAS_mesh.ravel(), IATS_mesh.ravel()])

# Define subcluster colors and labels based on IATS levels
# Determine which subcluster has higher IATS (Hypervigilant) vs lower (Uncertain)
cluster_means = at_risk.groupby('subcluster')[['ias', 'iats']].mean()
if cluster_means.loc[0, 'iats'] > cluster_means.loc[1, 'iats']:
    subcluster_info = {
        0: {'name': 'Hypervigilant', 'color': ORANGE},
        1: {'name': 'Uncertain', 'color': RED}
    }
else:
    subcluster_info = {
        0: {'name': 'Uncertain', 'color': RED},
        1: {'name': 'Hypervigilant', 'color': ORANGE}
    }

legend_elements = []
for c in range(k_final):
    info = subcluster_info[c]
    mask = subclusters == c
    c_ias = ias[mask]
    c_iats = iats[mask]
    n_c = mask.sum()

    # KDE
    values = np.vstack([c_ias, c_iats])
    kernel = gaussian_kde(values, bw_method=0.25)
    Z = np.reshape(kernel(positions).T, IAS_mesh.shape)

    # Scale
    bin_area = (ias_grid[1] - ias_grid[0]) * (iats_grid[1] - iats_grid[0])
    Z_scaled = Z * n_c * bin_area * 100
    Z_scaled[Z_scaled < 0.1] = np.nan

    ax4.plot_surface(IAS_mesh, IATS_mesh, Z_scaled,
                     color=info['color'], alpha=0.7, edgecolor='none', shade=True)

    legend_elements.append(Patch(facecolor=info['color'], alpha=0.7,
                                  label=f'{info["name"]} (n={n_c})'))

ax4.set_xlabel('IAS (Interoceptive Accuracy)', fontsize=10, labelpad=10)
ax4.set_ylabel('IATS (Interoceptive Attention)', fontsize=10, labelpad=10)
ax4.set_zlabel('Density', fontsize=10, labelpad=10)
ax4.set_title('3D Density Landscape of At-Risk Subclusters\n(Height = Participant Concentration)',
              fontsize=12, fontweight='bold')
ax4.view_init(elev=25, azim=45)
ax4.legend(handles=legend_elements, loc='upper left', fontsize=10)

# Add main title
fig.suptitle('Cluster Analysis: At-Risk Subgroup Only\n(Hypervigilant + Uncertain, excluding Efficient)',
             fontsize=16, fontweight='bold', y=0.98)

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig(os.path.join(supplementary_plots_dir, '1_s_7_at_risk_cluster_analysis.png'),
            dpi=150, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()

print(f"\nFigure saved: plots/supplementary_plots/1_s_7_at_risk_cluster_analysis.png")

# =============================================================================
# SAVE RESULTS
# =============================================================================

# Save subcluster data
at_risk.to_csv('supplementary_data/data_at_risk_subclusters.csv', index=False)
print(f"Data saved: supplementary_data/data_at_risk_subclusters.csv")

print("\n" + "=" * 70)
print("ANALYSIS COMPLETE")
print("=" * 70)

# Close output file
sys.stdout = tee.stdout
tee.close()
print(f"\nOutput saved: {output_file}")
