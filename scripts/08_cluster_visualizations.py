"""
================================================================================
CLUSTER VISUALIZATIONS
================================================================================

PART 1: 3D Cluster Density Surfaces
PART 2: Combined K-Selection Figure

Creates publication-quality visualizations for cluster analysis.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import gaussian_kde
from matplotlib.patches import Patch
from PIL import Image
import os

# Set working directory
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)
plots_dir = "C:/code/projects/intero_mod/plots"
sub_plots_dir = os.path.join(plots_dir, "sub_plots")
supplementary_plots_dir = os.path.join(plots_dir, "supplementary_plots")
os.makedirs(sub_plots_dir, exist_ok=True)
os.makedirs(supplementary_plots_dir, exist_ok=True)

# =============================================================================
# PART 1: 3D CLUSTER DENSITY SURFACES
# =============================================================================

print("=" * 70)
print("PART 1: 3D CLUSTER DENSITY VISUALIZATION")
print("=" * 70)

# Load data with cluster assignments
df = pd.read_csv('supplementary_data/data_with_k_means_clusters.csv')

print(f"Loaded {len(df)} participants")
print(f"\nCluster distribution:")
for cluster in sorted(df['cluster'].unique()):
    n = (df['cluster'] == cluster).sum()
    label = df[df['cluster'] == cluster]['cluster_ordered'].iloc[0]
    print(f"  Cluster {cluster}: {n} ({label})")

# Extract variables
ias = df['ias'].values
iats = df['iats'].values
clusters = df['cluster'].values

# Cluster configuration
cluster_info = {
    1: {'name': 'Hypervigilant', 'color': '#FF7F00'},
    2: {'name': 'Efficient', 'color': '#4DAF4A'},
    3: {'name': 'Uncertain', 'color': '#E41A1C'}
}

# Create grid for density estimation
grid_resolution = 80
kde_bandwidth = 0.2

ias_grid = np.linspace(ias.min() - 5, ias.max() + 5, grid_resolution)
iats_grid = np.linspace(iats.min() - 5, iats.max() + 5, grid_resolution)
IAS_mesh, IATS_mesh = np.meshgrid(ias_grid, iats_grid)
positions = np.vstack([IAS_mesh.ravel(), IATS_mesh.ravel()])

# Create 3D figure
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')

legend_elements = []
plot_order = [2, 1, 3]  # Efficient, Hypervigilant, Uncertain

for cluster_num in plot_order:
    info = cluster_info[cluster_num]
    name = info['name']
    color = info['color']

    mask = clusters == cluster_num
    cluster_ias = ias[mask]
    cluster_iats = iats[mask]
    n_cluster = mask.sum()

    print(f"\nProcessing {name} (n={n_cluster})...")

    # Compute KDE
    cluster_values = np.vstack([cluster_ias, cluster_iats])
    cluster_kernel = gaussian_kde(cluster_values, bw_method=kde_bandwidth)
    Z_cluster = np.reshape(cluster_kernel(positions).T, IAS_mesh.shape)

    # Scale to participant density
    bin_area = (ias_grid[1] - ias_grid[0]) * (iats_grid[1] - iats_grid[0])
    Z_scaled = Z_cluster * n_cluster * bin_area * 100
    Z_scaled[Z_scaled < 0.1] = np.nan

    # Plot surface
    ax.plot_surface(IAS_mesh, IATS_mesh, Z_scaled,
                    color=color, alpha=0.7, edgecolor='none', shade=True)

    legend_elements.append(Patch(facecolor=color, alpha=0.7,
                                  label=f'{name} (n={n_cluster})'))

# Styling
ax.set_xlabel('IAS (Interoceptive Accuracy)', fontsize=12, labelpad=12)
ax.set_ylabel('IATS (Interoceptive Attention)', fontsize=12, labelpad=12)
ax.set_zlabel('Participant Density', fontsize=12, labelpad=12)
ax.set_title('3D Density Landscape of Interoceptive Profiles',
             fontsize=14, fontweight='bold', pad=20)
ax.view_init(elev=25, azim=45)
ax.legend(handles=legend_elements, loc='upper left', fontsize=10)

plt.tight_layout()
cluster_3d_path = os.path.join(sub_plots_dir, "fig1_C_3d_density.png")
plt.savefig(cluster_3d_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
print(f"\n3D density figure saved: plots/sub_plots/fig1_C_3d_density.png")

# =============================================================================
# PART 2: COMBINED K-SELECTION FIGURE
# =============================================================================

print("\n" + "=" * 70)
print("PART 2: COMBINED K-SELECTION FIGURE")
print("=" * 70)

# Check if k-selection figure exists
k_selection_path = os.path.join(supplementary_plots_dir, "1_s_1_k_selection_methods.png")

if os.path.exists(k_selection_path):
    k_selection_img = Image.open(k_selection_path)
    cluster_3d_img = Image.open(cluster_3d_path)

    print(f"K-selection image size: {k_selection_img.size}")
    print(f"3D cluster image size: {cluster_3d_img.size}")

    # Get dimensions
    k_w, k_h = k_selection_img.size
    c3d_w, c3d_h = cluster_3d_img.size

    # Layout: 2 columns (50/50)
    col1_width = int(k_w * 0.50)
    col2_width = int(k_w * 0.50)
    final_width = col1_width + col2_width
    final_height = k_h

    # Crop from k_selection image
    support_crop = k_selection_img.crop((0, 0, int(k_w * 0.50), int(k_h * 0.50)))
    silhouette_crop = k_selection_img.crop((0, int(k_h * 0.50), int(k_w * 0.50), k_h))
    elbow_crop = k_selection_img.crop((int(k_w * 0.50), int(k_h * 0.50), k_w, k_h))

    # Resize
    support_resized = support_crop.resize((col1_width, int(final_height * 0.50)), Image.LANCZOS)
    silhouette_resized = silhouette_crop.resize((col1_width, int(final_height * 0.50)), Image.LANCZOS)

    elbow_height = int(final_height * 0.45)
    cluster_height = final_height - elbow_height
    elbow_resized = elbow_crop.resize((col2_width, elbow_height), Image.LANCZOS)
    cluster_3d_resized = cluster_3d_img.resize((col2_width, cluster_height), Image.LANCZOS)

    # Create canvas
    canvas = Image.new('RGB', (final_width, final_height), 'white')
    canvas.paste(support_resized, (0, 0))
    canvas.paste(silhouette_resized, (0, int(final_height * 0.50)))
    canvas.paste(elbow_resized, (col1_width, 0))
    canvas.paste(cluster_3d_resized, (col1_width, elbow_height))

    combined_path = os.path.join(supplementary_plots_dir, "1_s_2_k_selection_cluster.png")
    canvas.save(combined_path, dpi=(300, 300))
    print(f"\nCombined figure saved: plots/supplementary_plots/1_s_2_k_selection_cluster.png")
    print(f"Final size: {canvas.size}")
else:
    print(f"\nWarning: {k_selection_path} not found. Run R cluster analysis first.")

print("\nDone!")
