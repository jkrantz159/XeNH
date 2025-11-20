"""
Create visualization plots from success data

Generates KDE plots and parameter distribution visualizations
similar to the MATLAB VariableKDEs.m script.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path

# Set publication-quality defaults
plt.style.use('seaborn-v0_8-darkgrid')
plt.rcParams['figure.figsize'] = (10, 6)
plt.rcParams['font.size'] = 11
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 14

def load_success_data(filename='success_example.csv'):
    """Load success data from CSV file"""
    data_path = Path(__file__).parent.parent / 'results' / filename
    df = pd.read_csv(data_path)
    return df

def create_kde_plot(data, xlabel, title, filename, bins=50):
    """Create a kernel density estimate plot"""
    fig, ax = plt.subplots(figsize=(10, 6))

    # Histogram
    ax.hist(data, bins=bins, density=True, alpha=0.6, color='skyblue', edgecolor='black')

    # KDE
    kde = stats.gaussian_kde(data)
    x_range = np.linspace(data.min(), data.max(), 200)
    ax.plot(x_range, kde(x_range), 'k-', linewidth=2, label='KDE')

    ax.set_xlabel(xlabel, fontsize=12, fontweight='bold')
    ax.set_ylabel('Probability Density', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Save figure
    figures_dir = Path(__file__).parent.parent / 'figures'
    figures_dir.mkdir(exist_ok=True)
    save_path = figures_dir / filename
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✓ Saved: {filename}")

def create_parameter_plots(df):
    """Create all parameter distribution plots"""

    # 1. Processing Rate (eta)
    eta_data = df['eta'].values * 1e10
    create_kde_plot(eta_data, r'$\eta \times 10^{10}$ (yr$^{-1}$)',
                   'Processing Rate Distribution', 'processing_rate_kde.png')

    # 2. Inflection Point (beta)
    beta_data = df['beta'].values / 1e9
    create_kde_plot(beta_data, 'Time (Gyr)',
                   'Inflection Point Distribution', 'inflection_point_kde.png')

    # 3. Growth Rate (alpha) - log scale
    alpha_data = np.log10(df['alpha'].values)
    create_kde_plot(alpha_data, r'log$_{10}$($\alpha$)',
                   'Growth Rate Distribution', 'growth_rate_kde.png')

    # 4. Xe Capacity - log scale
    xe_data = np.log10(df['xe_cap'].values)
    create_kde_plot(xe_data, r'log$_{10}$(Xe Capacity) (atoms/g)',
                   'Xenon Carrying Capacity Distribution', 'xe_capacity_kde.png')

    # 5. N Capacity - log scale
    n_data = np.log10(df['n_cap'].values)
    create_kde_plot(n_data, r'log$_{10}$(N Capacity) (atoms/g)',
                   'Nitrogen Carrying Capacity Distribution', 'n_capacity_kde.png')

    # 6. Combined Recycling Plot (Xe vs N)
    fig, ax = plt.subplots(figsize=(10, 6))

    xe_kde = stats.gaussian_kde(xe_data)
    n_kde = stats.gaussian_kde(n_data)

    x_range = np.linspace(5, 20, 200)
    ax.plot(x_range, xe_kde(x_range), 'r-', linewidth=2, label='Xe')
    ax.plot(x_range, n_kde(x_range), 'g-', linewidth=2, label='N')

    ax.set_xlabel(r'log$_{10}$(Recycling Capacity)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Probability Density', fontsize=12, fontweight='bold')
    ax.set_title('Recycling Capacity Comparison', fontsize=14, fontweight='bold')
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    figures_dir = Path(__file__).parent.parent / 'figures'
    plt.tight_layout()
    plt.savefig(figures_dir / 'recycling_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("✓ Saved: recycling_comparison.png")

    # 7. 2D Parameter Space (Beta vs Alpha)
    fig, ax = plt.subplots(figsize=(10, 8))

    scatter = ax.scatter(df['beta'].values / 1e9, np.log10(df['alpha'].values),
                        c=df['eta'].values * 1e10, cmap='viridis',
                        s=100, alpha=0.6, edgecolors='black', linewidth=0.5)

    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label(r'$\eta \times 10^{10}$ (yr$^{-1}$)', fontsize=11)

    ax.set_xlabel('Inflection Point (Gyr)', fontsize=12, fontweight='bold')
    ax.set_ylabel(r'log$_{10}$($\alpha$)', fontsize=12, fontweight='bold')
    ax.set_title('Parameter Space: Inflection Point vs Growth Rate', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(figures_dir / 'parameter_space_2d.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("✓ Saved: parameter_space_2d.png")

    # 8. Summary Statistics Table Figure
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.axis('tight')
    ax.axis('off')

    # Calculate statistics
    stats_data = [
        ['Parameter', 'Mean', 'Std Dev', 'Min', 'Max'],
        ['α (×10⁻⁹)', f'{df["alpha"].mean()*1e9:.2f}', f'{df["alpha"].std()*1e9:.2f}',
         f'{df["alpha"].min()*1e9:.2f}', f'{df["alpha"].max()*1e9:.2f}'],
        ['β (Gyr)', f'{df["beta"].mean()/1e9:.2f}', f'{df["beta"].std()/1e9:.2f}',
         f'{df["beta"].min()/1e9:.2f}', f'{df["beta"].max()/1e9:.2f}'],
        ['η (×10⁻¹⁰)', f'{df["eta"].mean()*1e10:.2f}', f'{df["eta"].std()*1e10:.3f}',
         f'{df["eta"].min()*1e10:.2f}', f'{df["eta"].max()*1e10:.2f}'],
        ['log₁₀(Xe)', f'{np.log10(df["xe_cap"]).mean():.2f}', f'{np.log10(df["xe_cap"]).std():.2f}',
         f'{np.log10(df["xe_cap"]).min():.2f}', f'{np.log10(df["xe_cap"]).max():.2f}'],
        ['log₁₀(N)', f'{np.log10(df["n_cap"]).mean():.2f}', f'{np.log10(df["n_cap"]).std():.2f}',
         f'{np.log10(df["n_cap"]).min():.2f}', f'{np.log10(df["n_cap"]).max():.2f}'],
    ]

    table = ax.table(cellText=stats_data, cellLoc='center', loc='center',
                    colWidths=[0.25, 0.15, 0.15, 0.15, 0.15])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2)

    # Style header row
    for i in range(5):
        table[(0, i)].set_facecolor('#4CAF50')
        table[(0, i)].set_text_props(weight='bold', color='white')

    ax.set_title('Summary Statistics of Successful Parameters',
                fontsize=14, fontweight='bold', pad=20)

    plt.tight_layout()
    plt.savefig(figures_dir / 'summary_statistics.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("✓ Saved: summary_statistics.png")

if __name__ == "__main__":
    print("=" * 60)
    print("  Creating Visualizations")
    print("=" * 60)

    # Load data
    print("\nLoading success data...")
    df = load_success_data()
    print(f"✓ Loaded {len(df)} successful parameter sets\n")

    # Create plots
    print("Generating figures...")
    create_parameter_plots(df)

    print("\n" + "=" * 60)
    print("  Visualization Complete!")
    print("=" * 60)
    print(f"\nGenerated 8 figures in python/figures/")
    print("These demonstrate the expected output format and")
    print("parameter distributions from successful model runs.")
