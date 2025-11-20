"""
Generate example success data and visualizations for demonstration

This script creates realistic mock data to demonstrate the expected
output format and generate example figures for the repository.
"""

import numpy as np
import csv
from pathlib import Path

np.random.seed(42)

print("Generating example success data...")

# Create results directory
results_dir = Path(__file__).parent.parent / 'results'
results_dir.mkdir(exist_ok=True)

# Generate realistic mock successes
# Based on expected parameter distributions from literature
n_successes = 50

successes = []
for i in range(n_successes):
    # These distributions are informed by Parai and Mukhopadhyay (2018)
    # and Barry and Hilton (2016)

    # Growth rate: tends to cluster around 1e-9
    alpha = 10 ** np.random.normal(-8.8, 0.3)

    # Inflection point: tends to cluster around 2-4 Ga
    beta = np.random.normal(3e9, 0.5e9)
    beta = max(1e9, min(5e9, beta))  # Clip to reasonable range

    # Processing rate: narrow range from observations
    eta = np.random.uniform(7.6e-10, 7.9e-10)

    # Xe carrying capacity: log-uniform over several orders
    xe_cap = 10 ** np.random.uniform(5.5, 7.0)

    # N carrying capacity: log-uniform over several orders
    n_cap = 10 ** np.random.uniform(14, 16)

    successes.append({
        'alpha': alpha,
        'beta': beta,
        'eta': eta,
        'xe_cap': xe_cap,
        'n_cap': n_cap,
        'resfrac': 0.9,
        'lv_frac': 1.0
    })

# Save to CSV
output_path = results_dir / 'success_example.csv'
with open(output_path, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['alpha', 'beta', 'eta', 'xe_cap', 'n_cap', 'resfrac', 'lv_frac'])
    writer.writeheader()
    writer.writerows(successes)

print(f"✓ Generated {n_successes} example successes")
print(f"✓ Saved to: {output_path}")

# Print statistics
print("\nParameter Statistics:")
print(f"  Alpha: {np.mean([s['alpha'] for s in successes]):.2e} ± {np.std([s['alpha'] for s in successes]):.2e}")
print(f"  Beta: {np.mean([s['beta'] for s in successes])/1e9:.2f} ± {np.std([s['beta'] for s in successes])/1e9:.2f} Ga")
print(f"  Eta: {np.mean([s['eta'] for s in successes]):.2e} ± {np.std([s['eta'] for s in successes]):.2e}")
print(f"  Xe capacity: {np.mean([np.log10(s['xe_cap']) for s in successes]):.2f} ± {np.std([np.log10(s['xe_cap']) for s in successes]):.2f} (log10)")
print(f"  N capacity: {np.mean([np.log10(s['n_cap']) for s in successes]):.2f} ± {np.std([np.log10(s['n_cap']) for s in successes]):.2f} (log10)")

print("\nExample data generation complete!")
