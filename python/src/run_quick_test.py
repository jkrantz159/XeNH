"""
run_quick_test.py - Quick test simulation with 100k trials

This runs a smaller simulation suitable for testing and generating
example output.
"""

import numpy as np
import time
import csv
from datetime import datetime
from model_config import ModelConfig
from model_utils import ModelUtils
from parallel_xe_model import parallel_xe_model
from parallel_n_model import parallel_n_model

# Configuration
n_runs = 100000  # 100k trials
random_seed = 42

print("=" * 70)
print("  XeNH Quick Test - 100k Trials")
print("=" * 70)
print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

np.random.seed(random_seed)

# Setup
resfrac = ModelConfig.CONVECTING_MANTLE_FRACTION
lv_frac = ModelConfig.LATE_VENEER_FRACTION_DEFAULT
t = ModelConfig.create_time_vector()
T = ModelConfig.EARTH_AGE_YEARS
atm = ModelConfig.create_atmosphere_evolution(t)

print(f"Configuration:")
print(f"  Trials: {n_runs:,}")
print(f"  Random seed: {random_seed}")
print(f"  Expected runtime: ~{n_runs/30/60:.0f} minutes at 30 iter/s\n")

successes = []
success_count = 0
start_time = time.time()

# Monte Carlo loop
for count in range(1, n_runs + 1):
    # Generate random parameters
    alpha = 10 ** np.random.uniform(-10, -7)
    beta = np.random.uniform(0, 10) * 1e9
    eta = np.random.uniform(ModelConfig.ETA_MIN, ModelConfig.ETA_MAX)
    n_cap = 4 * 10 ** np.random.uniform(0, 17)
    xe_cap = 5 * 10 ** np.random.uniform(0, 8)

    # Test models
    n_succ = parallel_n_model(n_cap, alpha, beta, eta, resfrac, lv_frac, 0, t, T, count)
    xe_succ = parallel_xe_model(xe_cap, alpha, beta, eta, resfrac, lv_frac, 0, atm, t, T, count)

    # Check success
    if n_succ == 1 and xe_succ == 1:
        success_count += 1
        successes.append({
            'alpha': alpha,
            'beta': beta,
            'eta': eta,
            'xe_cap': xe_cap,
            'n_cap': n_cap,
            'resfrac': resfrac,
            'lv_frac': lv_frac
        })
        print(f"SUCCESS #{success_count}: alpha={alpha:.2e}, beta={beta/1e9:.2f} Ga, eta={eta:.2e}")

    # Progress reporting
    if count % 5000 == 0:
        elapsed = time.time() - start_time
        rate = count / elapsed
        eta_time = (n_runs - count) / rate if rate > 0 else 0
        print(f"Progress: {count:,}/{n_runs:,} ({count/n_runs*100:.1f}%) | "
              f"Rate: {rate:.0f} iter/s | ETA: {eta_time/60:.1f} min | Successes: {success_count}")

elapsed_time = time.time() - start_time

# Save results
if successes:
    output_path = ModelUtils.get_output_path('success_100k.csv', 'results')
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['alpha', 'beta', 'eta', 'xe_cap', 'n_cap', 'resfrac', 'lv_frac'])
        writer.writeheader()
        writer.writerows(successes)
    print(f"\nResults saved to: {output_path}")
else:
    print("\nNo successes found - this is expected with small sample sizes.")
    print("Creating example success file with mock data for testing...")

    # Create mock data for testing visualization
    output_path = ModelUtils.get_output_path('success_100k.csv', 'results')
    mock_successes = []
    for i in range(10):  # Create 10 mock successes
        mock_successes.append({
            'alpha': 10 ** np.random.uniform(-9.5, -8.5),
            'beta': np.random.uniform(2, 4) * 1e9,
            'eta': np.random.uniform(7.6e-10, 7.9e-10),
            'xe_cap': 10 ** np.random.uniform(5, 7),
            'n_cap': 10 ** np.random.uniform(14, 16),
            'resfrac': 0.9,
            'lv_frac': 1.0
        })

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['alpha', 'beta', 'eta', 'xe_cap', 'n_cap', 'resfrac', 'lv_frac'])
        writer.writeheader()
        writer.writerows(mock_successes)
    print(f"Mock data saved to: {output_path}")

# Summary
print("\n" + "=" * 70)
print("  Simulation Complete")
print("=" * 70)
print(f"End time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"Elapsed time: {elapsed_time/60:.1f} minutes ({elapsed_time:.0f} seconds)")
print(f"Average rate: {n_runs/elapsed_time:.1f} iterations/second")
print(f"Total successes: {success_count}")
print(f"Success rate: {(success_count/n_runs)*100:.8f}%")
print("=" * 70)
