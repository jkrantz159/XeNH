"""
Run REAL Monte Carlo simulation to find validated successes

This script runs a genuine Monte Carlo simulation that tests EVERY
parameter combination against actual success criteria for both Xe and N.

This will take several hours but produces REAL validated results.
"""

import numpy as np
import time
import csv
from datetime import datetime
from model_config import ModelConfig
from model_utils import ModelUtils
from parallel_xe_model import parallel_xe_model
from parallel_n_model import parallel_n_model

# Configuration for overnight run
N_RUNS = int(5e5)  # 500k trials (compromise between time and results)
RANDOM_SEED = 42

print("=" * 70)
print("  REAL Monte Carlo Simulation - Finding Validated Successes")
print("=" * 70)
print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
print("⚠️  This finds REAL validated successes (not synthetic data)")
print("   Parameters are tested against actual Xe and N criteria\n")

np.random.seed(RANDOM_SEED)

# Configuration
resfrac = ModelConfig.CONVECTING_MANTLE_FRACTION
lv_frac = ModelConfig.LATE_VENEER_FRACTION_DEFAULT

print(f"Configuration:")
print(f"  Trials: {N_RUNS:,}")
print(f"  Random seed: {RANDOM_SEED}")
print(f"  Expected time: ~{N_RUNS/30/3600:.1f} hours at 30 iter/s")
print(f"  Expected successes: Difficult to predict (very low rate)")
print()

# Create time vector
t = ModelConfig.create_time_vector()
T = ModelConfig.EARTH_AGE_YEARS
atm = ModelConfig.create_atmosphere_evolution(t)

print(f"Model setup:")
print(f"  Time steps: {len(t)}")
print(f"  Integration period: 0 to {T/1e9:.3f} Ga")
print(f"  Xe constraints: 130Xe=[{ModelConfig.XE130_MANTLE_MIN:.1e}, {ModelConfig.XE130_MANTLE_MAX:.1e}]")
print(f"                 128/130=[{ModelConfig.XE128_130_MANTLE_MIN}, {ModelConfig.XE128_130_MANTLE_MAX}]")
print(f"  N constraints:  14N=[{ModelConfig.N14_MANTLE_MIN_MOL:.1e}, {ModelConfig.N14_MANTLE_MAX_MOL:.1e}] mol/g")
print(f"                 15/14=[{ModelConfig.N15_14_RATIO_MANTLE_MIN}, {ModelConfig.N15_14_RATIO_MANTLE_MAX}]")
print()

# Results
successes = []
success_count = 0
start_time = time.time()

# Save progress periodically
SAVE_INTERVAL = 10000  # Save every 10k trials

print("Starting simulation...")
print("Progress reports every 5% of completion\n")

for count in range(1, N_RUNS + 1):
    # Generate random parameters
    alpha = 10 ** np.random.uniform(-10, -7)
    beta = np.random.uniform(0, 10) * 1e9
    eta = np.random.uniform(ModelConfig.ETA_MIN, ModelConfig.ETA_MAX)
    n_cap = 4 * 10 ** np.random.uniform(0, 17)
    xe_cap = 5 * 10 ** np.random.uniform(0, 8)

    # Test BOTH models against actual criteria
    n_succ = parallel_n_model(n_cap, alpha, beta, eta, resfrac, lv_frac, 0, t, T, count)
    xe_succ = parallel_xe_model(xe_cap, alpha, beta, eta, resfrac, lv_frac, 0, atm, t, T, count)

    # Only save if BOTH succeed
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
        print(f"\n✓ SUCCESS #{success_count} found at iteration {count}!")
        print(f"  alpha={alpha:.3e}, beta={beta/1e9:.2f} Ga, eta={eta:.3e}")
        print(f"  Xe_cap={xe_cap:.3e}, N_cap={n_cap:.3e}\n")

        # Save immediately when found
        output_path = ModelUtils.get_output_path('success_REAL.csv', 'results')
        with open(output_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['alpha', 'beta', 'eta', 'xe_cap', 'n_cap', 'resfrac', 'lv_frac'])
            writer.writeheader()
            writer.writerows(successes)

    # Progress reporting
    if count % (N_RUNS // 20) == 0:  # Every 5%
        elapsed = time.time() - start_time
        rate = count / elapsed if elapsed > 0 else 0
        eta_time = (N_RUNS - count) / rate if rate > 0 else 0
        print(f"Progress: {count:,}/{N_RUNS:,} ({count/N_RUNS*100:.1f}%) | "
              f"Rate: {rate:.1f} iter/s | ETA: {eta_time/3600:.2f} hrs | "
              f"Successes: {success_count} ({success_count/count*100:.6f}%)")

elapsed_time = time.time() - start_time

# Final save
if successes:
    output_path = ModelUtils.get_output_path('success_REAL.csv', 'results')
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['alpha', 'beta', 'eta', 'xe_cap', 'n_cap', 'resfrac', 'lv_frac'])
        writer.writeheader()
        writer.writerows(successes)
    print(f"\n✓ Final results saved to: {output_path}")
else:
    print("\n⚠️  No successes found in this run.")
    print("   This is possible with such strict constraints.")
    print("   Try increasing N_RUNS or widening parameter ranges.")

# Summary
print("\n" + "=" * 70)
print("  Simulation Complete")
print("=" * 70)
print(f"End time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"Elapsed time: {elapsed_time/3600:.2f} hours ({elapsed_time/60:.1f} minutes)")
print(f"Average rate: {N_RUNS/elapsed_time:.1f} iterations/second")
print(f"Total trials: {N_RUNS:,}")
print(f"Total successes: {success_count}")
print(f"Success rate: {(success_count/N_RUNS)*100:.6f}%")
print("=" * 70)
print("\nThese are REAL validated successes that satisfy both:")
print("  ✓ Xenon isotope constraints")
print("  ✓ Nitrogen isotope constraints")
