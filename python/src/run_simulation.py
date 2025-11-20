"""
run_simulation.py - Monte Carlo simulation runner for XeNH models

This script runs a Monte Carlo simulation to find parameter combinations
that simultaneously satisfy constraints from Xenon and Nitrogen isotope
systematics in Earth's mantle.
"""

import numpy as np
import time
import csv
from datetime import datetime
from model_config import ModelConfig
from model_utils import ModelUtils
from parallel_xe_model import parallel_xe_model
from parallel_n_model import parallel_n_model


def run_monte_carlo_simulation(n_runs=int(1e6), random_seed=None):
    """
    Run Monte Carlo simulation

    Parameters:
        n_runs (int): Number of Monte Carlo trials
        random_seed (int): Random seed for reproducibility

    Returns:
        list: Successful parameter combinations
    """
    print("=" * 70)
    print("  XeNH Monte Carlo Simulation - Python Version")
    print("=" * 70)
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

    # Set random seed if provided
    if random_seed is not None:
        np.random.seed(random_seed)
        print(f"Random seed: {random_seed}")

    # Configuration
    resfrac = ModelConfig.CONVECTING_MANTLE_FRACTION
    lv_frac = ModelConfig.LATE_VENEER_FRACTION_DEFAULT

    print(f"Configuration:")
    print(f"  Reservoir fraction: {resfrac * 100:.1f}%")
    print(f"  Late veneer fraction: {lv_frac:.1f}% Earth mass")
    print(f"  Monte Carlo trials: {n_runs:.0e}")
    print()

    # Create time vector and atmospheric evolution
    t = ModelConfig.create_time_vector()
    T = ModelConfig.EARTH_AGE_YEARS
    atm = ModelConfig.create_atmosphere_evolution(t)

    print(f"Time integration:")
    print(f"  Time steps: {len(t)}")
    print(f"  Integration period: 0 to {T/1e9:.3f} Ga\n")

    # Initialize results
    successes = []
    success_count = 0

    print("Starting Monte Carlo simulation...")
    print(f"Progress will be reported every 1% ({n_runs/100:.0e} iterations)\n")

    start_time = time.time()

    # Monte Carlo loop
    for count in range(1, n_runs + 1):
        # Generate random parameters
        alpha = 10 ** np.random.uniform(-10, -7)  # Log-uniform: 1E-10 to 1E-7
        beta = np.random.uniform(0, 10) * 1e9     # Uniform: 0 to 10 Gyr
        eta = np.random.uniform(ModelConfig.ETA_MIN, ModelConfig.ETA_MAX)
        n_cap = 4 * 10 ** np.random.uniform(0, 17)  # Log-uniform: 4 to 4E17
        xe_cap = 5 * 10 ** np.random.uniform(0, 8)  # Log-uniform: 5 to 5E8

        # Test Nitrogen model
        n_succ = parallel_n_model(n_cap, alpha, beta, eta, resfrac, lv_frac,
                                  0, t, T, count)

        # Test Xenon model
        xe_succ = parallel_xe_model(xe_cap, alpha, beta, eta, resfrac, lv_frac,
                                    0, atm, t, T, count)

        # Check if both models succeed
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
            print(f"SUCCESS #{success_count}: alpha={alpha:.2e}, beta={beta/1e9:.2f} Ga, "
                  f"eta={eta:.2e}, Xe={xe_cap:.2e}, N={n_cap:.2e}")

        # Progress reporting
        if count % (n_runs // 100) == 0:
            elapsed = time.time() - start_time
            progress = (count / n_runs) * 100
            rate = count / elapsed if elapsed > 0 else 0
            eta_time = (n_runs - count) / rate if rate > 0 else 0
            print(f"Progress: {progress:.1f}% ({count}/{n_runs}) | "
                  f"Rate: {rate:.0f} iter/s | ETA: {eta_time/60:.1f} min | "
                  f"Successes: {success_count}")

    elapsed_time = time.time() - start_time

    # Save results
    if successes:
        output_path = ModelUtils.get_output_path('success.csv', 'results')
        with open(output_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['alpha', 'beta', 'eta',
                                                   'xe_cap', 'n_cap',
                                                   'resfrac', 'lv_frac'])
            writer.writeheader()
            writer.writerows(successes)
        print(f"\nResults saved to: {output_path}")

    # Summary
    print("\n" + "=" * 70)
    print("  Simulation Complete")
    print("=" * 70)
    print(f"End time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Elapsed time: {elapsed_time/60:.2f} minutes")
    print(f"Total successes: {success_count}")
    print(f"Success rate: {(success_count/n_runs)*100:.6f}%")
    print("=" * 70)

    return successes


if __name__ == "__main__":
    # Run simulation with 1 million trials
    results = run_monte_carlo_simulation(n_runs=int(1e6), random_seed=42)

    if results:
        print(f"\nFound {len(results)} successful parameter combinations!")
    else:
        print("\nNo successful combinations found. Try increasing n_runs.")
