# Running XeNH Simulations on Google Colab

This guide explains how to run the XeNH Monte Carlo simulations using Google Colab's free cloud computing resources.

## ⚠️ Important: Branch Information

**Current Status:** The notebook is configured to use the `claude/google-colab-support-01CL1eMMoLYfXqHoiGmoNCdM` branch, which contains both the Python implementation and Colab notebook.

**After merging to main:** Once this branch is merged into the main branch, the notebook will be updated to use the main branch automatically. Until then, the notebook will continue to use the specified feature branch.

## Quick Start

### Option 1: Open from GitHub (Recommended)

1. Go to [Google Colab](https://colab.research.google.com/)
2. Click `File` → `Open notebook` → `GitHub` tab
3. Enter the repository URL: `jkrantz159/XeNH`
4. Select the `XeNH_Simulation_Colab.ipynb` notebook
5. Run all cells in order

### Option 2: Upload Notebook

1. Download `XeNH_Simulation_Colab.ipynb` from this repository
2. Go to [Google Colab](https://colab.research.google.com/)
3. Click `File` → `Upload notebook`
4. Upload the downloaded file
5. Run all cells in order

## Configuration

Before running the simulation, you can adjust these parameters in the **Configuration** section:

```python
NUM_TRIALS = 500000  # Number of Monte Carlo trials
RANDOM_SEED = 42     # For reproducibility
SAVE_EVERY = 5000    # Save frequency
```

### Recommended Trial Counts

| Trials | Estimated Time | Use Case |
|--------|----------------|----------|
| 100,000 | ~1.5 hours | Quick test run |
| 500,000 | ~4-5 hours | Moderate search |
| 1,000,000 | ~9 hours | Thorough search |

**Note:** Google Colab has runtime limits:
- Free tier: 12 hours max per session
- Colab Pro: 24 hours max per session

## What the Simulation Does

The notebook:

1. **Installs dependencies**: NumPy, SciPy, Matplotlib, Pandas
2. **Clones the repository**: Gets the latest Python model code
3. **Runs Monte Carlo simulation**: Tests random parameter combinations
4. **Validates against real constraints**: Both Xe AND N isotope criteria must be satisfied
5. **Saves validated results**: Only parameter sets that pass both models
6. **Creates visualizations**: Distribution plots for successful parameters
7. **Enables download**: Results CSV and figures

## Understanding Results

### Success Criteria

A parameter combination is considered a "success" only if **BOTH** models pass:

**Xenon Model:**
- Final ¹³⁰Xe: 4.3×10⁵ - 9.2×10⁵
- Final ¹²⁸Xe/¹³⁰Xe: 0.475 - 0.478

**Nitrogen Model:**
- Final ¹⁴N: 7.1×10¹⁹ - 9.8×10²¹ mol/g
- Final ¹⁵N/¹⁴N: 0.0036275 - 0.0036425

### Output Files

- `success_REAL_colab.csv`: Validated parameter combinations
- `parameter_distributions_colab.png`: Visualization of successful parameters

### Parameter Descriptions

| Parameter | Description | Range |
|-----------|-------------|-------|
| **η (eta)** | Processing rate | 5.0 - 20.0 |
| **α (alpha)** | Growth rate for downwelling | 1×10⁻¹⁰ - 1×10⁻⁷ |
| **β (beta)** | Inflection point (Gyr) | 0 - 10 |
| **xe_cap** | Xenon reservoir capacity | 1×10⁵ - 1×10⁶ |
| **n_cap** | Nitrogen reservoir capacity | 1×10²⁰ - 1×10²² |
| **resfrac** | Recycling fraction | 0.01 - 0.99 |
| **lv_frac** | Late veneer fraction | 1×10⁻⁴ - 1×10⁻² |

## Tips for Google Colab

### Keeping Session Alive

Google Colab may disconnect after periods of inactivity. To prevent this:

1. Keep the browser tab active
2. Consider using Colab Pro for longer sessions
3. The notebook saves progress, so you can resume if disconnected

### Monitoring Progress

The simulation prints updates every 5% of completion:

```
Progress:   5.0% (25,000/500,000) | Successes: 2 | Rate: 31.2 iter/s | Time remaining: 4.23 hrs
Progress:  10.0% (50,000/500,000) | Successes: 5 | Rate: 30.8 iter/s | Time remaining: 4.05 hrs
...
```

Each validated success is announced immediately:
```
✓ SUCCESS #3 found at iteration 73,421!
```

### Downloading Results

The final cell allows you to download:
- `success_REAL_colab.csv`: All validated parameter combinations
- `parameter_distributions_colab.png`: Visualization of distributions

## Troubleshooting

### "Runtime disconnected"

If your session disconnects:
1. Reconnect to the runtime
2. Re-run all cells
3. The simulation will start fresh (no automatic checkpointing yet)

### "No successes found"

Success rates can be very low (< 0.01%). Solutions:
- Increase `NUM_TRIALS` (try 1M trials)
- Check that parameter ranges make sense
- Verify the model constraints are reasonable

### Memory issues

If you encounter memory errors:
- Reduce `NUM_TRIALS`
- Use Colab Pro with more RAM
- The current setup should work fine with standard Colab

## Comparison: Colab vs Local

| Aspect | Google Colab | Local Python |
|--------|--------------|--------------|
| Setup | No setup needed | Requires Python install |
| Speed | ~30 iter/s | Depends on hardware |
| Runtime | 12 hours free / 24 hours Pro | Unlimited |
| Cost | Free or $9.99/month | Free (use own hardware) |
| Accessibility | Any device with browser | Requires local machine |

## Next Steps

After obtaining validated results:

1. **Analyze distributions**: Look for clustering in parameter space
2. **Scientific interpretation**: Relate successful parameters to geophysical processes
3. **Visualization**: Create additional plots using the CSV data
4. **Compare with literature**: Check if results align with published constraints
5. **Extended runs**: Try longer simulations for more comprehensive sampling

## Support

For issues with:
- **The model itself**: See main repository README
- **Google Colab**: Visit [Colab FAQ](https://research.google.com/colaboratory/faq.html)
- **This notebook**: Open an issue on the GitHub repository

## Citation

If you use this code in published research, please cite:

```
Krantz, J. et al. (2024). XeNH: Xenon and Nitrogen Geochemical Evolution Model.
GitHub repository: https://github.com/jkrantz159/XeNH
```

## License

This notebook is part of the XeNH project and is released under the MIT License.
