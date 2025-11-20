# XeNH Python Implementation

This is a Python port of the MATLAB XeNH code for modeling coupled H-N-Xe recycling in Earth's mantle.

## Requirements

```bash
pip install numpy scipy matplotlib pandas
```

## Quick Start

```python
# Run a 1 million trial simulation
cd python/src
python run_simulation.py
```

## ⚠️ Important Note About Example Data

The data in `results/success_example.csv` and figures in `figures/` are **SYNTHETIC** for demonstration purposes. They show expected output format but are NOT validated successes.

To generate real validated successes, run the simulation (see below).

## Files

- `model_config.py`: Configuration and constants
- `model_utils.py`: Utility functions
- `parallel_xe_model.py`: Xenon isotope model
- `parallel_n_model.py`: Nitrogen isotope model
- `run_simulation.py`: Main Monte Carlo simulation (finds REAL successes)

## Differences from MATLAB Version

- Uses NumPy instead of MATLAB arrays
- Single-threaded (no parallel processing yet)
- CSV output instead of txt
- Simplified visualization

## Performance

The Python version runs approximately 50-100 iterations per second on a single core.
For 1 million trials, expect ~3-5 hours of runtime.

## Future Improvements

- Add multiprocessing for parallel execution
- Add visualization module (matplotlib)
- Add Neon model
- Add KDE visualization
