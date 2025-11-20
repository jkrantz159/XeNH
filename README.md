# XeNH: Coupled H-N-Xe Recycling Research Code

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

MATLAB code for modeling coupled H-N-Xe recycling in Earth's mantle, tracing mantle regassing and the onset of subduction using observations of H, N, and Xe isotopes.

## Publication

This repository contains the code used to produce:

> **"Tracing Mantle Regassing and the Onset of Subduction Using Observations of H, N, and Xe"**
>
> John A. Krantz¹, Peter H. Barry², Stephen W. Parman¹
>
> 1. Department of Earth, Environmental, and Planetary Sciences, Brown University, Providence, RI
> 2. Woods Hole Oceanographic Institution, Woods Hole, MA
>
> Prepared for submission to *Earth and Planetary Science Letters*, May 2020

## Overview

This code implements a Monte Carlo simulation framework to explore parameter space for mantle degassing and volatile recycling models. The models use a box model approach with sigmoidal growth for subduction onset, testing against observational constraints from noble gas and nitrogen isotope systematics.

### Key Features

- **Monte Carlo parameter exploration**: Large-scale random sampling of parameter space
- **Multi-isotope constraints**: Simultaneous fitting to Xe, N, and H systematics
- **Parallel processing**: Efficient computation using MATLAB's Parallel Computing Toolbox
- **Configurable parameters**: Centralized configuration for all physical constants
- **Comprehensive visualization**: Kernel density estimation plots for parameter distributions

## Installation

### Requirements

- **MATLAB** R2018b or later (recommended: R2020a+)
- **Required Toolboxes**:
  - Parallel Computing Toolbox
  - Statistics and Machine Learning Toolbox

### Setup

1. Clone this repository:
   ```bash
   git clone https://github.com/jkrantz159/XeNH.git
   cd XeNH
   ```

2. Add the source directory to your MATLAB path:
   ```matlab
   addpath('src')
   ```

3. Verify installation by running the test script:
   ```matlab
   cd tests
   run_basic_tests
   ```

## Quick Start

### Running a Basic Simulation

```matlab
% Navigate to src directory
cd src

% Run a small test simulation (1 million trials)
% Edit ParallelCombinedRandomModel.m and set runs = 1E6
ParallelCombinedRandomModel

% Visualize results
VariableKDEs
```

### Understanding the Output

Results are saved to `results/success.txt` with the following columns:
1. **alpha**: Growth rate parameter (/yr)
2. **beta**: Sigmoid inflection point (years)
3. **eta**: Processing rate parameter (/yr)
4. **Xed**: Xe carrying capacity (atoms/gram)
5. **Nd**: N carrying capacity (atoms/gram)
6. **Resfrac**: Reservoir fraction
7. **LVfrac**: Late veneer fraction (%)

## Directory Structure

```
XeNH/
├── src/                    # MATLAB source code
│   ├── ModelConfig.m       # Configuration and constants
│   ├── ModelUtils.m        # Utility functions
│   ├── BaseGeochemicalModel.m  # Base model class
│   ├── ParallelXeModel.m   # Xenon isotope model
│   ├── ParallelNewNModel.m # Nitrogen isotope model
│   ├── ParallelNeModel.m   # Neon isotope model
│   ├── ParallelCombinedRandomModel.m  # Main Monte Carlo driver
│   ├── VariableKDEs.m      # Visualization script
│   ├── MassProcessingCalculation.m  # Post-processing analysis
│   └── parsave.m           # Parallel-safe file saving
├── data/                   # Input data (if needed)
├── results/                # Output CSV files
├── figures/                # Generated plots
├── tests/                  # Test scripts
├── docs/                   # Additional documentation
├── README.md               # This file
├── LICENSE                 # MIT License
└── CITATION.cff            # Citation information
```

## Usage Guide

### 1. Configuration

Edit `src/ModelConfig.m` to modify physical constants, parameter ranges, or success criteria. All model parameters are centralized in this file.

### 2. Running Simulations

#### Full Monte Carlo Simulation (100 million trials)
```matlab
cd src
ParallelCombinedRandomModel  % Default: 1E8 trials
```

**Note**: This may take several hours to days depending on your hardware and parallel pool size.

#### Test Simulation (smaller sample)
```matlab
% Edit ParallelCombinedRandomModel.m
% Change line: runs = 1E8;
% To: runs = 1E6;  % 1 million trials for testing
```

#### Individual Model Testing
```matlab
% Test Xenon model alone
t = ModelConfig.createTimeVector();
atm = ModelConfig.createAtmosphereEvolution(t);
deltaXe = @ModelConfig.ratioToDeltaXe;
XeSucc = ParallelXeModel(1e6, 1e-9, 3e9, 7.5e-10, 0.9, 1, ...
                         1, atm, t, 4.568e9, deltaXe, 1);

% Test Nitrogen model alone
deltaN = @ModelConfig.ratioToDeltaN;
NSucc = ParallelNewNModel(1e15, 1e-9, 3e9, 7.5e-10, 0.9, 1, ...
                          1, t, 4.568e9, deltaN, 1);
```

### 3. Analyzing Results

#### Visualize Parameter Distributions
```matlab
cd src
VariableKDEs  % Creates KDE plots of successful parameters
```

#### Calculate Total Mass Processed
```matlab
cd src
MassProcessingCalculation  % Analyzes mass flux through time
```

## Model Description

### Scientific Framework

The models are based on:
- **Box model approach**: Mantle reservoir with time-varying input/output
- **Exponential degassing**: Q(t) = Q_p × exp(η(T-t))
- **Sigmoidal downwelling growth**: Capacity/(1 + exp(-α(t-β)))
- **Late veneer initial conditions**: AVCC composition (1% Earth mass default)

### Parameters

| Parameter | Symbol | Range | Units | Description |
|-----------|--------|-------|-------|-------------|
| Growth rate | α | 10⁻¹⁰ to 10⁻⁷ | /yr | Controls downwelling growth rate |
| Inflection point | β | 0 to 10 | Gyr | Time of maximum downwelling growth |
| Processing rate | η | 7.5×10⁻¹⁰ to 8×10⁻¹⁰ | /yr | Mantle degassing rate parameter |
| Xe capacity | X_d | 5 to 5×10⁸ | atoms/g | Xenon downwelling concentration |
| N capacity | N_d | 4 to 4×10¹⁷ | atoms/g | Nitrogen downwelling concentration |

### Success Criteria

#### Xenon
- ¹³⁰Xe concentration: 4.3×10⁵ to 9.2×10⁵ atoms/gram
- ¹²⁸Xe/¹³⁰Xe ratio: 0.475 to 0.478

#### Nitrogen
- ¹⁴N concentration: 7.06×10¹⁹ to 9.78×10²¹ mol/g (converted to atoms/gram)
- ¹⁵N/¹⁴N ratio: 0.0036275 to 0.0036425 (δ¹⁵N ~ -5‰)

## Performance Optimization

### Parallel Processing

The code uses MATLAB's `parfor` for parallel execution. To optimize performance:

```matlab
% Start a parallel pool before running simulations
parpool('local', 4);  % Use 4 workers (adjust based on CPU cores)

% Run simulation
ParallelCombinedRandomModel

% Close pool when done
delete(gcp('nocreate'));
```

### Memory Management

For very large simulations (>10⁸ trials), results are streamed to disk to avoid memory issues.

## Troubleshooting

### Common Issues

1. **"Undefined function or variable 'ModelConfig'"**
   - Solution: Ensure `src/` is in your MATLAB path: `addpath('src')`

2. **"Parallel pool not available"**
   - Solution: The code will run serially if Parallel Computing Toolbox is unavailable, but will be slower

3. **"Out of memory" errors**
   - Solution: Reduce `runs` parameter or increase system RAM

4. **Figures not saving**
   - Solution: Ensure `figures/` directory exists and has write permissions

### Getting Help

- Open an issue on [GitHub](https://github.com/jkrantz159/XeNH/issues)

## Testing

Run the test suite to verify installation:

```matlab
cd tests
run_basic_tests     % Basic functionality tests
run_model_tests     % Individual model tests
```

## References

### Key Publications

- **Parai, R., & Mukhopadhyay, S. (2018)**. Xenon isotopic constraints on the history of volatile recycling into the mantle. *Nature*, 560(7717), 223-227.

- **Barry, P. H., & Hilton, D. R. (2016)**. Release of subducted sedimentary nitrogen throughout Earth's mantle. *Geochemistry, Geophysics, Geosystems*, 17(5), 1762-1773.

- **Marty, B. (2012)**. The origins and concentrations of water, carbon, nitrogen and noble gases on Earth. *Earth and Planetary Science Letters*, 313, 56-66.

- **Pepin, R. O. (2000)**. On the isotopic composition of primordial xenon in terrestrial planet atmospheres. *Space Science Reviews*, 92(3), 371-395.

### Data Sources

- **Porcelli, D., Ballentine, C. J., & Wieler, R. (2002)**. An overview of noble gas geochemistry and cosmochemistry. *Reviews in Mineralogy and Geochemistry*, 47(1), 1-19.

- **Sephton, M. A., et al. (2003)**. High molecular weight organic matter in martian meteorites. *Planetary and Space Science*, 51(6), 363-369.

- **Williams, C. D., & Mukhopadhyay, S. (2018)**. Capture of nebular gases during Earth's accretion is preserved in deep-mantle neon. *Nature*, 565(7737), 78-81.

- **Owen, T., et al. (2001)**. Protosolar nitrogen. *The Astrophysical Journal*, 553(1), L77.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this code in your research, please cite both the software and the associated paper.

See [CITATION.cff](CITATION.cff) for citation information in standard format.

## Contact

For questions about the code or methodology, please open an issue on GitHub.
