%% PARALLELCOMBINEDRANDOMMODEL - Monte Carlo simulation of coupled H-N-Xe recycling
%   This script runs a large-scale Monte Carlo simulation to find parameter
%   combinations that simultaneously satisfy constraints from Xenon and Nitrogen
%   isotope systematics in Earth's mantle.
%
%   Model Framework:
%       Based on Parai and Mukhopadhyay (2018) with nitrogen systematics from
%       Barry and Hilton (2016). Uses sigmoidal growth model for subduction
%       onset combined with exponential mantle degassing.
%
%   Parameters Explored:
%       - alpha: Growth rate (1E-10 to 1E-7 /yr)
%       - beta: Sigmoid inflection point (0 to 10 Gyr)
%       - eta: Processing rate (7.5E-10 to 8E-10 /yr)
%       - XeCap: Xe carrying capacity (5 to 5E8 atoms/gram)
%       - NCap: N carrying capacity (4 to 4E17 atoms/gram)
%
%   Fixed Parameters:
%       - Resfrac: Reservoir fraction = 0.9 (90% convecting mantle)
%       - LVfrac: Late veneer fraction = 1% of Earth mass (AVCC composition)
%
%   Output:
%       Successful parameter combinations saved to results/success.txt
%       Format: alpha, beta, eta, Xed, Nd, Resfrac, LVfrac
%
%   Usage:
%       Run this script directly. Modify 'runs' variable for different
%       sample sizes (default: 1E8 Monte Carlo trials).
%       Requires Parallel Computing Toolbox for parfor loops.
%
%   References:
%       - Parai and Mukhopadhyay (2018): Model framework
%       - Barry and Hilton (2016): Nitrogen isotope systematics
%       - Porcelli, Ballentine, and Wieler (2002): Atmospheric composition

%% Configuration and Setup
fprintf('=================================================================\n');
fprintf('  Parallel Combined Random Model - H-N-Xe Recycling Simulation\n');
fprintf('=================================================================\n');
fprintf('Start time: %s\n\n', datestr(now));

% Set random seed for reproducibility (optional - uncomment to enable)
% rng(12345, 'twister');

%% Model Parameters
Resfrac = ModelConfig.CONVECTING_MANTLE_FRACTION;  % 0.9 = 90% convecting mantle
LVfrac = ModelConfig.LATE_VENEER_FRACTION_DEFAULT;  % 1% Late veneer

% Number of Monte Carlo trials
runs = 1E8;  % 100 million trials (reduce for testing, e.g., 1E6)

fprintf('Configuration:\n');
fprintf('  Reservoir fraction: %.1f%%\n', Resfrac * 100);
fprintf('  Late veneer fraction: %.1f%% Earth mass\n', LVfrac);
fprintf('  Monte Carlo trials: %.0e\n', runs);
fprintf('  Parallel processing: %s\n', ...
        iif(isempty(gcp('nocreate')), 'Not initialized', 'Active'));
fprintf('\n');

%% Create Time Vector and Atmospheric Evolution
t = ModelConfig.createTimeVector();
T = ModelConfig.EARTH_AGE_YEARS;

% Atmospheric 128Xe/130Xe evolution (linear from fractionated to modern over 2.568 Ga)
atm = ModelConfig.createAtmosphereEvolution(t);

% Delta notation conversion functions
deltaXe = @ModelConfig.ratioToDeltaXe;
deltaN = @ModelConfig.ratioToDeltaN;

fprintf('Time integration:\n');
fprintf('  Time steps: %d\n', length(t));
fprintf('  Integration period: 0 to %.3f Ga\n\n', T / 1E9);

%% Monte Carlo Simulation
fprintf('Starting Monte Carlo simulation...\n');
fprintf('Progress will be reported every 1%% (%.0e iterations)\n\n', runs / 100);

% Initialize counters
successCount = 0;

% Parallel loop over random parameter combinations
parfor count = 1:runs
    % Generate random parameters from uniform distributions
    % Note: Variable names are more descriptive than original i,j,m,n,p

    % Growth rate alpha: 1E-10 to 1E-7 /yr (log-uniform distribution)
    alpha = 10 .^ (-10 + (-7 + 10) .* rand(1));

    % Inflection point beta: 0 to 10 Gyr (uniform distribution)
    beta = (10 .* rand(1)) * 1E9;

    % Processing rate eta: 7.5E-10 to 8E-10 /yr (uniform distribution)
    eta = ModelConfig.ETA_MIN + rand(1) * (ModelConfig.ETA_MAX - ModelConfig.ETA_MIN);

    % N carrying capacity: 4 to 4E17 atoms/gram (log-uniform distribution)
    NCap = 4 * 10 .^ (17 .* rand(1));

    % Xe carrying capacity: 5 to 5E8 atoms/gram (log-uniform distribution)
    XeCap = 5 * 10 .^ (8 .* rand(1));

    % Test Nitrogen model
    NSucc = ParallelNewNModel(NCap, alpha, beta, eta, Resfrac, LVfrac, ...
                              0, t, T, deltaN, count);

    % Test Xenon model
    XeSucc = ParallelXeModel(XeCap, alpha, beta, eta, Resfrac, LVfrac, ...
                             0, atm, t, T, deltaXe, count);

    % Check if both models succeed
    if NSucc == 1 && XeSucc == 1
        % Re-run models with plotting enabled to visualize successful case
        ParallelXeModel(XeCap, alpha, beta, eta, Resfrac, LVfrac, ...
                        1, atm, t, T, deltaXe, count);
        ParallelNewNModel(NCap, alpha, beta, eta, Resfrac, LVfrac, ...
                          1, t, T, deltaN, count);

        % Report success
        fprintf('SUCCESS #%d: alpha=%.2e, beta=%.2e Ga, eta=%.2e, Xe=%.2e, N=%.2e\n', ...
                count, alpha, beta / 1E9, eta, XeCap, NCap);

        % Save parameters to file (parallel-safe)
        parsave(alpha, beta, eta, XeCap, NCap, Resfrac, LVfrac);
    end

    % Progress reporting (every 1 million iterations)
    if mod(count, 1E6) == 0
        fprintf('Progress: %.1f%% (%d/%d iterations)\n', ...
                (count / runs) * 100, count, runs);
    end
end

%% Summary
fprintf('\n=================================================================\n');
fprintf('  Simulation Complete\n');
fprintf('=================================================================\n');
fprintf('End time: %s\n', datestr(now));
fprintf('\nResults saved to: results/success.txt\n');
fprintf('Review successful parameter combinations for further analysis.\n');

%% Helper function for inline if-else (MATLAB doesn't have ternary operator)
function result = iif(condition, trueVal, falseVal)
    if condition
        result = trueVal;
    else
        result = falseVal;
    end
end
