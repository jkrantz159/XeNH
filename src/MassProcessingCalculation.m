%% MASSPROCESSINGCALCULATION - Calculate total mass processed over Earth history
%   This script calculates the total mass of mantle processed over Earth's
%   history based on successful model parameters. It computes minimum and
%   maximum values based on the range of eta and reservoir fractions from
%   successful model runs.
%
%   Input:
%       Reads from: results/LowLVIndXeNHsuccess.txt
%       Expected columns: eta (col 1), other params, Resfrac (col 11)
%
%   Output:
%       Displays: TotalMassProcessed / ReservoirMass ratios
%       This ratio indicates how many times the reservoir has been processed
%
%   Scientific Context:
%       Understanding mantle processing rates is crucial for constraining
%       degassing history and volatile recycling over geological time.
%
%   References:
%       - Based on He flux measurements at mid-ocean ridges
%       - Processing rate model from Parai and Mukhopadhyay (2018)

%% Load successful model parameters
try
    filepath = ModelUtils.getOutputPath('LowLVIndXeNHsuccess.txt', 'results');
    if ~exist(filepath, 'file')
        error('MassProcessing:FileNotFound', ...
              'Success file not found: %s', filepath);
    end
    successes = readmatrix(filepath);
catch
    % Fallback for older MATLAB versions or if file is in current directory
    warning('MassProcessing:FallbackLoad', ...
            'Using fallback file loading method');
    successes = csvread('LowLVIndXeNHsuccess.txt');
end

%% Extract parameters from successful runs
eta = successes(:, 1);                  % Processing rate parameter (/yr)
Resfrac = successes(:, 11);             % Reservoir fraction

% Get min and max values for sensitivity analysis
eta_min = min(eta);
eta_max = max(eta);
Resfrac_min = min(Resfrac);
Resfrac_max = max(Resfrac);

fprintf('Parameter ranges from successful models:\n');
fprintf('  eta: %.4e to %.4e /yr\n', eta_min, eta_max);
fprintf('  Reservoir fraction: %.4f to %.4f\n', Resfrac_min, Resfrac_max);

%% Define time vector and constants
% Non-uniform time spacing: finer resolution in early Earth
t = [0:0.1E6:200E6, ...                 % Early Earth: 0-200 Ma, 0.1 Ma steps
     201E6:1E6:3300E6, ...              % Middle period: 200-3300 Ma, 1 Ma steps
     3305E6:5E6:4.565E9, ...            % Late period: 3300-4565 Ma, 5 Ma steps
     4.568E9];                          % Present day

% Physical constants
T = ModelConfig.EARTH_AGE_YEARS;        % 4.568 Ga
Qp = ModelConfig.PRESENT_DAY_PROCESSING_RATE;  % 6.1E17 g/yr
ME = ModelConfig.EARTH_MASS_GRAMS;      % 5.972E27 g

% Calculate reservoir masses
Mres_min = ME * Resfrac_min;
Mres_max = ME * Resfrac_max;

%% Calculate total mass processed - Minimum case (eta_min)
TotalMassProcessed_Min = 0;

for i = 2:length(t)
    % Calculate mass processed in this time step
    dM = ModelUtils.calculateMassProcessed(eta_min, t(i), t(i-1), T, Qp);
    TotalMassProcessed_Min = TotalMassProcessed_Min + dM;
end

%% Calculate total mass processed - Maximum case (eta_max)
TotalMassProcessed_Max = 0;

for i = 2:length(t)
    % Calculate mass processed in this time step
    dM = ModelUtils.calculateMassProcessed(eta_max, t(i), t(i-1), T, Qp);
    TotalMassProcessed_Max = TotalMassProcessed_Max + dM;
end

%% Display results
fprintf('\nMass processing results:\n');
fprintf('  Total mass processed (min): %.4e g\n', TotalMassProcessed_Min);
fprintf('  Total mass processed (max): %.4e g\n', TotalMassProcessed_Max);
fprintf('  Reservoir mass (min): %.4e g\n', Mres_min);
fprintf('  Reservoir mass (max): %.4e g\n', Mres_max);
fprintf('\nNumber of times reservoir processed:\n');
fprintf('  Minimum case: %.2f\n', TotalMassProcessed_Min / Mres_min);
fprintf('  Maximum case: %.2f\n', TotalMassProcessed_Max / Mres_max);

% Also display as direct output for backward compatibility
TotalMassProcessed_Min / Mres_min
TotalMassProcessed_Max / Mres_max
