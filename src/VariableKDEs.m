%% VARIABLEKDES - Kernel Density Estimation visualization of successful parameters
%   This script creates kernel density estimate (KDE) plots for successful
%   parameter combinations from Monte Carlo simulations. It visualizes the
%   distributions of model parameters that successfully match observational
%   constraints.
%
%   Input Files:
%       - results/success.txt: Xe-N combined success (LV=1%)
%       - results/IndXeNHsuccess.txt: Independent Xe-N-H success (LV=1%)
%       - results/LowLVIndXeNHsuccess.txt: Low late veneer case (LV=0.05%)
%
%   Output:
%       - KDE plots for processing rate (eta)
%       - KDE plots for inflection point (beta)
%       - KDE plots for growth rate (alpha)
%       - KDE plots for carrying capacities (Xe, N, H)
%
%   Figures:
%       Figure 1-4: Parameter distributions for different scenarios
%       Figure 11-14: Combined visualizations
%
%   Requirements:
%       - Statistics and Machine Learning Toolbox (for fitdist)
%
%   Usage:
%       Run this script after ParallelCombinedRandomModel.m has generated
%       success files. Modify file paths as needed for different datasets.

%% Load and Process Xe-N Combined Data (LV=1%)
try
    filepath = ModelUtils.getOutputPath('success.txt', 'results');
    if ~exist(filepath, 'file')
        warning('VariableKDEs:FileNotFound', ...
                'File not found: %s. Skipping this dataset.', filepath);
    else
        successes = readmatrix(filepath);

        % Extract parameters (columns: alpha, beta, eta, Xed, Nd, Resfrac, LVfrac)
        alpha = successes(:, 1);
        beta = successes(:, 2);
        eta = successes(:, 3);
        Xe = successes(:, 4);
        N = successes(:, 5);

        fprintf('Xe-N Combined Data (LV=1%%):\n');
        fprintf('  Success count: %d\n\n', size(successes, 1));

        % Plot Processing Rate (eta)
        plotKDE(14, eta * 1E10, 'Kernel', 0.003, 7.6:0.001:8, ...
                'Processing Rate KDE', '\eta \times 10^{10}', 'KDE', ...
                'k-', [6 8 0 50]);

        % Plot Inflection Point (beta)
        plotKDE(11, beta / 1E9, 'Kernel', 0.15, 0:0.1:10, ...
                'Inflection Point KDE', 'Time (Gyr)', 'KDE', 'k-', []);

        % Plot Growth Rate (alpha)
        plotKDE(12, log(alpha), 'Kernel', 0.2, -22:0.1:-15, ...
                'Growth Rate KDE', 'Log \alpha', 'KDE', 'k-', []);

        % Plot Carrying Capacities (Xe and N)
        figure(13);
        hold on;
        plotKDEOverlay(log10(Xe), 'Kernel', 0.2, 5:0.1:25, 'r-', 2);
        plotKDEOverlay(log10(N), 'Kernel', 0.2, 5:0.1:25, 'g-', 2);
        title('Recycling KDE', 'FontSize', 12, 'FontWeight', 'bold');
        xlabel('log_{10} Recycling', 'FontSize', 11);
        ylabel('KDE', 'FontSize', 11);
        legend('Xe', 'N', 'Location', 'best');
        grid on;
        hold off;

        % Save figures
        saveFigureSet([11, 12, 13, 14], 'XeN_Combined');
    end
catch ME
    warning('VariableKDEs:ProcessingError', ...
            'Error processing Xe-N combined data: %s', ME.message);
end

%% Load and Process Independent Xe-N-H Data (LV=1%)
try
    filepath = ModelUtils.getOutputPath('IndXeNHsuccess.txt', 'results');
    if ~exist(filepath, 'file')
        warning('VariableKDEs:FileNotFound', ...
                'File not found: %s. Skipping this dataset.', filepath);
    else
        successes = readmatrix(filepath);

        % Extract parameters (eta, Hd, Ha, Hb, Nd, Na, Nb, Xd, Xa, Xb)
        eta = successes(:, 1);
        Hd = successes(:, 2);
        Ha = successes(:, 3);
        Hb = successes(:, 4);
        Nd = successes(:, 5);
        Na = successes(:, 6);
        Nb = successes(:, 7);
        Xd = successes(:, 8);
        Xa = successes(:, 9);
        Xb = successes(:, 10);

        fprintf('Independent Xe-N-H Data (LV=1%%):\n');
        fprintf('  Success count: %d\n\n', size(successes, 1));

        % Plot Processing Rate
        plotKDE(4, eta * 1E10, 'Kernel', 0.003, 7.6:0.001:8, ...
                'Processing Rate KDE', '\eta \times 10^{10}', 'KDE', ...
                'k-', [6 8 0 50]);

        % Plot Inflection Points (Xe, N, H)
        figure(1);
        hold on;
        plotKDEOverlay(Xb / 1E9, 'Kernel', 0.15, 0:0.1:10, 'r-', 2);
        plotKDEOverlay(Nb / 1E9, 'Kernel', 0.15, 0:0.1:10, 'g-', 2);
        plotKDEOverlay(Hb / 1E9, 'Kernel', 0.15, 0:0.1:10, 'b-', 2);
        title('Inflection Point KDE', 'FontSize', 12, 'FontWeight', 'bold');
        xlabel('Time (Gyr)', 'FontSize', 11);
        ylabel('KDE', 'FontSize', 11);
        legend('Xe', 'N', 'H', 'Location', 'best');
        grid on;
        hold off;

        % Plot Growth Rates (Xe, N, H)
        figure(2);
        hold on;
        plotKDEOverlay(log(Xa), 'Kernel', 0.2, -25:0.1:-15, 'r-', 2);
        plotKDEOverlay(log(Na), 'Kernel', 0.2, -25:0.1:-15, 'g-', 2);
        plotKDEOverlay(log(Ha), 'Kernel', 0.2, -25:0.1:-15, 'b-', 2);
        title('Growth Rate KDE', 'FontSize', 12, 'FontWeight', 'bold');
        xlabel('Log \alpha', 'FontSize', 11);
        ylabel('KDE', 'FontSize', 11);
        legend('Xe', 'N', 'H', 'Location', 'best');
        grid on;
        hold off;

        % Plot Carrying Capacities (Xe, N, H)
        figure(3);
        hold on;
        plotKDEOverlay(log10(Xd), 'Kernel', 0.2, 5:0.1:25, 'r-', 2);
        plotKDEOverlay(log10(Nd), 'Kernel', 0.2, 5:0.1:25, 'g-', 2);
        plotKDEOverlay(log10(Hd), 'Kernel', 0.2, 5:0.1:25, 'b-', 2);
        title('Recycling KDE', 'FontSize', 12, 'FontWeight', 'bold');
        xlabel('log_{10} Recycling', 'FontSize', 11);
        ylabel('KDE', 'FontSize', 11);
        legend('Xe', 'N', 'H', 'Location', 'best');
        grid on;
        hold off;

        % Save figures
        saveFigureSet([1, 2, 3, 4], 'IndXeNH_LV1');
    end
catch ME
    warning('VariableKDEs:ProcessingError', ...
            'Error processing independent Xe-N-H data: %s', ME.message);
end

%% Load and Process Low LV Data (LV=0.05%)
try
    filepath = ModelUtils.getOutputPath('LowLVIndXeNHsuccess.txt', 'results');
    if ~exist(filepath, 'file')
        warning('VariableKDEs:FileNotFound', ...
                'File not found: %s. Skipping this dataset.', filepath);
    else
        successes = readmatrix(filepath);

        % Extract parameters
        eta = successes(:, 1);
        Hd = successes(:, 2);
        Ha = successes(:, 3);
        Hb = successes(:, 4);
        Nd = successes(:, 5);
        Na = successes(:, 6);
        Nb = successes(:, 7);
        Xd = successes(:, 8);
        Xa = successes(:, 9);
        Xb = successes(:, 10);

        fprintf('Low Late Veneer Data (LV=0.05%%):\n');
        fprintf('  Success count: %d\n\n', size(successes, 1));

        % Overlay on existing figures with dashed lines
        figure(4);
        hold on;
        plotKDEOverlay(eta * 1E10, 'Kernel', 0.003, 6.25:0.001:6.75, 'k--', 2);
        hold off;

        figure(1);
        hold on;
        plotKDEOverlay(Xb / 1E9, 'Kernel', 0.15, 0:0.1:10, 'r--', 2);
        plotKDEOverlay(Nb / 1E9, 'Kernel', 0.15, 0:0.1:10, 'g--', 2);
        plotKDEOverlay(Hb / 1E9, 'Kernel', 0.15, 0:0.1:10, 'b--', 2);
        hold off;

        figure(2);
        hold on;
        plotKDEOverlay(log(Xa), 'Kernel', 0.2, -22:0.1:-15, 'r--', 2);
        plotKDEOverlay(log(Na), 'Kernel', 0.2, -25:0.1:-15, 'g--', 2);
        plotKDEOverlay(log(Ha), 'Kernel', 0.2, -25:0.1:-15, 'b--', 2);
        hold off;

        figure(3);
        hold on;
        plotKDEOverlay(log10(Xd), 'Kernel', 0.2, 5:0.1:25, 'r--', 2);
        plotKDEOverlay(log10(Nd), 'Kernel', 0.2, 5:0.1:25, 'g--', 2);
        plotKDEOverlay(log10(Hd), 'Kernel', 0.2, 5:0.1:25, 'b--', 2);
        hold off;

        % Update legends to include LV comparison
        figure(1);
        legend('Xe (LV=1%)', 'N (LV=1%)', 'H (LV=1%)', ...
               'Xe (LV=0.05%)', 'N (LV=0.05%)', 'H (LV=0.05%)', ...
               'Location', 'best');

        % Save updated figures
        saveFigureSet([1, 2, 3, 4], 'IndXeNH_Comparison');
    end
catch ME
    warning('VariableKDEs:ProcessingError', ...
            'Error processing low LV data: %s', ME.message);
end

fprintf('\nVisualization complete. Figures saved to figures/ directory.\n');

%% Helper Functions

function plotKDE(figNum, data, distType, bandwidth, xRange, titleStr, xlabelStr, ylabelStr, lineSpec, axisLimits)
    %PLOTKDE Create a kernel density estimate plot
    figure(figNum);
    hold on;
    pd = fitdist(data, distType, 'BandWidth', bandwidth);
    y = pdf(pd, xRange);
    plot(xRange, y, lineSpec, 'LineWidth', 2);
    title(titleStr, 'FontSize', 12, 'FontWeight', 'bold');
    xlabel(xlabelStr, 'FontSize', 11);
    ylabel(ylabelStr, 'FontSize', 11);
    if ~isempty(axisLimits)
        axis(axisLimits);
    end
    grid on;
    hold off;
end

function plotKDEOverlay(data, distType, bandwidth, xRange, lineSpec, lineWidth)
    %PLOTKDEOVERLAY Add a KDE curve to the current figure
    pd = fitdist(data, distType, 'BandWidth', bandwidth);
    y = pdf(pd, xRange);
    plot(xRange, y, lineSpec, 'LineWidth', lineWidth);
end

function saveFigureSet(figNums, prefix)
    %SAVEFIGURESET Save a set of figures to the figures directory
    for figNum = figNums
        try
            filename = sprintf('%s_fig%d.png', prefix, figNum);
            filepath = ModelUtils.getOutputPath(filename, 'figures');
            saveas(figure(figNum), filepath);
            fprintf('Saved: %s\n', filename);
        catch ME
            warning('VariableKDEs:SaveError', ...
                    'Failed to save figure %d: %s', figNum, ME.message);
        end
    end
end
