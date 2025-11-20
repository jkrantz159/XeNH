function [NSucc] = ParallelNewNModel(NCap, alpha, beta, eta, Resfrac, LVfrac, PLOTCHECK, t, T, deltaN, count)
    %PARALLELNEWNMODEL Nitrogen isotope evolution model for Earth's mantle
    %   Models the evolution of 14N and 15N in Earth's mantle over
    %   geological time using a box model with sigmoidal downwelling growth.
    %
    %   Based on: Parai and Mukhopadhyay (2018) and Barry and Hilton (2016)
    %
    %   Parameters:
    %       NCap      - N carrying capacity/present day downwelling (atoms/gram)
    %                   Range: 4 to 4E17 atoms/gram
    %       alpha     - Growth rate parameter (/yr)
    %                   Range: 1E-10 to 1E-7 /yr
    %       beta      - Sigmoid inflection point (years)
    %                   Range: 0 to 10 Gyr
    %       eta       - Processing rate parameter (/yr)
    %                   Range: 7E-10 to 8E-10 /yr
    %       Resfrac   - Reservoir fraction (0.9 = 90% convecting mantle)
    %       LVfrac    - Late veneer fraction (% of Earth mass)
    %       PLOTCHECK - 1 to generate plots, 0 to skip plotting
    %       t         - Time vector (years)
    %       T         - Total age of Earth (years)
    %       deltaN    - Function handle to convert ratio to delta notation
    %       count     - Iteration counter for file naming
    %
    %   Returns:
    %       NSucc - 1 if model meets success criteria, 0 otherwise
    %
    %   Success Criteria:
    %       1. 14N concentration: 7.06E19 to 9.78E21 mol/g mantle * conversion
    %       2. 15N/14N ratio: 0.0036275 to 0.0036425 (delta-15N ~ -5 per mille)
    %
    %   Initial Conditions:
    %       - AVCC (Average Carbonaceous Chondrite) composition
    %       - N concentration from Sephton et al. (2003): 0.1519 wt%
    %       - 15N/14N ratio from Owen et al. (2001): 2.3E-3 (PSN)
    %       - Sediment 15N/14N = 0.003671565 (delta-15N = +5 per mille)
    %
    %   Example:
    %       t = ModelConfig.createTimeVector();
    %       deltaN = @ModelConfig.ratioToDeltaN;
    %       success = ParallelNewNModel(1e15, 1e-9, 3e9, 7.5e-10, 0.9, 1, ...
    %                                    0, t, 4.568e9, deltaN, 1);

    %% Validate inputs
    try
        ModelUtils.validateInputs(...
            'alpha', alpha, ModelConfig.ALPHA_MIN * 0.001, ModelConfig.ALPHA_MAX * 10, ...
            'beta', beta, 0, ModelConfig.BETA_MAX * 2, ...
            'eta', eta, ModelConfig.ETA_MIN * 0.1, ModelConfig.ETA_MAX * 10, ...
            'Resfrac', Resfrac, 0.1, 1.0, ...
            'LVfrac', LVfrac, 0, 10);
    catch ME
        warning('ParallelNewNModel:ValidationFailed', ...
                'Input validation failed: %s. Returning failure.', ME.message);
        NSucc = 0;
        return;
    end

    %% Initialize arrays and constants
    Nd = zeros(1, length(t));           % Downwelling N concentration
    N14M = zeros(1, length(t));         % 14N in mantle
    N15M = zeros(1, length(t));         % 15N in mantle

    % Physical constants
    Mres = ModelConfig.EARTH_MASS_GRAMS * Resfrac;  % Reservoir mass (g)
    ME = ModelConfig.EARTH_MASS_GRAMS;              % Earth mass (g)

    %% Calculate initial mantle composition
    % Late veneer of carbonaceous chondrites (AVCC: Murchison and Orgueil)
    % Concentration of N in CC = 0.1519 wt% (Sephton et al. 2003)
    g14M = (ModelConfig.N_CONCENTRATION_CC * ME * LVfrac / 100);  % g of N in mantle
    Initial14 = (g14M / ModelConfig.N_ATOMIC_MASS_14 * ModelConfig.N14_FRACTION * ...
                ModelConfig.AVOGADRO_NUMBER) / Mres;  % atoms/gram 14N
    Initial15 = Initial14 * ModelConfig.N15_14_RATIO_INITIAL;  % PSN ratio

    %% Calculate sigmoidal downwelling evolution
    Nd = ModelUtils.calculateSigmoidalDownwelling(NCap, alpha, beta, t);

    %% Time evolution loop - Box model integration
    for i = 2:length(t)
        % Calculate mass of mantle processed in this time step
        dM = ModelUtils.calculateMassProcessed(...
            eta, t(i), t(i-1), T, ModelConfig.PRESENT_DAY_PROCESSING_RATE);

        % Normalized mass processed (fraction of reservoir)
        dM_Mres = dM / Mres;

        % Get values from previous time step
        if i == 2
            N14Mlast = Initial14;
            N15Mlast = Initial15;
            Ndlast = Nd(i);
        else
            N14Mlast = N14M(i-1);
            N15Mlast = N15M(i-1);
            Ndlast = Nd(i-1);
        end

        % Update concentrations using box model equations
        % 14N: downwelling with constant composition
        N14M(i) = ModelUtils.updateConcentration(N14Mlast, dM_Mres, Ndlast, 1);

        % 15N: downwelling with sediment isotopic composition (delta-15N = +5)
        N15M(i) = ModelUtils.updateConcentration(...
            N15Mlast, dM_Mres, Ndlast, ModelConfig.N15_14_RATIO_SEDIMENT);
    end

    %% Check success criteria
    % Criterion 1: 14N concentration in present-day mantle
    % Convert from mol/g to atoms/gram reservoir
    N14_min = ModelConfig.N14_MANTLE_MIN_MOL * ModelConfig.N14_FRACTION * ...
              ModelConfig.AVOGADRO_NUMBER / Mres;
    N14_max = ModelConfig.N14_MANTLE_MAX_MOL * ModelConfig.N14_FRACTION * ...
              ModelConfig.AVOGADRO_NUMBER / Mres;
    Succ1 = ModelUtils.checkSuccessCriteria(N14M(end), N14_min, N14_max);

    % Criterion 2: 15N/14N ratio in present-day mantle
    ratio = N15M(end) / N14M(end);
    Succ2 = ModelUtils.checkSuccessCriteria(...
        ratio, ModelConfig.N15_14_RATIO_MANTLE_MIN, ModelConfig.N15_14_RATIO_MANTLE_MAX);

    % Overall success
    NSucc = Succ1 && Succ2;

    %% Generate plots if requested
    if PLOTCHECK == 1
        figHandle = ModelUtils.createStandardFigure(...
            1, 'Nitrogen Isotope Evolution', 'Time (Myr)', '\delta ^{15}N (‰)');
        plot(t / 1E6, deltaN(N15M ./ N14M), '-k', 'LineWidth', 1.5);

        % Save figure to figures directory
        filename = sprintf('N_evolution_%d.jpg', count);
        ModelUtils.saveFigureSafely(figHandle, filename, 'figures');
    end
end
