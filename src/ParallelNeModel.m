function [NeSucc] = ParallelNeModel(NeCap, alpha, beta, eta, Resfrac, LVfrac, PLOTCHECK, t, T, count)
    %PARALLELNEMODEL Neon isotope evolution model for Earth's mantle
    %   Models the evolution of 20Ne and 22Ne in Earth's mantle over
    %   geological time using a box model with sigmoidal downwelling growth.
    %
    %   Based on: Parai and Mukhopadhyay (2018)
    %
    %   Parameters:
    %       NeCap     - Ne carrying capacity/present day downwelling (atoms/gram)
    %                   Range: 0 to 5E8 atoms/gram
    %       alpha     - Growth rate parameter (/yr)
    %                   Range: 1E-10 to 1E-8 /yr
    %       beta      - Sigmoid inflection point (years)
    %                   Range: 0.08 to 10 Gyr
    %       eta       - Processing rate parameter (/yr)
    %       Resfrac   - Reservoir fraction (0.9 = 90% convecting mantle)
    %       LVfrac    - Late veneer fraction (% of Earth mass)
    %       PLOTCHECK - 1 to generate plots, 0 to skip plotting
    %       t         - Time vector (years)
    %       T         - Total age of Earth (years)
    %       count     - Iteration counter for file naming
    %
    %   Returns:
    %       NeSucc - 1 if model meets success criteria, 0 otherwise
    %
    %   Success Criteria:
    %       1. 22Ne concentration: (5.8 ± 3.2)E-15 mol/g * conversion
    %       2. 20Ne/22Ne ratio: 12 ± 0.5
    %
    %   Initial Conditions:
    %       - AVCC (Average Carbonaceous Chondrite) composition
    %       - 22Ne concentration from Marty (2012): 1.62E-12 mol/g
    %       - 20Ne/22Ne ratio from Williams and Mukhopadhyay (2018): 13.36 (PSN)
    %
    %   Example:
    %       t = ModelConfig.createTimeVector();
    %       success = ParallelNeModel(1e6, 1e-9, 3e9, 7.5e-10, 0.9, 1, ...
    %                                 0, t, 4.568e9, 1);

    %% Validate inputs
    try
        ModelUtils.validateInputs(...
            'alpha', alpha, ModelConfig.ALPHA_MIN * 0.001, ModelConfig.ALPHA_MAX * 10, ...
            'beta', beta, 0, ModelConfig.BETA_MAX * 2, ...
            'eta', eta, ModelConfig.ETA_MIN * 0.1, ModelConfig.ETA_MAX * 10, ...
            'Resfrac', Resfrac, 0.1, 1.0, ...
            'LVfrac', LVfrac, 0, 10);
    catch ME
        warning('ParallelNeModel:ValidationFailed', ...
                'Input validation failed: %s. Returning failure.', ME.message);
        NeSucc = 0;
        return;
    end

    %% Initialize arrays and constants
    Ned = zeros(1, length(t));          % Downwelling Ne concentration
    Ne20M = zeros(1, length(t));        % 20Ne in mantle
    Ne22M = zeros(1, length(t));        % 22Ne in mantle

    % Physical constants
    Mres = ModelConfig.EARTH_MASS_GRAMS * Resfrac;  % Reservoir mass (g)
    ME = ModelConfig.EARTH_MASS_GRAMS;              % Earth mass (g)

    %% Calculate initial mantle composition
    % Late veneer of carbonaceous chondrites
    % Concentration of 22Ne in CC = 1.62E-12 mol/g (Marty 2012)
    % 20/22 AVCC(PSN) = 13.36 (Williams and Mukhopadhyay 2018)
    mol22M = (ModelConfig.NE22_CONCENTRATION_CC * ME) * (LVfrac / 100);  % mol 22Ne
    Initial22 = mol22M * ModelConfig.AVOGADRO_NUMBER / Mres;  % atoms/gram
    Initial20 = Initial22 * ModelConfig.NE20_22_RATIO_AVCC_PSN;

    %% Calculate sigmoidal downwelling evolution
    Ned = ModelUtils.calculateSigmoidalDownwelling(NeCap, alpha, beta, t);

    %% Time evolution loop - Box model integration
    for i = 2:length(t)
        % Calculate mass of mantle processed in this time step
        dM = ModelUtils.calculateMassProcessed(...
            eta, t(i), t(i-1), T, ModelConfig.PRESENT_DAY_PROCESSING_RATE);

        % Normalized mass processed (fraction of reservoir)
        dM_Mres = dM / Mres;

        % Get values from previous time step
        if i == 2
            Ne20Mlast = Initial20;
            Ne22Mlast = Initial22;
            Nedlast = Ned(i);
        else
            Ne20Mlast = Ne20M(i-1);
            Ne22Mlast = Ne22M(i-1);
            Nedlast = Ned(i-1);
        end

        % Update concentrations using box model equations
        % 20Ne: downwelling with constant composition
        Ne20M(i) = ModelUtils.updateConcentration(Ne20Mlast, dM_Mres, Nedlast, 1);

        % 22Ne: downwelling with fractionation
        Ne22M(i) = ModelUtils.updateConcentration(...
            Ne22Mlast, dM_Mres, Nedlast, 1 / ModelConfig.NE130_FRACTIONATION);
    end

    %% Check success criteria
    % Criterion 1: 22Ne concentration in present-day mantle
    Ne22_min = (ModelConfig.NE22_MANTLE_CENTER - ModelConfig.NE22_MANTLE_UNCERTAINTY) * ...
               ModelConfig.AVOGADRO_NUMBER / Resfrac;
    Ne22_max = (ModelConfig.NE22_MANTLE_CENTER + ModelConfig.NE22_MANTLE_UNCERTAINTY) * ...
               ModelConfig.AVOGADRO_NUMBER / Resfrac;
    Succ1 = ModelUtils.checkSuccessCriteria(Ne22M(end), Ne22_min, Ne22_max);

    % Criterion 2: 20Ne/22Ne ratio in present-day mantle
    ratio = Ne20M(end) / Ne22M(end);
    Ne_ratio_min = ModelConfig.NE20_22_RATIO_MANTLE_CENTER - ModelConfig.NE20_22_RATIO_MANTLE_TOLERANCE;
    Ne_ratio_max = ModelConfig.NE20_22_RATIO_MANTLE_CENTER + ModelConfig.NE20_22_RATIO_MANTLE_TOLERANCE;
    Succ2 = ModelUtils.checkSuccessCriteria(ratio, Ne_ratio_min, Ne_ratio_max);

    % Overall success
    NeSucc = Succ1 && Succ2;

    %% Generate plots if requested
    if PLOTCHECK >= 1
        figHandle = ModelUtils.createStandardFigure(...
            2, 'Neon Isotope Evolution', 'Time (Myr)', '^{20}Ne/^{22}Ne');

        if PLOTCHECK == 1
            lineColor = '-k';
        elseif PLOTCHECK == 3
            lineColor = '-g';
        else % PLOTCHECK == 4
            lineColor = '-b';
        end

        plot(t / 1E6, Ne20M ./ Ne22M, lineColor, 'LineWidth', 1.5);

        % Save figure if PLOTCHECK == 1
        if PLOTCHECK == 1
            filename = sprintf('Ne_evolution_%d.jpg', count);
            ModelUtils.saveFigureSafely(figHandle, filename, 'figures');
        end
    end
end
