function [XeSucc] = ParallelXeModel(XeCap, alpha, beta, eta, Resfrac, LVfrac, PLOTCHECK, atm, t, T, deltaXe, count)
    %PARALLELXEMODEL Xenon isotope evolution model for Earth's mantle
    %   Models the evolution of 128Xe and 130Xe in Earth's mantle over
    %   geological time using a box model with sigmoidal downwelling growth.
    %
    %   Based on: Parai and Mukhopadhyay (2018)
    %
    %   Parameters:
    %       XeCap     - Xe carrying capacity/present day downwelling (atoms/gram)
    %                   Range: 5 to 5E8 atoms/gram
    %       alpha     - Growth rate parameter (/yr)
    %                   Range: 1E-10 to 1E-7 /yr
    %       beta      - Sigmoid inflection point (years)
    %                   Range: 0 to 10 Gyr
    %       eta       - Processing rate parameter (/yr)
    %                   Range: 7E-10 to 8E-10 /yr
    %       Resfrac   - Reservoir fraction (0.9 = 90% convecting mantle)
    %       LVfrac    - Late veneer fraction (% of Earth mass)
    %       PLOTCHECK - 1 to generate plots, 0 to skip plotting
    %       atm       - Atmospheric 128Xe/130Xe evolution vector
    %       t         - Time vector (years)
    %       T         - Total age of Earth (years)
    %       deltaXe   - Function handle to convert ratio to delta notation
    %       count     - Iteration counter for file naming
    %
    %   Returns:
    %       XeSucc - 1 if model meets success criteria, 0 otherwise
    %
    %   Success Criteria:
    %       1. 130Xe concentration: 4.3E5 to 9.2E5 atoms/gram
    %       2. 128Xe/130Xe ratio: 0.475 to 0.478
    %
    %   Initial Conditions:
    %       - AVCC (Average Carbonaceous Chondrite) composition
    %       - 130Xe concentration from Marty (2012): 5.38E-14 mol/g
    %       - 128Xe/130Xe ratio from Pepin (2000): 0.5073
    %
    %   Example:
    %       t = ModelConfig.createTimeVector();
    %       atm = ModelConfig.createAtmosphereEvolution(t);
    %       deltaXe = @ModelConfig.ratioToDeltaXe;
    %       success = ParallelXeModel(1e6, 1e-9, 3e9, 7.5e-10, 0.9, 1, ...
    %                                 0, atm, t, 4.568e9, deltaXe, 1);

    %% Validate inputs
    try
        ModelUtils.validateInputs(...
            'alpha', alpha, ModelConfig.ALPHA_MIN * 0.001, ModelConfig.ALPHA_MAX * 10, ...
            'beta', beta, 0, ModelConfig.BETA_MAX * 2, ...
            'eta', eta, ModelConfig.ETA_MIN * 0.1, ModelConfig.ETA_MAX * 10, ...
            'Resfrac', Resfrac, 0.1, 1.0, ...
            'LVfrac', LVfrac, 0, 10);
    catch ME
        warning('ParallelXeModel:ValidationFailed', ...
                'Input validation failed: %s. Returning failure.', ME.message);
        XeSucc = 0;
        return;
    end

    %% Initialize arrays and constants
    Xed = zeros(1, length(t));          % Downwelling Xe concentration
    Xe130M = zeros(1, length(t));       % 130Xe in mantle
    Xe128M = zeros(1, length(t));       % 128Xe in mantle

    % Physical constants
    Mres = ModelConfig.EARTH_MASS_GRAMS * Resfrac;  % Reservoir mass (g)
    ME = ModelConfig.EARTH_MASS_GRAMS;              % Earth mass (g)

    %% Calculate initial mantle composition
    % Late veneer of carbonaceous chondrites
    % Concentration of 130Xe in CC = 5.38E-14 mol/g (Marty 2012)
    % 128/130 AVCC = 0.5073 (Pepin 2000)
    mol130M = (ModelConfig.XE130_CONCENTRATION_CC * ME) * (LVfrac / 100);  % mol 130Xe
    Initial130 = mol130M * ModelConfig.AVOGADRO_NUMBER / Mres;  % atoms/gram
    Initial128 = Initial130 * ModelConfig.XE128_130_RATIO_AVCC;

    %% Calculate sigmoidal downwelling evolution
    Xed = ModelUtils.calculateSigmoidalDownwelling(XeCap, alpha, beta, t);

    %% Time evolution loop - Box model integration
    for i = 2:length(t)
        % Calculate mass of mantle processed in this time step
        dM = ModelUtils.calculateMassProcessed(...
            eta, t(i), t(i-1), T, ModelConfig.PRESENT_DAY_PROCESSING_RATE);

        % Normalized mass processed (fraction of reservoir)
        dM_Mres = dM / Mres;

        % Get values from previous time step
        if i == 2
            Xe130Mlast = Initial130;
            Xe128Mlast = Initial128;
            atmlast = atm(1);
            Xedlast = Xed(i);
        else
            Xe130Mlast = Xe130M(i-1);
            Xe128Mlast = Xe128M(i-1);
            atmlast = atm(i-1);
            Xedlast = Xed(i-1);
        end

        % Update concentrations using box model equations
        % 130Xe: downwelling with constant isotopic composition
        Xe130M(i) = ModelUtils.updateConcentration(Xe130Mlast, dM_Mres, Xedlast, 1);

        % 128Xe: downwelling with time-varying atmospheric composition
        Xe128M(i) = ModelUtils.updateConcentration(Xe128Mlast, dM_Mres, Xedlast, atmlast);
    end

    %% Check success criteria
    % Criterion 1: 130Xe concentration in present-day mantle
    Succ1 = ModelUtils.checkSuccessCriteria(...
        Xe130M(end), ModelConfig.XE130_MANTLE_MIN, ModelConfig.XE130_MANTLE_MAX);

    % Criterion 2: 128Xe/130Xe ratio in present-day mantle
    ratio = Xe128M(end) / Xe130M(end);
    Succ2 = ModelUtils.checkSuccessCriteria(...
        ratio, ModelConfig.XE128_130_MANTLE_MIN, ModelConfig.XE128_130_MANTLE_MAX);

    % Overall success
    XeSucc = Succ1 && Succ2;

    %% Generate plots if requested
    if PLOTCHECK == 1
        figHandle = ModelUtils.createStandardFigure(...
            3, 'Xenon Isotope Evolution', 'Time (Myr)', '\delta ^{128}Xe (‰)');
        plot(t / 1E6, deltaXe(Xe128M ./ Xe130M), '-k', 'LineWidth', 1.5);

        % Save figure to figures directory
        filename = sprintf('Xe_evolution_%d.jpg', count);
        ModelUtils.saveFigureSafely(figHandle, filename, 'figures');
    end
end
