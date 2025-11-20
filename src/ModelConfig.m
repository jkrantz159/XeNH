classdef ModelConfig
    %MODELCONFIG Configuration and constants for XeNH mantle recycling models
    %   This class contains all physical constants, initial conditions, and
    %   model parameters used in the H-N-Xe mantle recycling simulations.
    %
    %   References:
    %   - Marty (2012): Noble gas concentrations in chondrites
    %   - Pepin (2000): Isotopic composition of primordial Xe
    %   - Parai and Mukhopadhyay (2018): Model framework
    %   - Barry and Hilton (2016): Nitrogen isotope systematics
    %   - Porcelli, Ballentine, and Wieler (2002): Earth atmosphere composition
    %   - Sephton et al. (2003): Nitrogen content in carbonaceous chondrites
    %   - Williams and Mukhopadhyay (2018): Neon isotope ratios
    %   - Owen et al. (2001): Nitrogen isotope ratios

    properties (Constant)
        %% Physical Constants
        EARTH_MASS_GRAMS = 5.972E27;            % Mass of Earth (g)
        EARTH_AGE_YEARS = 4.568E9;              % Age of Earth (years)
        AVOGADRO_NUMBER = 6.022E23;             % Avogadro's number (atoms/mol)

        %% Mantle Properties
        CONVECTING_MANTLE_FRACTION = 0.9;       % Fraction of mantle that is convecting
        CONVECTING_MANTLE_MASS = 3.6E27;        % Mass of convecting mantle (g) = EARTH_MASS * 0.9

        %% Mantle Processing Rates
        PRESENT_DAY_PROCESSING_RATE = 6.1E17;   % Present day mantle processing rate (g/yr)
                                                % Based on He flux at mid-ocean ridges

        %% Late Veneer
        LATE_VENEER_FRACTION_DEFAULT = 1.0;     % Late veneer as percentage of Earth mass (%)
                                                % Default 1%, consider also 0.5%, 0.1%

        %% Xenon Initial Conditions (AVCC - Average Carbonaceous Chondrite)
        % From Marty (2012) and Pepin (2000)
        XE130_CONCENTRATION_CC = 5.38E-14;      % Concentration of 130Xe in CC (mol/g)
        XE128_130_RATIO_AVCC = 0.5073;          % 128Xe/130Xe ratio in AVCC (Pepin 2000)

        %% Xenon Present-Day Atmosphere (Porcelli et al. 2002)
        XE128_132_RATIO_ATM = 0.0714;           % Atmospheric 128Xe/132Xe ratio
        XE130_132_RATIO_ATM = 0.1514;           % Atmospheric 130Xe/132Xe ratio
        XE128_130_RATIO_ATM_TODAY = 0.4716;     % Modern atmospheric 128Xe/130Xe = 0.0714/0.1514
        XE128_130_RATIO_ATM_INITIAL = 0.5178;   % Initial atmospheric 128Xe/130Xe (at t=0)
                                                % Corresponds to -39 per mille fractionation
        ATMOSPHERE_EVOLUTION_TIME = 2.568E9;    % Time for atmosphere evolution to modern (years)

        %% Xenon Success Criteria (Present-Day Mantle)
        XE130_MANTLE_MIN = 4.3E5;               % Minimum 130Xe concentration in mantle (atoms/g)
        XE130_MANTLE_MAX = 9.2E5;               % Maximum 130Xe concentration in mantle (atoms/g)
        XE128_130_MANTLE_MIN = 0.475;           % Minimum 128Xe/130Xe ratio in mantle
        XE128_130_MANTLE_MAX = 0.478;           % Maximum 128Xe/130Xe ratio in mantle

        %% Nitrogen Initial Conditions
        % From Sephton et al. (2003)
        N_CONCENTRATION_CC = 0.001519;          % Nitrogen concentration in CC (weight fraction, 0.1519 wt%)
        N14_FRACTION = 0.99636;                 % Fraction of nitrogen that is 14N
        N15_14_RATIO_INITIAL = 2.3E-3;          % Initial 15N/14N ratio in PSN (Owen et al. 2001)
        N15_14_RATIO_SEDIMENT = 0.003671565;    % 15N/14N ratio in sediments (from delta-15N = +5 per mille)
        N_ATOMIC_MASS_14 = 14;                  % Atomic mass of 14N

        %% Nitrogen Success Criteria
        N14_MANTLE_MIN_MOL = 7.06E19;           % Minimum 14N in mantle (mol/g mantle)
        N14_MANTLE_MAX_MOL = 9.78E21;           % Maximum 14N in mantle (mol/g mantle)
        N15_14_RATIO_MANTLE_MIN = 0.0036275;    % Minimum 15N/14N ratio in mantle
        N15_14_RATIO_MANTLE_MAX = 0.0036425;    % Maximum 15N/14N ratio in mantle

        %% Neon Initial Conditions
        % From Marty (2012) and Williams and Mukhopadhyay (2018)
        NE22_CONCENTRATION_CC = 1.62E-12;       % Concentration of 22Ne in CC (mol/g)
        NE20_22_RATIO_AVCC_PSN = 13.36;         % 20Ne/22Ne ratio in AVCC (PSN)

        %% Neon Success Criteria
        NE22_MANTLE_CENTER = 5.8E-15;           % Central value for 22Ne in mantle (mol/g)
        NE22_MANTLE_UNCERTAINTY = 3.2E-15;      % Uncertainty in 22Ne concentration (mol/g)
        NE20_22_RATIO_MANTLE_CENTER = 12;       % Central value for 20Ne/22Ne ratio
        NE20_22_RATIO_MANTLE_TOLERANCE = 0.5;   % Tolerance for 20Ne/22Ne ratio
        NE130_FRACTIONATION = 9.80;             % Fractionation factor for Ne isotopes

        %% Model Parameters - Default Ranges for Monte Carlo
        % Growth rate (alpha) range
        ALPHA_MIN = 1E-10;                      % Minimum growth rate (/Gyr)
        ALPHA_MAX = 1E-7;                       % Maximum growth rate (/Gyr)

        % Sigmoid inflection point (beta) range
        BETA_MIN = 0;                           % Minimum inflection point (years)
        BETA_MAX = 10E9;                        % Maximum inflection point (years)

        % Processing rate parameter (eta) range
        ETA_MIN = 7.5E-10;                      % Minimum eta (/yr)
        ETA_MAX = 8.0E-10;                      % Maximum eta (/yr)

        % Carrying capacity ranges
        XE_CAPACITY_MIN = 5;                    % Minimum Xe carrying capacity (atoms/gram)
        XE_CAPACITY_MAX = 5E8;                  % Maximum Xe carrying capacity (atoms/gram)
        N_CAPACITY_MIN = 4;                     % Minimum N carrying capacity (atoms/gram)
        N_CAPACITY_MAX = 4E17;                  % Maximum N carrying capacity (atoms/gram)

        %% Time Vector Parameters
        TIME_START = 0;                         % Start time (years)
        TIME_END = 4.568E9;                     % End time (years) = Earth age
        TIME_STEP_EARLY = 0.1E6;                % Time step for early Earth (years)
        TIME_STEP_MIDDLE = 1E6;                 % Time step for middle period (years)
        TIME_STEP_LATE = 5E6;                   % Time step for late period (years)
        TIME_EARLY_END = 200E6;                 % End of early period (years)
        TIME_MIDDLE_START = 201E6;              % Start of middle period (years)
        TIME_MIDDLE_END = 3300E6;               % End of middle period (years)
        TIME_LATE_START = 3305E6;               % Start of late period (years)
    end

    methods (Static)
        function t = createTimeVector()
            %CREATETIMEVECTOR Create the time vector for model integration
            %   Creates a non-uniform time vector with finer resolution in
            %   early Earth history and coarser resolution later.
            %
            %   Returns:
            %       t - Time vector from 0 to 4.568 Ga (years)

            t = [ModelConfig.TIME_START : ...
                 ModelConfig.TIME_STEP_EARLY : ...
                 ModelConfig.TIME_EARLY_END, ...
                 ModelConfig.TIME_MIDDLE_START : ...
                 ModelConfig.TIME_STEP_MIDDLE : ...
                 ModelConfig.TIME_MIDDLE_END, ...
                 ModelConfig.TIME_LATE_START : ...
                 ModelConfig.TIME_STEP_LATE : ...
                 ModelConfig.TIME_END, ...
                 ModelConfig.TIME_END + 3E6];
        end

        function atm = createAtmosphereEvolution(t)
            %CREATEATMOSPHEREEVOLUTION Create 128Xe/130Xe atmospheric evolution
            %   Models the evolution of atmospheric 128Xe/130Xe ratio from
            %   initial fractionated value to modern value over 2.568 Ga.
            %
            %   Parameters:
            %       t - Time vector (years)
            %
            %   Returns:
            %       atm - Atmospheric 128Xe/130Xe ratio at each time step

            atm = zeros(1, length(t));
            for i = 1:length(t)
                if t(i) < ModelConfig.ATMOSPHERE_EVOLUTION_TIME
                    % Linear evolution from initial to modern
                    atm(i) = ModelConfig.XE128_130_RATIO_ATM_INITIAL - ...
                             (t(i) / ModelConfig.ATMOSPHERE_EVOLUTION_TIME) * ...
                             (ModelConfig.XE128_130_RATIO_ATM_INITIAL - ...
                              ModelConfig.XE128_130_RATIO_ATM_TODAY);
                else
                    % Modern atmospheric value
                    atm(i) = ModelConfig.XE128_130_RATIO_ATM_TODAY;
                end
            end
        end

        function ratio = deltaXeToRatio(deltaXe)
            %DELTAXETORATIO Convert delta notation to 128Xe/130Xe ratio
            %   deltaXe: Delta notation in per mille
            %   Returns: 128Xe/130Xe ratio
            ratio = (deltaXe / 1000 + 1) * ModelConfig.XE128_130_RATIO_ATM_TODAY;
        end

        function deltaXe = ratioToDeltaXe(ratio)
            %RATIOTODELTAXE Convert 128Xe/130Xe ratio to delta notation
            %   ratio: 128Xe/130Xe ratio
            %   Returns: Delta notation in per mille
            deltaXe = ((ratio / ModelConfig.XE128_130_RATIO_ATM_TODAY) - 1) * 1000;
        end

        function ratio = deltaNToRatio(deltaN)
            %DELTANTORATIO Convert delta notation to 15N/14N ratio
            %   deltaN: Delta notation in per mille
            %   Returns: 15N/14N ratio
            airRatio = (1 - ModelConfig.N14_FRACTION) / ModelConfig.N14_FRACTION;
            ratio = (deltaN / 1000 + 1) * airRatio;
        end

        function deltaN = ratioToDeltaN(ratio)
            %RATIODELTAN Convert 15N/14N ratio to delta notation
            %   ratio: 15N/14N ratio
            %   Returns: Delta notation in per mille
            airRatio = (1 - ModelConfig.N14_FRACTION) / ModelConfig.N14_FRACTION;
            deltaN = ((ratio / airRatio) - 1) * 1000;
        end
    end
end
