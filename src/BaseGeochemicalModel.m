classdef BaseGeochemicalModel
    %BASEGEOCHEMICALMODEL Base class for geochemical box models
    %   Provides common functionality for modeling mantle degassing and
    %   isotope evolution through Earth history using box model approach.
    %
    %   This class implements the framework from Parai and Mukhopadhyay (2018)
    %   for modeling volatile recycling in Earth's mantle.

    methods (Static)
        function [success, isotope1, isotope2] = runBoxModel(params)
            %RUNBOXMODEL Execute a two-isotope box model simulation
            %   Generic implementation of mantle degassing and recycling model
            %
            %   Parameters (in params struct):
            %       capacity    - Carrying capacity for downwelling
            %       alpha       - Growth rate parameter (/yr)
            %       beta        - Sigmoid inflection point (years)
            %       eta         - Processing rate parameter (/yr)
            %       Resfrac     - Reservoir fraction (0.9 = 90% convecting mantle)
            %       LVfrac      - Late veneer fraction (% of Earth mass)
            %       t           - Time vector (years)
            %       T           - Total age (years)
            %       initial1    - Initial concentration of isotope 1
            %       initial2    - Initial concentration of isotope 2
            %       downwellingRatio - Ratio for isotope 2 in downwelling
            %       criteria1_min, criteria1_max - Success criteria for isotope 1
            %       criteria2_min, criteria2_max - Success criteria for ratio
            %       downwellingModifier - Optional modifier for downwelling (e.g., atm for Xe)
            %
            %   Returns:
            %       success  - 1 if model meets success criteria, 0 otherwise
            %       isotope1 - Concentration history of isotope 1
            %       isotope2 - Concentration history of isotope 2

            % Validate inputs
            ModelUtils.validateInputs(...
                'alpha', params.alpha, ModelConfig.ALPHA_MIN * 1e-3, ModelConfig.ALPHA_MAX * 10, ...
                'beta', params.beta, 0, ModelConfig.BETA_MAX * 2, ...
                'eta', params.eta, ModelConfig.ETA_MIN * 0.1, ModelConfig.ETA_MAX * 10, ...
                'Resfrac', params.Resfrac, 0.1, 1.0);

            % Calculate reservoir mass
            Mres = ModelConfig.EARTH_MASS_GRAMS * params.Resfrac;

            % Initialize arrays
            isotope1 = zeros(1, length(params.t));
            isotope2 = zeros(1, length(params.t));

            % Calculate downwelling using sigmoidal model
            downwelling = ModelUtils.calculateSigmoidalDownwelling(...
                params.capacity, params.alpha, params.beta, params.t);

            % Time evolution loop
            for i = 2:length(params.t)
                % Calculate mass processed in this time step
                dM = ModelUtils.calculateMassProcessed(...
                    params.eta, params.t(i), params.t(i-1), params.T, ...
                    ModelConfig.PRESENT_DAY_PROCESSING_RATE);

                % Normalized mass processed
                dM_Mres = dM / Mres;

                % Get previous concentrations
                if i == 2
                    isotope1_prev = params.initial1;
                    isotope2_prev = params.initial2;
                    downwelling_prev = downwelling(i);
                    if isfield(params, 'downwellingModifier')
                        modifier_prev = params.downwellingModifier(1);
                    else
                        modifier_prev = 1;
                    end
                else
                    isotope1_prev = isotope1(i-1);
                    isotope2_prev = isotope2(i-1);
                    downwelling_prev = downwelling(i-1);
                    if isfield(params, 'downwellingModifier')
                        modifier_prev = params.downwellingModifier(i-1);
                    else
                        modifier_prev = 1;
                    end
                end

                % Update concentrations using box model equations
                isotope1(i) = ModelUtils.updateConcentration(...
                    isotope1_prev, dM_Mres, downwelling_prev, 1);

                isotope2(i) = ModelUtils.updateConcentration(...
                    isotope2_prev, dM_Mres, downwelling_prev, ...
                    params.downwellingRatio * modifier_prev);
            end

            % Check success criteria
            % Criterion 1: Absolute concentration of isotope 1
            succ1 = ModelUtils.checkSuccessCriteria(...
                isotope1(end), params.criteria1_min, params.criteria1_max);

            % Criterion 2: Isotope ratio
            ratio = isotope2(end) / isotope1(end);
            succ2 = ModelUtils.checkSuccessCriteria(...
                ratio, params.criteria2_min, params.criteria2_max);

            % Overall success
            success = succ1 && succ2;
        end
    end
end
