classdef ModelUtils
    %MODELUTILS Utility functions for XeNH mantle recycling models
    %   This class contains common calculations used across different
    %   geochemical models to reduce code duplication.

    methods (Static)
        function dM = calculateMassProcessed(eta, t_curr, t_prev, T, Qp)
            %CALCULATEMASSPROCESSED Calculate mass of mantle processed in time step
            %   Uses exponential degassing model based on He flux measurements
            %
            %   Parameters:
            %       eta    - Processing rate parameter (/yr)
            %       t_curr - Current time (years)
            %       t_prev - Previous time (years)
            %       T      - Total age of Earth (years)
            %       Qp     - Present day processing rate (g/yr)
            %
            %   Returns:
            %       dM - Mass of mantle processed (g)
            %
            %   Model: Q(t) = Qp * exp(eta * (T - t))
            %   Integrated mass processed between t_prev and t_curr

            dM = (Qp / eta) * (exp(eta * (T - t_prev)) - exp(eta * (T - t_curr)));
        end

        function downwelling = calculateSigmoidalDownwelling(capacity, alpha, beta, t)
            %CALCULATESIGMOIDALDOWNWELLING Calculate downwelling using sigmoidal model
            %   Models the growth of subduction/downwelling over Earth history
            %
            %   Parameters:
            %       capacity - Carrying capacity (atoms/gram or similar)
            %       alpha    - Growth rate parameter (/yr)
            %       beta     - Inflection point (years)
            %       t        - Time vector (years)
            %
            %   Returns:
            %       downwelling - Downwelling flux at each time (same units as capacity)

            downwelling = capacity ./ (1 + exp(-alpha .* (t - beta)));
        end

        function [concentration] = updateConcentration(prevConc, dM_Mres, downwelling, downwellingRatio)
            %UPDATECONCENTRATION Update isotope concentration using box model
            %   Standard box model equation for mixing of downwelling material
            %   into mantle reservoir
            %
            %   Parameters:
            %       prevConc          - Previous concentration in mantle
            %       dM_Mres           - Normalized mass processed (dM/Mres)
            %       downwelling       - Downwelling concentration
            %       downwellingRatio  - Isotopic ratio of downwelling (optional, default 1)
            %
            %   Returns:
            %       concentration - Updated concentration
            %
            %   Model: C(t+dt) = C(t) + (dM/Mres) * (Cd - C(t))
            %   where Cd = downwelling * downwellingRatio

            if nargin < 4
                downwellingRatio = 1;
            end

            concentration = prevConc + dM_Mres * (downwelling * downwellingRatio - prevConc);
        end

        function success = checkSuccessCriteria(value, minVal, maxVal)
            %CHECKSUCCESSCRITERIA Check if value falls within acceptable range
            %
            %   Parameters:
            %       value  - Value to check
            %       minVal - Minimum acceptable value
            %       maxVal - Maximum acceptable value
            %
            %   Returns:
            %       success - 1 if value in range, 0 otherwise

            success = (value >= minVal) && (value <= maxVal);
        end

        function validateInputs(varargin)
            %VALIDATEINPUTS Validate input parameters for models
            %   Checks that parameters are scalar, positive, and within reasonable ranges
            %
            %   Usage: ModelUtils.validateInputs('alpha', alpha, 1e-12, 1e-6, ...)
            %   Parameters come in groups of 4: name, value, min, max

            numArgs = length(varargin);
            if mod(numArgs, 4) ~= 0
                error('ModelUtils:validateInputs', ...
                      'Arguments must come in groups of 4: name, value, min, max');
            end

            for i = 1:4:numArgs
                name = varargin{i};
                value = varargin{i+1};
                minVal = varargin{i+2};
                maxVal = varargin{i+3};

                if ~isscalar(value)
                    error('ModelUtils:validateInputs', ...
                          '%s must be a scalar value', name);
                end

                if ~isnumeric(value) || ~isfinite(value)
                    error('ModelUtils:validateInputs', ...
                          '%s must be a finite numeric value', name);
                end

                if value < minVal || value > maxVal
                    error('ModelUtils:validateInputs', ...
                          '%s = %g is outside valid range [%g, %g]', ...
                          name, value, minVal, maxVal);
                end
            end
        end

        function ensureDirectoryExists(dirPath)
            %ENSUREDIRECTORYEXISTS Create directory if it doesn't exist
            %
            %   Parameters:
            %       dirPath - Path to directory

            if ~exist(dirPath, 'dir')
                mkdir(dirPath);
            end
        end

        function filePath = getOutputPath(filename, subdir)
            %GETOUTPUTPATH Construct output file path
            %   Creates platform-independent path for output files
            %
            %   Parameters:
            %       filename - Name of output file
            %       subdir   - Subdirectory (e.g., 'results', 'figures')
            %
            %   Returns:
            %       filePath - Full path to output file

            if nargin < 2
                subdir = 'results';
            end

            % Get project root (assumes we're in src/ directory)
            projectRoot = fileparts(fileparts(mfilename('fullpath')));
            outputDir = fullfile(projectRoot, subdir);

            % Ensure directory exists
            ModelUtils.ensureDirectoryExists(outputDir);

            filePath = fullfile(outputDir, filename);
        end

        function reportProgress(current, total, message)
            %REPORTPROGRESS Display progress information
            %
            %   Parameters:
            %       current - Current iteration
            %       total   - Total iterations
            %       message - Optional message to display

            if nargin < 3
                message = 'Progress';
            end

            if mod(current, max(1, floor(total/100))) == 0
                percentComplete = (current / total) * 100;
                fprintf('%s: %.1f%% (%d/%d)\n', message, percentComplete, current, total);
            end
        end

        function figHandle = createStandardFigure(figNum, titleStr, xlabelStr, ylabelStr)
            %CREATESTANDARDFIGURE Create figure with standard formatting
            %
            %   Parameters:
            %       figNum     - Figure number
            %       titleStr   - Figure title
            %       xlabelStr  - X-axis label
            %       ylabelStr  - Y-axis label
            %
            %   Returns:
            %       figHandle - Handle to created figure

            figHandle = figure(figNum);
            hold on;
            if nargin >= 2 && ~isempty(titleStr)
                title(titleStr, 'FontSize', 12, 'FontWeight', 'bold');
            end
            if nargin >= 3 && ~isempty(xlabelStr)
                xlabel(xlabelStr, 'FontSize', 11);
            end
            if nargin >= 4 && ~isempty(ylabelStr)
                ylabel(ylabelStr, 'FontSize', 11);
            end
            grid on;
        end

        function saveFigureSafely(figHandle, filename, subdir)
            %SAVEFIGURESAFELY Save figure to file with error handling
            %
            %   Parameters:
            %       figHandle - Handle to figure
            %       filename  - Output filename
            %       subdir    - Subdirectory (default 'figures')

            if nargin < 3
                subdir = 'figures';
            end

            try
                filepath = ModelUtils.getOutputPath(filename, subdir);
                saveas(figHandle, filepath);
                fprintf('Figure saved: %s\n', filepath);
            catch ME
                warning('ModelUtils:saveFigureSafely', ...
                        'Failed to save figure: %s', ME.message);
            end
        end
    end
end
