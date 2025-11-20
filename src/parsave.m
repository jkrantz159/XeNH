function parsave(alpha, beta, eta, Xed, Nd, Resfrac, LVfrac, varargin)
    %PARSAVE Save successful model parameters to file (parallel-safe)
    %   This function is designed to be called from parallel workers to save
    %   successful model run parameters without file conflicts.
    %
    %   Parameters:
    %       alpha   - Growth rate parameter (/yr)
    %       beta    - Sigmoid inflection point (years)
    %       eta     - Processing rate parameter (/yr)
    %       Xed     - Xe carrying capacity/downwelling (atoms/gram)
    %       Nd      - N carrying capacity/downwelling (atoms/gram)
    %       Resfrac - Reservoir fraction (e.g., 0.9 for 90% convecting mantle)
    %       LVfrac  - Late veneer fraction (% of Earth mass)
    %       varargin - Optional: 'filename', custom_filename
    %
    %   Output:
    %       Appends parameters as CSV to results/success.txt
    %       Format: alpha,beta,eta,Xed,Nd,Resfrac,LVfrac
    %
    %   Example:
    %       parsave(1e-9, 3e9, 7.5e-10, 1e6, 1e15, 0.9, 1.0)
    %       parsave(1e-9, 3e9, 7.5e-10, 1e6, 1e15, 0.9, 1.0, 'filename', 'custom_success.txt')

    % Parse optional arguments
    p = inputParser;
    addParameter(p, 'filename', 'success.txt', @ischar);
    parse(p, varargin{:});
    filename = p.Results.filename;

    % Construct full output path
    try
        filepath = ModelUtils.getOutputPath(filename, 'results');
    catch
        % Fallback if ModelUtils not available
        filepath = filename;
    end

    % Validate inputs
    if ~all(isfinite([alpha, beta, eta, Xed, Nd, Resfrac, LVfrac]))
        warning('parsave:InvalidInput', ...
                'One or more parameters are not finite. Skipping save.');
        return;
    end

    % Open file and write data with error handling
    try
        fid = fopen(filepath, 'a+');
        if fid == -1
            error('parsave:FileOpenError', 'Cannot open file: %s', filepath);
        end

        % Write data in CSV format
        fprintf(fid, '%E,%E,%E,%E,%E,%E,%E\n', ...
                alpha, beta, eta, Xed, Nd, Resfrac, LVfrac);

        % Close file
        fclose(fid);

    catch ME
        % Ensure file is closed even if error occurs
        if exist('fid', 'var') && fid ~= -1
            fclose(fid);
        end
        warning('parsave:WriteError', ...
                'Failed to save parameters: %s', ME.message);
    end
end
