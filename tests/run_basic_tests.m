%% RUN_BASIC_TESTS - Basic functionality tests for XeNH code
%   This script tests the core functionality of the XeNH models to verify
%   proper installation and basic operation.
%
%   Tests:
%       1. ModelConfig constants are accessible
%       2. ModelUtils functions work correctly
%       3. Time vector generation
%       4. Atmospheric evolution calculation
%       5. Individual model functions execute without errors
%
%   Usage:
%       cd tests
%       run_basic_tests

fprintf('====================================================\n');
fprintf('  XeNH Basic Functionality Tests\n');
fprintf('====================================================\n\n');

% Add src to path if not already there
if ~contains(path, 'src')
    addpath(fullfile(fileparts(pwd), 'src'));
end

testsPassed = 0;
testsFailed = 0;

%% Test 1: ModelConfig Constants
fprintf('Test 1: ModelConfig constants... ');
try
    assert(ModelConfig.EARTH_MASS_GRAMS == 5.972E27, 'Earth mass incorrect');
    assert(ModelConfig.EARTH_AGE_YEARS == 4.568E9, 'Earth age incorrect');
    fprintf('PASSED\n');
    testsPassed = testsPassed + 1;
catch ME
    fprintf('FAILED: %s\n', ME.message);
    testsFailed = testsFailed + 1;
end

%% Test 2: Time Vector Generation
fprintf('Test 2: Time vector generation... ');
try
    t = ModelConfig.createTimeVector();
    assert(length(t) > 0, 'Time vector is empty');
    assert(t(1) == 0, 'Time vector should start at 0');
    assert(t(end) >= 4.568E9, 'Time vector should extend to present');
    fprintf('PASSED\n');
    testsPassed = testsPassed + 1;
catch ME
    fprintf('FAILED: %s\n', ME.message);
    testsFailed = testsFailed + 1;
end

%% Test 3: Atmospheric Evolution
fprintf('Test 3: Atmospheric evolution... ');
try
    t = ModelConfig.createTimeVector();
    atm = ModelConfig.createAtmosphereEvolution(t);
    assert(length(atm) == length(t), 'Atmosphere vector length mismatch');
    assert(atm(end) == ModelConfig.XE128_130_RATIO_ATM_TODAY, 'Modern atmosphere incorrect');
    fprintf('PASSED\n');
    testsPassed = testsPassed + 1;
catch ME
    fprintf('FAILED: %s\n', ME.message);
    testsFailed = testsFailed + 1;
end

%% Test 4: ModelUtils Functions
fprintf('Test 4: ModelUtils.calculateMassProcessed... ');
try
    eta = 7.5E-10;
    t_curr = 1E9;
    t_prev = 0;
    T = 4.568E9;
    Qp = 6.1E17;
    dM = ModelUtils.calculateMassProcessed(eta, t_curr, t_prev, T, Qp);
    assert(dM > 0, 'Mass processed should be positive');
    assert(isfinite(dM), 'Mass processed should be finite');
    fprintf('PASSED\n');
    testsPassed = testsPassed + 1;
catch ME
    fprintf('FAILED: %s\n', ME.message);
    testsFailed = testsFailed + 1;
end

%% Test 5: Sigmoidal Downwelling
fprintf('Test 5: ModelUtils.calculateSigmoidalDownwelling... ');
try
    t = linspace(0, 4.568E9, 100);
    capacity = 1E6;
    alpha = 1E-9;
    beta = 2E9;
    downwelling = ModelUtils.calculateSigmoidalDownwelling(capacity, alpha, beta, t);
    assert(all(downwelling >= 0), 'Downwelling should be non-negative');
    assert(all(downwelling <= capacity), 'Downwelling should not exceed capacity');
    fprintf('PASSED\n');
    testsPassed = testsPassed + 1;
catch ME
    fprintf('FAILED: %s\n', ME.message);
    testsFailed = testsFailed + 1;
end

%% Test 6: ParallelXeModel Execution
fprintf('Test 6: ParallelXeModel execution... ');
try
    t = ModelConfig.createTimeVector();
    atm = ModelConfig.createAtmosphereEvolution(t);
    deltaXe = @ModelConfig.ratioToDeltaXe;

    % Use parameters that should fail (to test execution, not success)
    XeSucc = ParallelXeModel(1e3, 1e-9, 3e9, 7.5e-10, 0.9, 1, ...
                             0, atm, t, 4.568e9, deltaXe, 1);
    assert(XeSucc == 0 || XeSucc == 1, 'Return value should be 0 or 1');
    fprintf('PASSED\n');
    testsPassed = testsPassed + 1;
catch ME
    fprintf('FAILED: %s\n', ME.message);
    testsFailed = testsFailed + 1;
end

%% Test 7: ParallelNewNModel Execution
fprintf('Test 7: ParallelNewNModel execution... ');
try
    t = ModelConfig.createTimeVector();
    deltaN = @ModelConfig.ratioToDeltaN;

    % Use parameters that should fail (to test execution, not success)
    NSucc = ParallelNewNModel(1e10, 1e-9, 3e9, 7.5e-10, 0.9, 1, ...
                              0, t, 4.568e9, deltaN, 1);
    assert(NSucc == 0 || NSucc == 1, 'Return value should be 0 or 1');
    fprintf('PASSED\n');
    testsPassed = testsPassed + 1;
catch ME
    fprintf('FAILED: %s\n', ME.message);
    testsFailed = testsFailed + 1;
end

%% Test 8: Delta Notation Conversions
fprintf('Test 8: Delta notation conversions... ');
try
    % Test Xe conversion
    ratio = 0.4716;
    delta = ModelConfig.ratioToDeltaXe(ratio);
    assert(abs(delta) < 1, 'Delta should be close to 0 for modern ratio');

    % Test N conversion
    ratio = 0.003647;
    delta = ModelConfig.ratioToDeltaN(ratio);
    assert(isfinite(delta), 'Delta should be finite');

    fprintf('PASSED\n');
    testsPassed = testsPassed + 1;
catch ME
    fprintf('FAILED: %s\n', ME.message);
    testsFailed = testsFailed + 1;
end

%% Summary
fprintf('\n====================================================\n');
fprintf('  Test Summary\n');
fprintf('====================================================\n');
fprintf('Tests passed: %d\n', testsPassed);
fprintf('Tests failed: %d\n', testsFailed);
fprintf('Success rate: %.1f%%\n', (testsPassed / (testsPassed + testsFailed)) * 100);

if testsFailed == 0
    fprintf('\nAll tests PASSED! Installation verified.\n');
else
    fprintf('\nSome tests FAILED. Please check installation and dependencies.\n');
end
fprintf('====================================================\n');
