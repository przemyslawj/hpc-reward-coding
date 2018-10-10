function tests = testCalculateVelocity
tests = functiontests(localfunctions);
end

function testCalculatedVelocity(testCase)
timestamp = (0:100:700)';
smooth_trans_x = [0 0 2 4 4 8 12 12]';
smooth_trans_y = zeros(size(smooth_trans_x));
exp_dist = [0 0 2 2 0 4 4 0]';
T = table(timestamp, smooth_trans_x, smooth_trans_y, 'VariableNames', ...
    {'timestamp', 'smooth_trans_x', 'smooth_trans_y'});
result = calculateVelocity(T);
testCase.verifyEqual(result.dist, exp_dist);
testCase.verifyEqual(result.velocity, exp_dist / 0.1);
end
