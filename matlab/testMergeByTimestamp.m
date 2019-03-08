function tests = testMergeByTimestamp
tests = functiontests(localfunctions);
end

function testMergedByTimestamp(testCase)
trace = [0:5; 5:10; 10:15];
traceTimestamp = [-5 0:5:20];
T = table((0:2:20)', (10:20)', (20:30)', (10:-1:0)', (10:-1:0)', ones(11,1), ... 
          'VariableNames', {'timestamp', 'smooth_trans_x', 'smooth_trans_y', ...
          'dist_reward0', 'dist_reward1', 'inside_roi'});

resultTable = mergeByTimestamp(T, trace, traceTimestamp);
testCase.verifyEqual(resultTable{:,'timestamp'}, traceTimestamp(2:end)');
testCase.verifyEqual(resultTable{:,'trace'}, trace(:,2:end)');
testCase.verifyEqual(resultTable{1, 'smooth_trans_x'}, T{1, 'smooth_trans_x'});
testCase.verifyEqual(resultTable{1, 'smooth_trans_y'}, T{1, 'smooth_trans_y'});
testCase.verifyEqual(resultTable{2, 'smooth_trans_x'}, 12);
testCase.verifyEqual(resultTable{2, 'smooth_trans_y'}, 22);
testCase.verifyEqual(resultTable{5, 'smooth_trans_x'}, 20);
testCase.verifyEqual(resultTable{4, 'smooth_trans_x'}, 17);
end
