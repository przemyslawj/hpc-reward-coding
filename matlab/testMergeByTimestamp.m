function tests = testMergeByTimestamp
tests = functiontests(localfunctions);
end

function testMergedByTimestamp(testCase)
trace = [1:5; 6:10; 11:15];
traceTimestamp = 0:5:20;
T = table((0:2:20)', (10:20)', (20:30)', (10:-1:0)', (10:-1:0)', ones(11,1), ... 
          'VariableNames', {'timestamp', 'smooth_trans_x', 'smooth_trans_y', ...
          'dist_reward0', 'dist_reward1', 'inside_roi'});

resultTable = mergeByTimestamp(T, trace, traceTimestamp);
testCase.verifyEqual(resultTable{:,'timestamp'}, traceTimestamp');
testCase.verifyEqual(resultTable{:,'trace'}, trace');
testCase.verifyEqual(resultTable{1, 'smooth_trans_x'}, T{1, 'smooth_trans_x'});
testCase.verifyEqual(resultTable{1, 'smooth_trans_y'}, T{1, 'smooth_trans_y'});
testCase.verifyEqual(resultTable{2, 'smooth_trans_x'}, 12.5);
testCase.verifyEqual(resultTable{2, 'smooth_trans_y'}, 22.5);
testCase.verifyEqual(resultTable{5, 'smooth_trans_x'}, 20);
testCase.verifyEqual(resultTable{4, 'smooth_trans_x'}, 17.5);
end
