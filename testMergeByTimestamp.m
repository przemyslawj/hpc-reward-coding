function tests = testMergeByTimestamp
tests = functiontests(localfunctions);
end

function testMergedByTimestamp(testCase)
trace = [1:5; 6:10; 11:15];
traceTimestamp = 0:5:20;
T = table((0:2:20)', (10:20)', (20:30)', (10:-1:0)', (10:-1:0)', ... 
          'VariableNames', {'timestamp', 'trans_x', 'trans_y', 'dist_reward0', 'dist_reward1'});

resultTable = mergeByTimestamp(T, trace, traceTimestamp);
testCase.verifyEqual(resultTable{:,'timestamp'}, traceTimestamp');
testCase.verifyEqual(resultTable{:,'trace'}, trace');
testCase.verifyEqual(resultTable{1, 'trans_x'}, T{1, 'trans_x'});
testCase.verifyEqual(resultTable{1, 'trans_y'}, T{1, 'trans_y'});
testCase.verifyEqual(resultTable{2, 'trans_x'}, 12.5);
testCase.verifyEqual(resultTable{2, 'trans_y'}, 22.5);
testCase.verifyEqual(resultTable{5, 'trans_x'}, 20);
testCase.verifyEqual(resultTable{4, 'trans_x'}, 17.5);
end