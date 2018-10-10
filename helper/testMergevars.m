function tests = testMergevars
tests = functiontests(localfunctions);
end

function testMergevarsSuccess(testCase)
T = table((1:10)', (11:20)', (0:9)', (21:30)', (31:40)', ... 
          'VariableNames', {'A', 'trace_1', 'trace_2', 'B', 'trace_3'});

resultTable = mergevars(T, 'trace');
testCase.verifyEqual(size(resultTable), [10, 3]);
testCase.verifyEqual(resultTable{:,'A'}, T{:,'A'});
testCase.verifyEqual(resultTable{:,'B'}, T{:,'B'});
resultTraces = resultTable{:, 'trace'};
testCase.verifyEqual(size(resultTraces), [10, 3]);

testCase.verifyEqual(resultTraces(:, 1), T{:, 'trace_1'});
testCase.verifyEqual(resultTraces(:, 2), T{:, 'trace_2'});
testCase.verifyEqual(resultTraces(:, 3), T{:, 'trace_3'});
end
