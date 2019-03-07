function tests = testIsAtReward()
tests = functiontests(localfunctions);
end

function testAtRewardInTheMiddleOfTrial(testCase)
velocities = [9 9 0 1 1 0 1 0 0 9 6 5];
distances = [10 9 6 7 5 4 5 4 4 7 9 6];
expResult =  [0 0 0 0 1 1 1 1 1 0 0 0];

result = isAtReward(velocities, distances, 6, 2);
testCase.verifyEqual(result, expResult);
end

function testAtRewardInTheEndOfTrial(testCase)
velocities = [9 9 0 1 1 0 1 0 0];
distances = [10 9 6 7 5 4 5 4 4];
expResult =  [0 0 0 0 1 1 1 1 1];

result = isAtReward(velocities, distances, 6, 2);
testCase.verifyEqual(result, expResult);
end
