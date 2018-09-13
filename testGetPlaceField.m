function tests = testGetPlaceField
tests = functiontests(localfunctions);
end

function testCorrectPlaceField(testCase)
binSize = 10;
pos_x = [-1 0 0 1 2 3 3 2] * binSize;
pos_y = [0 0 0 0 0 0 0 0];
trace = [0 1 1 2 1 20 22 3];

expectedPlaceField = zeros(200 / binSize, 200 / binSize);
expectedPlaceField(10, 9) = 0;
expectedPlaceField(10, 10) = 1;
expectedPlaceField(10, 11) = 2;
expectedPlaceField(10, 12) = 2;
expectedPlaceField(10, 13) = 21;

[placeField, PCI] = getPlaceField(pos_x, pos_y, trace, 10);

testCase.verifyEqual(placeField, expectedPlaceField);
testCase.verifyGreaterThan(PCI, 1.05);
testCase.verifyLessThan(PCI, 1.25);

end

function testLowPCI(testCase)
binSize = 10;
pos_x = (1:10) * binSize;
pos_y = zeros(size(pos_x));
trace = ones(size(pos_x)) * 10;
trace(1) = 0;

[~, PCI] = getPlaceField(pos_x, pos_y, trace, 10);
testCase.verifyLessThan(PCI, 0.2);
end
