function [ result ] = isAtReward(velocity, rewardDist, maxDistReward, minTimestampsAtReward)

RUNNING_VELOCITY_THRESH = 4;
if ~exist('maxDistReward','var') || isempty(maxDistReward)
    maxDistReward = 6;
end    
if ~exist('minTimestampsAtReward', 'var') || isempty(minTimestampsAtReward)
    minTimestampsAtReward = 10;
end

atReward = (velocity <= RUNNING_VELOCITY_THRESH) & ...
        (rewardDist <= maxDistReward);

% satisfy minimum count of the consecutive timestamps at reward
startsAtReward = find(diff(atReward) == 1) + 1;
endsAtReward = find(diff(atReward) == -1);
if numel(endsAtReward) < numel(startsAtReward)
    endsAtReward = [endsAtReward, numel(atReward)];
end

longerStartsAtRewardIndex = find(endsAtReward - startsAtReward >= minTimestampsAtReward);
result = zeros(size(atReward));
for i = 1:numel(longerStartsAtRewardIndex)
    index = longerStartsAtRewardIndex(i);
    result(startsAtReward(index):endsAtReward(index)) = 1;
end

end

