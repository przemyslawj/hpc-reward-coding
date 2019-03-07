function [dataTable] = addAlignedTrialTimestamps(dataTable)
% Adds to the table columns with timestamps aligned to have time 0 at time
% when arriving to reward. Added columns: timestampFromReward_0
% timestampFromReward_1

[grouping, groupingName] = grp2idx(dataTable.trial_id);

timestamps = zeros(size(dataTable,1), 2);

for trialIndex = 1:numel(groupingName)
    trialTable = dataTable(grouping == trialIndex,:);
    arrivedAtRewardIndex = find(trialTable.arrivedAtReward > 0);
    for i = 1:numel(arrivedAtRewardIndex)
        rewardIndex = trialTable.arrivedAtReward(arrivedAtRewardIndex(i));
        rewardTimestamp = int32(trialTable.timestamp(arrivedAtRewardIndex(i)));
        shiftedTimestamps = (int32(trialTable.timestamp) - rewardTimestamp);
        timestamps(grouping == trialIndex, rewardIndex) = shiftedTimestamps;
    end
end

dataTable.timestampFromReward_0 = timestamps(:,1);
dataTable.timestampFromReward_1 = timestamps(:,2);

end
