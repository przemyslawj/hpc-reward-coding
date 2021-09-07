function [] = plotCellTraces(dataTable, cellIndex)
%PLOTCELLTRACE Plots traces for each trial and the mean trace.
% Timeline on x axis is centered around two times mice arrive at reward.

nrewards = max(dataTable.arrivedAtReward);

[grouping, groupingName] = grp2idx(dataTable.trial_id);

hold on;
subplot(nrewards,1,1);

minTimestamps = zeros(nrewards, 1);
maxTimestamps = zeros(nrewards, 1);

for i=1:nrewards
    minTimestamps(i) = min(dataTable{:, ['timestampFromReward_' num2str(i-1)]});
    maxTimestamps(i) = max(dataTable{:, ['timestampFromReward_' num2str(i-1)]});
end
timestampStep = dataTable.timestamp(2) - dataTable.timestamp(1);

maxYVal = max(abs(dataTable.trace(:,cellIndex)));

%% Align trial traces and calculate mean trace at each reward
alignedTrialTraces = nan(nrewards, numel(groupingName), (max(maxTimestamps)-min(minTimestamps))/timestampStep);
alignedTrialEvents = nan(size(alignedTrialTraces));
for trialIndex = 1:numel(groupingName)
    trialTable = dataTable(grouping == trialIndex,:);
    arrivedAtRewardIndex = find(trialTable.arrivedAtReward > 0);
    for i = 1:numel(arrivedAtRewardIndex)
        rewardIndex = trialTable.arrivedAtReward(arrivedAtRewardIndex(i));
        timestamps = dataTable{grouping == trialIndex, ['timestampFromReward_' num2str(rewardIndex-1)]};

        timestampDiff = abs(minTimestamps(rewardIndex) - timestamps(1));
        offset = timestampDiff / timestampStep;
        alignedTrialTraces(rewardIndex, trialIndex, (1 + offset):(offset + numel(timestamps))) = ...
            trialTable.trace(:, cellIndex)';
        alignedTrialEvents(rewardIndex, trialIndex, (1 + offset):(offset + numel(timestamps))) = ...
            trialTable.events(:, cellIndex)';
    end
end
globalTimestamps = zeros(nrewards, (max(maxTimestamps)-min(minTimestamps))/timestampStep);
for i = 1:nrewards
    globalTimestamps(i,:) = linspace(minTimestamps(i),...
            double(timestampStep)*size(globalTimestamps,2), ...
            size(globalTimestamps, 2));
end

meanRewardTrace = mean(alignedTrialTraces, 2, 'omitnan');

%% Plot traces and vertical bars marking reward times
for i = 1:nrewards
    subplot(2,1,i);
    hold on;
    for trialIndex = 1:numel(groupingName)
       alignedTrace = alignedTrialTraces(i, trialIndex,:);
       plot(globalTimestamps(i,:), alignedTrace(:));
       eventTimes = find(alignedTrialEvents(i, trialIndex,:) == 1);
       plot(globalTimestamps(i, eventTimes), ...
            reshape(alignedTrace(eventTimes), [], 1), ...
            'r*');
    end
     
    plot(globalTimestamps(i,:), meanRewardTrace(i,:), 'black', 'LineWidth', 1.5);
    plot([0 0], [-10 0.8*int32(maxYVal)], 'r');
    xlim([minTimestamps(i) maxTimestamps(i)]);
    xlabel('Time from reward (ms)');
    ylabel('dF/F');
end

hold off;

end
