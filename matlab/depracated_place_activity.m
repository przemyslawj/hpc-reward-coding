% Analyses place activity. Requires running loadCnmfeTrial.m first.

%% Load cheeseboard map and reward locations
locationsFile = [ datedRootDir filesep 'locations.csv' ];
locationsTable = readtable(locationsFile);

cheeseboardMapFile = [ rootDir filesep 'cheeseboard_map.csv' ];
cheeseboardMapTable = readtable(cheeseboardMapFile);
rewardLocationsTable = locationsTable(...
        strcmp(locationsTable.Animal, animal) == 1 ...
        & strcmp(locationsTable.Valence, 'Positive') == 1, :);
negativeLocationsTable = locationsTable(...
        strcmp(locationsTable.Animal, animal) == 1 ...
        & strcmp(locationsTable.Valence, 'Negative') == 1, :);

rewardedWellsMap = [];
if ~isempty(rewardLocationsTable)
    rewardedWellsMap = cheeseboardMapTable(...
        (cheeseboardMapTable.Row_X == rewardLocationsTable.Well_row(1) ...
         & cheeseboardMapTable.Row_Y == rewardLocationsTable.Well_col(1))...
        | ...
         (cheeseboardMapTable.Row_X == rewardLocationsTable.Well_row(2) ...
         & cheeseboardMapTable.Row_Y == rewardLocationsTable.Well_col(2)), ...
        :);
end

negativeWellsMap = [];
if ~isempty(negativeLocationsTable)
    negativeWellsMap = cheeseboardMapTable(...
            cheeseboardMapTable.Row_X == negativeLocationsTable.Well_row(1) ...
            & cheeseboardMapTable.Row_Y == negativeLocationsTable.Well_col(1));
end
%% Calculate place fields
% TODO: evaluate different velocity thresholds
RUNNING_VELOCITY_THRESH = 2;

ncells = size(sessionData.trace, 2);
PCIs = zeros(size(1, ncells));

runningData = allData(allData.velocity > RUNNING_VELOCITY_THRESH, :);
binSize = 5;
cheeseboardSize = 100 / binSize;
for i = 1:ncells
    [ placeField, PCI ] = getPlaceField(...
        runningData.smooth_trans_x, runningData.smooth_trans_y,...
        runningData.trace(:, i), binSize);
    PCIs(i) = PCI;
    
    if PCI > 0.2
        figure;
        image(placeField, 'CDataMapping','scaled'), colorbar
        hold on;
        % draw cheeseboard
        cheeseboardCentre = (cheeseboardSize + 1) / 2;
        drawCircle(cheeseboardSize / 2, cheeseboardCentre, cheeseboardCentre);

        if ~isempty(rewardedWellsMap)
            drawCircle(0.5, ...
                    rewardedWellsMap.trans_x(1) / binSize,...
                    rewardedWellsMap.trans_y(1) / binSize,...
                    'Color', 'w');
            drawCircle(0.5, ...
                    rewardedWellsMap.trans_x(2) / binSize,...
                    rewardedWellsMap.trans_y(2) / binSize,...
                    'Color', 'w');
        end
        if ~isempty(negativeWellsMap)
             drawCircle(2, ...
                    negativeWellsMap.trans_x(1) / binSize,...
                    negativeWellsMap.trans_y(1) / binSize,...
                    'Color', 'r');
        end
        
        hold off;
        waitforbuttonpress
    end
end

%% Reward responsiveness
% Avg activity as function of reward distance
% bin reward distances
maxDist = 100;
rewDistBinSize = 10;
rew0 = allData.dist_reward0 / rewDistBinSize;
rew1 = allData.dist_reward1 / rewDistBinSize;

for i = 1:ncells
    if PCIs(i) > 0.2
        figure;
        plot(allData.dist_reward0, allData.trace(:,i), '*');
        hold on;
        plot(allData.dist_reward1, allData.trace(:,i), '*');
        xlabel('Dist to reward');
        hold off;
        waitforbuttonpress
    end
end

allData = addAlignedTrialTimestamps(allData);
%% 
% Calculates mean trace and mean events count before and after the reward
beforeRewardMs = 2000;
afterRewardMs = 2000;

arrivedAtRewardTimestamps = allData.timestamp(allData.arrivedAtReward > 0);
beforeIndecies = [];
afterIndecies = [];
for t_index = 1:numel(arrivedAtRewardTimestamps)
    beforeIndecies = [ beforeIndecies; ...
        find(allData.timestamp <= arrivedAtRewardTimestamps(t_index) ...
            & allData.timestamp >= arrivedAtRewardTimestamps(t_index) - beforeRewardMs)];
    afterIndecies = [ afterIndecies; ...
        find(allData.timestamp >= arrivedAtRewardTimestamps(t_index) ...
            & allData.timestamp <= arrivedAtRewardTimestamps(t_index) + afterRewardMs)];
end

beforeStats = grpstats(allData(beforeIndecies, {'trace', 'events', 'trial_id'}), ...
        'trial_id', {'mean'});
trialMeanTraceBefore = table2array(beforeStats(:, {'mean_trace'}));
trialMeanEventsBefore = table2array(beforeStats(:, {'mean_events'}));

afterStats = grpstats(allData(afterIndecies, {'trace', 'events', 'trial_id'}), ...
        'trial_id', {'mean'});
trialMeanTraceAfter = table2array(afterStats(:, {'mean_trace'}));
trialMeanEventsAfter = table2array(afterStats(:, {'mean_events'}));


% TODO: 
% - use interquartile interval for estimation of STD in find peaks
% - calculate place fields for events
% - calculate reward responsiveness of cells
% - save the place fields images
