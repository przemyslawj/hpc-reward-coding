addpath(genpath(pwd)) 

rootDir = '/mnt/DATA/Prez/conditioning/';
caimg_analysis_rootdir = '/mnt/DATA/Prez/conditioning/results_all_caimg/';
animal = 'F';
sessions = 1:9;

EVENT_THRESH_NUM_STD = 4;


%% Load joint ca img data
%caimg_analysis_dir = '/media/prez/DATA/Prez/ca_img/sample_data/results_single_session3/M1/Session1/sorted';
caimg_analysis_dir = [caimg_analysis_rootdir filesep animal filesep 'jointExtraction/sorted'];
%sortedCellActivityFile = [caimg_analysis_dir filesep 'intermediateAnnotationResult.mat'];
sortedCellActivityFile = [caimg_analysis_dir filesep 'PCAICAsorted.mat'];
% loads variables traces, valid, filters
load(sortedCellActivityFile);

%h5file = [caimg_analysis_dir filesep '..' filesep 'preprocessed' filesep 'preprocessedMovie.h5'];
h5file = [caimg_analysis_dir filesep '..' filesep 'alignment' filesep 'jointMovie.h5'];

info = h5info(h5file, '/timestamps');
ts = h5read(h5file, '/timestamps', [1 1], info.Dataspace.Size);

info = h5info(h5file, '/sessionLengths');
sessionLengths = h5read(h5file, '/sessionLengths', [1 1], info.Dataspace.Size);
timestampsBySession = mat2cell(ts, sessionLengths, 1);

%% Load behavioural data for each session
sessionsInfo = readtable([caimg_analysis_rootdir filesep 'session_info.csv']);

allData = [];

for session = sessions
    sessionName = ['Session' num2str(session)];
    sessionMeta = sessionsInfo(strcmp(sessionsInfo.SessionName, sessionName), :);
    dateStr = datestr(sessionMeta.Date, 'yyyy-mm-dd');
    datedRootDir = [ rootDir filesep dateStr ];
    trackingDir = [ datedRootDir filesep 'movie' filesep 'tracking' ];
    tracesBySession = mat2cell(traces, size(traces, 1), sessionLengths);
    sessionTimestamps = timestampsBySession{session}';
    freq = (sessionTimestamps(end) - sessionTimestamps(2)) / numel(sessionTimestamps);
    
    trackingFile = [ dateStr '_' animal '_' 'trial_' num2str(sessionMeta.Trial) '_positions.csv' ];

    trackingFilepath = [trackingDir filesep trackingFile]
    opts = detectImportOptions(trackingFilepath);
    index = find(cellfun(@(x) strcmp(x, 'inside_roi'), opts.VariableNames, 'UniformOutput', 1));
    opts.VariableTypes(index) = { 'logical' };

    trialPositions = readtable(trackingFilepath, opts);

    sessionTraces = tracesBySession{session};
    if size(sessionTraces,1) == size(valid,1)
        sessionTraces = sessionTraces(valid == 1, :);
    end

    sessionData = mergeByTimestamp(trialPositions, sessionTraces, sessionTimestamps);
    taskStartedIndecies = find(sessionData.smooth_trans_x > -100 | sessionData.smooth_trans_y > -100);
    sessionData = sessionData(taskStartedIndecies,:);
    sessionData = calculateVelocity(sessionData);
    sessionData.atReward0 = isAtReward(sessionData.velocity, sessionData.dist_reward0);
    sessionData.atReward1 = isAtReward(sessionData.velocity, sessionData.dist_reward1);
    
    n = size(sessionData.atReward0, 1);
    sessionData.arrivedAtReward = zeros(n, 1);
    sessionData.arrivedAtReward(find(sessionData.atReward0, 1, 'first')) = 1;
    sessionData.arrivedAtReward(find(sessionData.atReward1, 1, 'first')) = 2;
    trial_id = [ dateStr '_' num2str(sessionMeta.Trial) ];
    sessionData.trial_id = mat2cell(repmat(trial_id, n, 1), ones(n, 1), numel(trial_id));
    sessionData.date = repmat(sessionMeta.Date, n, 1);
    sessionData.trial = repmat(sessionMeta.Trial, n, 1);
    fe = findEvents(sessionData.trace, EVENT_THRESH_NUM_STD, freq);
    sessionData.events = fe;
    
    if isempty(allData)
        allData = sessionData;
    else
        allData = [allData; sessionData];
    end
end

trial_data_path = [caimg_analysis_rootdir filesep 'traces_and_positions.csv'];
writetable(allData, trial_data_path);


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
    beforeIndecies = [ beforeIndecies ...
        find(allData.timestamp <= arrivedAtRewardTimestamps(t_index) ...
            & allData.timestamp >= arrivedAtRewardTimestamps(t_index) - beforeRewardMs)];
    afterIndecies = [ afterIndecies ...
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

%% 
% Evaluate with anova difference in traces and events before and after the reward
for i = 1:ncells
    [p_trace, stats] = anova1([trialMeanTraceBefore(:, i), trialMeanTraceAfter(:, i)], ...
            {'before', 'atReward'}, 'off');
        
    % TODO: test other test than anova, since variance will be
    % zero in case of 0 events.
    [p_events, stats] = anova1([trialMeanEventsBefore(:, i), trialMeanEventsAfter(:, i)], ...
        {'before', 'atReward'}, 'off');
    if p_trace < 0.3
        i
        anova1([trialMeanTraceBefore(:, i), trialMeanTraceAfter(:, i)], ...
            {'before', 'atReward'}, 'on');
        figure;
        plotCellTraces(allData, i);
    end
end

% TODO: 
% - use interquartile interval for estimation of STD in find peaks
% - calculate place fields for events
% - calculate reward responsiveness of cells
% - save the place fields images
