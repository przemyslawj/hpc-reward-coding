addpath(genpath('/mnt/DATA/Prez/code/cheeseboard_analysis/matlab/')) 

rootDir = '/mnt/DATA/Prez/cheeseboard/2019-08/2019-08-habituation/';
dateStr = '2019-08-28';
datedRootDir = fullfile(rootDir, dateStr);
animal = 'E-BL';

EVENT_THRESH_NUM_STD = 5;
is_habituation = strfind(rootDir, 'habituation') > 0;

%% Load ca img data
%caimg_analysis_dir = [datedRootDir filesep 'cnmfe' filesep animal];
caimg_analysis_dir = [datedRootDir filesep 'trials' filesep 'mv_caimg' filesep animal]
load([caimg_analysis_dir filesep 'ms.mat']);

sessionDirs = dir(fullfile(datedRootDir, 'trials', 'mv_caimg', animal));
sessionDirs = sessionDirs(cellfun(...
    @startsWith, {sessionDirs.name}, repmat({'Session'},1,numel(sessionDirs))));

%% Align timestamps in case of missing frames
if iscell(ms.sessionLengths)
    sessionLengths = cell2mat(ms.sessionLengths)';
else
    sessionLengths = ms.sessionLengths';
end
if numel(ms.time) > size(ms.RawTraces, 1)
    sessionStartsIndex = find(diff(ms.time)<0) + 1;
    sessionStartsIndex = [1; sessionStartsIndex];
    sessionEndsIndex = sessionStartsIndex - 1 + sessionLengths;
    
    time2 = [];
    for i = 1:numel(sessionStartsIndex)
        time2 = [time2; ms.time(sessionStartsIndex(i):sessionEndsIndex(i))];
    end
    ms.time = time2;
end
ts = ms.time;

timestampsBySession = mat2cell(ts, sessionLengths, 1);

%% Load behavioural data for each session

allData = [];

for session_i = 1:numel(sessionDirs)
    trackingDir = fullfile(datedRootDir, 'trials', 'movie', 'tracking');
    tracesBySession = mat2cell(ms.RawTraces, sessionLengths);
    sessionTimestamps = timestampsBySession{session_i}';
    freq = (sessionTimestamps(end) - sessionTimestamps(2)) / numel(sessionTimestamps);
    sessionNo = str2num(replace(sessionDirs(session_i).name, 'Session', ''));
    
    trackingFile = [ dateStr '_' animal '_' 'trial_' num2str(sessionNo) '_positions.csv' ];

    trackingFilepath = [trackingDir filesep trackingFile]
    opts = detectImportOptions(trackingFilepath);
    index = find(cellfun(@(x) strcmp(x, 'inside_roi'), opts.VariableNames, 'UniformOutput', 1));
    opts.VariableTypes(index) = { 'logical' };

    trialPositions = readtable(trackingFilepath, opts);

    sessionTraces = (tracesBySession{session_i})';

    sessionData = mergeByTimestamp(trialPositions, sessionTraces, sessionTimestamps);
    taskStartedIndecies = find(sessionData.smooth_trans_x > -100 | sessionData.smooth_trans_y > -100);
    sessionData = sessionData(taskStartedIndecies,:);
    sessionData = calculateVelocity(sessionData);
    n = size(sessionData, 1);
    if ismember('dist_reward0', sessionData.Properties.VariableNames)
        sessionData.atReward0 = isAtReward(sessionData.velocity, sessionData.dist_reward0);
        sessionData.arrivedAtReward = zeros(n, 1);
        sessionData.arrivedAtReward(find(sessionData.atReward0, 1, 'first')) = 1;
    end
    if ismember('dist_reward1', sessionData.Properties.VariableNames)
        sessionData.atReward1 = isAtReward(sessionData.velocity, sessionData.dist_reward1);
        sessionData.arrivedAtReward(find(sessionData.atReward1, 1, 'first')) = 2;
    end
    
    trial_id = [ dateStr '_' num2str(sessionNo) ];
    sessionData.trial_id = mat2cell(repmat(trial_id, n, 1), ones(n, 1), numel(trial_id));
    sessionData.date = repmat(dateStr, n, 1);
    sessionData.trial = repmat(sessionNo, n, 1);
    fe = findEvents(sessionData.trace, EVENT_THRESH_NUM_STD, freq);
    sessionData.events = fe;
    
    if isempty(allData)
        allData = sessionData;
    else
        allData = [allData; sessionData];
    end
end

trial_data_path = [caimg_analysis_dir filesep 'traces_and_positions.csv'];
writetable(allData, trial_data_path);

