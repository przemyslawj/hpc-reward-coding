%% Merges calcium traces with position tracking producing .csv
% Requres following variables in the env:
% rootDir, dateStr, animal

EVENT_THRESH_NUM_STD = 5;
addpath(genpath('matlab'))

datedRootDir = fullfile(rootDir, dateStr);

%% Load ca img data
caimg_analysis_dir = [datedRootDir filesep 'caiman' filesep animal]
caimg_result_file = [caimg_analysis_dir filesep 'filtered' filesep 'ms.mat'];
if exist(caimg_result_file, 'file')
    load(caimg_result_file);
else
    exit
end

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

timestampsBySession = mat2cell(ts', sessionLengths', 1);

%% Load session info
json_info_fpath= [caimg_analysis_dir filesep 'session_info.json'];
session_info_txt = fileread(json_info_fpath);
session_info = jsondecode(session_info_txt);

%% Load behavioural data for each session

allData = [];
trackingVars = {'dist', 'inside_roi', 'smooth_trans_x', 'smooth_trans_y',...
                'velocity', 'dist_reward0', 'dist_reward1', ...
                'atReward0', 'atReward1', 'arrivedAtReward'};

for session_i = 1:numel(session_info.session_fpaths)
    session_fpath_parts = split(session_info.session_fpaths{session_i}, filesep);
    exp_title = session_fpath_parts{end-3};

    traceBySession = mat2cell(ms.RawTraces, sessionLengths);
    deconvTraceBySession = mat2cell(ms.DeconvTraces, sessionLengths);
    sessionTimestamps = timestampsBySession{session_i}';
    sessionNo = str2num(replace(session_fpath_parts{end}, 'Session', ''));

    if strcmp(exp_title, 'homecage')
        sessionData = table();
        sessionData.timestamp = sessionTimestamps';
        sessionData.trace = traceBySession{session_i};
        sessionData.deconvTrace = deconvTraceBySession{session_i};
    else
        trackingDir = fullfile(datedRootDir, 'trial', 'movie', 'tracking');
        trackingFile = [ dateStr '_' animal '_' 'trial_' num2str(sessionNo) '_positions.csv' ];
        trackingFilepath = [trackingDir filesep trackingFile]
        opts = detectImportOptions(trackingFilepath);
        index = find(cellfun(@(x) strcmp(x, 'inside_roi'), opts.VariableNames, 'UniformOutput', 1));
        opts.VariableTypes(index) = { 'logical' };

        trialPositions = readtable(trackingFilepath, opts);

        sessionTraces = (traceBySession{session_i})';

        sessionData = mergeByTimestamp(trialPositions, sessionTraces, sessionTimestamps);
        taskStartedIndecies = find(sessionData.smooth_trans_x > -100 | sessionData.smooth_trans_y > -100);

        sessionData = sessionData(taskStartedIndecies,:);
        sessionData = calculateVelocity(sessionData);
        deconvTraceData = mergeByTimestamp(trialPositions, (deconvTraceBySession{session_i})', sessionTimestamps);
        sessionData.deconvTrace = deconvTraceData.trace(taskStartedIndecies,:);

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
    end

    % Populate unknown tracking variables with zeroes
    n = size(sessionData, 1);
    otherTrackingVars = setdiff(trackingVars, sessionData.Properties.VariableNames);
    zeroValsTable = array2table(zeros(n, numel(otherTrackingVars)),...
        'VariableNames', otherTrackingVars);
    sessionData = [sessionData zeroValsTable];

    trial_id = [ animal '_' dateStr '_' exp_title '_' num2str(sessionNo) ];
    sessionData.trial_id = mat2cell(repmat(trial_id, n, 1), ones(n, 1), numel(trial_id));
    sessionData.date = repmat(dateStr, n, 1);
    sessionData.trial = repmat(sessionNo, n, 1);
    sessionData.exp_title = repmat({exp_title}, n, 1);
    sessionData.animal = repmat({animal}, n, 1);

    freq = (sessionTimestamps(end) - sessionTimestamps(2)) / numel(sessionTimestamps);
    fe = findEvents(sessionData.trace, EVENT_THRESH_NUM_STD, freq);
    sessionData.events = fe;

    if isempty(allData)
        allData = sessionData;
    else
        allData = [allData; sessionData];
    end
end

trial_data_path = [caimg_analysis_dir filesep 'filtered' filesep 'traces_and_positions.csv']
writetable(allData, trial_data_path);

%% Save cell Id mappings
cell_mapping=table();
cell_mapping.cell_no = (1:numel(ms.cellId))';
cell_mapping.cell_id = ms.cellId';
writetable(cell_mapping, [caimg_analysis_dir filesep 'filtered' filesep 'cell_mapping.csv']);

