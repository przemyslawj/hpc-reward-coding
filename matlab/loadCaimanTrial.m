%% Merges calcium traces with position tracking producing .csv
% Requres following variables in the env:
% rootDir, dateStr, animal

addpath(genpath('matlab'))

datedRootDir = fullfile(rootDir, dateStr);

%% Load ca img data
caimg_analysis_dir = [datedRootDir filesep 'caiman' filesep animal]
caimg_result_file = [caimg_analysis_dir filesep 'filtered' filesep 'ms.mat'];
if exist(caimg_result_file, 'file')
    load(caimg_result_file);
else
	sprintf('Input file does not exist: %s', caimg_result_file)
    return
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

size_diff = size(ms.RawTraces, 1) - numel(ts);
if size_diff > 0
	avg_timestamp = floor(mean(diff(ts)));
	last_ts = ts(numel(ts)) + avg_timestamp;
	missing_ts = linspace(last_ts, size_diff*avg_timestamp, size_diff);
	ts = [ts missing_ts];
	warning('Too few timestamps, approximating timestamp based on frame rate')
end
timestampsBySession = mat2cell(ts', sessionLengths', 1);

%% Load session info
json_info_fpath= [caimg_analysis_dir filesep 'session_info.json'];
session_info_txt = fileread(json_info_fpath);
session_info = jsondecode(session_info_txt);

%% Load behavioural data for each session

allData = [];
trackingVars = {'inside_roi', 'smooth_trans_x', 'smooth_trans_y',...
                'velocity', 'dist_reward0', 'dist_reward1', ...
                'atReward0', 'atReward1', 'arrivedAtReward', ...
				'smooth_heading_angle', 'is_headdip'};

for session_i = 1:numel(session_info.session_fpaths)
    session_fpath_parts = split(session_info.session_fpaths{session_i}, filesep);

    traceBySession = mat2cell(ms.RawTraces, sessionLengths);
    deconvTraceBySession = mat2cell(ms.DeconvTraces, sessionLengths);
    sessionTimestamps = timestampsBySession{session_i}';

    %%
    if v3 == 1
        sessionNo = str2num(replace(session_fpath_parts{end}, 'Session', ''));
        exp_title = session_fpath_parts{end-3};
        if strcmp(exp_title, 'test')
            exp_title = 'beforetest';
        end
        trackingFilepath = getTrackingFilepathV3(datedRootDir, dateStr, ...
            animal, exp_title, sessionNo);
    else
        ts_file_parts = split(session_info.timestamp_files{session_i}, filesep);
        sessionNo = str2num(replace(ts_file_parts{end-3}, 'Session', ''));
        trackingFilepath = getTrackingFilepathV4(datedRootDir, ...
            dateStr, animal, ts_file_parts);
        exp_title = ts_file_parts{end-5};
    end
    if 	~exist(trackingFilepath, 'file')
        warning('No tracking file for session: %s', session_info.session_fpaths{session_i})
    end
    if strcmp(exp_title, 'homecage') || ...
           ~exist(trackingFilepath, 'file')
        sessionData = table();
        sessionData.timestamp = sessionTimestamps';
        sessionData.trace = traceBySession{session_i};
        sessionData.deconvTrace = deconvTraceBySession{session_i};
    else
        opts = detectImportOptions(trackingFilepath);
        index = find(cellfun(@(x) strcmp(x, 'inside_roi'), opts.VariableNames, 'UniformOutput', 1));
        opts.VariableTypes(index) = { 'logical' };
        index = find(cellfun(@(x) strcmp(x, 'is_headdip'), opts.VariableNames, 'UniformOutput', 1));
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
	sprintf('Zeroing tracking variables: ')
	otherTrackingVars
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

%%
function trackingFilepath = getTrackingFilepathV3(datedRootDir,...
        dateStr, animal, exp_title, sessionNo)
	trackingDir = fullfile(datedRootDir, exp_title, 'movie', 'tracking');
    filenamepattern = [ dateStr '*' '_' animal '_' 'trial_' num2str(sessionNo) '_positions.csv' ];
    trackfiles = dir(fullfile(trackingDir, filenamepattern));
    if numel(trackfiles) > 0
        trackingFilepath = fullfile(trackingDir, trackfiles(1).name)
    else
        trackingFilepath = '/non-existant-path/';
    end
end

function trackingFilepath = getTrackingFilepathV4(datedRootDir, ...
        dateStr, animal, ts_file_parts)
    trackingDir = fullfile(datedRootDir, ts_file_parts{(end-5):(end-2)},...
        'BehavCam', 'tracking');
    filenamepattern = [ dateStr '*' '_' animal '_trial_*_positions.csv' ];
    trackfiles = dir(fullfile(trackingDir, filenamepattern));
    if numel(trackfiles) > 0
        trackingFilepath = fullfile(trackingDir, trackfiles(1).name)
    else 
        trackingFilepath = '/non-existant-path/';
    end
end
