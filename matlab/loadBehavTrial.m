addpath(genpath(pwd)) 

rootDir = '/mnt/DATA/Prez/cheeseboard/2019-02-learning/';
caimg_analysis_rootdir = '/mnt/DATA/Prez/cheeseboard/2019-02-learning/joined_caimg';
animal = '1BR';
sessions = 1:81;

EVENT_THRESH_NUM_STD = 4;


%% Load joint ca img data
%caimg_analysis_dir = '/media/prez/DATA/Prez/ca_img/sample_data/results_single_session3/M1/Session1/sorted';
caimg_analysis_dir = [caimg_analysis_rootdir filesep animal filesep 'jointExtraction/sorted'];
%sortedCellActivityFile = [caimg_analysis_dir filesep 'intermediateAnnotationResult.mat'];
sortedCellActivityFile = [caimg_analysis_dir filesep 'PCAICAsorted.mat'];
% loads variables traces, valid, filters
load(sortedCellActivityFile);
ncells = size(traces, 1);

%h5file = [caimg_analysis_dir filesep '..' filesep 'preprocessed' filesep 'preprocessedMovie.h5'];
h5file = [caimg_analysis_dir filesep '..' filesep 'alignment' filesep 'jointMovie.h5'];

info = h5info(h5file, '/timestamps');
ts = h5read(h5file, '/timestamps', [1 1], info.Dataspace.Size);

info = h5info(h5file, '/sessionLengths');
sessionLengths = h5read(h5file, '/sessionLengths', [1 1], info.Dataspace.Size);
timestampsBySession = mat2cell(ts, sessionLengths, 1);

%% Load behavioural data for each session
session_info_filepath = [caimg_analysis_rootdir filesep 'session_info_template.csv'];
sessionsInfo = readtable(session_info_filepath);

allData = [];

for session = sessions
    sessionName = ['Session' num2str(session)];
    sessionMeta = sessionsInfo(strcmp(sessionsInfo.SessionName, sessionName), :);
    dateStr = datestr(sessionMeta.Date, 'yyyy-mm-dd');
    is_test_str = '';
    is_test = 0;
    if strcmp(sessionMeta.is_test{1}, 'TRUE')
        is_test_str = '_test';
        is_test = 1;
    end
    datedRootDir = [ rootDir filesep dateStr is_test_str];
    trackingDir = [ datedRootDir filesep 'movie' filesep 'tracking' ];
    tracesBySession = mat2cell(traces, size(traces, 1), sessionLengths);
    sessionTimestamps = timestampsBySession{session}';
    freq = 1000 / median(diff(sessionTimestamps));
    
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
    taskStartedIndecies = find(sessionData.smooth_trans_x > 0 | sessionData.smooth_trans_y > 0);
    sessionData = sessionData(taskStartedIndecies,:);
    sessionData = calculateVelocity(sessionData);
    
    n = size(sessionData, 1);

    trial_id = [ dateStr '_' num2str(sessionMeta.Trial) ];
    if is_test
        trial_id = [trial_id '_test' ];
    end
    sessionData.date = repmat(sessionMeta.Date, n, 1);
    sessionData.trial = repmat(sessionMeta.Trial, n, 1);
    sessionData.trial_id = mat2cell(repmat(trial_id, n, 1), ones(n, 1), numel(trial_id));
    sessionData.is_test = repmat(is_test, n, 1);
    [eventsVec, normTrace] = findEvents(sessionData.trace, EVENT_THRESH_NUM_STD, freq);
    sessionData.events = eventsVec;
    sessionData.normTrace = normTrace;
    
    if isempty(allData)
        allData = sessionData;
    else
        allData = [allData; sessionData];
    end
end

trial_data_path = [caimg_analysis_rootdir filesep 'traces_and_positions.csv'];
writetable(allData, trial_data_path);

%% Events stats
for i = 1:ncells
    nevents = sum(allData.events(:,i));
    disp(['Cell ' num2str(i) ' events: ' num2str(nevents)])
end

%% Calculate place fields
% TODO: evaluate different velocity thresholds
RUNNING_VELOCITY_THRESH = -1;

ncells = size(sessionData.trace, 2);
PCIs = zeros(size(1, ncells));

runningData = allData(allData.velocity > RUNNING_VELOCITY_THRESH, :);
binSize = 10;
for i = 1:ncells
    [ placeField, PCI, occupancyMap, eventMap ] = getPlaceField(...
        runningData.smooth_trans_x, runningData.smooth_trans_y,...
        runningData.events(:, i), binSize);
    PCIs(i) = PCI;
    
    if PCI > 1.2
        nevents = sum(allData.events(:,i));
        figure('Name', ['Cell ' num2str(i) ', nevents: ' num2str(nevents) ...
                        ', PCI: ' num2str(PCI, 2)]);
        %subplot(1, 3, 1);
        %image(occupancyMap, 'CDataMapping','scaled'), colorbar
        %hold on;
        %subplot(1, 3, 2);
        %image(eventMap, 'CDataMapping','scaled'), colorbar
        %subplot(1, 3, 3);
        image(placeField, 'CDataMapping','scaled'), colorbar
     
        
        hold off;
        waitforbuttonpress
    end
end


%% 

% TODO: 
% - use interquartile interval for estimation of STD in find peaks
% - calculate place fields for events
% - calculate reward responsiveness of cells
% - save the place fields images
