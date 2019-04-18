rootDir = '/mnt/DATA/Prez/cheeseboard/2019-02-learning/';
caimg_analysis_rootdir = [rootDir filesep 'joined_caimg'];
animal = '1BR';
session = 35;

caimg_analysis_dir = [caimg_analysis_rootdir filesep animal filesep 'jointExtraction' filesep 'sorted'];
sessionName = ['Session' num2str(session)];

% Load h5 file timestamps
h5file = [caimg_analysis_dir filesep '..' filesep 'alignment' filesep 'jointMovie.h5'];
info = h5info(h5file, '/sessionLengths');
sessionLengths = h5read(h5file, '/sessionLengths', [1 1], info.Dataspace.Size);

info = h5info(h5file, '/timestamps');
ts = h5read(h5file, '/timestamps', [1 1], info.Dataspace.Size);
timestampsBySession = mat2cell(ts, sessionLengths, 1);
caimg_timestamps = timestampsBySession{session};

% Load traces
sortedCellActivityFile = [caimg_analysis_dir filesep 'PCAICAsorted.mat'];
load(sortedCellActivityFile);
ncells = size(traces, 1);
tracesBySession = mat2cell(traces, size(traces, 1), sessionLengths);

% Load tracking positions
session_info_filepath = [caimg_analysis_rootdir filesep 'session_info_template.csv'];
sessionsInfo = readtable(session_info_filepath);
sessionMeta = sessionsInfo(strcmp(sessionsInfo.SessionName, sessionName), :);

% Load video
dateStr = datestr(sessionMeta.Date, 'yyyy-mm-dd');
tracking_vid_filename = [dateStr '_' animal '_trial_' num2str(sessionMeta.Trial) '.avi'];
movie_dir = [rootDir filesep dateStr filesep 'movie'];
tracking_vid_filepath = [movie_dir filesep tracking_vid_filename];
tracking_vid = readAvi(tracking_vid_filepath);

% Load tracking
trackingFile = [ dateStr '_' animal '_' 'trial_' num2str(sessionMeta.Trial) '_positions.csv' ];
trackingFilepath = [movie_dir filesep 'tracking' filesep trackingFile];
opts = detectImportOptions(trackingFilepath);
index = find(cellfun(@(x) strcmp(x, 'inside_roi'), opts.VariableNames, 'UniformOutput', 1));
opts.VariableTypes(index) = { 'logical' };

trialPositions = readtable(trackingFilepath, opts); 

addpath('~/code/miniscope_dat_analysis/misc/')
cumLenghts = double([1 cumsum(sessionLengths)]);
processedCamovieMat = loadMovie(h5file, cumLenghts(session), cumLenghts(session+1));

h5file = [caimg_analysis_rootdir filesep animal filesep sessionName filesep 'raw' filesep 'rawMovie.h5'];
rawCamovieMat = loadMovie(h5file);
X=permute(rawCamovieMat, [3, 1, 2]);
rawCamovieMat=permute(downsample(X, 2), [2, 3, 1]);

playMovies(0, 1000 * 100, 10.0, tracking_vid_filename, ...
    tracesBySession{session},...
    tracking_vid, trialPositions.timestamp,...
    rawCamovieMat, caimg_timestamps,...
    processedCamovieMat, caimg_timestamps);



