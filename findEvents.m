function [cellPeaks] = findEvents( sessionTraces, numStdsForThresh, freq)
% Returns cell array with event times
% Event times are found using findpeak function. Only peaks higher than
% numStdsForThresh * std(signal) are found.
    minpeakdistance = 0.5 * freq;
    
    ncells = size(sessionTraces, 2);
    cellPeaks = zeros(size(sessionTraces), 'logical');
    for i=1:ncells
        trace = sessionTraces(:,i);
        threshold = std(trace(:)) * numStdsForThresh;

        [~,peakTimes] = findpeaks(trace,'minpeakheight', threshold,...
                'minpeakdistance',minpeakdistance);
        cellPeaks(peakTimes,i) = 1;
    end
    
end