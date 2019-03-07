function [eventVec, normApt] = findEvents( sessionTraces, numStdsForThresh, freq)
% Returns cell array with event times
% Event times are found using findpeak function. Only peaks higher than
% numStdsForThresh * std(signal) are found.
    minpeakdistance = 0.5 * freq;
    
    ncells = size(sessionTraces, 2);
    eventVec = zeros(size(sessionTraces), 'logical');
    normApt = zeros(size(eventVec));
    for i=1:ncells
        trace = sessionTraces(:,i);
        norm_trace = trace / std(trace(:));

        [~, peakTimes] = findpeaks(norm_trace,'minpeakheight', numStdsForThresh,...
                'minpeakdistance',minpeakdistance);
        eventVec(peakTimes,i) = 1;
        normApt(:,i) = norm_trace;
    end
    
end
