function [eventVec, normApt, thresholds] = findEvents( sessionTraces, numStdsForThresh, freq)
% Returns cell array with event times
% Event times are found using findpeak function. Only peaks higher than
% numStdsForThresh * std(signal) are found.
    
    ncells = size(sessionTraces, 2);
    eventVec = zeros(size(sessionTraces), 'logical');
    normApt = zeros(size(eventVec));
    thresholds = zeros(ncells,1);
    for i=1:ncells
        trace = sessionTraces(:,i);
        trace = trace - mean(trace);
        std_val = iqr(trace) / 1.349;
        thresholds(i) = std_val * numStdsForThresh;
        norm_trace = trace / std_val;

        cell_event_times = findCellEvents(norm_trace, numStdsForThresh, freq);
        eventVec(cell_event_times(:,1),i) = 1;
        normApt(:,i) = norm_trace;
    end

end

function [cell_event_times] = findCellEvents(norm_trace, std_thr, freq)
    minpeakdistance = 2 * freq;
    [~, peakTimes] = findpeaks(norm_trace,...
        'minpeakheight', std_thr,...
        'minpeakdistance', minpeakdistance,...
        'minpeakprominence', 3,...
        'minpeakwidth', 1);

    thresholded = norm_trace > (std_thr * 0.75);
    cell_event_times = zeros(numel(peakTimes), 3);
    events_to_remove = [];
    prev_peak_i = -1;
    for peak_i = 1:numel(peakTimes)
        left_i = peakTimes(peak_i);
        while thresholded(left_i) > 0 && left_i > 1
            left_i = left_i - 1;
        end
        right_i = peakTimes(peak_i) + 1;
        while thresholded(right_i) > 0 && right_i < size(norm_trace,1)
            right_i = right_i + 1;
        end
        cell_event_times(peak_i,:) = [peakTimes(peak_i), left_i, right_i];

        % Remove smaller event if the previous event overlaps
        if prev_peak_i > 0 && cell_event_times(prev_peak_i, 3) > left_i
            if norm_trace(peakTimes(prev_peak_i)) > ...
                    norm_trace(peakTimes(peak_i))
                to_remove = peak_i;
            else
                to_remove = prev_peak_i;
                prev_peak_i = peak_i;
            end
            events_to_remove = [events_to_remove, to_remove];
        else
            prev_peak_i = peak_i;
        end
    end
    cell_event_times(events_to_remove,:) = [];
    peakTimes(events_to_remove) = [];
end

