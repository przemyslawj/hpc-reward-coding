function [resultTable] = mergeByTimestamp(behavData, trace, traceTimestamps)
% Merges trace into table with behaviour data using timestamps.
%   behaviourData - table with timestamp column
%   trace - vector with dF/F values
%   traceTimestamp - vector with timestamp values for each trace value

ntimestamps = numel(traceTimestamps);
assert(size(trace, 2) == ntimestamps);


trans_x = zeros(ntimestamps, 1);
trans_y = zeros(ntimestamps, 1);
dist_reward0 = zeros(ntimestamps, 1);
dist_reward1 = zeros(ntimestamps, 1);
inside_roi = zeros(ntimestamps, 1);

function avg = weightedAvg(datatable, indecies, colName, weights)
    values = datatable{indecies, colName};
    avg = double(weights) * values / sum(weights);
end

traceIndex = 1;
behavIndex = 1;
while traceIndex <= numel(traceTimestamps)
    while (behavIndex <= size(behavData, 1)) && ...
            (behavData{behavIndex, 'timestamp'} < traceTimestamps(traceIndex))
        behavIndex = behavIndex + 1;
    end
    
    if behavIndex >= size(behavData, 1)
        behavIndex = size(behavData, 1);
    end
    
    prevBehavIndex = behavIndex;
    avgWeights = [1 1];
    if behavData{behavIndex, 'timestamp'} > traceTimestamps(traceIndex)
        prevBehavIndex = behavIndex - 1;
        avgWeights = [...
            behavData{behavIndex, 'timestamp'} - traceTimestamps(traceIndex), ...
            traceTimestamps(traceIndex) - behavData{prevBehavIndex, 'timestamp'}
        ];
    end
    
    smooth_trans_x(traceIndex) = weightedAvg(behavData,...
        [prevBehavIndex behavIndex], 'smooth_trans_x', avgWeights);
    smooth_trans_y(traceIndex) = weightedAvg(behavData,...
        [prevBehavIndex behavIndex], 'smooth_trans_y', avgWeights);
    dist_reward0(traceIndex) = weightedAvg(behavData,...
        [prevBehavIndex behavIndex], 'dist_reward0', avgWeights);
    dist_reward1(traceIndex) = weightedAvg(behavData,...
        [prevBehavIndex behavIndex], 'dist_reward1', avgWeights);
    inside_roi(traceIndex) = any([behavData{prevBehavIndex, 'inside_roi'}(1),...
        behavData{behavIndex, 'inside_roi'}(1)]);
        
    traceIndex = traceIndex + 1;
end

resultTable = table(traceTimestamps', smooth_trans_x', smooth_trans_y', inside_roi, ...
    dist_reward0, dist_reward1, trace', 'VariableNames', ...
    {'timestamp', 'smooth_trans_x', 'smooth_trans_y', 'inside_roi', ...
     'dist_reward0', 'dist_reward1', 'trace'});

end

