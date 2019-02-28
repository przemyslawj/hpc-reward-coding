function [resultTable] = mergeByTimestamp(behavData, trace, traceTimestamps)
% Merges trace into table with behaviour data using timestamps.
%   behaviourData - table with timestamp column
%   trace - vector with dF/F values
%   traceTimestamp - vector with timestamp values for each trace value

ntimestamps = numel(traceTimestamps);
assert(size(trace, 2) == ntimestamps);


function avg = weightedAvg(datatable, indecies, colName, weights)
    values = datatable{indecies, colName};
    avg = double(weights) * values / sum(weights);
end

function val = closerVal(datatable, indecies, colName, time_diff)
    [~, I] = min(time_diff);
    values = datatable{indecies, colName};
    val = values(I);
end

traceIndex = 1;
behavIndex = 1;

avgedVariables = intersect({'smooth_trans_x', 'smooth_trans_y', ...
                           'dist_reward0', 'dist_reward1', 'inside_roi'}, ...
                          behavData.Properties.VariableNames);
avgedVals = zeros(ntimestamps, numel(avgedVariables)); 

while traceIndex <= ntimestamps
    while (behavIndex <= size(behavData, 1)) && ...
            (behavData{behavIndex, 'timestamp'} < traceTimestamps(traceIndex))
        behavIndex = behavIndex + 1;
    end
    
    if behavIndex >= size(behavData, 1)
        behavIndex = size(behavData, 1);
    end
    
    prevBehavIndex = behavIndex;
    timestampDiff = [1 1];
    if behavData{behavIndex, 'timestamp'} > traceTimestamps(traceIndex)
        prevBehavIndex = behavIndex - 1;
        timestampDiff = [...
            behavData{behavIndex, 'timestamp'} - traceTimestamps(traceIndex), ...
            traceTimestamps(traceIndex) - behavData{prevBehavIndex, 'timestamp'}
        ];
    end
    
    for i = 1:numel(avgedVariables)
        avgedVals(traceIndex, i) = closerVal(behavData,...
                [prevBehavIndex behavIndex], avgedVariables{i}, timestampDiff);
    end
        
    traceIndex = traceIndex + 1;
end

resultTable = table(traceTimestamps', trace', 'VariableNames', ...
        {'timestamp', 'trace'});
for i = 1:numel(avgedVariables)
    resultTable{:,avgedVariables{i}} = avgedVals(:,i);
end

end

