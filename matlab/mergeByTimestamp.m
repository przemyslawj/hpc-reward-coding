function [resultTable] = mergeByTimestamp(behavData, trace, traceTimestamps)
% Merges trace into table with behaviour data using timestamps.
%   behaviourData - table with timestamp column
%   trace - vector with dF/F values
%   traceTimestamp - vector with timestamp values for each trace value

ntimestamps = numel(traceTimestamps);
assert(size(trace, 2) == ntimestamps);


function avg = weightedAvg(datatable, indecies, colName, weights)
    if ismember(colName, datatable.Properties.VariableNames)
        values = datatable{indecies, colName};
        if values(1) <= 0
            avg = values(2);
        else
            avg = double(weights) * values / sum(weights);
        end
    else
        avg = 0;
    end
end

function val = closerVal(datatable, indecies, colName, time_diff)
    if ismember(colName, datatable.Properties.VariableNames)
        [~, I] = min(time_diff);
        values = datatable{indecies, colName};
        val = values(I);
    else
        val = 0;
    end
end


avgedVariables = intersect({'smooth_trans_x', 'smooth_trans_y', ...
                            'smooth_heading_angle',...
                            'dist_reward0', 'dist_reward1'}, ...
                          behavData.Properties.VariableNames);
chooseCloserVariables = {'inside_roi', 'is_headdip'};

traceIndex = 1;
behavIndex = 1;
% Skip caimg trace before the tracking positions start
while traceIndex <= ntimestamps &&...
        behavData{behavIndex, 'timestamp'} > traceTimestamps(traceIndex)
    traceIndex = traceIndex + 1;
end
skippedTraceIndecies = traceIndex - 1;

avgedVals = zeros(ntimestamps - traceIndex + 1, ...
    numel(avgedVariables) + numel(chooseCloserVariables));

savedTraceValueIndex = 1;
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
        avgedVals(savedTraceValueIndex, i) = weightedAvg(behavData,...
                [prevBehavIndex behavIndex], avgedVariables{i}, timestampDiff);
    end
    for i = 1:numel(chooseCloserVariables)
        avgedVals(savedTraceValueIndex, numel(avgedVariables) + i) = closerVal(behavData,...
                [prevBehavIndex behavIndex], chooseCloserVariables{i}, timestampDiff);
    end

    traceIndex = traceIndex + 1;
    savedTraceValueIndex = savedTraceValueIndex + 1;
end

resultTable = table(traceTimestamps(skippedTraceIndecies+1:end)', ...
                    trace(:,skippedTraceIndecies+1:end)', ...
                    'VariableNames', {'timestamp', 'trace'});
for i = 1:numel(avgedVariables)
    resultTable{:,avgedVariables{i}} = avgedVals(:,i);
end
for i = 1:numel(chooseCloserVariables)
    j = numel(avgedVariables) + i;
    resultTable{:,chooseCloserVariables{i}} = avgedVals(:,j);
end

end

