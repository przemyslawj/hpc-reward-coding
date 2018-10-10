function [U] = mergevars(T, varPrefix)
%MERGEVARS Merges variables within the table sharing the common prefix into a
% single variable
varIndecies = startsWith(T.Properties.VariableNames, varPrefix);
U = T(:, ~varIndecies);

U.varPrefix = T{:, varIndecies};
U.Properties.VariableNames = [ T.Properties.VariableNames(~varIndecies), ...
    varPrefix ];
end

