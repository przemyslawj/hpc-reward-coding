function [ field, PCI ] = getPlaceField(posX, posY, trace, binSize)
% Calculates spatial heatmap of cell activity and Place Cell Information
% defined as bits of information per 1% change in fluorescence. The metric
% translated for calcium imaging from the original definition in Skiggs et 
% al 1993.
% 
% TODO: consider calculating MI measure which might be more meaningful:
% https://www.biorxiv.org/content/biorxiv/early/2017/09/15/189084.full.pdf
%

if ~exist('binSize','var') || isempty(binSize)
    binSize = 5;
end    

minX = 0;
minY = 0;
posX = posX - minX;
posY = posY - minY;

maxX = 100 - minX;
maxY = 100 - minY;

binnedPosX = max(int32(round(posX / binSize)), 1);
binnedPosY = max(int32(round(posY / binSize)), 1);

trace = trace - min(trace);

totalActivityMap = zeros(ceil(maxY / binSize), ceil(maxX / binSize));
occupancyMap = zeros(size(totalActivityMap)) - 1;

for i = 1:numel(trace)
    x = binnedPosX(i);
    y = binnedPosY(i);

    if occupancyMap(x, y) < 0
        occupancyMap(x, y) = 0;
    end
    occupancyMap(x, y) = occupancyMap(x, y) + 1;
    totalActivityMap(x, y) = totalActivityMap(x, y) + trace(i);
end

PCI = 0;
mfr = mean(trace);
for x = 1:size(totalActivityMap, 1)
    for y = 1:size(totalActivityMap, 2)
        occupancyProb = max(0, occupancyMap(x, y)) / numel(trace);
        fr = totalActivityMap(x, y) / occupancyMap(x, y);
        if occupancyProb > 0 && fr > 0
            PCI = PCI + occupancyProb * fr * log2(fr / mfr);
        end
    end
end
if mfr > 0
    PCI = PCI / mfr;
end

field = (totalActivityMap ./ occupancyMap)';
end
