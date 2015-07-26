%% function to get indices (returned as 'strInd') of elements in a cell array of strings
% ('cellArray') that match with 'str'
% Dependencies: NONE
% Used by: GetType, GetUnits, ProcessBSNEs, StressFluxBSNE,
% VelocityProfileAnalysis, TimeIntervalUnion, find_indices

function strInd = CellStrFind(cellArray,str)

%s = strfind(cellArray,str); %get cellArray elements with string
%strInd = find(~cellfun(@isempty,s)); %get index of matching element
strInd = find(strcmp(cellArray,str));