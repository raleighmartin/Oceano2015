function [d_10, d_50, d_90, Q3_1mm, Q3_2mm, Q3_4mm, SPAN3, U3, ...
    Q3_SPHT0pt9, Q3_Symm0pt9, Q3_bl0pt9, SPHT3_bar, Symm3_bar, bl3_bar, ...
    Sizeclass_lower,Sizeclass_upper,retained,passing,...
    SPHT3,Symm3,bl3,particlesdetected] = gsimport(filename)
%function [] = gsimport(filename)

%%%%%%%%%%% SUMMARY DATA %%%%%%%%%%%%%
%% Initialize variables.
delimiter = '\t';
startRow = 2;
endRow = 15;
formatSpec = '%*s%f%*s%*s%*s%*s%*s%*s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
textscan(fileID, '%[^\n\r]', startRow(1)-1, 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    textscan(fileID, '%[^\n\r]', startRow(block)-1, 'ReturnOnError', false);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
    dataArray{1} = [dataArray{1};dataArrayBlock{1}];
end

%% Close the text file.
fclose(fileID);

%% Allocate imported array to column variable names
data = dataArray{:, 1};
d_10 = data(1); %10th percentile (by volume) sieve diameter (mm)
d_50 = data(2); %50th percentile (by volume) sieve diameter (mm)
d_90 = data(3); %90th percentile (by volume) sieve diameter (mm)
Q3_1mm = data(4); %percentile (by volume) corresponding to 1 mm sieve diameter
Q3_2mm = data(5); %percentile (by volume) corresponding to 2 mm sieve diameter
Q3_4mm = data(6); %percentile (by volume) corresponding to 4 mm sieve diameter
SPAN3 = data(7); %a measure of range of particle sizes
U3 = data(8); %"non-uniformity" = d_60/d_10
Q3_SPHT0pt9 = data(9);  %proportion (by volume) of particles with sphericity < 0.9
Q3_Symm0pt9 = data(10); %proportion (by volume) of particles with symmetry < 0.9
Q3_bl0pt9 = data(11);  %proportion (by volume) of particles with b/l ratio < 0.9
SPHT3_bar = data(12); %volume-weighted mean sphericity
Symm3_bar = data(13); %volume-weighted mean symmetry
bl3_bar = data(14); %volume-weighted mean b/l

%%%%%%%%%%% FULL DATA %%%%%%%%%%%%%%%%
%% Initialize variables.
delimiter = '\t';
startRow = 17;
endRow = inf;
formatSpec = '%f%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
textscan(fileID, '%[^\n\r]', startRow(1)-1, 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    textscan(fileID, '%[^\n\r]', startRow(block)-1, 'ReturnOnError', false);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Allocate imported array to column variable names
Sizeclass_lower = dataArray{:, 1};
Sizeclass_upper = dataArray{:, 2};
retained = dataArray{:, 3};
passing = dataArray{:, 4};
SPHT3 = dataArray{:, 5};
Symm3 = dataArray{:, 6};
bl3 = dataArray{:, 7};
particlesdetected = dataArray{:, 8};
