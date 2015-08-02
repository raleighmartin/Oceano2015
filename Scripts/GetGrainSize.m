%% function to take info about grain size collections
% and use this to crunch through outputs of Camsizer analyses
% Dependencies: 
% Used by: Processing_Master

function GrainSizeArray = GetGrainSize(GrainSizeMetadata,folder_GrainSize)

%% initialize structured array of grain size data based on metadata structured array
GrainSizeArray = GrainSizeMetadata;

% get number of samples
N_Samples = length(GrainSizeMetadata);

%% go through each flux calculation time interval
for i = 1:N_Samples;
    
    %% Get file location
    filePath = strcat(folder_GrainSize,GrainSizeMetadata(i).Filename,'_xc_min_.xle')
    
    %process file
    [d_10_mm, d_50_mm, d_90_mm, Q3_1mm, Q3_2mm, Q3_4mm, SPAN3, U3, ...
    Q3_SPHT0pt9, Q3_Symm0pt9, Q3_bl0pt9, SPHT3_bar, Symm3_bar, bl3_bar, ...
    Sizeclass_lower_mm,Sizeclass_upper_mm,retained,passing,...
    SPHT3,Symm3,bl3,particlesdetected] = gsimport(filePath);

    %compute mid-point for sizeclass bins
    Sizeclass_mid_mm = sqrt(Sizeclass_lower_mm.*Sizeclass_upper_mm);

    %generate structured array for grain-size distribution
    gsd = struct('Sizeclass_lower_mm',num2cell(Sizeclass_lower_mm),...
        'Sizeclass_upper_mm',num2cell(Sizeclass_upper_mm),...
        'Sizeclass_mid_mm',num2cell(Sizeclass_mid_mm),...
        'retained',num2cell(retained),...
        'passing',num2cell(passing),...
        'SPHT3',num2cell(SPHT3),...
        'Symm3',num2cell(Symm3),...
        'bl3',num2cell(bl3),...
        'particlesdetected',num2cell(particlesdetected));

    %add elements to structured array
    GrainSizeArray(i).d_10_mm = d_10_mm;
    GrainSizeArray(i).d_50_mm = d_50_mm;
    GrainSizeArray(i).d_90_mm = d_90_mm;
    GrainSizeArray(i).Q3_1mm = Q3_1mm;
    GrainSizeArray(i).Q3_2mm = Q3_2mm;
    GrainSizeArray(i).Q3_4mm = Q3_4mm;
    GrainSizeArray(i).SPAN3 = SPAN3;
    GrainSizeArray(i).U3 = U3;
    GrainSizeArray(i).Q3_SPHT0pt9 = Q3_SPHT0pt9;
    GrainSizeArray(i).Q3_Symm0pt9 = Q3_Symm0pt9;
    GrainSizeArray(i).Q3_bl0pt9 = Q3_bl0pt9;
    GrainSizeArray(i).SPHT3_bar = SPHT3_bar;
    GrainSizeArray(i).Symm3_bar = Symm3_bar;
    GrainSizeArray(i).bl3_bar = bl3_bar;
    GrainSizeArray(i).gsd = gsd;
end