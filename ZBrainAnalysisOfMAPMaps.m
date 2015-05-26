% This function will quantify the ammount of 'significant delta median'
% from a MAP-Map signal in each 3D ROI in the Z-Brain database. These are
% then ranked, treating positive (activation) and negative (supression)
% signals separately. 

% We then loop through the Z-Brain anatomy labels, and ask, in each ROI,
% what is the average signal overlapping with the activated/supressed
% voxels. We then normalize this signal for the background/surrounding, by
% dividing by the mean signal in a 50 voxel ring surrounding the ROI.
% Numbers > 1 indicate enrichement of the label signal overlapping with the
% activity pattern. This can help generate hypotheses as to the relevant
% cell types.

%%%%%% Input arguments

% CalculateLabels
%%%% 1 (default), overlap with antomy labels will be calulcated
%%%% 0 - this calculation will be skipped

%%%%%% OUTPUT

% First will be a stack (*SignalInEachROIImage.tif), with each regions
% drawn with signal proportional to its intensity and scaled to the max of
% the data set. Green = activation, magneta = supression. ROIs are split at
% the middle of the stack, to reveal any strong left/rigth differences.
% These stacks are then made into a 2D image of the Z and X maximum
% projections (*SignalInEachROIProjection.tif)
%
% Two tables are output, one for each the activated and supressed regions,
% ranked for highest average signal within the ROI. The value for signal is
% a number between 0 and 65535, where 65535 would be saturated signal in
% the entire ROI. Also in these tables will be a ranking of the different
% transgenes/labels in the anatomy database for those who show the most
% signal within the activation signals in that ROI. This number represent
% the average brighness of the signal wihtin the activated region of that
% ROI divided by the average brigheness in the 50 pixels surrounding that
% ROI. Therefore, numbers above 1 indicate enrichment of that label in the
% active regions, whereas numbers less than 1 indicate that the label is
% brighter in the adjacent areas. Finally we will output a .mat file called
% '*SignalInEachROI.mat', which contains the names of the anatomical
% regions (in their normal alphabetical order), and the positive and
% negative signal within each anatomical region.

% Input files.

% This function uses the pre-downsampled versions of the Z-Brain contained
% in files 'AnatomyLabelDatabaseDownsampled.hdf5', and
% 'MaskDatabaseDownsampled.mat', and a '*SignificantDeltaMedians.tif'
% MAP-Map image output from the 'MakeTheMAPMap.m' function. The funciton
% starts with a promt to direct towards the MAPMap, and the other files if
% they are not in Matlab's path.
%
%
function ZBrainAnalysisOfMAPMaps(CalculateLabels)


if nargin ~=1
CalculateLabels = 1; % set to true if you want to scroll through the anatomy label database to look for overlapping signals
end

CalculateLabels = logical(CalculateLabels);

XYRez = 1.65; % the resolution of MAP-Maps
ZRez = 3.45;

% point to the delta medians file for analysis
[DeltaMedianFilename, DeltaMedianPath] = uigetfile('*Medians.tif', 'Select the significant delta median file for analysis');


try % if we are already in the correct directory, load in the mask database files. The 'AnatomyLabelDatabaseDownsampled.hdf5' file also needs to be in this directory
    
    load('MaskDatabaseDownsampled.mat');
    MaskDatabasePath = pwd;
    
catch
    MaskDatabasePath = uigetdir(pwd,'Point to the folder with the ZBrain files'); %'MaskDatabaseDownsampled.mat' and 'AnatomyLabelDatabaseDownsampled.hdf5' need to be in this folder
    cd(MaskDatabasePath);
    addpath(MaskDatabasePath);
    load('MaskDatabaseDownsampled.mat');
end

% remove the eyes since Nacre fish dont show real signal there.
for j = length(MaskDatabaseNames):-1:1
    if strcmp('Ganglia - Eyes', MaskDatabaseNames{j})
        MaskDatabaseNames(j) = [];
        MaskDatabaseDownsamp(:,j) = [];
        MaskDatabase50OutsideDownsamp(:,j) = [];
        MaskDatabaseOutlinesDownsamp(:,j) = [];
    end
end

% remove any commas from mask names, as this messes up the csv output file

for j = 1:length(MaskDatabaseNames)
        MaskDatabaseNames{j} = strrep(MaskDatabaseNames{j},',', '-');
end

cd(DeltaMedianPath);

info = imfinfo(DeltaMedianFilename);

Zs = length(info);
width = info(1).Width;
height = info(1).Height;

DeltaMedianImage = uint16(zeros(height, width, Zs, 3));

for z = 1:Zs
    DeltaMedianImage(:,:,z,:) = imread(DeltaMedianFilename, z);
end

%
% To find the most relevant ROIs we determine the 'significant delta
% median' signal within each roi. we remove the top ROI and re-run the
% analysis with that signal missing. This will avoid redundantly calling
% two ROIs that overlap with a single focus of signal we also loop through
% the transgene/label database and find which label shows the greatest
% signal within the calculated ROI region. This is weighted based on the
% deltaMedian signal.

DeltaMedianImagePositive = double(DeltaMedianImage(:,:,:,2));
DeltaMedianImageNegative = double(DeltaMedianImage(:,:,:,1));

clear DeltaMedianImage

%
nMasks = size(MaskDatabaseDownsamp, 2);
PosResults = cell(nMasks, 2);
NegResults = cell(nMasks, 2);

SignalInMasks = zeros(nMasks, 2); % here we will get the summed delta median value divided by the number of pixels within that ROI - AKA the mean signal within the ROI. First column is for positive (acitvation) signals. Second column is for negative (supressive) signals.

MaskImagePositive = zeros(height, width, Zs);
MaskImageNegative = zeros(height, width, Zs);

wait = waitbar(0, 'Calculating the signal in each mask');
% loop through and calculate signal in each mask. Also make an image which
% draws the mask, scaled to the ammount of signal within that mask. This is
% done on left vs right independently, as some signals are assymetric
% across the brain.

leftInds = false(height, width, Zs);
leftInds(:, 1:round(width/2), :) = 1;
rightInds = logical(true(height, width, Zs) - leftInds);

PosRight = zeros(height, width, Zs);
PosLeft = zeros(height, width, Zs);
NegRight = zeros(height, width, Zs);
NegLeft = zeros(height, width, Zs);


for i = 1:nMasks
    
    mask = reshape(full(MaskDatabaseDownsamp(:,i)), [height, width, Zs]);
    maskRight = mask;
    maskRight(leftInds) = 0;
    
    maskLeft = mask;
    maskLeft(rightInds) = 0;
    
    SignalInMasks(i, 1) = mean(DeltaMedianImagePositive(mask)); % get mean signal within logical 1s of the mask
    SignalInMasks(i, 2) = mean(DeltaMedianImageNegative(mask));
    
    PosRight(maskRight) = mean(DeltaMedianImagePositive(maskRight));
    PosLeft(maskLeft) = mean(DeltaMedianImagePositive(maskLeft));
    
    MaskImagePositive = max(cat(4, MaskImagePositive, PosRight, PosLeft), [], 4);
    
    NegRight(maskRight) = mean(DeltaMedianImageNegative(maskRight));
    NegLeft(maskLeft) = mean(DeltaMedianImageNegative(maskLeft));
    
    MaskImageNegative = max(cat(4, MaskImageNegative, NegRight, NegLeft), [], 4);
    waitbar(i/nMasks)
end
close(wait)



maxScaledMaskValue = max(max(MaskImagePositive(:)), max(MaskImageNegative(:)));
MaskImageNegative = MaskImageNegative./maxScaledMaskValue;
MaskImagePositive = MaskImagePositive./maxScaledMaskValue;

ScaledMaskStackName = strcat(strrep(DeltaMedianFilename, 'SignificantDeltaMedians.tif', ''), 'SignalInEachROIImage.tif');


for z = 1:Zs
    while true % avoid the annoying write error in windows, keep tryin to write until successful
        try
            if z == 1;
                imwrite(cat(3, squeeze(MaskImageNegative(:,:,z)), squeeze(MaskImagePositive(:,:,z)), squeeze(MaskImageNegative(:,:,z))), ScaledMaskStackName, 'writemode', 'overwrite')
            else
                imwrite(cat(3, squeeze(MaskImageNegative(:,:,z)), squeeze(MaskImagePositive(:,:,z)), squeeze(MaskImageNegative(:,:,z))), ScaledMaskStackName, 'writemode', 'append')
            end
            break
        catch
        end
    end
end
% take the maximum projection to make a 2D image
MaskImagePosZProject = max(MaskImagePositive,[], 3);
MaskImagePosXProject = squeeze(max(MaskImagePositive,[], 2));
MaskImagePosXProject = imresize(MaskImagePosXProject, [height, Zs.*(ZRez/XYRez)]);

MaskImagePosProject = double(cat(2, MaskImagePosZProject, MaskImagePosXProject));


MaskImageNegZProject = max(MaskImageNegative,[], 3);
MaskImageNegXProject = squeeze(max(MaskImageNegative,[], 2));
MaskImageNegXProject = imresize(MaskImageNegXProject, [height, Zs.*(ZRez/XYRez)]);

MaskImageNegProject = double(cat(2, MaskImageNegZProject, MaskImageNegXProject));

% normalize to the maximum value

maxProject = max(max(MaskImagePosProject(:)), max(MaskImageNegProject(:)));
MaskImagePosProject = MaskImagePosProject./maxProject;
MaskImageNegProject = MaskImageNegProject./maxProject;

ProjectionName = strcat(strrep(DeltaMedianFilename, 'SignificantDeltaMedians.tif', ''), 'SignalInEachROIProjection.tif');

imwrite(cat(3,MaskImageNegProject, MaskImagePosProject, MaskImageNegProject), ProjectionName, 'writemode', 'overwrite');

ROIFileStackName = strcat(strrep(DeltaMedianFilename, 'SignificantDeltaMedians.tif', ''), 'SignalInEachROI.mat');

save(ROIFileStackName, 'SignalInMasks', 'MaskDatabaseNames')

% rank the ROIs according to signal
if CalculateLabels
[SortedPositiveSignalInMasks, SortedPositiveSignalInMasksIX] = sort(SignalInMasks(:,1), 'descend');
[SortedNegativeSignalInMasks, SortedNegativeSignalInMasksIX] = sort(SignalInMasks(:,2), 'descend');

% get the names associated with each ROI

for i = 1:nMasks
    ind = SortedPositiveSignalInMasksIX(i);
    PosResults(i,1) = MaskDatabaseNames(ind);
end

for i = 1:nMasks
    ind = SortedNegativeSignalInMasksIX(i);
    NegResults(i,1) = MaskDatabaseNames(ind);
end

% in the second column we place the average signal per pixel in that region
for i = 1:nMasks
    PosResults(i,2) = {num2str(SortedPositiveSignalInMasks(i))};
    NegResults(i,2) = {num2str(SortedNegativeSignalInMasks(i))};
end
%
% loop through all the transgenes and rank them for how much signal they
% have within each ROI. look only at the mask pixels that are also hot in
% the deltaMedian image so that the signal only coming from areas within
% the ROI that show significant pERK signal. We then divide by the average
% value of the pixels in the 50 pixels surronding the ROI. This will give a
% value of 1 if the signal is the same inside and outside (or not-enriched
% in the activated area), greater than 1 if the signal is enriched, or less
% than one of the signal is greater ouside the region

% start by making new masks, where the region is only that covered by the
% map-map signal.
nMasksWPosSignal = find(SortedPositiveSignalInMasks == 0, 1, 'first')-1;
if isempty(nMasksWPosSignal) % if all masks have signal
    nMasksWPosSignal = nMasks;
end
PosMasks = false(height, width, Zs, nMasksWPosSignal);
Masks50OutsidePos = false(height, width, Zs, nMasksWPosSignal);


for n = 1:nMasksWPosSignal
    PosMasks(:,:,:,n) = reshape(full(MaskDatabaseDownsamp(:, SortedPositiveSignalInMasksIX(n))), [height, width, Zs]);
    PosMasks(:,:,:,n) = PosMasks(:,:,:,n).*(logical(DeltaMedianImagePositive)); % ,multiply to only look at mask that overlaps with activated region
    Masks50OutsidePos(:,:,:,n) = reshape(full(MaskDatabase50OutsideDownsamp(:, SortedPositiveSignalInMasksIX(n))), [height, width, Zs]);
end

nMasksWNegSignal = find(SortedNegativeSignalInMasks == 0, 1, 'first')-1;
if isempty(nMasksWNegSignal)
    nMasksWNegSignal = nMasks;
end
Masks50OutsideNeg = false(height, width, Zs, nMasksWNegSignal);


NegMasks = false(height, width, Zs, nMasksWNegSignal);
for n = 1:nMasksWNegSignal
    NegMasks(:,:,:,n) = reshape(full(MaskDatabaseDownsamp(:, SortedNegativeSignalInMasksIX(n))), [height, width, Zs]);
    NegMasks(:,:,:,n) = NegMasks(:,:,:,n).*(logical(DeltaMedianImageNegative)); % ,multiply to only look at mask that overlaps with activated region
    Masks50OutsideNeg(:,:,:,n) = reshape(full(MaskDatabase50OutsideDownsamp(:, SortedNegativeSignalInMasksIX(n))), [height, width, Zs]);
    
end

% now we loop through each transgene, and then find the average signal
% within the Mask divided by the average signal in the 50 Outside mask.

cd(MaskDatabasePath)
TransgeneDatabaseInfo = h5info('AnatomyLabelDatabaseDownsampled.hdf5');
maxMasks = max(nMasksWNegSignal, nMasksWPosSignal);
nTransgenes = length(TransgeneDatabaseInfo.Datasets);

PosTransgeneSignalWithinROI = zeros(nMasksWPosSignal, nTransgenes);
NegTransgeneSignalWithinROI = zeros(nMasksWNegSignal, nTransgenes);
TransgeneNames = struct();
wait = waitbar(0, 'sorting through the transgenes');
for i = 1:nTransgenes
    TransgeneNameFull = TransgeneDatabaseInfo.Datasets(i).Name;
    endChar = strfind(TransgeneNameFull, 'dpf') - 3; % make a list of the transgene names_ dont include the 'Xdpf_MeanOfXFish.tif' part of the string
     if isempty(endChar) % check to see if dpf is missing from file name, if so, just use whole name
        TransgeneNames(i).Name = TransgeneNameFull;
     else
         TransgeneNames(i).Name = TransgeneNameFull(1:endChar);
    end 
     
  
    Transgene = h5read('AnatomyLabelDatabaseDownsampled.hdf5', strcat('/', TransgeneNameFull));
    
    for n = 1:maxMasks;
        
        if n <= nMasksWPosSignal
            
            VoxelsWithinROI = Transgene(PosMasks(:,:,:,n)); % get the vales of the voxels within the ROI that overlap with the significant signal
            VoxelsOutsideROI = Transgene(Masks50OutsidePos(:,:,:,n)); % get the values of the voxels in a ring surrounding the ROI for normalization
            PosTransgeneSignalWithinROI(n, i) = mean(VoxelsWithinROI)./ mean(VoxelsOutsideROI); % calculate the enrichement ratio
            
        end
        
        if n <= nMasksWNegSignal
            
            VoxelsWithinROI = Transgene(NegMasks(:,:,:,n)); % get the vales of the voxels within the ROI that overlap with the significant signal
            VoxelsOutsideROI = Transgene(Masks50OutsideNeg(:,:,:,n)); % get the values of the voxels in a ring surrounding the ROI for normalization
            NegTransgeneSignalWithinROI(n, i) = mean(VoxelsWithinROI)./ mean(VoxelsOutsideROI); % calculate the enrichement ratio
        end
        
        
        
    end
    waitbar(i/nTransgenes)
end
close(wait)
%



% sort to rank the transgene signals within the ROI:

[SortedPositiveTransgeneValues, SortedPositiveTransgeneIX] = sort(PosTransgeneSignalWithinROI,2, 'descend');
[SortedNegativeTransgeneValues, SortedNegativeTransgeneIX] = sort(NegTransgeneSignalWithinROI,2, 'descend');


for j = 1:nMasksWPosSignal
    for m = 1:5
        PosResults(j,2+m*2-1) = {TransgeneNames(SortedPositiveTransgeneIX(j,m)).Name};
        PosResults(j,2+m*2) = {num2str(SortedPositiveTransgeneValues(j,m))};
    end
end

for j = 1:nMasksWNegSignal
    for m = 1:5
        NegResults(j,2+m*2-1) = {TransgeneNames(SortedNegativeTransgeneIX(j,m)).Name};
        NegResults(j,2+m*2) = {num2str(SortedNegativeTransgeneValues(j,m))};
    end
end


cd(DeltaMedianPath);
fid = fopen(strcat(strrep(DeltaMedianFilename, 'SignificantDeltaMedians.tif', ''), 'PositiveSignalResults.csv'),'w');
numColumns = size(PosResults,2);
numRows = size(PosResults,1);
fprintf(fid,'%s\n','ROI name, Signal in ROI, Top Label, Signal, 2nd Label, Signal, 3rd Label, Signal, 4th Label, Signal, 5th Label, Signal');


for j = 1:numRows
    if SortedPositiveSignalInMasks(j) > 0
        for i = 1:numColumns-1
            fprintf(fid,'%s,',PosResults{j,i});
        end
        fprintf(fid,'%s\n',PosResults{j,numColumns});
    end
end
fclose(fid);

fid = fopen(strcat(strrep(DeltaMedianFilename, 'SignificantDeltaMedians.tif', ''), 'NegativeSignalResults.csv'),'w');
fprintf(fid,'%s\n','ROI name, Signal in ROI, Top Label, Signal, 2nd Label, Signal, 3rd Label, Signal, 4th Label, Signal, 5th Label, Signal');



for j = 1:numRows
    if SortedNegativeSignalInMasks(j) > 0
        for i = 1:numColumns-1
            
            fprintf(fid,'%s,',NegResults{j,i});
        end
        fprintf(fid,'%s\n',NegResults{j,numColumns});
    end
end
fclose(fid);
end


beep on; beep
%
close all
clear all

end