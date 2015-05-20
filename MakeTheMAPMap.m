% This function will create the MAP-Map whole-brain activity map from
% registered pERK/tERK data by performing differential intensity analysis
% at each voxel.
%
% Morphed pERK data from two groups needs to be in two folders, and we are
% propted for the folders when the program starts. Each parwise comparison
% between all the groups is made. The input images need to be single
% channel 8-bit tiff stacks. Any registered stacks registered to any
% reference brain should work, but to keep consistent with the formatting
% of the Z-Brain and subsequent analysis, images should be registered into
% the reference brain ('Ref20131120pt14pl2.nrrd'), downsampled to
% 300x679x80 (x,y,z), and pre-smoothed with a 2D gaussian filter (sigma =
% 2). I use an Imagej macro script 'PrepareStacksForMAPMapping.ijm' to do
% this, which converts the '.nrrd' files output from CMTK into Tiffs,
% downsample and smooths appropriately.
%
% We look for which files are pERK stains, and which are total ERK stains
% based on the naming, where channel 1 has the string '1_warp' within it,
% while channel 2 is '2_warp', etc. This needs to be set appropriately
% before running the script.

% Pairwise comparisons are made between all of the input folders. the two
% groups analyzed (GroupA and GroupB) are combined and then ranked. A
% Mann-Whitney U statistic Z score is calculated at each pixel. To then set an FDR threshold,
% pseudogroups are assembled by randomization such that they have an equal
% number of fish from GroupA and GroupB, randomly assigned, and a Z score
% is calculated as before. This mixing procedure is repeated 500 times,
% creating a distribution of control Zs at each voxel. From this
% distribution we then set a FDR threshold 0.00005 (0.005% probability of
% finding a value greater than this threshold). This is done assuming a
% normal distribution. The GroupA vs GroupB Z score is then thresholded at
% this value of Z (all pixels with lesser value are set to 0. The image is
% normalized such that 10 = 255 in an 8bit tiff, and the image stacks are
% written. I also set a requirement that there is a minimal pERK staining
% intensity at the pixel to get rid of most non-brain pixels (average
% across all brains must be 5/255), and there also needs to be a mininmum
% of 5 warped fish in each group at the relevant pixel, set based on the
% fact that the only '0' pixels in the stack are those missing after
% warping - ie. not real image values.

% The main file output from this analysis is the
% '*SignificantDeltaMedian.tif' file, which we refer to as the MAP-Map .
% This stack depicts the delta Median value, at each pixel foud to be above
% the FDR threshold.
% Significant pixels are assigned a non-0 value. Green signal indicates
% higher signals in the first Group (for example, GroupX in a
% 'GroupX_over_GroupY' comparison), magneta signals indicate higher signals
% in the second group. These signals are mapped linearly between 0 and
% 65535 in a 16-bit stack, where 65535 represents a delta pERK level of
% 0.5.  We also write the thresholded Z-score comparison ('*ZScores.tif'),
% which represents is mapped 0-255 = ZScore 0-10 in an 8bit stack, As a
% control we write the final mixed group analysis from which the FDR
% threshold is calculated in order to give a sense of what false signals
% look like (there should be almost none). Finally we also write stacks
% depicting the median pERK level at each voxel, and the standard deviation
% at each voxel, for each group.
%
% Written by Owen Randlett (owen.randlett@gmail.com)
%
% This function uses two non-standard matlab functions, which I gratefully acknowldege the authors of:
% 
% uipickfiles.m by Douglas Schwarz 
%
% http://www.mathworks.com/matlabcentral/fileexchange/10867-uipickfiles--uigetfile-on-steroids
% http://www.mathworks.com/matlabcentral/fileexchange/view_license?file_info_id=10867
%
% parfor_progress.m by Jeremy Scheff 
% 
% http://www.mathworks.com/matlabcentral/fileexchange/32101-progress-monitor--progress-bar--that-works-with-parfor/content/parfor_progress.m
% http://www.mathworks.com/matlabcentral/fileexchange/view_license?file_info_id=32101

%

%-------Variables that can be tweaked
%
% ncores -- number of cores on your computer (i.e., num workers). over
% which to parallelize parfor loops. Make sure not to set too hight that
% you run out of RAM!
%
%
% ERK and pERK labels -- '*N_warp*' where N = post-warp channel number
%
% nPermutes -- the number of times we make pseudogroup comparisons for the
% FDR calculation
%
% FDRThresh -- The FRD thresholding of the Z-score. A value of 0.00005
% gives very clean control comparisons.
%
% UsingERK -- logical, set to 1 to run the total-ERK based normalization to
% calculate pERK levels (pERK/tERK). If you are not using this normalization (for
% example if you dont have the tERK stain, or want to quantify differences
% in another channel, set to 0).

%---------------------------------------------

function MakeTheMAPMap
close all
clear all
ncores = 4;
ERKLabel = '*1_warp*';
pERKLabel = '*2_warp*';
nPermutes = 500;
FDRThresh = 0.00005; %
UsingERK = 1;

disp(strcat('Running with the following settings: ncores = ', num2str(ncores), ' ERKLabel = ', num2str(ERKLabel), ' pERKLabel = ', num2str(pERKLabel), ' nPermutes = ', num2str(nPermutes), ' FDRThresh = ', num2str(FDRThresh), 'UsingERK = ', num2str(UsingERK)))


prompt = 'Select all the groups';

Groups.dirs = uipickfiles('Prompt', prompt);

NumGroups = length(Groups.dirs);

cd(Groups.dirs{1});
cd('..');

prompt = 'select output directory';

dirOut = cell2mat(uipickfiles('Prompt', prompt));

if length(dirOut) < 1
    error('Invalide dirOut')
end



ParObject = parpool; % open up the parallel workers



for i = 1:NumGroups;
    
    [~, name, ext] = fileparts(Groups.dirs{i}); % automatically get the group nape from the folder name
    
    Groups.names{i} = strcat(name, ext);
    
    disp(['Group ',num2str(i),' was selected as: ', Groups.names{i}]);
    
end

%setup a progress bar
NumComparisons = NumGroups.*(NumGroups-1)./2;

% now we need to define all the group comparisons
for m = 1:NumGroups
    for n = m+1:NumGroups % these loops define all the pariwise comparisons
        tic
        GrA = Groups.names{m};
        dirA = Groups.dirs{m};
        GrB = Groups.names{n};
        dirB = Groups.dirs{n};
        ChurnThroughGroups(GrA, dirA, GrB, dirB, dirOut, UsingERK, ERKLabel, pERKLabel, ncores, nPermutes, FDRThresh) % process each two group comparison
        toc
    end
end

delete(ParObject); %shut down the paralell pool
clear all

end


function ChurnThroughGroups(GrA, dirA, GrB, dirB, dirOut, UsingERK, ERKLabel, pERKLabel, ncores, nPermutes, FDRThresh)

cd(dirA)
if UsingERK == 1;
    fnames = dir(ERKLabel); %Get the file names: this relies on the CMTK program to add a "warp" after the channel numnber to the warped file images
else
    fnames = dir('*.tif*');
end
info = imfinfo(fnames(1).name); % get the image paramaters from the info file

height = info(1).Height; % assign the image dimesnsions
width = info(1).Width;
Zs = numel(info);

if info(1).BitsPerSample ~=8
    error('Input images not of correct bit depth. Should be 8-bit images') % if not using 8-bit images, the thresholds will be different
end


cd(dirA)

if UsingERK == 1;
    
    if length(dir(ERKLabel)) ~= length(dir(pERKLabel)) % check to make sure that there are the same number of pERK and ERK stacks
        warndlg('The number of stacks in group A arent equal');
    end
end

if UsingERK == 1;
    
    fnamesERKA = dir(ERKLabel);
    fnamespERKA = dir(pERKLabel);
    
else
    fnamespERKA = dir(pERKLabel);
    fnamesERKA = [];
end

nGroupA = length(fnamespERKA); % get the number of fish from the number of images found


cd(dirB) % do the same for groupB

if UsingERK == 1;
    if length(dir(ERKLabel)) ~= length(dir(pERKLabel)) % check to make sure that there are the same number of pERK and ERK stacks
        warndlg('The number of stacks in group B arent equal');
    end
end

if UsingERK == 1;
    fnamesERKB = dir(ERKLabel);
    fnamespERKB = dir(pERKLabel);
else
    fnamespERKB = dir(pERKLabel);
    fnamesERKB = [];
end

nGroupB = length(fnamespERKB);



% preallocate images


Amedian = zeros(height, width, Zs);
AStd = zeros(height, width, Zs);


Bmedian= zeros(height, width, Zs);
BStd = zeros(height, width, Zs);

mwuZScoresImage = zeros(height, width, Zs);
mwuMixedZScoresImageSingleIteration = zeros(height, width, Zs); % when I finally write the file, I dont want the one that is the mean that has been iterated over many times so I will use this single iteration to write as an example image of the FDR thresholding

%Set up the groups for the mixed group analysis

RandMixedGrA = zeros(floor(nGroupA./2)*2,nPermutes);
RandMixedGrB = zeros(floor(nGroupB./2)*2,nPermutes);

for nperm = 1:nPermutes % pre-assign all of the mixed groups. This is done so that the computation is only done once, and so that we have the same exact mixed group comparisons to write as the control image
    RandMixedGrA(:,nperm) = datasample(1:nGroupA, floor(nGroupA./2)*2, 'Replace', false); % pick random indicies for the psuedo group
    RandMixedGrB(:,nperm) = datasample(nGroupA+1:nGroupB+nGroupA, floor(nGroupB./2)*2, 'Replace', false);
end

cd(dirOut);
parfor_progress(Zs);

parfor (z = 1:Zs, ncores) % loop through each z plane one at a time to build up the final Z-score image stack
    
    BpERK = zeros(height, width, nGroupB);
    ApERK = zeros(height, width, nGroupA);
    
    if UsingERK == 1
        BERK = zeros(height, width, nGroupB);
        AERK = zeros(height, width, nGroupA);
    end
    
    cd(dirA);
    
    for j = 1:nGroupA % we loop through load all images from one group into an image stack
        
        infopERKA = imfinfo(fnamespERKA(j).name);
        
        % now we will calculate the pERK signal over the total signal. This
        % will be our matrix.
        
        
        ApERK(:,:,j) = imread(fnamespERKA(j).name, z, 'Info', infopERKA);
        if UsingERK == 1
            infoERKA = imfinfo(fnamesERKA(j).name);
            AERK(:,:,j) = imread(fnamesERKA(j).name, z, 'Info', infoERKA);
        end
    end
    
    
    if UsingERK == 1
        A = ApERK./AERK; %calculate pERK signal over total signal
    elseif UsingERK == 0
        A = ApERK;
    end
    
    
  
    
    cd (dirB);
    
    for j = 1:nGroupB % we loop through load all images from one group into an image stack
        
        infopERKB = imfinfo(fnamespERKB(j).name);
        
        % now we will calculate the pERK signal over the total signal. This
        % will be our matrix.
        
        
        BpERK(:,:,j) = imread(fnamespERKB(j).name, z, 'Info', infopERKB);
        if UsingERK == 1
            infoERKB = imfinfo(fnamesERKB(j).name);
            BERK(:,:,j) = imread(fnamesERKB(j).name, z, 'Info', infoERKB);
        end
    end
    
    
    
    if UsingERK == 1
        B = BpERK./BERK; %calculate pERK signal over total signal
    elseif UsingERK == 0
        B = BpERK;
    end
    
    
    Amedian(:,:,z) = nanmedian(A, 3);
    Bmedian(:,:,z) = nanmedian(B, 3);
    
    AStd(:,:,z) = nanstd(A, 0, 3);
    BStd(:,:,z) = nanstd(B, 0, 3);
    
    meanpERK = mean(cat(3, ApERK, BpERK), 3);
    
    
    % Now we want to calculate a ranksum like statistic for the groups
    % rather than caculating a T-statistic
    %
    
    %Combine all into one big matrix to rank.
    ABCombined = cat(3,A,B);
    %Determine which indices are from which group
    FishInGroupA = 1:size(A,3);
    FishInGroupB = size(A,3)+1:size(A,3)+size(B,3);
    
    ranks = tiedrank(permute(ABCombined, [3,2,1])); % rank the fish. The permute calls are required because tiedrank only works along the first dimension of a matrix
    ranks = permute(ranks, [3,2,1]);
    
    %Fisrt we do a kind of boostrapping we randomly mixed data to find
    %the false discovery distribution. By mixing the groups with equal
    %numbers of fish from each group, we should destroy any of the
    %group by group variablility that we are trying to quantify, and
    %thus determine the noise level in our analysis. This is done
    %'nPermutes' number of times
    
    ZScoresMixedGrA = zeros(height, width, nPermutes);
    ZScoresMixedGrB = zeros(height, width, nPermutes);
    
    for perm = 1:nPermutes
        
        % get the vector for the indiecies of the fish in each mixed group
        RandMixedGrAThisPerm = squeeze(RandMixedGrA(:,perm))';
        RandMixedGrBThisPerm = squeeze(RandMixedGrB(:,perm))';
        
        FishInMixedGrA = [RandMixedGrAThisPerm(1:end/2),RandMixedGrBThisPerm(1:end/2)];
        FishInMixedGrB = [RandMixedGrAThisPerm(end/2+1:end), RandMixedGrBThisPerm(end/2+1:end)];
        
        % sum the ranks from each group
        SumRanksMixedGrA = nansum(ranks(:,:,FishInMixedGrA),3);
        SumRanksMixedGrB = nansum(ranks(:,:,FishInMixedGrB),3);
        
        % determine the number of fish contributing to each pixel. This is
        % not always the same, because if for example part of the forebrain
        % was not imaged in one fish, after warping there are 0-pixel
        % values at this missing region. These are NaN'd out above, but we
        % need to keep track of the Ns for the ranksum and Z-score
        % statistics
        
        NsInMixedGrA = length(FishInMixedGrA) - sum(isnan(ranks(:,:,FishInMixedGrA)),3);
        NsInMixedGrB = length(FishInMixedGrB) - sum(isnan(ranks(:,:,FishInMixedGrB)),3);
        
        %Calculate the Mann-Whitney U-statistic
        %(http://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U)
        
        UMixedGrA = SumRanksMixedGrA - (NsInMixedGrA.*(NsInMixedGrA + 1))./2;
        UMixedGrB = SumRanksMixedGrB - (NsInMixedGrB.*(NsInMixedGrB + 1))./2;
        
        %Calculate the Z-score. The disribution of Us can be assumed to be
        %normal, and from inspection of my data this assumption appears
        %valid
        
        %The value of U is the mean between the two groups
        MUmixed = (UMixedGrA + UMixedGrB)./2;
        SigmaUmixed = sqrt(NsInMixedGrA.*NsInMixedGrB.*(NsInMixedGrA + NsInMixedGrB + 1)./12);
        ZScoresMixedGrATemp = (UMixedGrA - MUmixed)./SigmaUmixed;
        ZScoresMixedGrBTemp = (UMixedGrB - MUmixed)./SigmaUmixed;
        
        ZScoresMixedGrA(:,:,perm) = ZScoresMixedGrATemp;
        ZScoresMixedGrB(:,:,perm) = ZScoresMixedGrBTemp;
        
    end
    % now we set the FDR threshold. The distrubition of Zs appears
    % normally distributed at each pixel, so I will use norminv to find
    % the critical value for Z at each pixel, with an frd rate set
    % above (ex 0.00005).
    
    ZScoreThresh = norminv((1-FDRThresh), 0, std(ZScoresMixedGrA, 0, 3));
    
    
    % now we do the actual comparison between groups. Fist sum the
    % ranks;
    SumRanksGrA = nansum(ranks(:,:,FishInGroupA),3);
    SumRanksGrB = nansum(ranks(:,:,FishInGroupB),3);
    
    %Determine the Ns based on the number of non-NaN values at each
    %pixel
    NsInGrA = size(A,3) - sum(isnan(ranks(:,:,FishInGroupA)),3);
    NsInGrB = size(B,3) - sum(isnan(ranks(:,:,FishInGroupB)),3);
    
    % get the mann-whiteny Us
    UGrA = SumRanksGrA - (NsInGrA.*(NsInGrA + 1))./2;
    UGrB = SumRanksGrB - (NsInGrB.*(NsInGrB + 1))./2;
    %Calculate the Z-score
    MU = (UGrA + UGrB)./2;
    SigmaU = sqrt(NsInGrA.*NsInGrB.*(NsInGrA + NsInGrB + 1)./12);
    ZScoresGrA = (UGrA - MU)./SigmaU;
    
    %throw away any pixel where the mean pERK value is less than 5 (as
    %the brain is highly stained, these are mostly non-brain pixels)
    ZScoresGrA(meanpERK < 5) = NaN;
    %throw away any pixel in which we dont have at least 5 fish
    ZScoresGrA(NsInGrA < 5) = NaN;
    ZScoresGrA(NsInGrB < 5) = NaN;
    
    %Threshold the ZScores based on the FDR threshold we calculated.
    %Note this is done separately at each pixel.
    ZScoresGrA(abs(ZScoresGrA)<ZScoreThresh) = NaN;
    
    %Assign these values in the image stack we will wite later
    mwuZScoresImage(:,:,z) = ZScoresGrA;
    
    % Do the same thing for the last mixed-group comparison as an
    % example of the false discoveries.
    ZScoresMixedSigleIterationForWriting = ZScoresMixedGrA(:,:,end); % we will save the last iteration and write to file to show what the control FDR rate images look like.
    ZScoresMixedSigleIterationForWriting(meanpERK < 5) = NaN;
    ZScoresMixedSigleIterationForWriting(NsInGrA < 5) = NaN;
    ZScoresMixedSigleIterationForWriting(NsInGrB < 5) = NaN;
    ZScoresMixedSigleIterationForWriting(abs(ZScoresMixedSigleIterationForWriting)<ZScoreThresh) = NaN;
    mwuMixedZScoresImageSingleIteration(:,:,z) = ZScoresMixedSigleIterationForWriting;
    
    cd(dirOut);
    parfor_progress; %update the lovely parfor progess bar.
end

cd(dirOut);
parfor_progress(0);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(strcat('Completed the comparison of: ', GrA, '-vs-', GrB)) % uodate where we are in the analysis

%now we threshold the images such that a Z-statistic of 10 is saturated. We
%also split into the positive and negative components such that these can
%be written in differnt colours.
mwuZScoresImage = mwuZScoresImage./10;
mwuZScoresImage(mwuZScoresImage > 1) = 1;
mwuZScoresImage(mwuZScoresImage < -1) = -1;
mwuZScoresImagePositive = mwuZScoresImage;
mwuZScoresImagePositive(mwuZScoresImagePositive < 0) = 0;
mwuZScoresImageNegative = mwuZScoresImage.*-1;
mwuZScoresImageNegative(mwuZScoresImageNegative < 0) = 0;


mwuMixedZScoresImageSingleIteration = mwuMixedZScoresImageSingleIteration./10;

mwuMixedZScoresImageSingleIteration(mwuMixedZScoresImageSingleIteration > 1) = 1;
mwuMixedZScoresImageSingleIteration(mwuMixedZScoresImageSingleIteration < -1) = -1;
mwuMixedZScoresImageSingleIterationPositive = mwuMixedZScoresImageSingleIteration;
mwuMixedZScoresImageSingleIterationPositive(mwuMixedZScoresImageSingleIterationPositive < 0) = 0;
mwuMixedZScoresImageSingleIterationNegative = mwuMixedZScoresImageSingleIteration.*-1;
mwuMixedZScoresImageSingleIterationNegative(mwuMixedZScoresImageSingleIterationNegative < 0) = 0;

% Compute the delta Median images thresholded with the ZScores

DeltaMedian = (Amedian - Bmedian)./(Amedian + Bmedian);

DeltaMedianPos = zeros(size(DeltaMedian));
DeltaMedianNeg = zeros(size(DeltaMedian));

DeltaMedianPos(mwuZScoresImagePositive >0) = abs(DeltaMedian(mwuZScoresImagePositive >0));
DeltaMedianNeg(mwuZScoresImageNegative >0) = abs(DeltaMedian(mwuZScoresImageNegative >0));



DeltaMedianPos = uint16(DeltaMedianPos./0.5.*65535);
DeltaMedianNeg = uint16(DeltaMedianNeg./0.5.*65535);



cd(dirOut);

%Write the tiffs. The try/catch loop avoids the bug in windows when maltab
%tries to write to the tiff that windows explorer is accessing.

%Currently I write two files for the comparison of GrA over GrB (where what
%is brighter in GrA is green and brighter in GrB is Magenta), and also
%write GrB over GrA (Where GrB is green). This is redundant, but probably
%useful as I am not defining which is control vs experimental, and I like
%to have experimental as green.



if length(dir(strcat(GrA, '_over_', GrB, '_ZScores.tif'))) < 1 % check to make sure we havent already written this file
    tiffFile = strcat(GrA, '_over_', GrB, '_ZScores.tif');
    for ind = 1:Zs
        while true
            try
                imwrite(cat(3, mwuZScoresImageNegative(:,:,ind), mwuZScoresImagePositive(:,:,ind), mwuZScoresImageNegative(:,:,ind)), tiffFile,'WriteMode','append', 'Compression', 'none');
                break
            catch
                disp('Encountered a Tiff write error')
            end
        end
    end
else
    warning(strcat('Warning! Already written:  ', GrA, '_over_', GrB, '_ZScores.tif, Skipping file write'))
end

if length(dir(strcat(GrB, '_over_', GrA, '_ZScores.tif'))) < 1 % check to make sure we havent already written this file
    tiffFile = strcat(GrB, '_over_', GrA, '_ZScores.tif');
    for ind = 1:Zs
        while true
            try
                imwrite(cat(3, mwuZScoresImagePositive(:,:,ind), mwuZScoresImageNegative(:,:,ind), mwuZScoresImagePositive(:,:,ind)), tiffFile,'WriteMode','append', 'Compression', 'none');
                break
            catch
                disp('Encountered a Tiff write error')
            end
        end
    end
else
    warning(strcat('Warning! Already written:', GrB, '_over_', GrA, '_ZScores.tif, Skipping file write'))
end

% write the mixeed control images

if length(dir(strcat('MixedControl', GrA, '_over_', GrB, '_ZScores.tif'))) < 1 % check to make sure we havent already written this file
    tiffFile = strcat('MixedControl', GrA, '_over_', GrB, '_ZScores.tif');
    for ind = 1:Zs
        while true
            try
                imwrite(cat(3, mwuMixedZScoresImageSingleIterationNegative(:,:,ind), mwuMixedZScoresImageSingleIterationPositive(:,:,ind), mwuMixedZScoresImageSingleIterationNegative(:,:,ind)), tiffFile,'WriteMode','append', 'Compression', 'none');
                break
            catch
                disp('Encountered a Tiff write error')
            end
        end
    end
else
    warning(strcat('Warning! Already written:', 'MixedControl', GrA, '_over_', GrB, '_ZScores.tif, Skipping file write'))
end

if length(dir(strcat('MixedControl', GrB, '_over_', GrA, '_ZScores.tif'))) < 1 % check to make sure we havent already written this file
    tiffFile = strcat('MixedControl', GrB, '_over_', GrA, '_ZScores.tif');
    for ind = 1:Zs
        while true
            try
                imwrite(cat(3, mwuMixedZScoresImageSingleIterationPositive(:,:,ind), mwuMixedZScoresImageSingleIterationNegative(:,:,ind), mwuMixedZScoresImageSingleIterationPositive(:,:,ind)), tiffFile,'WriteMode','append', 'Compression', 'none');
                break
            catch
                disp('Encountered a Tiff write error')
            end
        end
    end
else
    warning(strcat('Warning! Already written:', 'MixedControl', GrB, '_over_', GrA, '_ZScores.tif, Skipping file write'))
end


% Write the median stack to file.


if UsingERK == 0; % if raw averages these should be 8 bit, if pERK/ERK+pERK analysis, they will already be from 0 to 1.
    Amedian = uint8(Amedian);
    Bmedian = uint8(Bmedian);
end

if length(dir(strcat(GrA, '_medianStack.tif'))) < 1 % check to make sure we havent already written this file
    
    tiffFile = strcat(GrA, '_medianStack.tif');
    for ind = 1:Zs
        while true
            try
                imwrite(Amedian(:,:,ind), tiffFile,'WriteMode','append', 'Compression', 'none');
                break
            catch
                disp('Encountered a Tiff write error')
            end
        end
        
    end
else
    %warning(strcat('Warning! Already written:',GrA, '_medianStack.tif,
    %Skipping file write'))
end


if length(dir(strcat(GrB, '_medianStack.tif')))<1
    
    tiffFile = strcat(GrB, '_medianStack.tif');
    for ind = 1:Zs
        while true
            try
                imwrite(Bmedian(:,:,ind), tiffFile,'WriteMode','append', 'Compression', 'none');
                break
            catch
                disp('Encountered a Tiff write error')
            end
        end
        
        
    end
else
    %warning(strcat('Warning! Already written:',GrB, '_medianStack.tif,
    %Skipping file write'))
    
    
end

% Write the std stack to file.


if UsingERK == 0; % if raw averages these should be 8 bit, if pERK/ERK+pERK analysis, they will already be from 0 to 1.
    AStd = uint8(AStd);
    BStd = uint8(BStd);
elseif UsingERK == 1;
    AStd = uint16(AStd.*65535);
    BStd = uint16(BStd.*65535);
end

if length(dir(strcat(GrA, '_StdStack.tif'))) < 1 % check to make sure we havent already written this file
    
    tiffFile = strcat(GrA, '_StdStack.tif');
    for ind = 1:Zs
        while true
            try
                imwrite(AStd(:,:,ind), tiffFile,'WriteMode','append', 'Compression', 'none');
                break
            catch
                disp('Encountered a Tiff write error')
            end
        end
        
    end
else
    %warning(strcat('Warning! Already written:',GrA, '_medianStack.tif,
    %Skipping file write'))
end


if length(dir(strcat(GrB, '_StdStack.tif')))<1
    
    tiffFile = strcat(GrB, '_StdStack.tif');
    for ind = 1:Zs
        while true
            try
                imwrite(BStd(:,:,ind), tiffFile,'WriteMode','append', 'Compression', 'none');
                break
            catch
                disp('Encountered a Tiff write error')
            end
        end
        
        
    end
else
    %warning(strcat('Warning! Already written:',GrB, '_medianStack.tif,
    %Skipping file write'))
    
    
end


TiffName = strcat(GrA, '_over_', GrB, '_SignificantDeltaMedians.tif');

if length(dir(TiffName)) < 1 % check to make sure we havent already written this file
    for ind = 1:Zs
        while true
            try
                imwrite(cat(3, DeltaMedianNeg(:,:,ind), DeltaMedianPos(:,:,ind), DeltaMedianNeg(:,:,ind)), TiffName,'WriteMode','append', 'Compression', 'none');
                break
            catch
                disp('Encountered a Tiff write error')
            end
        end
    end
else
    warning(strcat('Warning! Already written:', TiffName, 'Skipping file write'))
end

TiffName = strcat(GrB, '_over_', GrA, '_SignificantDeltaMedians.tif');

if length(dir(TiffName)) < 1 % check to make sure we havent already written this file
    for ind = 1:Zs
        while true
            try
                imwrite(cat(3, DeltaMedianPos(:,:,ind), DeltaMedianNeg(:,:,ind), DeltaMedianPos(:,:,ind)), TiffName,'WriteMode','append', 'Compression', 'none');
                break
            catch
                disp('Encountered a Tiff write error')
            end
        end
    end
else
    warning(strcat('Warning! Already written:', TiffName, 'Skipping file write'))
end





end





