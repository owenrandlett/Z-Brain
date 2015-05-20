

function ZBrainViewer(ViewerMode, InitialGreyStackNo, SectDim)

%Optional Input Arguments = ZebrafishBrainVeiwer7(ViewerMode, InitialGreyStackNo, SectDim)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function will allow the user to explore the larval zebrafish anatomy
% database. Instructions for interacting with the data are displayed when
% the program is started. Requires at least the 2013a version of Matlab.

% Confocal data (mean across fish) is contained within the file
% 'AnatomyLabelDatabase.hdf5'. The masks outlining the anatomical regions
% are contained within the file 'MaskDatabase.mat'. If the function
% 'ZBrainViewer.m' and the files are in Matlab's path, the images should
% open when the function is started. If not, you will need to point to the
% folder containing the files.
%
%
% Tested and written using Matlab R2014a running on OSx 10.9.3 by Owen
% Randlett (owen.randlett@gmail.com). October 17, 2014.
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input Arguments = ZBrainViewer(ViewerMode, InitialGreyStackNo, SectDim)
%
%%%%% ViewerMode = 0  --- (default). Anatomy database browsing mode only
%
% no additional image stacks are asked for
%
%
%
%%%%% ViewerMode = 1  --- Comparison to other Tiff stack
%
% This mode will overlay the anatomy database with another tiff stack. This
% can be used, for example, to compare the anatomy database to a new
% transgene or label, or to a MAP-Map. The label must be prepared as either
% a signle channel or three channel (RGB/composite) tiff stack (can be done
% by saving as a Tiff in FIJI/ImageJ). The stack must have been already
% registered into the reference brain, so as to be the correct resolution.
% Intensityies will be scaled to the max values in the satck.
%
%%%%%% InitialGreyStackNo - The number of the stack to load into the grey
%%%%%% channel at the beginning of the program
%
% InitialGreyStackNo = 0 -- load in the Elavl3-H2BRFP stack (default)
%
% InitialGreyStackNo >= 1 -- this stack number is loaded into the grey channel,
% ordered by the stack order in the hdf5 file. If InitialGreyStackNo >
% number of stacks, one is selected randomly
%
%%%%%% SectDim - the dimension/plane of sectioning. Note that ROI creation
%%%%%% and 'click to define regions' currently only work in the default
%%%%%% (coronal) view
%
% SectDim = 0 -- coronal (default)
% SectDim = 1 -- transverse
% SectDom = 2 -- saggital
%%%%%%
% Exporting data
%
% The 'AnatomyLabelDatabase.hdf5' file is in hdf5 format, and therefore
% should be directly readable by some imaging software packages (for
% example using the HDF5 plugin for FIJI). Or these can be read into Matlab
% using 'h5read' function, and then exported in whichever format you like.
% The mean-stacks in the HDF5 file are 3D image stacks, with an image size of
% 1406 x 621 x 138 (y, x, z), and voxel physical size of 0.798 um (x/y) and
% 2um (z).
%
% The 'MaskDatabse' variable in the MaskDatabase.mat file contains the
% regional mask definitions. These are stored as 2D sparse logical arrays,
% where each column is a regional mask. These can be reconstructed using
% reshape(full(MaskDatabase(x,:)), [height, width, Zs]), where height, width, Zs
% is the native the image size (1406, 621, 138). 'MaskDatabaseNames' lists
% the names of the regions. Once reconstructed these could similarly be
% exported as image stacks.
%
% 'MaskDatabaseDownsampled.mat' and 'AnatomyLabelDatabase.mat' are
% downsampled versions of the same files, at a resolution that matches the
% MAP-Maps: 679x300x80 (x, y, z). These are used by the
% 'ZBrainAnalysisOfMAPMaps.m' function.

%% Parse Inputs
if nargin == 0
    ViewerMode = 0; % if no value is passed, just load the anatomy database, do not ask for a MAP-Map
    SectDim = 0;
    InitialGreyStackNo = 0;
    disp('running "Default Viewermode", loading Elavl3:H2BRFP Stack')
elseif nargin == 1
    InitialGreyStackNo = 0; % load tERK stack if variable is not passed
    if ViewerMode == 0; % if no value is passed, just load the anatomy database, do not ask for a MAP-Map
        disp('running "Default Viewermode", loading Elavl3:H2BRFP Stack')
    elseif ViewerMode == 1
         disp('running "Comparison To new stack Viewermode"')       
    else
        warning('Invalid input!!! overriding and running "Default Viewermode", loading Elavl3:H2BRFP Stack')
    end
    SectDim = 0;
elseif nargin == 2
    
    if ViewerMode == 0;
        disp('running "Default Viewermode", loading specified stack')
    elseif ViewerMode == 1
        disp('running "MAP-Map Comparison Viewermode", loading specified Stack')
    elseif ViewerMode == 2
        disp('running "Comparison To Other Tiff Stack Viewermode"')
    else
        warning('Invalid input!!! overriding and running "Default Viewermode"')
    end
    SectDim = 0;
    
elseif nargin == 3
    
    if ViewerMode == 0;
        disp('running "Default Viewermode", loading specified stack')    
    elseif ViewerMode == 1
        disp('running "Comparison To Other Tiff Stack Viewermode"')
    else
        warning('Invalid input!!! overriding and running "Default Viewermode"')
    end
    
    if SectDim == 0
        disp('Coronal View')
    elseif SectDim == 1
        disp('Transverse View')
    elseif SectDim == 2
        disp('Saggital View')
    else
        warning('Invalid sectioning input, running Coronal View')
        SectDim = 0;
    end
    
elseif nargin > 3
    ViewerMode = 0;
    warning('Invalid input!!! overriding and running "Default Viewermode", loading tERK Stack')
end

% check version
Ver = version('-release');
if str2double(Ver(1:4)) < 2014
    beep on;
    beep
    warning('You are running an old version of Matlab - upgrade to avoid errors related to how MATLAB handles HDF5 files');
end

%% Intialize variables
global imFig
global sliderFig
global rect
global DisplaySlice
global PrintSlice
global GetMaskName
global ShiftingTime
global currentPos
global lastPos
global xSubtractor
global ySubtractor
global ClearMasks
global firstPos
global MakeNewRegionMask
global FinishNewMask
global SaveMask
global ClearMaskSlice
global WriteMask
global ChangeDim
global PrintStack;
global FirstSlice;
global upSlice;
global downSlice;
global ZoomIn;
global ZoomOut;

ChangeDim = 0;
rect = [0,0,0,0];
PrintSlice = 0;
GetMaskName = 0;
ShiftingTime = 0;
currentPos = 0;
lastPos = [0,0];
MakeNewRegionMask = 0;
FinishNewMask = 0;
SaveMask = 0;
ClearMaskSlice = 0;
MakingNewMask = 0;
ZRez = 2;  % the voxel size of the reference brain
XYRez = 0.789;
PrintStack = 0;
upSlice = 0;
downSlice = 0;
ZoomIn = 0;
ZoomOut = 0;


%

try
    StacksDir = which('AnatomyLabelDatabase.hdf5');
    StacksDir = strrep(StacksDir, 'AnatomyLabelDatabase.hdf5', '');
    cd(StacksDir)
    
    
catch
    
    StacksDir = uigetdir('*.hdf5', 'Select Directory With the Z-Brain Files (AnatomyLabelDatabase.hdf5 and MaskDatabase.mat)');
    cd(StacksDir);
    addpath(StacksDir);
    
end

% get the HDF5 file info
hdf5fileInfo = h5info('AnatomyLabelDatabase.hdf5');
filesLabels = hdf5fileInfo.Datasets;

ImageSize = filesLabels(1).Dataspace.Size;
height = ImageSize(1);
width = ImageSize(2);
Zs = ImageSize(3);


% load in the mask database '.mat' file, containing the anatomy ROIs
load('MaskDatabase.mat');


% set up the string for the dropdown menu to load in the stacks
for i = 1:length(filesLabels)+1
    if i == 1
        popupnamesLabels = '-';
    else
        FullName = filesLabels(i-1).Name;
        EndPos = strfind(filesLabels(i-1).Name, '_');
        if numel(EndPos) > 0
            TruncName = FullName(1:EndPos-1);
        else
            TruncName = FullName;
        end
        popupnamesLabels = strcat(popupnamesLabels, '|', TruncName);
        
        if strfind(filesLabels(i-1).Name, 'Elavl3-H2BRFP') == 1
            H2BStack = i;
        end
    end
end




popupnamesMasks = cell(1, length(MaskDatabaseNames) + 1);

popupnamesMasks(1, 2:end) = MaskDatabaseNames;

if ViewerMode == 1; % load in the new stack, either to RGB or grey channel. upsize appropriately
    
    [TiffName, TiffPath] = uigetfile('*.tif', 'Select the new stack');
    
    cd(TiffPath)
    info = imfinfo(TiffName);
    ZsTiff = length(info);
    heightTiff = info(1).Height;
    widthTiff = info(1).Width;
    nChan = length(info(1).BitsPerSample);
    
    if nChan ~= 3; % check if not RGB image
        
        
        greyStack = zeros(heightTiff, widthTiff, ZsTiff);
        
        wait = waitbar(0, 'Loading new Stack', 'Color', [1 1 1]);
        set(findobj(wait,'type','patch'), ...
            'edgecolor','k','facecolor','k')
        
        for Sect = 1:ZsTiff
            TempSlice = double(imread(TiffName, Sect, 'info', info));
            if numel(info(1).BitsPerSample) == 1 % if it is a 1 channel image
                greyStack(:,:,Sect) = TempSlice;
            elseif numel(info(1).BitsPerSample) > 1 % someone tries to load in a 3 channel image or something, sum the channels to make a greyscale image
                greyStack(:,:,Sect) = squeeze(sum(TempSlice, 3));
                if Sect == 1
                    h = warndlg('loading in a mutlichannel image, it will be converted to one channel');
                    warning('loading in a mutlichannel image, it will be converted to one channel')
                elseif Sect==ZsTiff
                    close(h)
                end
            end
            waitbar(Sect/ZsTiff)
        end
        close(wait)
        
        if ZsTiff ~=Zs || heightTiff ~= height || widthTiff ~= width
            wait = waitbar(0, 'RE-Sizing stack to fit Anatomy Database', 'Color', [1 1 1]);
            set(findobj(wait,'type','patch'), ...
                'edgecolor','k','facecolor','k')
            
            nLoops = ZsTiff + width;
            
            
            StackTemp = zeros(height, width, ZsTiff);
            
            
            for Sect = 1:ZsTiff
                StackTemp(:,:,Sect) = imresize(greyStack(:,:,Sect), [height, width], 'method', 'nearest');
                
                waitbar(Sect/nLoops)
            end
            
            greyStack = zeros(height, width, Zs);
            
            
            
            for w = 1:width
                greyStack(:,w,:) = imresize(squeeze(StackTemp(:, w, :)), [height, Zs], 'method', 'nearest');
                waitbar((ZsTiff + w)/nLoops)
            end
            
            
            greyStack = uint16(65535.*greyStack./max(greyStack(:)));
            close(wait);
        end
    end
    
    if nChan == 3
        
        rStackOriginal= zeros(heightTiff, widthTiff, ZsTiff);
        gStackOriginal= zeros(heightTiff, widthTiff, ZsTiff);
        bStackOriginal= zeros(heightTiff, widthTiff, ZsTiff);
        
        wait = waitbar(0, 'Loading MAP-Map', 'Color', [0 1 0]);
        set(findobj(wait,'type','patch'), ...
            'edgecolor','m','facecolor','m')
        
        for Sect = 1:ZsTiff
            TempSlice = imread(TiffName, Sect, 'info', info);
            rStackOriginal(:,:,Sect) = squeeze(TempSlice(:,:,1));
            gStackOriginal(:,:,Sect) = squeeze(TempSlice(:,:,2));
            bStackOriginal(:,:,Sect) = squeeze(TempSlice(:,:,3));
            waitbar(Sect/ZsTiff)
        end
        
        StackMax = max(max([rStackOriginal(:), gStackOriginal(:), bStackOriginal(:)]));
        
        rStackOriginal = 65535.*rStackOriginal./StackMax;
        gStackOriginal = 65535.*gStackOriginal./StackMax;
        bStackOriginal = 65535.*bStackOriginal./StackMax;
        
        close(wait)
        
        wait = waitbar(0, 'UP-Sizing new stack to fit Anatomy Database', 'Color', [0 1 0]);
        set(findobj(wait,'type','patch'), ...
            'edgecolor','m','facecolor','m')
        nLoops = ZsTiff + width;
        
        rStackTemp = zeros(height, width, ZsTiff);
        gStackTemp = zeros(height, width, ZsTiff);
        bStackTemp = zeros(height, width, ZsTiff);
        
        for Sect = 1:ZsTiff
            rStackTemp(:,:,Sect) = imresize(rStackOriginal(:,:,Sect), [height, width], 'method', 'nearest');
            gStackTemp(:,:,Sect) = imresize(gStackOriginal(:,:,Sect), [height, width], 'method', 'nearest');
            bStackTemp(:,:,Sect) = imresize(bStackOriginal(:,:,Sect), [height, width], 'method', 'nearest');
            waitbar(Sect/nLoops)
        end
        
        rStack = zeros(height, width, Zs);
        gStack = zeros(height, width, Zs);
        bStack = zeros(height, width, Zs);
        
        
        for w = 1:width
            rStack(:,w,:) = imresize(squeeze(rStackTemp(:, w, :)), [height, Zs], 'method', 'nearest');
            gStack(:,w,:) = imresize(squeeze(gStackTemp(:, w, :)), [height, Zs], 'method', 'nearest');
            bStack(:,w, :) = imresize(squeeze(bStackTemp(:, w, :)), [height, Zs], 'method', 'nearest');
            waitbar((ZsTiff + w)/nLoops)
        end
        
        rStack = uint16(rStack);
        gStack = uint16(gStack);
        bStack = uint16(bStack);
        

        
        close(wait)
        
        clear pERKSlice rStackOriginal gStackOriginal bStackOriginal rStackTemp gStackTemp bStackTemp
    end
end

if ~exist('InitialGreyStackNo', 'var') || InitialGreyStackNo == 0
    greyVal = H2BStack;
elseif InitialGreyStackNo > length(filesLabels) || InitialGreyStackNo < 0
    greyVal = randi(length(filesLabels) + 1);
else
    greyVal = InitialGreyStackNo + 1;
end

rVal = 1;
gVal = 1;
bVal = 1;
cVal = 1;
mVal = 1;
yVal = 1;

greyValMask = 1;
rValMask = 1;
gValMask = 1;
bValMask = 1;
cValMask = 1;
mValMask = 1;
yValMask = 1;

greyMaskStack = false(height, width, Zs);
rMaskStack = greyMaskStack;
gMaskStack = greyMaskStack;
bMaskStack = greyMaskStack;
cMaskStack = greyMaskStack;
mMaskStack = greyMaskStack;
yMaskStack = greyMaskStack;

%

%% Define the initial image paramaters

% initial image min/max scalings
greyMin = 1000;
greyMax = 25000;
rMin = 1000;
rMax = 25000;
gMin = 1000;
gMax = 25000;
bMin = 1000;
bMax = 25000;
cMin = 1000;
cMax = 25000;
mMin = 1000;
mMax = 25000;
yMin = 1000;
yMax = 25000;

% initial image opacities (the intenstiy values will be multiplied by this
% number (0 to 1) to reduce intensity or make image more transparent.
greyvalTransp = 1;
rvalTransp = 1;
gvalTransp = 1;
bvalTransp = 1;
cvalTransp = 1;
mvalTransp = 1;
yvalTransp = 1;

% initial image gammas
greyvalGamma = 1;
rvalGamma = 1;
gvalGamma = 1;
bvalGamma = 1;
cvalGamma = 1;
mvalGamma = 1;
yvalGamma = 1;

% opacity value of the mask outlines
greyvalMaskTransp = 0.7;
rvalMaskTransp = 0.7;
gvalMaskTransp = 0.7;
bvalMaskTransp = 0.7;
cvalMaskTransp = 0.7;
mvalMaskTransp = 0.7;
yvalMaskTransp = 0.7;


if SectDim == 0;
    ImHeight = height;
    ImWidth = width;
elseif SectDim == 1;
    ImHeight = width;
    ImWidth = Zs;
elseif SectDim == 2
    ImHeight = height;
    ImWidth = Zs;
end

greySlice = zeros(ImHeight, ImWidth);
rSlice = zeros(ImHeight, ImWidth);
gSlice = zeros(ImHeight, ImWidth);
bSlice = zeros(ImHeight, ImWidth);
cSlice = zeros(ImHeight, ImWidth);
mSlice = zeros(ImHeight, ImWidth);
ySlice = zeros(ImHeight, ImWidth);

greySliceMask = zeros(ImHeight, ImWidth);
rSliceMask = zeros(ImHeight, ImWidth);
gSliceMask = zeros(ImHeight, ImWidth);
bSliceMask = zeros(ImHeight, ImWidth);
cSliceMask = zeros(ImHeight, ImWidth);
mSliceMask = zeros(ImHeight, ImWidth);
ySliceMask = zeros(ImHeight, ImWidth);

nIter = 1;

%% Set up control figure with sliders

% this could be done in a more pretty way with loops, but for now just
% define each panel separately

sliderFig = figure;
set(sliderFig, 'menubar',...
    'none','Units', 'normalized',...
    'OuterPosition', [0, 0.05, 1,0.2],...
    'Color', [0 0 0], 'Numbertitle',...
    'off', 'Name', 'Control Window');

GammaText = uicontrol('Style','text',...
    'Units', 'normalized',...
    'Position',[0.61 0.93 0.1 0.065],...
    'String','Gamma Adjustment');

minText = uicontrol('Style','text',...
    'Units', 'normalized',...
    'Position',[0.2 0.93 0.2 0.065],...
    'String','Min Displayed Intensity');

maxText = uicontrol('Style','text',...
    'Units', 'normalized',...
    'Position',[0.405 0.93 0.2 0.065],...
    'String','Max Displayed Intensity');

TranspText = uicontrol('Style','text',...
    'Units', 'normalized',...
    'Position',[0.715 0.93 0.07 0.065],...
    'String','Opacity');

MaskText = uicontrol('Style','text',...
    'Units', 'normalized',...
    'Position',[0.788 0.93 0.211 0.065],...
    'String','Brain Region Masks');


MaskSeparator = uicontrol('Style','text',...
    'Units', 'normalized',...
    'Position',[0.787 0.2 0.002 1],...
    'BackgroundColor', [1 1 1]);
%
if SectDim == 0;
    TotalSlices = Zs;
elseif SectDim == 1;
    TotalSlices = height;
elseif SectDim == 2;
    TotalSlices = width;
end

Sect = round(TotalSlices./1.5);

sliderObj = uicontrol('Style', 'slider',...
    'Min',1,'Max',TotalSlices,'Value',Sect,...
    'Units', 'normalized',...
    'String', 'Z Slider',...
    'SliderStep', [1,5]/(Zs-1),...
    'BackgroundColor', [0.5, 0.5, 0.5],...
    'Position', [0, 0, 1, 0.22]);



if ViewerMode ~= 1 || nChan ~=1
    greyObj = uicontrol('Style', 'popup',...
        'String', popupnamesLabels,...
        'Units', 'normalized',...
        'Value',greyVal,...
        'Position', [0 0.82 0.15 0.1]);
elseif ViewerMode == 1 && nChan == 1;
    greyObj = uicontrol('Style', 'popup',...
        'String', TiffName,...
        'Units', 'normalized',...
        'Position', [0 0.62 0.15 0.1]);
end

greyText = uicontrol('Style','text',...
    'Units', 'normalized',...
    'Position',[0.15 0.85 0.046 0.07],...
    'String','Grey');

greyGammaSliderObj = uicontrol('Style', 'slider',...
    'Min',0,'Max',2,'Value',1,...
    'Units', 'normalized',...
    'SliderStep', [0.05, 0.2],...
    'BackgroundColor', [0.5, 0.5, 0.5],...
    'Position',[0.61 0.82 0.1 0.1]);

greyMinSliderObj = uicontrol('Style', 'slider',...
    'Min',0,'Max',65535,'Value',greyMin,...
    'Units', 'normalized',...
    'SliderStep', [250 5000]./65535,...
    'BackgroundColor', [0.5, 0.5, 0.5],...
    'Position',[0.2 0.82 0.2 0.1]);

greyMaxSliderObj = uicontrol('Style', 'slider',...
    'Min',0,'Max',65535,'Value',greyMax,...
    'Units', 'normalized',...
    'SliderStep', [250 5000]./65535,...
    'BackgroundColor', [0.5, 0.5, 0.5],...
    'Position',[0.405 0.82 0.2 0.1]);

greyTransparencyObj = uicontrol('Style', 'slider',...
    'Min', 0, 'Max', 1, 'Value', 1,...
    'SliderStep', [0.1, 1],...
    'Units', 'normalized',...
    'BackgroundColor', [0.5, 0.5, 0.5],...
    'Position',[0.715 0.82 0.07 0.1]);

greyMaskObj = uicontrol('Style', 'popup',...
    'String', popupnamesMasks,...
    'Units', 'normalized',...
    'Value',greyValMask,...
    'Position', [0.788 0.82 0.21 0.1]);


if ViewerMode ~= 1 || nChan == 1
    rObj = uicontrol('Style', 'popup',...
        'String', popupnamesLabels,...
        'Units', 'normalized',...
        'Value',rVal,...
        'Position', [0 0.72 0.15 0.1]);
    
elseif ViewerMode == 1 && nChan == 3;
    rObj = uicontrol('Style', 'popup',...
        'String', strcat(TiffName, '-RedCh'),...
        'Units', 'normalized',...
        'Value',rVal,...
        'Position', [0 0.72 0.15 0.1]);
end

rText = uicontrol('Style','text',...
    'Units', 'normalized',...
    'Position',[0.15 0.75 0.046 0.07],...
    'BackgroundColor', [1 0 0],...
    'String','Red');

rGammaSliderObj = uicontrol('Style', 'slider',...
    'Min',0,'Max',2,'Value',1,...
    'Units', 'normalized',...
    'SliderStep', [0.05, 0.2],...
    'BackgroundColor', [1 0 0],...
    'Position',[0.61 0.72 0.1 0.1]);


rMinSliderObj = uicontrol('Style', 'slider',...
    'Min',0,'Max',65535,'Value',rMin,...
    'Units', 'normalized',...
    'SliderStep', [250 5000]./65535,...
    'BackgroundColor', [1 0 0],...
    'Position',[0.2 0.72 0.2 0.1]);

rMaxSliderObj = uicontrol('Style', 'slider',...
    'Min',0,'Max',65535,'Value',rMax,...
    'Units', 'normalized',...
    'SliderStep', [250 5000]./65535,...
    'BackgroundColor', [1 0 0],...
    'Position',[0.405 0.72 0.2 0.1]);

rTransparencyObj = uicontrol('Style', 'slider',...
    'Min', 0, 'Max', 1, 'Value', 1,...
    'SliderStep', [0.1, 1],...
    'Units', 'normalized',...
    'BackgroundColor', [1 0 0],...
    'Position',[0.715 0.72 0.07 0.1]);

rMaskObj = uicontrol('Style', 'popup',...
    'String', popupnamesMasks,...
    'Units', 'normalized',...
    'Value',rValMask,...
    'Position', [0.788 0.72 0.21 0.1]);


if ViewerMode ~= 1 || nChan == 1
    gObj = uicontrol('Style', 'popup',...
        'String', popupnamesLabels,...
        'Units', 'normalized',...
        'Value',gVal,...
        'Position', [0 0.62 0.15 0.1]);
elseif ViewerMode == 1 && nChan == 3;
    gObj = uicontrol('Style', 'popup',...
        'String', strcat(TiffName, '-GrCh'),...
        'Units', 'normalized',...
        'Position', [0 0.62 0.15 0.1]);
end

gText = uicontrol('Style','text',...
    'Units', 'normalized',...
    'BackgroundColor', [0 1 0],...
    'Position',[0.15 0.65 0.046 0.07],...
    'String','Green');

gGammaSliderObj = uicontrol('Style', 'slider',...
    'Min',0,'Max',2,'Value',1,...
    'Units', 'normalized',...
    'SliderStep', [0.05, 0.2],...
    'BackgroundColor', [0 1 0],...
    'Position',[0.61 0.62 0.1 0.1]);

gMinSliderObj = uicontrol('Style', 'slider',...
    'Min',0,'Max',65535,'Value',gMin,...
    'Units', 'normalized',...
    'BackgroundColor', [0 1 0],...
    'SliderStep', [250 5000]./65535,...
    'Position',[0.2 0.62 0.2 0.1]);

gMaxSliderObj = uicontrol('Style', 'slider',...
    'Min',0,'Max',65535,'Value',gMax,...
    'Units', 'normalized',...
    'BackgroundColor', [0 1 0],...
    'SliderStep', [250 5000]./65535,...
    'Position',[0.405 0.62 0.2 0.1]);

gTransparencyObj = uicontrol('Style', 'slider',...
    'Min', 0, 'Max', 1, 'Value', 1,...
    'SliderStep', [0.1, 1],...
    'Units', 'normalized',...
    'BackgroundColor', [0 1 0],...
    'Position',[0.715 0.62 0.07 0.1]);

gMaskObj = uicontrol('Style', 'popup',...
    'String', popupnamesMasks,...
    'Units', 'normalized',...
    'Value',gValMask,...
    'Position', [0.788 0.62 0.21 0.1]);



if ViewerMode ~= 1 || nChan == 1
bObj = uicontrol('Style', 'popup',...
    'String', popupnamesLabels,...
    'Units', 'normalized',...
    'Value',bVal,...
    'Position', [0 0.52 0.15 0.1]);
elseif ViewerMode == 1 && nChan == 3;
    bObj = uicontrol('Style', 'popup',...
    'String', strcat(TiffName, '-BlCh'),...
    'Units', 'normalized',...
    'Value',bVal,...
    'Position', [0 0.52 0.15 0.1]);
end

bText = uicontrol('Style','text',...
    'Units', 'normalized',...
    'BackgroundColor', [0 0 1],...
    'Position',[0.15 0.55 0.046 0.07],...
    'String','Blue');

bGammaSliderObj = uicontrol('Style', 'slider',...
    'Min',0,'Max',2,'Value',1,...
    'Units', 'normalized',...
    'SliderStep', [0.05, 0.2],...
    'BackgroundColor', [0 0 1],...
    'Position',[0.61 0.52 0.1 0.1]);


bMinSliderObj = uicontrol('Style', 'slider',...
    'Min',0,'Max',65535,'Value',bMin,...
    'Units', 'normalized',...
    'BackgroundColor', [0 0 1],...
    'SliderStep', [250 5000]./65535,...
    'Position',[0.2 0.52 0.2 0.1]);

bMaxSliderObj = uicontrol('Style', 'slider',...
    'Min',0,'Max',65535,'Value',bMax,...
    'Units', 'normalized',...
    'BackgroundColor', [0 0 1],...
    'SliderStep', [250 5000]./65535,...
    'Position',[0.405 0.52 0.2 0.1]);

bTransparencyObj = uicontrol('Style', 'slider',...
    'Min', 0, 'Max', 1, 'Value', 1,...
    'SliderStep', [0.1, 1],...
    'Units', 'normalized',...
    'BackgroundColor', [0 0 1],...
    'Position',[0.715 0.52 0.07 0.1]);


bMaskObj = uicontrol('Style', 'popup',...
    'String', popupnamesMasks,...
    'Units', 'normalized',...
    'Value',bValMask,...
    'Position', [0.788 0.52 0.21 0.1]);



cObj = uicontrol('Style', 'popup',...
    'String', popupnamesLabels,...
    'Units', 'normalized',...
    'Value',cVal,...
    'Position', [0 0.42 0.15 0.1]);

cText = uicontrol('Style','text',...
    'Units', 'normalized',...
    'BackgroundColor', [0 1 1],...
    'Position',[0.15 0.45 0.046 0.07],...
    'String','Cyan');

cGammaSliderObj = uicontrol('Style', 'slider',...
    'Min',0,'Max',2,'Value',1,...
    'Units', 'normalized',...
    'SliderStep', [0.05, 0.2],...
    'BackgroundColor', [0 1 1],...
    'Position',[0.61 0.42 0.1 0.1]);

cMinSliderObj = uicontrol('Style', 'slider',...
    'Min',0,'Max',65535,'Value',cMin,...
    'Units', 'normalized',...
    'BackgroundColor', [0 1 1],...
    'SliderStep', [250 5000]./65535,...
    'Position',[0.2 0.42 0.2 0.1]);

cMaxSliderObj = uicontrol('Style', 'slider',...
    'Min',0,'Max',65535,'Value',cMax,...
    'Units', 'normalized',...
    'BackgroundColor', [0 1 1],...
    'SliderStep', [250 5000]./65535,...
    'Position',[0.405 0.42 0.2 0.1]);

cTransparencyObj = uicontrol('Style', 'slider',...
    'Min', 0, 'Max', 1, 'Value', 1,...
    'SliderStep', [0.1, 1],...
    'Units', 'normalized',...
    'BackgroundColor', [0 1 1],...
    'Position',[0.715 0.42 0.07 0.1]);

cMaskObj = uicontrol('Style', 'popup',...
    'String', popupnamesMasks,...
    'Units', 'normalized',...
    'Value',cValMask,...
    'Position', [0.788 0.42 0.21 0.1]);


    mObj = uicontrol('Style', 'popup',...
        'String', popupnamesLabels,...
        'Units', 'normalized',...
        'Position', [0 0.32 0.15 0.1]);


mText = uicontrol('Style','text',...
    'Units', 'normalized',...
    'BackgroundColor', [1 0 1],...
    'Position',[0.15 0.35 0.046 0.07],...
    'String','Magenta');

mGammaSliderObj = uicontrol('Style', 'slider',...
    'Min',0,'Max',2,'Value',1,...
    'Units', 'normalized',...
    'SliderStep', [0.05, 0.2],...
    'BackgroundColor', [1 0 1],...
    'Position',[0.61 0.32 0.1 0.1]);

mMinSliderObj = uicontrol('Style', 'slider',...
    'Min',0,'Max',65535,'Value',mMin,...
    'Units', 'normalized',...
    'BackgroundColor', [1 0 1],...
    'SliderStep', [250 5000]./65535,...
    'Position',[0.2 0.32 0.2 0.1]);

mMaxSliderObj = uicontrol('Style', 'slider',...
    'Min',0,'Max',65535,'Value',mMax,...
    'Units', 'normalized',...
    'BackgroundColor', [1 0 1],...
    'SliderStep', [250 5000]./65535,...
    'Position',[0.405 0.32 0.2 0.1]);


mTransparencyObj = uicontrol('Style', 'slider',...
    'Min', 0, 'Max', 1, 'Value', 1,...
    'SliderStep', [0.1, 1],...
    'Units', 'normalized',...
    'BackgroundColor', [1 0 1],...
    'Position',[0.715 0.32 0.07 0.1]);

mMaskObj = uicontrol('Style', 'popup',...
    'String', popupnamesMasks,...
    'Units', 'normalized',...
    'Value',mValMask,...
    'Position', [0.788 0.32 0.21 0.1]);

yObj = uicontrol('Style', 'popup',...
    'String', popupnamesLabels,...
    'Units', 'normalized',...
    'Value',yVal,...
    'Position', [0 0.22 0.15 0.1]);

yText = uicontrol('Style','text',...
    'Units', 'normalized',...
    'BackgroundColor', [1 1 0],...
    'Position',[0.15 0.25 0.046 0.07],...
    'String','Yellow');

yGammaSliderObj = uicontrol('Style', 'slider',...
    'Min',0,'Max',2,'Value',1,...
    'Units', 'normalized',...
    'SliderStep', [0.05, 0.2],...
    'BackgroundColor', [1 1 0],...
    'Position',[0.61 0.22 0.1 0.1]);

yMinSliderObj = uicontrol('Style', 'slider',...
    'Min',0,'Max',65535,'Value',yMin,...
    'Units', 'normalized',...
    'SliderStep', [250 5000]./65535,...
    'BackgroundColor', [1 1 0],...
    'Position',[0.2 0.22 0.2 0.1]);

yMaxSliderObj = uicontrol('Style', 'slider',...
    'Min',0,'Max',65535,'Value',yMax,...
    'Units', 'normalized',...
    'BackgroundColor', [1 1 0],...
    'SliderStep', [250 5000]./65535,...
    'Position',[0.405 0.22 0.2 0.1]);

yTransparencyObj = uicontrol('Style', 'slider',...
    'Min', 0, 'Max', 1, 'Value', 1,...
    'SliderStep', [0.1, 1],...
    'Units', 'normalized',...
    'BackgroundColor', [1 1 0],...
    'Position',[0.715 0.22 0.07 0.1]);

yMaskObj = uicontrol('Style', 'popup',...
    'String', popupnamesMasks,...
    'Units', 'normalized',...
    'Value',yValMask,...
    'Position', [0.788 0.22 0.21 0.1]);

%% Make the image display window
imFig = figure;

set(imFig, 'menubar', 'none',...
    'units', 'normalized',...
    'outerposition', [0 0.3 1 0.70],...
    'Color', [0 0 0],...
    'Numbertitle', 'off',...
    'Name', 'Image Window',...
    'doublebuffer','off',...
    'resize', 'off',...
    'KeyPressFcn', @IMFigCallback,...
    'WindowButtonDownFcn', @imFigButtonDown,...
    'WindowButtonUpFcn', @imFigButtonUp);

set(imFig, 'units', 'pixels');
iptsetpref('ImshowBorder', 'tight');
figSizeVec = get(imFig, 'position');
figWidth = floor(figSizeVec(3));
figHeight = floor(figSizeVec(4));
CurrentImage = imshow(zeros(figHeight, figWidth)); % for some reason this call shifts the figure, put it back up
set(imFig, 'units', 'normalized');
set(imFig, 'outerposition', [0 0.3 1 0.70])
set(CurrentImage,  'erasemode', 'none');
set(imFig, 'units', 'pixels');


axis manual

%% Continuous loop that looks for updates to sliders and dropdown menus, and updates the dispaly image accordingly


DisplaySlice = 1;
while true
    tic
    if ChangeDim == 1; % If 'v' is pressed, change sectioning dimension
        DisplaySlice = 1;
        
        if SectDim == 2;
            SectDim = 0;
        else
            SectDim = SectDim + 1;
        end
        
        if SectDim == 0;
            ImHeight = height;
            ImWidth = width;
        elseif SectDim == 1;
            ImHeight = width;
            ImWidth = Zs;
        elseif SectDim == 2
            ImHeight = height;
            ImWidth = Zs;
        end
        
        if SectDim == 0;
            TotalSlices = Zs;
        elseif SectDim == 1;
            TotalSlices = height;
        elseif SectDim == 2;
            TotalSlices = width;
        end
        figure(sliderFig)
        Sect = round(TotalSlices/3);
        sliderObj = uicontrol('Style', 'slider',...
            'Min',1,'Max',TotalSlices,'Value',Sect,...
            'Units', 'normalized',...
            'String', 'Z Slider',...
            'SliderStep', [1,5]/(Zs-1),...
            'BackgroundColor', [0.5, 0.5, 0.5],...
            'Position', [0, 0, 1, 0.22]);
        
        ChangeDim = 0;
    end
    if MakeNewRegionMask == 1   % if 'r' has been selected, enter the ROI creation mode
        if    SectDim == 0 % only allow this in the coronal imaging dimension for now
            if ~exist('newMaskStack', 'var')
                newMaskStack = false(height, width, Zs);
                NewMaskText = text(0.2,0.95, 'Draw the region in this slice using the "roipoly" tools. The program is frozen on this slice until you draw something', 'Units', 'Normalized', 'Color', [1,1,1], 'FontUnits', 'Normalized', 'FontSize', 0.03);
                NewMaskText2 = text(0.2,0.92, 'When done, move to the next slice you want to draw and press "r" again to draw', 'Units', 'Normalized', 'Color', [1,1,1], 'FontUnits', 'Normalized', 'FontSize', 0.03);
                NewMaskText3 = text(0.2,0.89, 'If you are unhappy with a slice you have drawn, press "x" to clear it', 'Units', 'Normalized', 'Color', [1,1,1], 'FontUnits', 'Normalized', 'FontSize', 0.03);
                NewMaskText4 = text(0.2,0.86, 'When you have drawn all of the relevant slices, press "i" to interpolate any missing slices, press "s" to complete the region and give it a name', 'Units', 'Normalized', 'Color', [1,1,1], 'FontUnits', 'Normalized', 'FontSize', 0.03);
                NewMaskText5 = text(0.2,0.83, 'When you have drawn all of the new regions you want, press "w" to write the updated database to disk', 'Units', 'Normalized', 'Color', [1,1,1], 'FontUnits', 'Normalized', 'FontSize', 0.03);
                
            end
            
            figure(imFig);
            %try
            [~, ~, ROIslice, ~, ~] = roipoly;
            
            if sum(rect) == 0
                if SlicePadHeight > 0
                    ROIslice(1:SlicePadHeight, :) = []; % remove the padded parts
                end
                if SlicePadWidth > 0
                    ROIslice(:, 1:SlicePadWidth) = [];
                end
            end
            
            if sum(rect) > 0; % add here what to do when
                %warning('Regions have to be drawn when zoomed out!');
                
                ROIslice(1:PaddingHeight, :) = [];
                ROIslice(end-PaddingHeight+1:end, :) = [];
                ROIslice(:, 1:PaddingWidth) = [];
                ROIslice(:, end-PaddingWidth+1:end) = [];
                
                ROIsliceRectFill = imresize(ROIslice, [ShiftCorrectedRect(4), ShiftCorrectedRect(3)], 'method', 'nearest');
                
                %ROIslice = false(SliceHeight, SliceWidth);
                ROIslice = false(width, height);
                ROIslice(ShiftCorrectedRect(2):ShiftCorrectedRect(2)+ShiftCorrectedRect(4)-1, ShiftCorrectedRect(1):ShiftCorrectedRect(1)+ShiftCorrectedRect(3)-1) = logical(ROIsliceRectFill);
                
                
            end
            
            
            newMaskStack(:,:,Sect) = newMaskStack(:,:,Sect) + imresize(rot90(ROIslice, 3), [height, width]);
            MakeNewRegionMask = 0;
            DisplaySlice = 1;
            %catch % catch any ROI input errors
            %    warning('ROI Input error, try again')
            %end
            
            
            try
                delete([NewMaskText, NewMaskText2, NewMaskText3, NewMaskText4, NewMaskText5, NewMaskText5]);
            catch
            end
        else
            beep
            warning('Only can draw regions when in coronal view')
            MakeNewRegionMask = 0;
        end
        
    end
    
    if FinishNewMask == 1
        if exist('newMaskStack', 'var') % check to make sure we actually started making a new mask, instead of pressing "i" by accident
            InterpText = text(0.2,0.90, 'InterpolatingSlices....', 'Units', 'Normalized', 'Color', [1,1,1], 'FontUnits', 'Normalized', 'FontSize', 0.07);
            
            DrawnZs = logical(squeeze(squeeze(sum(sum(newMaskStack))))); % Zs that have drwan masks in them
            DrawnInd = find(DrawnZs == 1);
            
            for nInterp = 1:length(DrawnInd) - 1;
                if DrawnInd(nInterp+1) - DrawnInd(nInterp) > 1
                    newMaskStack(:,:,DrawnInd(nInterp)+1:DrawnInd(nInterp+1)-1)= interp_shape_embedded(squeeze(newMaskStack(:,:,DrawnInd(nInterp+1))),squeeze(newMaskStack(:,:,DrawnInd(nInterp))), DrawnInd(nInterp+1) - DrawnInd(nInterp) - 1);
                end
            end
            delete(InterpText)
        end
        
        FinishNewMask = 0;
        DisplaySlice = 1;
        
    end
    
    if SaveMask == 1
        if exist('newMaskStack', 'var') % check to make sure we actually started making a new mask, instead of pressing "s" by accident
            
            DrawnZs = logical(squeeze(squeeze(sum(sum(newMaskStack))))); % Zs that have drwan masks in them
            DrawnInd = find(DrawnZs == 1);
            MissingInds = find(diff(DrawnInd) > 1);
            
            if numel(MissingInds) > 0
                FinishNewMask = 1;
                warning('You forgot to interpolate the missing slices... doing that now before creating the mask')
            else
                
                [newMaskRegionNames] = inputdlg({sprintf('What is the major region? :\n\nTelencephalon\nDiencephalon\nMesencephalon\nRhombencephalon\nSpinal Cord\nGanglia\n '), 'What is the descriptive name?'});
                newMaskName = [newMaskRegionNames{1}, ' - ', newMaskRegionNames{2}];
                newMaskSparse = sparse(reshape(newMaskStack, [height* width*Zs, 1]));
                
                MaskDatabaseNames(end+1) = {newMaskName};
                MaskDatabase(:, end+1) = newMaskSparse;
                MaskDatabaseOutlines(:, end+1) = sparse(reshape(bwdist(newMaskStack) ==1, [height* width*Zs, 1]));
                
                [ResortedNames, I] = sort(MaskDatabaseNames);
                MaskDatabaseNames = ResortedNames;
                MaskDatabase = MaskDatabase(:, I);
                MaskDatabaseOutlines = MaskDatabaseOutlines(:,I);
                
                cd(StacksDir)
                save(strcat('NewAnatomyRegion--', newMaskName, '.mat'), 'newMaskName', 'newMaskSparse') % write the new mask to file
                clear newMaskStack newMaskSparse
                
                WriteMask = 0;
                DisplaySlice = 1;
                NewMaskStackWritten = 1;
            end
            
            
        else
            SaveMask = 0;
        end
        
        
    end
    
    if WriteMask == 1 % command was called to save the new mask
        if exist('NewMaskStackWritten', 'var')
            msg = msgbox('Creating and saving the updated mask database... viewer will close when complete');
            movefile('MaskDatabase.mat',strcat('MaskDatabase_OldVersion_From_', DateCreated, '.mat'))
            DateCreated = datestr(now,'yymmdd_HHMMSS');
            save('MaskDatabase.mat', 'MaskDatabase', 'MaskDatabaseOutlines', 'MaskDatabaseNames', 'DateCreated', 'height', 'width', 'Zs')
            close(msg)
            close all
        end
        WriteMask = 0;
    end
    
    if ClearMaskSlice == 1
        newMaskStack(:,:,Sect) = false(height, width);
        ClearMaskSlice = 0;
        DisplaySlice = 1;
    end
    
    if ClearMasks == 1 % populate the masks upon location selection with mouse click
        
        % reset all of the mask objects
        figure(sliderFig)
        set(greyMaskObj, 'value', 1)
        set(rMaskObj, 'value', 1)
        set(gMaskObj, 'value', 1)
        set(bMaskObj, 'value', 1)
        set(cMaskObj, 'value', 1)
        set(mMaskObj, 'value', 1)
        set(yMaskObj, 'value', 1)
        
        
        
        
        ClearMasks = 0;
    end
    
    if GetMaskName == 1 && SectDim == 0
        
        try
            figure(imFig)
            ClickText = text(0.2,0.9, 'click on the region of interest', 'Units', 'Normalized', 'Color', [1,1,1], 'FontUnits', 'Normalized', 'FontSize', 0.05);
            
            
            [gInputX, gInputY] = ginput(1);
            delete(ClickText)
            
            
            if sum(rect) > 0
                
                maskX = (ShiftCorrectedRect(1) + ceil((gInputX - PaddingWidth)/CropMag));
                maskY = (ShiftCorrectedRect(2) + ceil((gInputY - PaddingHeight)/CropMag));
                
                maskIDSlice = false(width, height);
                maskIDSlice(maskY, maskX) = 1;
                
            else
                
                maskX = round((gInputX - SlicePadWidth));
                maskY = round(gInputY);
                
                maskIDSlice = false(SliceHeight, SliceWidth);
                
                maskIDSlice(maskY, maskX) = 1;
                
                maskIDSlice = imresize(maskIDSlice, [width, height], 'nearest');
                
                
            end
            
            maskIDSlice = rot90(maskIDSlice, 3);
            maskIDMatrix = false(height, width, Zs);
            
            maskIDMatrix(:,:,Sect) = maskIDSlice;
            maskIDInd = reshape(maskIDMatrix, [height*width*Zs, 1]);
            
            maskIDVec = full(MaskDatabase(maskIDInd, :));
            
            maskIDs = find(maskIDVec == 1);
            
            %FoundText = text(0.2,0.9, ['Found ', num2str(length(maskIDs)),
            %' anatomical masks in this region ... Loading them now...'],
            %'Units', 'Normalized', 'Backgroundcolor', [0,0,0], 'Color',
            %[1,1,1], 'FontUnits', 'Normalized', 'FontSize', 0.05);
            
            maskIDSize = sum(MaskDatabase(:,maskIDs), 1);
            
            [~, IX] = sort(maskIDSize);
            
            switch length(maskIDs)
                case 0
                    warning('Input error, please select again')
                    figure(imFig)
                    errorText = text(0.2,0.82, 'OOPS! Nothing found here, press "m" to try again....', 'Units', 'Normalized', 'BackgroundColor', [0,0,0], 'Color', [1,1,1], 'FontUnits', 'Normalized', 'FontSize', 0.05);
                    pause(1)
                    delete(errorText)
                case 1
                    set(greyMaskObj, 'value',  maskIDs(IX(1))+1)
                case 2
                    set(greyMaskObj, 'value', maskIDs(IX(1))+1)
                    set(rMaskObj, 'value', maskIDs(IX(2))+1)
                case 3
                    set(greyMaskObj, 'value', maskIDs(IX(1))+1)
                    set(rMaskObj, 'value', maskIDs(IX(2))+1)
                    set(gMaskObj, 'value', maskIDs(IX(3))+1)
                case 4
                    set(greyMaskObj, 'value', maskIDs(IX(1))+1)
                    set(rMaskObj, 'value', maskIDs(IX(2))+1)
                    set(gMaskObj, 'value', maskIDs(IX(3))+1)
                    set(bMaskObj, 'value', maskIDs(IX(4))+1)
                case 5
                    set(greyMaskObj, 'value', maskIDs(IX(1))+1)
                    set(rMaskObj, 'value', maskIDs(IX(2))+1)
                    set(gMaskObj, 'value', maskIDs(IX(3))+1)
                    set(bMaskObj, 'value', maskIDs(IX(4))+1)
                    set(cMaskObj, 'value', maskIDs(IX(5))+1)
                case 6
                    set(greyMaskObj, 'value', maskIDs(IX(1))+1)
                    set(rMaskObj, 'value', maskIDs(IX(2))+1)
                    set(gMaskObj, 'value', maskIDs(IX(3))+1)
                    set(bMaskObj, 'value', maskIDs(IX(4))+1)
                    set(cMaskObj, 'value', maskIDs(IX(5))+1)
                    set(mMaskObj, 'value', maskIDs(IX(6))+1)
                    
                otherwise
                    set(greyMaskObj, 'value', maskIDs(IX(1))+1)
                    set(rMaskObj, 'value', maskIDs(IX(2))+1)
                    set(gMaskObj, 'value', maskIDs(IX(3))+1)
                    set(bMaskObj, 'value', maskIDs(IX(4))+1)
                    set(cMaskObj, 'value', maskIDs(IX(5))+1)
                    set(mMaskObj, 'value', maskIDs(IX(6))+1)
                    set(yMaskObj, 'value', maskIDs(IX(7))+1)
                    
            end
            
            
            DisplaySlice = 1;
            GetMaskName = 0;
            %delete(FoundText)
        catch
            warning('Input error, please select again')
            figure(imFig)
            ErrorText = text(0.2,0.9, 'OOPS! Error, press "m" to try again...', 'Units', 'Normalized', 'BackgroundColor', [0,0,0], 'Color', [1,1,1], 'FontUnits', 'Normalized', 'FontSize', 0.05);
            pause(3)
            delete(ErrorText)
            GetMaskName = 0;
            
        end
        
    end
    
    if ShiftingTime == 1
        if numel(rect) == 4
            MoveX = round((currentPos(1) - firstPos(1))./CropMag) - xSubtractor;
            MoveY = round((currentPos(2) - firstPos(2))./CropMag) - ySubtractor;
            
            xSubtractor = xSubtractor + MoveX;
            ySubtractor = ySubtractor + MoveY;
            rect(1) = rect(1) - MoveX;
            rect(2) = rect(2) + MoveY;
            
            DisplaySlice = 1;
            
        end
    end
    
    if upSlice == 1; % update the sections if left/right arrow key pressed
        set(sliderObj, 'Value', Sect + 1)
        upSlice = 0;
    end
    if downSlice == 1; % update the sections if left/right arrow key pressed
        set(sliderObj, 'Value', Sect - 1)
        downSlice = 0;
    end
    
    if ZoomIn == 1; % Zoom in 2x
        ZoomF = 1.3;
        DisplaySlice = 1;
        if sum(rect) == 0;
            rect = [(figWidth - figWidth/ZoomF)/2,...
                (figHeight - figHeight/ ZoomF)/2,...
                figWidth/ ZoomF, figHeight/ ZoomF];
        end
        rect(1) = rect(1) + (rect(3) - rect(3)/ZoomF)/2;
        rect(2) = rect(2) + (rect(4) - rect(4)/ZoomF)/2;
        rect(3) = rect(3)/ ZoomF;
        rect(4) = rect(4)/ ZoomF;
        
        rect = round(rect);
        ZoomIn = 0;
    end
    
    if ZoomOut == 1; % Zoom in 2x
        ZoomF = 1.3;
        DisplaySlice = 1;
        if sum(rect) == 0;
            rect = [0,0,0,0];
        end
        rect(1) = rect(1) - (rect(3) - rect(3)/ZoomF)/2;
        rect(2) = rect(2) - (rect(4) - rect(4)/ZoomF)/2;
        rect(3) = rect(3)* ZoomF;
        rect(4) = rect(4)* ZoomF;
        
        rect = round(rect);
        ZoomOut = 0;
    end
    
    
    try % try to get new object values
        newSect = round(get(sliderObj, 'Value'));
        
        newvalgrey = get(greyObj, 'Value');
        newvalr = get(rObj, 'Value');
        newvalg = get(gObj, 'Value');
        newvalb = get(bObj, 'Value');
        newvalc = get(cObj, 'Value');
        newvalm = get(mObj, 'Value');
        newvaly = get(yObj, 'Value');
        
        newvalgreyMin = get(greyMinSliderObj, 'Value');
        newvalgreyMax = get(greyMaxSliderObj, 'Value');
        newvalrMin = get(rMinSliderObj, 'Value');
        newvalrMax = get(rMaxSliderObj, 'Value');
        newvalgMin = get(gMinSliderObj, 'Value');
        newvalgMax = get(gMaxSliderObj, 'Value');
        newvalbMin = get(bMinSliderObj, 'Value');
        newvalbMax = get(bMaxSliderObj, 'Value');
        newvalcMin = get(cMinSliderObj, 'Value');
        newvalcMax = get(cMaxSliderObj, 'Value');
        newvalmMin = get(mMinSliderObj, 'Value');
        newvalmMax = get(mMaxSliderObj, 'Value');
        newvalyMin = get(yMinSliderObj, 'Value');
        newvalyMax = get(yMaxSliderObj, 'Value');
        
        newvalgreyTransp = get(greyTransparencyObj, 'Value');
        newvalrTransp = get(rTransparencyObj, 'Value');
        newvalgTransp = get(gTransparencyObj, 'Value');
        newvalbTransp = get(bTransparencyObj, 'Value');
        newvalcTransp = get(cTransparencyObj, 'Value');
        newvalmTransp = get(mTransparencyObj, 'Value');
        newvalyTransp = get(yTransparencyObj, 'Value');
        
        newvalgreyGamma = get(greyGammaSliderObj, 'Value');
        newvalrGamma = get(rGammaSliderObj, 'Value');
        newvalgGamma = get(gGammaSliderObj, 'Value');
        newvalbGamma = get(bGammaSliderObj, 'Value');
        newvalcGamma = get(cGammaSliderObj, 'Value');
        newvalmGamma = get(mGammaSliderObj, 'Value');
        newvalyGamma = get(yGammaSliderObj, 'Value');
        
        newvalgreyMask = get(greyMaskObj, 'Value');
        newvalrMask = get(rMaskObj, 'Value');
        newvalgMask = get(gMaskObj, 'Value');
        newvalbMask = get(bMaskObj, 'Value');
        newvalcMask = get(cMaskObj, 'Value');
        newvalmMask = get(mMaskObj, 'Value');
        newvalyMask = get(yMaskObj, 'Value');
        
        
        
        
    catch
        close all
        break
    end
    
    
    % check to see if the stack index for any colour has changed, and if
    % so, load in the new stack.
    try
        if newvalgrey ~=greyVal || nIter == 1  % on first round load in the tERK stain
            if ViewerMode ~= 1 || nChan ~=1
                cd(StacksDir);
                greyVal = newvalgrey;
                if greyVal ~= 1;
                    wait = waitbar(0, 'Loading Stack', 'Color', [0.5 0.5 0.5]);
                    set(findobj(wait,'type','patch'), ...
                        'edgecolor','k','facecolor','k')
                    
                    waitbar(0.1)
                    greyStack = h5read('AnatomyLabelDatabase.hdf5', strcat('/', filesLabels(greyVal - 1).Name));
                    
                    waitbar(0.75)
                    
                    % check to make sure the stack is 16bit, and has an
                    % appropriate max value
                    maxStack = max(greyStack(:));
                    if maxStack < 65530 % if the stack isnt already maxed out to uint16 range, force it to do so
                        minStack = min(greyStack(greyStack > 0));
                        greyStack = uint16(65535.*(greyStack - minStack)./(maxStack - minStack));
                    end
                    close(wait)
                else
                    clear greyStack
                    greySlice = zeros(ImHeight, ImWidth);
                end
                DisplaySlice = 1;
            end
        end
        
        if newvalr ~=rVal
            cd(StacksDir);
            rVal = newvalr;
            if rVal ~= 1;
                wait = waitbar(0, 'Loading Stack', 'Color', [0.5 0.5 0.5]);
                set(findobj(wait,'type','patch'), ...
                    'edgecolor','r','facecolor','r')
                
                waitbar(0.1)
                rStack = h5read('AnatomyLabelDatabase.hdf5', strcat('/', filesLabels(rVal - 1).Name));
                
                waitbar(0.75)
                
                
                % check to make sure the stack is 16bit, and has an
                % appropriate max value
                maxStack = max(rStack(:));
                if maxStack < 65530 % if the stack isnt already maxed out to uint16 range, force it to do so
                    minStack = min(rStack(rStack > 0));
                    rStack = uint16(65535.*(rStack - minStack)./(maxStack - minStack));
                end
                close(wait)
            else
                clear rStack
                rSlice = zeros(ImHeight, ImWidth);
            end
            DisplaySlice = 1;
        end
        
        if newvalg ~=gVal
            cd(StacksDir);
            gVal = newvalg;
            if gVal ~= 1;
                wait = waitbar(0, 'Loading Stack', 'Color', [0.5 0.5 0.5]);
                set(findobj(wait,'type','patch'), ...
                    'edgecolor','g','facecolor','g')
                
                waitbar(0.1)
                gStack = h5read('AnatomyLabelDatabase.hdf5', strcat('/', filesLabels(gVal - 1).Name));
                
                waitbar(0.75)
                
                
                
                % check to make sure the stack is 16bit, and has an
                % appropriate max value
                maxStack = max(gStack(:));
                if maxStack < 65530 % if the stack isnt already maxed out to uint16 range, force it to do so
                    minStack = min(gStack(gStack > 0));
                    gStack = uint16(65535.*(gStack - minStack)./(maxStack - minStack));
                end
                close(wait)
            else
                clear gStack
                gSlice = zeros(ImHeight, ImWidth);
            end
            DisplaySlice = 1;
        end
        
        if newvalb ~=bVal
            cd(StacksDir);
            bVal = newvalb;
            if bVal ~= 1;
                wait = waitbar(0, 'Loading Stack', 'Color', [0.5 0.5 0.5]);
                set(findobj(wait,'type','patch'), ...
                    'edgecolor','b','facecolor','b')
                
                waitbar(0.1)
                bStack = h5read('AnatomyLabelDatabase.hdf5', strcat('/', filesLabels(bVal - 1).Name));
                
                waitbar(0.75)
                
                
                % check to make sure the stack is 16bit, and has an
                % appropriate max value
                maxStack = max(bStack(:));
                if maxStack < 65530 % if the stack isnt already maxed out to uint16 range, force it to do so
                    minStack = min(bStack(bStack > 0));
                    bStack = uint16(65535.*(bStack - minStack)./(maxStack - minStack));
                end
                close(wait)
            else
                clear bStack
                bSlice = zeros(ImHeight, ImWidth);
            end
            DisplaySlice = 1;
        end
        
        if newvalc ~=cVal
            cd(StacksDir);
            cVal = newvalc;
            if cVal ~= 1;
                wait = waitbar(0, 'Loading Stack', 'Color', [0.5 0.5 0.5]);
                set(findobj(wait,'type','patch'), ...
                    'edgecolor','c','facecolor','c')
                
                waitbar(0.1)
                cStack = h5read('AnatomyLabelDatabase.hdf5', strcat('/', filesLabels(cVal - 1).Name));
                
                waitbar(0.75)
                
                
                % check to make sure the stack is 16bit, and has an
                % appropriate max value
                maxStack = max(cStack(:));
                if maxStack < 65530 % if the stack isnt already maxed out to uint16 range, force it to do so
                    minStack = min(cStack(cStack > 0));
                    cStack = uint16(65535.*(cStack - minStack)./(maxStack - minStack));
                end
                close(wait)
            else
                clear cStack
                cSlice = zeros(ImHeight, ImWidth);
            end
            DisplaySlice = 1;
        end
        
        if newvalm ~=mVal
            cd(StacksDir);
            mVal = newvalm;
            if mVal ~= 1;
                wait = waitbar(0, 'Loading Stack', 'Color', [0.5 0.5 0.5]);
                set(findobj(wait,'type','patch'), ...
                    'edgecolor','m','facecolor','m')
                
                waitbar(0.1)
                mStack = h5read('AnatomyLabelDatabase.hdf5', strcat('/', filesLabels(mVal - 1).Name));
                
                waitbar(0.75)
                
                
                % check to make sure the stack is 16bit, and has an
                % appropriate max value
                maxStack = max(mStack(:));
                if maxStack < 65530 % if the stack isnt already maxed out to uint16 range, force it to do so
                    minStack = min(mStack(mStack > 0));
                    mStack = uint16(65535.*(mStack - minStack)./(maxStack - minStack));
                end
                close(wait)
            else
                clear mStack
                mSlice = zeros(ImHeight, ImWidth);
            end
            DisplaySlice = 1;
        end
        
        if newvaly ~=yVal
            cd(StacksDir);
            yVal = newvaly;
            if yVal ~= 1;
                wait = waitbar(0, 'Loading Stack', 'Color', [0.5 0.5 0.5]);
                set(findobj(wait,'type','patch'), ...
                    'edgecolor','y','facecolor','y')
                
                waitbar(0.1)
                yStack = h5read('AnatomyLabelDatabase.hdf5', strcat('/', filesLabels(yVal - 1).Name));
                
                waitbar(0.75)
                
                % check to make sure the stack is 16bit, and has an
                % appropriate max value
                maxStack = max(yStack(:));
                if maxStack < 65530 % if the stack isnt already maxed out to uint16 range, force it to do so
                    minStack = min(yStack(yStack > 0));
                    yStack = uint16(65535.*(yStack - minStack)./(maxStack - minStack));
                end
                close(wait)
            else
                clear yStack
                ySlice = zeros(ImHeight, ImWidth);
            end
            DisplaySlice = 1;
        end
        
    catch
    end
    
    % update Min, Max, Transp, Gamma values
    try
        if newvalgreyMin ~= greyMin ||...
                newvalrMin ~= rMin ||...
                newvalgMin ~= gMin ||...
                newvalbMin ~= bMin ||...
                newvalcMin ~= cMin ||...
                newvalmMin ~= mMin ||...
                newvalyMin ~= yMin
            
            greyMin = newvalgreyMin;
            rMin = newvalrMin;
            gMin = newvalgMin;
            bMin = newvalbMin;
            cMin = newvalcMin;
            mMin = newvalmMin;
            yMin = newvalyMin;
            DisplaySlice = 1;
            
        end
        
        if newvalgreyMax ~= greyMax ||...
                newvalrMax ~= rMax ||...
                newvalgMax ~= gMax ||...
                newvalbMax ~= bMax ||...
                newvalcMax ~= cMax ||...
                newvalmMax ~= mMax ||...
                newvalyMax ~= yMax
            
            greyMax = newvalgreyMax;
            rMax = newvalrMax;
            gMax = newvalgMax;
            bMax = newvalbMax;
            cMax = newvalcMax;
            mMax = newvalmMax;
            yMax = newvalyMax;
            
            
            DisplaySlice = 1;
            
        end
        
        if newvalgreyTransp ~= greyvalTransp ||...
                newvalrTransp ~= rvalTransp ||...
                newvalgTransp ~= gvalTransp ||...
                newvalbTransp ~= bvalTransp ||...
                newvalcTransp ~= cvalTransp ||...
                newvalmTransp ~= mvalTransp ||...
                newvalyTransp ~= yvalTransp
            
            greyvalTransp = newvalgreyTransp;
            rvalTransp = newvalrTransp;
            gvalTransp = newvalgTransp;
            bvalTransp = newvalbTransp;
            cvalTransp = newvalcTransp;
            mvalTransp = newvalmTransp;
            yvalTransp = newvalyTransp;
            DisplaySlice = 1;
            
        end
        
        if newvalgreyGamma ~= greyvalGamma ||...
                newvalrGamma ~= rvalGamma ||...
                newvalgGamma ~= gvalGamma ||...
                newvalbGamma ~= bvalGamma ||...
                newvalcGamma ~= cvalGamma ||...
                newvalmGamma ~= mvalGamma ||...
                newvalyGamma ~= yvalGamma
            
            greyvalGamma = newvalgreyGamma;
            rvalGamma = newvalrGamma;
            gvalGamma = newvalgGamma;
            bvalGamma = newvalbGamma;
            cvalGamma = newvalcGamma;
            mvalGamma = newvalmGamma;
            yvalGamma = newvalyGamma;
            DisplaySlice = 1;
        end
        
    catch
    end
    
    % update the regional masks if necessary
    try
        if newvalgreyMask ~=greyValMask % if a new mask is selected
            greyValMask = newvalgreyMask;
            greyMaskStack = false(height, width, Zs);
            
            if greyValMask ~= 1;
                greyMaskStack = reshape(full(MaskDatabaseOutlines(:,greyValMask-1)), [height, width, Zs]); % pull the correct mask 'outline' image from the MaskDatabaseOutlines file
                DisplaySlice = 1;
            else
                DisplaySlice = 1;
            end
        end
        
        if newvalrMask ~=rValMask % if a new mask is selected
            rValMask = newvalrMask;
            rMaskStack = false(height, width, Zs);
            
            if rValMask ~= 1;
                rMaskStack = reshape(full(MaskDatabaseOutlines(:,rValMask-1)), [height, width, Zs]); % pull the correct mask 'outline' image from the MaskDatabaseOutlines file
                DisplaySlice = 1;
            else
                DisplaySlice = 1;
            end
        end
        if newvalgMask ~=gValMask % if a new mask is selected
            gValMask = newvalgMask;
            gMaskStack = false(height, width, Zs);
            
            if gValMask ~= 1;
                gMaskStack = reshape(full(MaskDatabaseOutlines(:,gValMask-1)), [height, width, Zs]); % pull the correct mask 'outline' image from the MaskDatabaseOutlines file
                DisplaySlice = 1;
            else
                DisplaySlice = 1;
            end
        end
        
        if newvalbMask ~=bValMask % if a new mask is selected
            bValMask = newvalbMask;
            bMaskStack = false(height, width, Zs);
            
            if bValMask ~= 1;
                bMaskStack = reshape(full(MaskDatabaseOutlines(:,bValMask-1)), [height, width, Zs]); % pull the correct mask 'outline' image from the MaskDatabaseOutlines file
                DisplaySlice = 1;
            else
                DisplaySlice = 1;
            end
        end
        
        if newvalcMask ~=cValMask % if a new mask is selected
            cValMask = newvalcMask;
            cMaskStack = false(height, width, Zs);
            
            if cValMask ~= 1;
                cMaskStack = reshape(full(MaskDatabaseOutlines(:,cValMask-1)), [height, width, Zs]); % pull the correct mask 'outline' image from the MaskDatabaseOutlines file
                DisplaySlice = 1;
            else
                DisplaySlice = 1;
            end
        end
        
        if newvalmMask ~=mValMask % if a new mask is selected
            mValMask = newvalmMask;
            mMaskStack = false(height, width, Zs);
            
            if mValMask ~= 1;
                mMaskStack = reshape(full(MaskDatabaseOutlines(:,mValMask-1)), [height, width, Zs]); % pull the correct mask 'outline' image from the MaskDatabaseOutlines file
                DisplaySlice = 1;
            else
                DisplaySlice = 1;
            end
        end
        
        if newvalyMask ~=yValMask % if a new mask is selected
            yValMask = newvalyMask;
            yMaskStack = false(height, width, Zs);
            
            if yValMask ~= 1;
                yMaskStack = reshape(full(MaskDatabaseOutlines(:,yValMask-1)), [height, width, Zs]); % pull the correct mask 'outline' image from the MaskDatabaseOutlines file
                DisplaySlice = 1;
            else
                DisplaySlice = 1;
            end
        end
        
        
        if newSect ~=Sect
            Sect = newSect;
            DisplaySlice = 1;
            
            % once we move in Z for the first time, delete the instruction
            % box
            try delete([InstText1, InstText2, InstText3, InstText4, InstText5, InstText6, InstText7, InstText8, InstText9, InstTextLabels, InstTextAdjust, InstTextMasks]), catch, end
            
        end
    catch
    end
    
    % 
    if PrintStack == 1;
        
        OldSlice = Sect;
        if ~exist('StackExportDir', 'var')
        StackExportDir = uigetdir(pwd, 'select the directory to write stack to');
        oldDir = pwd;
        end
        if FirstSlice == 1 % start at first slice
            cd(StackExportDir);
            StackExportName = strcat('StackGrab_', datestr(now,'yymmdd_HHMMSS'), '.tif');
            WriteSliceNo = 1;
            Sect = WriteSliceNo;
            FirstSlice = 0;
            DisplaySlice = 1;
        elseif WriteSliceNo < TotalSlices % loop through all slices
            WriteSliceNo = WriteSliceNo+1;
            Sect = WriteSliceNo;
            DisplaySlice = 1;
        else  
            Sect = OldSlice;  % go back to where we started
            DisplaySlice = 1;
            PrintStack = 0;
            cd(oldDir);
        end
    end
    
    % compute and display the new slice when updated
    
        
    if DisplaySlice == 1 % Coronal View
        
        if SectDim == 0;
            try
                greySlice = greyvalTransp.*((double(greyStack(:,:,Sect)) - greyMin)./(greyMax - greyMin)).^greyvalGamma;
            catch
                greySlice = zeros(ImHeight, ImWidth);
            end
            
            try
                rSlice = rvalTransp.*((double(rStack(:,:,Sect)) - rMin)./(rMax - rMin)).^rvalGamma;
            catch
                rSlice = zeros(ImHeight, ImWidth);
            end
            
            try
                gSlice = gvalTransp.*((double(gStack(:,:,Sect)) - gMin)./(gMax - gMin)).^gvalGamma;
            catch
                gSlice = zeros(ImHeight, ImWidth);
            end
            
            try
                bSlice = bvalTransp.*((double(bStack(:,:,Sect)) - bMin)./(bMax - bMin)).^bvalGamma;
            catch
                bSlice = zeros(ImHeight, ImWidth);
            end
            
            try
                cSlice = cvalTransp.*((double(cStack(:,:,Sect)) - cMin)./(cMax - cMin)).^cvalGamma;
            catch
                cSlice = zeros(ImHeight, ImWidth);
            end
            
            try
                mSlice = mvalTransp.*((double(mStack(:,:,Sect)) - mMin)./(mMax - mMin)).^mvalGamma;
            catch
                mSlice = zeros(ImHeight, ImWidth);
            end
            
            try
                ySlice = yvalTransp.*((double(yStack(:,:,Sect)) - yMin)./(yMax - yMin)).^yvalGamma;
            catch
                ySlice = zeros(ImHeight, ImWidth);
            end
            
            greySliceMask = greyvalMaskTransp.* greyMaskStack(:,:,Sect);
            try
                greySliceMask = greySliceMask + bwdist(newMaskStack(:,:,Sect)) == 1;
            catch
            end
            
            rSliceMask = rvalMaskTransp.* rMaskStack(:,:,Sect);
            gSliceMask = gvalMaskTransp.* gMaskStack(:,:,Sect);
            bSliceMask = bvalMaskTransp.* bMaskStack(:,:,Sect);
            cSliceMask = cvalMaskTransp.* cMaskStack(:,:,Sect);
            mSliceMask = mvalMaskTransp.* mMaskStack(:,:,Sect);
            ySliceMask = yvalMaskTransp.* yMaskStack(:,:,Sect);
            
        elseif SectDim == 1 % Transverse View
            
            try
                greySlice = greyvalTransp.*((double(greyStack(Sect,:,:)) - greyMin)./(greyMax - greyMin)).^greyvalGamma;
            catch
                greySlice = zeros(ImHeight, ImWidth);
            end
            
            try
                rSlice = rvalTransp.*((double(rStack(Sect,:,:)) - rMin)./(rMax - rMin)).^rvalGamma;
            catch
                rSlice = zeros(ImHeight, ImWidth);
            end
            
            try
                gSlice = gvalTransp.*((double(gStack(Sect,:,:)) - gMin)./(gMax - gMin)).^gvalGamma;
            catch
                gSlice = zeros(ImHeight, ImWidth);
            end
            
            try
                bSlice = bvalTransp.*((double(bStack(Sect,:,:)) - bMin)./(bMax - bMin)).^bvalGamma;
            catch
                bSlice = zeros(ImHeight, ImWidth);
            end
            
            try
                cSlice = cvalTransp.*((double(cStack(Sect,:,:)) - cMin)./(cMax - cMin)).^cvalGamma;
            catch
                cSlice = zeros(ImHeight, ImWidth);
            end
            
            try
                mSlice = mvalTransp.*((double(mStack(Sect,:,:)) - mMin)./(mMax - mMin)).^mvalGamma;
            catch
                mSlice = zeros(ImHeight, ImWidth);
            end
            
            try
                ySlice = yvalTransp.*((double(yStack(Sect,:,:)) - yMin)./(yMax - yMin)).^yvalGamma;
            catch
                ySlice = zeros(ImHeight, ImWidth);
            end
            
            greySliceMask = greyvalMaskTransp.* squeeze(greyMaskStack(Sect,:,:));
            try
                greySliceMask = greySliceMask + bwdist(newMaskStack(Sect,:,:)) == 1;
            catch
            end
            
            rSliceMask = rvalMaskTransp.* squeeze(rMaskStack(Sect,:,:));
            gSliceMask = gvalMaskTransp.* squeeze(gMaskStack(Sect,:,:));
            bSliceMask = bvalMaskTransp.* squeeze(bMaskStack(Sect,:,:));
            cSliceMask = cvalMaskTransp.* squeeze(cMaskStack(Sect,:,:));
            mSliceMask = mvalMaskTransp.* squeeze(mMaskStack(Sect,:,:));
            ySliceMask = yvalMaskTransp.* squeeze(yMaskStack(Sect,:,:));
            
            greySlice = squeeze(greySlice);
            rSlice = squeeze(rSlice);
            gSlice = squeeze(gSlice);
            bSlice = squeeze(bSlice);
            cSlice = squeeze(cSlice);
            mSlice = squeeze(mSlice);
            ySlice = squeeze(ySlice);
            
            
        elseif SectDim == 2 % Saggital View
            
            try
                greySlice = greyvalTransp.*((double(greyStack(:,Sect,:)) - greyMin)./(greyMax - greyMin)).^greyvalGamma;
            catch
                greySlice = zeros(ImHeight, ImWidth);
            end
            
            try
                rSlice = rvalTransp.*((double(rStack(:,Sect,:)) - rMin)./(rMax - rMin)).^rvalGamma;
            catch
                rSlice = zeros(ImHeight, ImWidth);
            end
            
            try
                gSlice = gvalTransp.*((double(gStack(:,Sect,:)) - gMin)./(gMax - gMin)).^gvalGamma;
            catch
                gSlice = zeros(ImHeight, ImWidth);
            end
            
            try
                bSlice = bvalTransp.*((double(bStack(:,Sect,:)) - bMin)./(bMax - bMin)).^bvalGamma;
            catch
                bSlice = zeros(ImHeight, ImWidth);
            end
            
            try
                cSlice = cvalTransp.*((double(cStack(:,Sect,:)) - cMin)./(cMax - cMin)).^cvalGamma;
            catch
                cSlice = zeros(ImHeight, ImWidth);
            end
            
            try
                mSlice = mvalTransp.*((double(mStack(:,Sect,:)) - mMin)./(mMax - mMin)).^mvalGamma;
            catch
                mSlice = zeros(ImHeight, ImWidth);
            end
            
            try
                ySlice = yvalTransp.*((double(yStack(:,Sect,:)) - yMin)./(yMax - yMin)).^yvalGamma;
            catch
                ySlice = zeros(ImHeight, ImWidth);
            end
            
            greySliceMask = greyvalMaskTransp.* squeeze(greyMaskStack(:,Sect,:));
            try
                greySliceMask = greySliceMask + bwdist(newMaskStack(:,Sect,:)) == 1;
            catch
            end
            
            rSliceMask = rvalMaskTransp.* squeeze(rMaskStack(:,Sect,:));
            gSliceMask = gvalMaskTransp.* squeeze(gMaskStack(:,Sect,:));
            bSliceMask = bvalMaskTransp.* squeeze(bMaskStack(:,Sect,:));
            cSliceMask = cvalMaskTransp.* squeeze(cMaskStack(:,Sect,:));
            mSliceMask = mvalMaskTransp.* squeeze(mMaskStack(:,Sect,:));
            ySliceMask = yvalMaskTransp.* squeeze(yMaskStack(:,Sect,:));
            
            greySlice = squeeze(greySlice);
            rSlice = squeeze(rSlice);
            gSlice = squeeze(gSlice);
            bSlice = squeeze(bSlice);
            cSlice = squeeze(cSlice);
            mSlice = squeeze(mSlice);
            ySlice = squeeze(ySlice);
            
            
        end
        
        
        
        greySlice = greySlice + greySliceMask;
        rSlice = rSlice + rSliceMask;
        gSlice = gSlice + gSliceMask;
        bSlice = bSlice + bSliceMask;
        cSlice = cSlice + cSliceMask;
        mSlice = mSlice + mSliceMask;
        ySlice = ySlice + ySliceMask;
        
        greySlice(greySlice < 0) = 0;
        rSlice(rSlice < 0) = 0;
        gSlice(gSlice < 0) = 0;
        bSlice(bSlice < 0) = 0;
        cSlice(cSlice < 0) = 0;
        mSlice(mSlice < 0) = 0;
        ySlice(ySlice < 0) = 0;
        
        
        Slice = cat(3, rSlice + mSlice + ySlice + greySlice, gSlice + cSlice + ySlice + greySlice, bSlice + mSlice + cSlice + greySlice);
        Slice = uint8(real(Slice).*255);
        

        
        if SectDim == 1 || SectDim == 2
            Slice = imresize(Slice, [ImHeight, ImWidth.*(ZRez/XYRez)]);
        end
        
        if PrintSlice == 1;
            if ~exist('dirOut', 'var')
                dirOut = uigetdir(pwd, 'select directory to ouput slice to');
                oldDir = pwd;
            end
            cd(dirOut)
            try
                imwrite(rot90(Slice), strcat('AnatomySliceGrab_', datestr(now,'yymmdd_HHMMSS'), '.tif'))
            catch % catch rot90 error in old versions
                imwrite(Slice, strcat('AnatomySliceGrab_', datestr(now,'yymmdd_HHMMSS'), '.tif'))
            end
            cd(oldDir);
            PrintSlice = 0;
        end
        
        if PrintStack == 1;
            try
            imwrite(rot90(Slice), StackExportName, 'writemode', 'append', 'compression', 'none');    
            catch % catch rot90 error in old versions
             imwrite(Slice, StackExportName, 'writemode', 'append', 'compression', 'none');    
            end
        end
            

        try
            Slice = rot90(Slice);
        catch % for versions before 2014a that cant rot90 in 3D
            tmp1 = rot90(squeeze(Slice(:,:,1)));
            tmp2 = rot90(squeeze(Slice(:,:,2)));
            tmp3 = rot90(squeeze(Slice(:,:,3)));
            Slice = cat(3, tmp1, tmp2, tmp3);
        end
        
        if SectDim ~=0 && sum(rect) ~=0;
            rect = [0,0,0,0];
            warning('reset zoom - can only zoom in coronal view')
        end
        
        figure(imFig);
        
        if sum(rect) == 0;
            % figure out whether to resize the image by setting either the
            % height of the width to the size of the figure
            
            RatioHeight = size(Slice, 1)/figHeight;
            RatioWidth = size(Slice, 2)/figWidth;
            
            if RatioHeight >= RatioWidth
                CurrenImageSizeRatio = RatioHeight;
                Slice = imresize(Slice, [figHeight, NaN], 'method', 'bilinear');
                [SliceHeight, SliceWidth, ~] = size(Slice);
                SlicePadHeight = 0;
                SlicePadWidth = figWidth - SliceWidth;
                Slice = padarray(Slice, [SlicePadHeight, SlicePadWidth], 0, 'pre');
                
                
            else
                CurrenImageSizeRatio = RatioWidth;
                Slice = imresize(Slice, [NaN, figWidth], 'method', 'bilinear');
                [SliceHeight, SliceWidth, ~] = size(Slice);
                SlicePadHeight = figHeight - SliceHeight;
                SlicePadWidth = 0;
                Slice = padarray(Slice, [SlicePadHeight, SlicePadWidth], 0, 'post');
            end
            % update the image display
           
            set(CurrentImage,'CData',Slice);
            
            
        elseif sum(rect) > 0; % now we need to crop the image, crop it, blow it up using imresize, then use padarray to center it in the figure
            try
                ShiftCorrectedRect = rect;
                ShiftCorrectedRect(1) = ShiftCorrectedRect(1)-SlicePadWidth;
                ShiftCorrectedRect(2) = ShiftCorrectedRect(2)-SlicePadHeight;
                
                if ShiftCorrectedRect(1) < 0; % deal with cropping at edges of image, and cropping ouside the real image
                    ShiftCorrectedRect(3) = ShiftCorrectedRect(3) - ShiftCorrectedRect(1);
                    ShiftCorrectedRect(1) = 0;
                end
                
                if ShiftCorrectedRect(2) < 0; % deal with cropping at edges of image, and cropping ouside the real image
                    ShiftCorrectedRect(4) = ShiftCorrectedRect(4) - ShiftCorrectedRect(2);
                    ShiftCorrectedRect(2) = 0;
                end
                
                % resize the ShiftCorrectedRect such that it is in the
                % original dimensions of Slice
                
                RatioHeight = width/figHeight;
                RatioWidth = height/figWidth;
                
                if RatioHeight >= RatioWidth
                    ShiftCorrectedRect = round(ShiftCorrectedRect).*RatioHeight;
                else
                    ShiftCorrectedRect = round(ShiftCorrectedRect).*RatioHeight;
                end
                
                ShiftCorrectedRect = ceil(ShiftCorrectedRect);
                
                SliceCrop = imcrop(Slice, ShiftCorrectedRect);
                
                %re-size the cropped slice to fill the image window
                
                if ShiftCorrectedRect(3)/figWidth > ShiftCorrectedRect(4)/figHeight
                    SliceCrop = imresize(SliceCrop, [NaN, figWidth]);
                    CropMag = figWidth/rect(3);
                else
                    SliceCrop = imresize(SliceCrop, [figHeight, NaN]);
                    CropMag = figHeight/rect(4);
                end
                
                PaddingHeight = floor((figHeight - size(SliceCrop, 1))/2) + 1;
                PaddingWidth = floor((figWidth - size(SliceCrop, 2))/2) + 1;
                Slice = padarray(SliceCrop, [PaddingHeight, PaddingWidth], 'both');
                set(CurrentImage,'CData',Slice);
                
            catch
                warning('Zoom Error! Resetting Zoom')
                rect = 0;
            end
        end
        
        

        DisplaySlice = 0;

    end
    
    
    
    
    
    
    if nIter == 1  % display instructions when first started
        figure(imFig)
        InstText1 = text(0.01,0.95,'to zoom use up/down arrows, or to zoom to a specific region press "z", define zoom box then double click', 'Units', 'Normalized', 'BackgroundColor', [0.2,0.2,0.2], 'Color', [1,1,1], 'FontUnits', 'Normalized', 'FontSize', 0.02);
        InstText2 = text(0.01,0.92, 'press "z" again to zoom back out',  'Units', 'Normalized', 'BackgroundColor', [0.2,0.2,0.2], 'Color', [1,1,1], 'FontUnits', 'Normalized', 'FontSize', 0.02);
        InstText3 = text(0.01,0.89, 'click and drag to pan while zoomed in',  'Units', 'Normalized', 'BackgroundColor', [0.2,0.2,0.2], 'Color', [1,1,1], 'FontUnits', 'Normalized', 'FontSize', 0.02);
        InstText4 = text(0.01,0.83,'"m" to enter click-to-find anatomical regions mode. Click and area and the regions that overlap here will be displayed',  'Units', 'Normalized', 'BackgroundColor', [0.2,0.2,0.2], 'Color', [1,1,1], 'FontUnits', 'Normalized', 'FontSize', 0.02);
        InstText5 = text(0.01,0.80,'"c" to clear all loaded anatomical regions away',  'Units', 'Normalized', 'BackgroundColor', [0.2,0.2,0.2], 'Color', [1,1,1], 'FontUnits', 'Normalized', 'FontSize', 0.02);
        InstText6 = text(0.01,0.74,'"p" to print the current slice to file, "o" to print the entire stack to file',  'Units', 'Normalized', 'BackgroundColor', [0.2,0.2,0.2], 'Color', [1,1,1], 'FontUnits', 'Normalized', 'FontSize', 0.02);
        InstText7 = text(0.01,0.68,'"r" to draw and define a new region/mask',  'Units', 'Normalized', 'BackgroundColor', [0.2,0.2,0.2], 'Color', [1,1,1], 'FontUnits', 'Normalized', 'FontSize', 0.02);
        InstText8 = text(0.01,0.58,'"v" to change view-mode (Coronal -> Transverse -> Saggital). WARNING - "z", "m", "r" functions may only work in default (Coronal) view mode',  'Units', 'Normalized', 'BackgroundColor', [0.2,0.2,0.2], 'Color', [1,1,1], 'FontUnits', 'Normalized', 'FontSize', 0.02);
        InstText9 = text(0.01,0.48,'Left/Right arrows move one plane step',  'Units', 'Normalized', 'BackgroundColor', [0.2,0.2,0.2], 'Color', [1,1,1], 'FontUnits', 'Normalized', 'FontSize', 0.02);
        InstTextLabels = text(0.01,0.03,'Pick from the menus below to load in labels',  'Units', 'Normalized', 'BackgroundColor', [0.2,0.2,0.2], 'Color', [1,1,1], 'FontUnits', 'Normalized', 'FontSize', 0.02);
        InstTextAdjust = text(0.3,0.03,'Adjust the image paramaters with the sliders, the slider along the bottom controls the Z-plane',  'Units', 'Normalized', 'BackgroundColor', [0.2,0.2,0.2], 'Color', [1,1,1], 'FontUnits', 'Normalized', 'FontSize', 0.02);
        InstTextMasks = text(0.8,0.03,'Pick below to view an annotated anatomical region',  'Units', 'Normalized', 'BackgroundColor', [0.2,0.2,0.2], 'Color', [1,1,1], 'FontUnits', 'Normalized', 'FontSize', 0.02);
        
    end
    
    while toc < 0.008 % run at ~120hz. More frequent than this is not necessary due to  for screen refresh rate
        pause(0.002)
    end
    
    nIter = nIter + 1;
end

end


%% callback functions on figures and other subfunctions


function IMFigCallback(handle, event)
global rect
global DisplaySlice
global PrintSlice
global GetMaskName
global ClearMasks
global MakeNewRegionMask
global FinishNewMask
global SaveMask
global ClearMaskSlice
global WriteMask
global ChangeDim;
global PrintStack;
global FirstSlice;
global upSlice;
global downSlice;
global ZoomIn;
global ZoomOut;

figure(handle)


if event.Character == 'z'
    if sum(rect) >0
        rect = 0;
        DisplaySlice = 1;
    else
        
        ZoomTextInstructions = text(0.2,0.9, 'define zoom region box, then double click', 'Units', 'Normalized', 'Color', [1,1,1], 'FontUnits', 'Normalized', 'FontSize', 0.05);
        [~, rect] = imcrop(handle);
        rect = round(rect);
        DisplaySlice = 1;
        delete(ZoomTextInstructions);
    end
    
end



if event.Character == 'p'
    PrintSlice = 1;
    DisplaySlice = 1;
end

if event.Character == 'o'

    PrintStack = 1;
    FirstSlice = 1;

end


if event.Character == 'm'
    ClearMasks = 1;
    GetMaskName = 1;
    
end

if event.Character == 'c'
    ClearMasks = 1;
end

if event.Character == 'r'
    MakeNewRegionMask = 1;
end

if event.Character == 'i' % interpolate the missing slices of the mask
    FinishNewMask = 1;
end

if event.Character == 's' % save the new mask
    SaveMask = 1;
end

if event.Character == 'w' % write the new database to file
    WriteMask = 1;
end

if event.Character == 'x'
    ClearMaskSlice = 1;
end

if event.Character == 'v'
    ChangeDim = 1;
    
end

if strcmp(event.Key, 'rightarrow')    % trigger a move up one slice
    upSlice = 1;
end

if strcmp(event.Key, 'leftarrow')   % trigger a move down one slice
    downSlice = 1; 
end

if strcmp(event.Key, 'uparrow')    % trigger a zoom in

    ZoomIn = 1;
end

if strcmp(event.Key, 'downarrow')    % trigger a zoom in

    ZoomOut = 1;
end

end



function imFigButtonDown(handle, event)
global rect
global xSubtractor
global ySubtractor
global firstPos
global currentPos


xSubtractor = 0;
ySubtractor = 0;
firstPos = get(handle, 'CurrentPoint');
currentPos = firstPos;

if sum(rect) > 0
    set(handle, 'WindowButtonMotionFcn', @dragWindowPos)
end


end

function imFigButtonUp(handle, event)

set(handle, 'WindowButtonMotionFcn', '')
end

function dragWindowPos (handle, event)
global currentPos
global ShiftingTime

currentPos = get(handle, 'CurrentPoint');
ShiftingTime = 1;

end

function out = interp_shape_embedded(top,bottom,num)

% from
% http://stackoverflow.com/questions/18084698/interpolating-between-two-planes-in-3d-space
% credit to Frederick, http://stackoverflow.com/users/2007836/frederick

if nargin<2;
    error('not enough args');
end
if nargin<3;
    num = 1;
end
if ~num>0 && round(num)== num;
    error('number of slices to be interpolated must be integer >0');
end

top = signed_bwdist(top); % see local function below
bottom = signed_bwdist(bottom);

r = size(top,1);
c = size(top,2);
t = num+2;

[x y z] = ndgrid(1:r,1:c,[1 t]); % existing data
[xi yi zi] = ndgrid(1:r,1:c,1:t); % including new slice

out = interpn(x,y,z,cat(3,bottom,top),xi,yi,zi);
out = out(:,:,2:end-1)>=0;
end

function im = signed_bwdist(im)
im = -bwdist(bwperim(im)).*~im + bwdist(bwperim(im)).*im;
end


