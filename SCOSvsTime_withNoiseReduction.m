%  ---------------------------------------------------------------------------------------------------------
%  [ timeVec,  , rawSpeckleContrast , rawSpeckleVar, corrSpeckleVar , corrSpeckleContrast, imMeanVec , info] = PlotSCOSvsTime(recordName,windowSize,plotFlag,maskInput)
%  GUI mode:   - Choose the recording *folder*
%              - Choose widow size on which to do the std ( it will be used in stdfilter() function in order to calc the local std)
%              - First frame of the recording will appear, and then the
%                user should create a circle with the ROI (Region of Interest)
%                New figure will be created with speckleNoiseVec and speckleNoiseVecRel,
%                which is speckleNoiseVec devided by the the mean intensity inside the mask)
%                In the second time this function will run for the same
%                recording , the ROI is already saved , so ne need to
%                choose it again.
%                If the recording is .tiff files, FR (FrameRate) parameter should appeare in the name of the folder.
%                Example: murecord_FR20Hz 
%               
%  Command Mode: Same as GUI mode , but recordName and windowSize variables must be specified. 
%                plotFlag - [optional] create figure with 5 graphs, defualt
%                = true
%                
%                maskInput - true (then all the image is taken as mask) or
%                boolaen map (same size as record two first dimentions
%  ---------------------------------------------------------------------------------------------------------

function [ timeVec, rawSpeckleContrast , rawSpeckleVar, corrSpeckleVar , corrSpeckleContrast, imMeanVec , info] = ...
    SCOSvsTime_withNoiseReduction(recordName,backgroundName,windowSize,plotFlag,maskInput)
if nargin <3
    plotFlag = true;
end
%% Constants
timePeriodForP2P = 2; % [s]

%% Check input parameters
if nargin == 0 % GUI mode
    plotFlag = 1;
    [recordName] = uigetdir();
    if recordName == 0; return; end % if 'Cancel' was pressed
    if numel(dir([recordName, '\*.avi' ])) > 1 
        [recordRawName, recordDir] = uigetfile([recordName '\*.avi']);
        if recordRawName == 0; return; end % if 'Cancel' was pressed
        recordName = fullfile(recordDir, recordRawName);
    elseif ( numel(dir([recordName, '\*.avi' ])) + numel(dir([recordName, '\*.tiff' ])) + numel(dir([recordName, '\*.tif' ])) +  numel(dir([recordName, '\*.mat' ])) ) < 1 
        errordlg(['No .avi or .tiff/.tif or .mat files found in ' recordName ])
        error(['No .avi or .tiff/.tif or .mat files found in ' recordName ]);
    elseif numel(dir([recordName, '\*.avi' ])) == 1 && ( numel(dir([recordName, '\*.tiff' ])) + numel(dir([recordName, '\*.tif' ])) ) > 1 
        % if in folder apear both .avi and .tiff files -> assume that .avi is the recording
        d = dir([recordName, '\*.avi' ]);
        recordName = fullfile( recordName , d(1).name );
    end
    
    maxWindowSize = 50; minWindowSize = 3;
    answer = inputdlg('Window Size','',[1 25],{'9'});
    windowSize = str2double(answer{1});
    if isnan(windowSize) || windowSize > maxWindowSize || windowSize < minWindowSize
        errordlg(['Window Size must be a number between ' num2str(minWindowSize)  ' and ' num2str( num2str(maxWindowSize) ) ]);
        error(['Window Size must be a number between ' num2str(minWindowSize)  ' and ' num2str( num2str(maxWindowSize) ) ])
    end
    
    clear answer
end

%% Create Mask
upFolders = strsplit(recordName,filesep);
rawName = strrep( strjoin(upFolders(end-2:end),'; '), '_',' ');

if exist(recordName,'file') == 7 % it's a folder
    recSavePrefix = [ recordName filesep ];
else % it's a file
    recSavePrefix = [ recordName(1:find(recordName=='.',1,'last')-1) '_' ];
end

[ first_frame ] = ReadRecord(recordName,1); 

maskFile = [recSavePrefix 'Mask.mat'];
if exist('maskInput','var')
    if isequal(maskInput, 1)
        mask = true(size(first_frame));
    elseif islogical(maskInput)
        if ~isequal( size(maskInput), size(first_frame) ) 
            error('wrong maskInput');
        end
        mask = maskInput;
    else
        error('wrong maskInput type');
    end
    loadExistingFile_flag = 0;
else    
    if ~exist(maskFile,'file')
        % TBD : add automatic circle recognition
        loadExistingFile_flag = 0;    
    else
        if nargin == 0 % GUI mode
            answer = questdlg('Mask file already exist, do you want to define it again?','','Yes','No','No');
            loadExistingFile_flag = strcmp(answer,'No');
        else      % Command line mode -> load the saved mask 
            loadExistingFile_flag = 1;
        end
    end
end

if ~exist('mask','var')    
    if ~loadExistingFile_flag
        % ask the user to mark a circle on the image
        f = figure('position',[50,50,1200,800]); imagesc(first_frame); colormap gray; colorbar
        title(rawName,'interpreter','none');
        answer = questdlg('One Channel or multichannel?', 'What is the test type', ...
            'One Channel','MultiChannel','One Channel');
        multiChFlag = strcmp(answer,'MultiChannel');      
        
        mask = false(size(first_frame,1),size(first_frame,2));
        drawOneMore = true;
        while drawOneMore
            circ = drawcircle('Color','r','FaceAlpha',0.2);
            c.Center = circ.Center; c.Radius = circ.Radius;
            [x,y] = meshgrid(1:size(first_frame,2),1:size(first_frame,1));
            mask((x-circ.Center(1)).^2 + (y-circ.Center(2)).^2 < circ.Radius^2 ) = true;
            if multiChFlag
                answer = questdlg('One more channel?', '','Yes','No','Yes');
                drawOneMore = strcmp(answer,'Yes');
            else
                drawOneMore = false;
            end
        end
        save(maskFile,'mask','circ');
        close(f)
    else 
        load(maskFile)
    %     figure; imshowpair(first_frame,mask);  
    end
else
    
end
%% Get Background
if ~exist(backgroundName,'file')
    error([backgroundName,' does not exist!']);
end

if exist(backgroundName,'file') == 7 % it's a folder
    if exist( [ backgroundName '\meanIm.mat'],'file')
        bgS = load([ backgroundName '\meanIm.mat']);
        background = bgS.meanIm;
    else
        background = mean(ReadRecord(backgroundName),3);
        meanIm = background;
        save([ backgroundName '\meanIm.mat'], 'meanIm');
        clear recMean 
    end
elseif endsWith(backgroundName,'.mat')
    bgS = load( backgroundName );
    background = bgS.meanIm;    
end

% background = zeros(size(first_frame));

%% Check info
nOfFrames = GetNumOfFrames(recordName);
info = GetRecordInfo(recordName);
info_background = GetRecordInfo(backgroundName);
im1 = ReadRecord(recordName,1);
if ~isequal(size(background),size(im1))
    error('The background should be the same picture size as the record. Record size is [%d,%d], but background size is [%d,%d]',size(im1,1), size(im1,2), size(background,1), size(background,1));
end

fields = {'expT','Gain'};
for fi = 1:numel(fields)
    param = fields{fi};
    if info_background.name.(param) ~= info_background.name.(param)
        error('Background.%s=%g   Record.%s=%g',param, info_background.name.(param), param, info_background.name.(param));
    end
end  

%% Get ReadoutNoise
dData = load('./cameraData/Basler_1440GS_Vika01/readNoiseVsGain.mat');  % detector Data
readoutN   = interp1(dData.gainArr,dData.readoutNoise, info.Gain ,'spline');
actualGain = 10^(info.Gain+detectorData.g0); 
H = fspecial('average', [1 1]*windowSize);

%% cut smaller ROI
roi_half_size = ceil((c.Radius+windowSize)) ;
roi = [ round(c.Center(2)) + (-roi_half_size:roi_half_size) ;
        round(c.Center(1)) + (-roi_half_size:roi_half_size) ];
cut_background  = background(  roi(1,:)  , roi(2,:));
cut_mask        = mask(roi(1,:)  , roi(2,:));

%% Calc spatialNoise 
numFramesForSPNoise = 500;
spRec = ReadRecord(recordName,numFramesForSPNoise);
cut_spRec = spRec(  roi(1,:)  , roi(2,:));
spIm = mean(cut_spRec,3);
spVar = stdfilt( spIm ,true(windowSize)).^2;
[fitI_A,fitI_B] = FitMeanIm(cut_spRec);

%% Calc Specle Contrast
disp(['Calculating SCOS on "' recordName '" ... ']);


% init loop vars
[rawSpeckleVar , rawSpeckleContrast , corrSpeckleVar , corrSpeckleContrast , imMeanVec] = InitNaN([nOfFrames 1]);
for i=1:nOfFrames
    im = ReadRecord(recordName,1,i);
    cut_im   = im(  roi(1,:)  , roi(2,:)) - cut_background;

    stdIm = stdfilt(cut_im,true(windowSize));
%     meanIm = imfilter(cut_im, true(windowSize)/windowSize^2,'conv','same'); % TBD!!! check if valid
    meanFrame = mean(cut_im(cut_mask));
    fittedI = fitI_A*meanFrame + fitI_B ;
    
    meanImSquare = fittedI.^2;
    Kraw = (stdIm.^2)./meanImSquare ;
    rawSpeckleContrast(i) = mean(Kraw(cut_mask));
    corrSpeckleContrast(i) = mean(Kraw(cut_mask) - actualGain./fittedI(cut_mask) - ( readoutN^2 + 1/12 + spVar(cut_mask))./meanImSquare(cut_mask) );
    
    imMeanVec(i) = meanFrame; 
end
IMean = mean(imMeanVec);

%% Create Time vector
if isfield(info,'cam') && isfield(info.cam,'AcquisitionFrameRate')
    frameRate = info.cam.AcquisitionFrameRate;
else    
    if ~isfield(info.name,'FR') || isnan(info.name.FR)
        error('Frame Rate must be part of the recording name as "FR"');
    end
    frameRate = info.name.FR; 
end

timeVec = (0:(nOfFrames-1))'*(1/frameRate) ;   % FR = FrameRate
p2p_time = timeVec<timePeriodForP2P;

%% Plot
stdStr = sprintf('Std%dx%d',windowSize,windowSize);

infoFields = fieldnames(info.name);
if isfield(info.name,'SDS')
    titleStr =  [ infoFields{1} ' SDS=' num2str(info.name.SDS)  '; exp=' num2str(info.name.expT)  'ms; Gain='  num2str(info.name.Gain) 'dB' ];
else
    titleStr =  [ infoFields{1} '; exp=' num2str(info.name.expT)  'ms; Gain='  num2str(info.name.Gain) 'dB' ];
end
    
if plotFlag
    fig = figure('name',['SCOS ' recordName]);
    subplot(3,1,1);
        plot(timeVec,corrSpeckleContrast)
        ylabel('Corrected Contrast (var/I^2)')
%         title({ rawName , ['p2p = ' ]})
        title(titleStr)
    subplot(3,1,2);
        plot(timeVec,rawSpeckleContrast)
        ylabel('Raw Contrast (var/I^2)')
    subplot(3,1,3);
        plot(timeVec,imMeanVec);
        ylabel('I [DU]')
end

for plot_i=1:3
    subplot(3,1,plot_i);
    xlabel('Time [s]')
end
%% Save
if exist('fig','var')
    savefig(fig,[recSavePrefix 'Local' stdStr '_plot.fig']);
    save([recSavePrefix 'Local' stdStr '.mat'],'timeVec', 'corrSpeckleContrast' , 'rawSpeckleContrast','rawSpeckleVar','corrSpeckleVar', 'imMeanVec', 'info', 'recordName','windowSize');
end
 
