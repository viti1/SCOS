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

function [ timeVec, rawSpeckleContrast , corrSpeckleContrast, meanVec , info] = ...
    SCOSvsTimeMultichannel(recordName,windowSize,plotFlag,resetMask)
if nargin <3
    plotFlag = true;
end
%% Constants
timePeriodForP2P = 2; % [s]

%% Check input parameters
if nargin == 0 % GUI mode
    plotFlag = 1;
    resetMask = 0;
    if exist('.\lastRec.mat','file')
        lastF = load('.\lastRec.mat');        
    else
        lastF.recordName = '';
    end
    [recordName] = uigetdir(fileparts(lastF.recordName));

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
    save('.\lastRec.mat','recordName')
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
maskFile = [recSavePrefix 'Mask.mat'];

% delete(maskFile)
if ~exist(maskFile,'file') || resetMask
    im = mean(ReadRecord(recordName, min(30,GetNumOfFrames(recordName))),3);
    expectedRadius = 125;
    [ masks, totMask, channels , figIm , h_circles] = autoFindRIOMultichannel(im,expectedRadius);
    set(figIm,'Position',[100,50,1200,800])
    if numel(masks)<7
        answer = questdlg('Draw more channels?', '','Yes','No','Discard Existing','Yes');
        drawOneMore = strcmp(answer,'Yes') || strcmp(answer,'Discard Existing') ;
        [x,y] = meshgrid(1:size(im,2),1:size(im,1));

        if strcmp(answer,'Discard Existing')
            masks = {};
            totMask = false(size(totMask));
            channels.Centers = [];
            channels.Radii = [];
            delete(h_circles);
        end
                    
        k = numel(masks)+1; 
        while drawOneMore
            circ = drawcircle('Color','r','FaceAlpha',0.2);
            channels.Centers(k,:) = circ.Center;
            channels.Radii(k) = circ.Radius;

            masks{k} = false(size(im,1),size(im,2));
            masks{k}((x-circ.Center(1)).^2 + (y-circ.Center(2)).^2 < circ.Radius^2 ) = true; 
            totMask = masks{k} | totMask ;
            answer = questdlg('One more channel?', '','Yes','No','Yes');
            drawOneMore = strcmp(answer,'Yes');
            k=k+1;
        end

        if exist('imfig','var'); close(imfig); end  
    end


    save(maskFile,'channels','masks','totMask');
else
    M = load(maskFile);
    if isfield(M,'mask')
        masks{1} = M.mask;
        totMask = M.mask;
        channels.Centers = M.circ.Center;
        channels.Radii  = M.circ.Radius;
    else
        load(maskFile);
    end
end

%% Read Record
disp(['Reading Record "' recordName '" ... '])
info = GetRecordInfo(recordName);

%% Set params
% detectorFolder = [fileparts(fileparts(mfilename('fullpath'))) '\Records\NoiseAndBackground\Basler_1440GS_Vika01\Mono8'];
% dData = load([detectorFolder '\ReadNoise\vsGain\readNoiseVsGain.mat']);  % detector Data
% nOfBits = 8;
% maxCapacity = 10.5e3;% [e]
% actualGain = ConvertGain(info.name.Gain,nOfBits,maxCapacity);
% readoutN   = interp1(dData.gainArr,dData.totNoise, info.name.Gain ,'spline');

%%  Calc Specle Contrast
disp(['Calculation SCOS on "' recordName '" ... '])
nOfFrames = GetNumOfFrames(recordName);
nOfChannels = numel(masks);

[ rawSpeckleContrast , corrSpeckleContrast , meanVec ] = InitNaN([nOfFrames 1],nOfChannels);
tic
for i=1:nOfFrames
%     fprintf('%d ',i);
    if mod(i,50)==0; fprintf('%d\t',i); end
    im = ReadRecord(recordName,1,i);
    stdIm = stdfilt(im,true(windowSize));
    meanIm = imfilter(im, true(windowSize)/windowSize^2,'conv','same'); % TBD!!! add Xiaojun algorithm
    meanImSquare = meanIm.^2;
    Kraw = (stdIm.^2)./meanImSquare ;
    for k = 1:nOfChannels
        rawSpeckleContrast{k}(i) = mean(Kraw(masks{k}));
        % corrSpeckleContrast{k}(i) = mean(Kraw(masks{k}) - actualGain./meanIm(masks{k}) - ( readoutN^2 + 1/12)./meanImSquare(masks{k}) );
        meanVec{k}(i) = mean(im(masks{k}));
    end
end
fprintf('\n');
toc
% TBD!! add SNR graph for each channel for corrected and not - corrected
%% Create Time vector
if ~isfield(info.name,'FR') || isnan(info.name.FR)
    error('Frame Rate must be part of the recording name as "FR"');
end
frameRate = info.name.FR; 

timeVec = (0:(nOfFrames-1))'*(1/frameRate) ;   % FR = FrameRate
p2p_time = timeVec<timePeriodForP2P;

%% Plot
stdStr = sprintf('Std%dx%d',windowSize,windowSize);

% infoFields = fieldnames(info.name);
% titleStr =  [ infoFields{1} '; exp=' num2str(info.name.expT)  'ms; Gain='  num2str(info.name.Gain) 'dB' ];

    
if plotFlag
    fig1 = figure('name',['SCOS ' recordName],'Position',[10 80 2025 970]);
    for k = 1:nOfChannels
        subplot(nOfChannels,3,3*k-2);
            plot(timeVec,meanVec{k});
            title(sprintf('Channel %d - mean I (<I>=%.0fDU)',k,mean(meanVec{k})));
        subplot(nOfChannels,3,3*k-1);
            plot(timeVec,rawSpeckleContrast{k})
            ylabel('Kraw^2')
            title(sprintf('Channel %d - Raw (<I>=%.0fDU)',k,mean(meanVec{k})));
        % subplot(nOfChannels,3,3*k);
        %     plot(timeVec,corrSpeckleContrast{k})
        %     ylabel('Kf^2')
        %     title(sprintf('Channel %d - Corrected partly ',k));
    end
end

for plot_i=1:nOfChannels*3
    subplot(nOfChannels,3,plot_i);
    xlim([0 timeVec(end)]);
    xlabel('Time [s]')
end

%TBD add SNR DATA!
%% Save
if exist('fig','var')
    savefig(fig1,[recSavePrefix 'Local' stdStr '_plot.fig']);
    save([recSavePrefix 'Local' stdStr '.mat'],'timeVec', 'corrSpeckleContrast' , 'rawSpeckleContrast','meanVec', 'info', 'recordName','windowSize');
end
 