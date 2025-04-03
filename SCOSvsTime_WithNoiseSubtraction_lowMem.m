%  ---------------------------------------------------------------------------------------------------------
%  [ timeVec,  , rawSpeckleContrast , rawSpeckleVar, corrSpeckleVar , corrSpeckleContrast, meanVec , info] = PlotSCOSvsTime(recName,windowSize,plotFlag,maskInput)
%  GUI mode:   - Choose the recording folder
%              - Choose widow size ( it will be used in stdfilter() function in order to calc the local std)
%              - Choose dark recording folder
%
%                First frame of the recording will appear, on witch the
%                user should draw a circle with the ROI (Region of Interest)
%                
%                In the second time this function will run for the same
%                recording , the ROI is already saved , so ne need to
%                choose it again.
%                Main and dark recording folder names should be in the following format :
%                <Name>_Gain<X>dB_expT<>ms_FR<X>Hz_BL<>DU 
%                Where expT -> exposure time in ms, FR -> Frame Rate, BL -> Black Level
%                If the dark recording is located in the same folder as the Main one, and starts with "background"
%                it is automatically recognized.
%
%  Command Mode: Same as GUI mode , but recName and windowSize variables must be specified. 
%                plotFlag - [optional] defualt= true                
%                maskInput - can be "true" (then all the image is taken as mask) or
%                a boolaen map the same size as the record.
%  ---------------------------------------------------------------------------------------------------------

function  [ timeVec, rawSpeckleContrast , corrSpeckleContrast, meanVec , info] = ...
    SCOSvsTime_WithNoiseSubtraction_Ver2(recName,backgroundName,windowSize,plotFlag,maskInput)
if nargin <3
    plotFlag = true;
end
addpath('.\baseFunc')
%% Constants
timePeriodForP2P = 2; % [s]

%% Check input parameters
if nargin == 0 % GUI mode
    plotFlag = 1;
    if exist('.\lastRec.mat','file')
        lastF = load('.\lastRec.mat');        
    else
        lastF.recName = [ fileparts(pwd) '\Records' ];
    end
    
    [recName] = uigetdir(fileparts(lastF.recName));
    if recName == 0; return; end % if 'Cancel' was pressed
    if numel(dir([recName, '\*.avi' ])) > 1 
        [recordRawName, recordDir] = uigetfile([recName '\*.avi']);
        if recordRawName == 0; return; end % if 'Cancel' was pressed
        recName = fullfile(recordDir, recordRawName);
    elseif ( numel(dir([recName, '\*.avi' ])) + numel(dir([recName, '\*.tiff' ])) + numel(dir([recName, '\*.tif' ])) +  numel(dir([recName, '\*.mat' ])) ) < 1 
        errordlg(['No .avi or .tiff/.tif or .mat files found in ' recName ])
        error(['No .avi or .tiff/.tif or .mat files found in ' recName ]);
    elseif numel(dir([recName, '\*.avi' ])) == 1 && ( numel(dir([recName, '\*.tiff' ])) + numel(dir([recName, '\*.tif' ])) ) > 1 
        % if in folder apear both .avi and .tiff files -> assume that .avi is the recording
        d = dir([recName, '\*.avi' ]);
        recName = fullfile( recName , d(1).name );
    end
    
    save('.\lastRec.mat','recName')
    
    maxWindowSize = 50; minWindowSize = 3;
    answer = inputdlg('Window Size','',[1 25],{'7'});
    windowSize = str2double(answer{1});
    if isnan(windowSize) || windowSize > maxWindowSize || windowSize < minWindowSize
        errordlg(['Window Size must be a number between ' num2str(minWindowSize)  ' and ' num2str( num2str(maxWindowSize) ) ]);
        error(['Window Size must be a number between ' num2str(minWindowSize)  ' and ' num2str( num2str(maxWindowSize) ) ])
    end
    clear answer
end

isRecordFile = exist(recName,'file') == 2;
if nargin == 0  || isempty(backgroundName)
    if exist([recName , '_dark'],'dir')
        backgroundName = [recName , '_dark'];
    else
        dir_Background = [ dir([fileparts(recName) , '\DarkIm*']) dir([fileparts(recName) , '\background*']) dir([fileparts(recName) , '\BG_*'])   ] ;
    
        if  isempty(dir_Background) || numel(dir_Background) > 1 
            if isRecordFile
                backgroundName = uigetfile( fileparts(recName) ,'Please Select the background');
            else
                backgroundName = uigetdir( fileparts(recName) ,'Please Select the background'); 
            end
        else
            backgroundName = fullfile(fileparts(recName),dir_Background(1).name);
        end
        if isequal(backgroundName,0)
            disp('Aborting...')
            return;
        end
    end
end

%% Create Mask
upFolders = strsplit(recName,filesep);
rawName = strrep( strjoin(upFolders(end-2:end),'; '), '_',' ');

if exist(recName,'file') == 7 % it's a folder
    recSavePrefix = [ recName filesep ];
else % it's a file
    recSavePrefix = [ recName(1:find(recName=='.',1,'last')-1) '_' ];
end

[ mean_frame ] = mean( ReadRecord(recName,20),3); 

maskFile = [recSavePrefix 'Mask.mat'];
if exist('maskInput','var')
    if isequal(maskInput, true)
        masks = {true(size(mean_frame))};
        totMask = masks{1};
        disp('Mask is the whole image')
    elseif iscell(maskInput)
        masks = maskInput;
        for k=1:numel(masks)
            if ~isequal( size(masks{k}), size(mean_frame) ) 
                error('wrong maskInput size, should be the same as the recording');
            end
        end
        totMask = masks{1};
        for k=2:numel(masks)
            totMask = totMask | masks{k};
        end
        disp('Input mask is a cell array')
    elseif islogical(maskInput)
        if ~isequal( size(maskInput), size(mean_frame) ) 
            error('wrong maskInput');
        end
        masks = {maskInput};
        totMask = masks{1};
        disp('Input single mask')
    else
        error('wrong maskInput type');
    end
    loadExistingFile_flag = 0;
else    
    if ~exist(maskFile,'file')
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

% -- get ROI if needed and cut the margins 
ws2 = ceil(windowSize/2); % for margins marking as false  
if ~exist('masks','var')
    if ~loadExistingFile_flag
        % [channels, masks, totMask, figIm] = CreateMask(recName);
        [ totMask , circ , figMask] = GetROI(mean(ReadRecord(recName,20),3),windowSize);
                
        masks{1} = totMask;       
        channels.Centers = circ.Center;
        channels.Radii = circ.Radius;
        save(maskFile,'masks','channels','totMask');
        savefig(figMask,[recName '\maskIm.fig'])
        % close(figIm)
    else 
        M = load(maskFile);
        if isfield(M,'mask')
            masks{1} = false(size(mask));
            masks{1}(ws2+1:end-ws2,ws2+1:end-ws2) = mask(ws2+1:end-ws2);
            totMask = M.mask > 0;
            channels.Centers = M.circ.Center;
            channels.Radii  = M.circ.Radius;
        else
            load(maskFile); %#ok<LOAD>
            for k=1:numel(masks)
                masks{k}( [ 1:ws2 (end-ws2+1):end ],:) = false; %#ok<AGROW>
                masks{k}( : , [ 1:ws2 (end-ws2+1):end ]) = false; %#ok<AGROW>
            end
            totMask( [ 1:ws2 (end-ws2+1):end ], : ) = false;
            totMask( : , [ 1:ws2 (end-ws2+1):end ]) = false;            
        end
    end
else
    for k=1:numel(masks)
        masks{k}( [ 1:ws2 (end-ws2+1):end ],:) = false; %#ok<AGROW>
        masks{k}( : , [ 1:ws2 (end-ws2+1):end ]) = false; %#ok<AGROW>
    end
    totMask( [ 1:ws2 (end-ws2+1):end ], : ) = false;
    totMask( : , [ 1:ws2 (end-ws2+1):end ]) = false;
end

%% Check info
upFolders = strsplit(recName,filesep);
shortRecName = strjoin(upFolders(end-2:end));

nOfFrames = GetNumOfFrames(recName);
info = GetRecordInfo(recName);
im1 = ReadRecord(recName,1);

if backgroundName~=0
    info_background = GetRecordInfo(backgroundName);
    fields = {'expT','Gain','BL'};
    for fi = 1:numel(fields)
        param = fields{fi};
        if info_background.name.(param) ~= info.name.(param)
            error('Background.%s=%g   Record.%s=%g',param, info_background.name.(param), param, info_background.name.(param));
        end
    end 
end

if isfield(info,'cam') && isfield(info.cam,'AcquisitionFrameRate')
    frameRate = info.cam.AcquisitionFrameRate;
else    
    if ~isfield(info.name,'FR') || isnan(info.name.FR)
        error('Frame Rate must be part of the recording name as "FR"');
    end
    frameRate = info.name.FR; 
end
%% Get Background and Background Noise
start_calib_time = tic;
disp('Load Background');
if ~exist(backgroundName,'file')
    error([backgroundName,' does not exist!']);
end

if exist(backgroundName,'file') == 7 % it's a folder
    if exist( [ backgroundName '\meanIm.mat'],'file')
        bgS = load([ backgroundName '\meanIm.mat']);
        if isfield(bgS,'recVar')
            background = bgS.recMean - info_background.name.BL;
            darkVar = bgS.recVar;
        else
            delete([ backgroundName '\meanIm.mat']);
            [ background , darkVar ] = ReadRecordVarAndMean( backgroundName );
            background =  background - info_background.name.BL;
        end
        clear bgS
    else
        [ background , darkVar ] = ReadRecordVarAndMean( backgroundName );
        background =  background - info_background.name.BL;

        if abs(mean2(background)) > 3
            my_imagesc(background); title('Background');
            warning('Suspicious level of the background %gDU !', round(mean2(background),2));
        end        
    end    
elseif endsWith(backgroundName,'.mat')
    bgS = load( backgroundName );
    fields = fieldnames(bgS);
    if ismember(fields, 'recMean')
        background = bgS.recMean - info_background.name.BL;
    elseif startsWith(fields{1}, 'Video')
        darkRec = bgS.(fields{1}); 
        bgS.recMean = mean(darkRec,3);
        background = bgS.recMean - info_background.name.BL;
        darkVar = std(darkRec,0,3).^2;
        bgS.recVar   = darkVar;
        save(backgroundName,'-struct','bgS')
    else
        error('wrong fields')
    end
    clear bgS
end

darkVarPerWindow = imboxfilt(darkVar,windowSize) ;

if ~isequal(size(background),size(im1))
    error('The background should be the same picture size as the record. Record size is [%d,%d], but background size is [%d,%d]',size(im1,1), size(im1,2), size(background,1), size(background,1));
end
%% Get G[DU/e]
nOfBits = info.nBits;   
actualGain = GetActualGain(info);

%% Calc spatialNoise
disp('Calculate Spatial Noise')
if ~isfield(info.name , 'BL' )    
    BlackLevel = 0;
else
    BlackLevel = info.name.BL;
end

if isRecordFile
    meanImFile = [fileparts(recName)  '\meanIm.mat'];
else
    meanImFile = [recName  '\meanIm.mat'];
end

numFramesForSPNoise = 600;

if exist(meanImFile,'file') 
    M = load(meanImFile);
    spIm  = M.recMean;
    clear M
else
    spIm = ReadRecordVarAndMean(recName,numFramesForSPNoise);
end
spVar = stdfilt( spIm ,true(windowSize)).^2;
%% Decrease Image Size
[y,x] = find(totMask) ;
roi_lims  = [ min(y)-windowSize  , max(y)+windowSize
              min(x)-windowSize  , max(x)+windowSize ];

roi_lims( roi_lims < 1 ) = 1;        
if roi_lims(1,2) > info.imageSize(1)
    roi_lims(1,2) = info.imageSize(1);
end
if roi_lims(2,2) > info.imageSize(2)
    roi_lims(2,2) = info.imageSize(2);
end
roi.y = roi_lims(1,1):roi_lims(1,2);
roi.x = roi_lims(2,1):roi_lims(2,2);

masks_cut = cell(size(masks));
for ch = 1:numel(masks)    
    masks_cut{ch} = masks{ch}(roi.y  , roi.x); 
end
   
% fitI_A_cut =  fitI_A(roi.y  , roi.x);
% fitI_B_cut =  fitI_B(roi.y  , roi.x);
% bpMap_cut = bpMap(roi.y  , roi.x);

disp('End Calib')
toc(start_calib_time)
%% Calc Specle Contrast
disp(['Calculating SCOS on "' recName '" ... ']);
disp(['Mono' num2str(nOfBits)]);
nOfChannels = numel(masks);
frameNames = dir([recName '\*.tiff']);
[~,sort_ind] = sort([frameNames.datenum]);
frameNames = frameNames(sort_ind);

% init loop vars
[ rawSpeckleContrast , corrSpeckleContrast , meanVec] =InitNaN([nOfFrames 1],nOfChannels);

if isRecordFile
    rec = ReadRecord(recName);
    im1 = double(rec(:,:,1));
else
	im1 = double(imread([recName,filesep,frameNames(1).name])) ;
end
devide_by_16 = nOfBits == 12  && all(mod(im1(:),16) == 0);
devide_by_64 = nOfBits == 10  && all(mod(im1(:),64) == 0);

start_scos = tic;
               
for i=1:nOfFrames
    if i == 50
        time50frames = toc(start_scos);
        fprintf('\n Estimated Time = %g min (%d frames)\n',round(time50frames/50*nOfFrames/60,2), nOfFrames)
    end
    if mod(i,200) == 0 
        fprintf('%d\t',i);
        if mod(i,2000) == 0
            fprintf('\n');
        end
    end
    if isRecordFile 
        im_raw = double(rec(:,:,i));
    else
        frameNames(i).name
        im_raw = double(imread([recName,filesep,frameNames(i).name])) ;  
        if devide_by_16
            im_raw = im_raw/16;
        elseif devide_by_64
            im_raw = im_raw/64;
        end
        im_raw = im_raw - BlackLevel;
    end
    
    im = im_raw - background;
    im_cut = im(roi.y,roi.x);
    stdIm = stdfilt(im_cut,true(windowSize));
    meanIm = imboxfilt(im_cut,windowSize); % TBD remove
    
    for ch = 1:nOfChannels
        meanFrame = mean(im_cut(masks_cut{ch}));
%         fittedI = fitI_A_cut*meanFrame + fitI_B_cut ;         
%         fittedISquare = fittedI.^2;
        fittedI = meanIm;
        fittedISquare = meanIm.^2;

        rawSpeckleContrast{ch}(i) = mean((stdIm(masks_cut{ch}).^2 ./ fittedISquare(masks_cut{ch})));
        corrSpeckleContrast{ch}(i) = mean( ( stdIm(masks_cut{ch}).^2 - actualGain.*fittedI(masks_cut{ch})  - spVar(masks_cut{ch}) - 1/12 - darkVarPerWindow(masks_cut{ch}))./fittedISquare(masks_cut{ch}) ); 
        meanVec{ch}(i) = meanFrame;
        if i==1
            fprintf('<I>=%.3gDU , K_raw = %.5g , Ks=%.5g , Kr=%.5g, Ksp=%.5g, Kq=%.5g, Kf=%.5g\n',meanFrame,rawSpeckleContrast{ch}(i), ...
               mean(actualGain.*fittedI(masks_cut{ch})./fittedISquare(masks_cut{ch})),mean(darkVarPerWindow(masks_cut{ch})./fittedISquare(masks_cut{ch})),...
               mean(spVar(masks_cut{ch})./fittedISquare(masks_cut{ch})),mean(1./(12*fittedISquare(masks_cut{ch}))),corrSpeckleContrast{ch}(i));
        end
    end 
end
fprintf('\n');
%% Create Time vector
% timeVec = (0:(nOfFrames-1))'*(1/frameRate) ;   % FR = FrameRate
timeVec = ([frameNames.datenum]' - frameNames(1).datenum )*24*60*60;  % time vec in seconds
p2p_time = timeVec<timePeriodForP2P;
%% Save
stdStr = sprintf('Std%dx%d',windowSize,windowSize);
if exist([recSavePrefix 'Local' stdStr '.mat'],'file'); delete([recSavePrefix 'Local' stdStr '.mat']); end % just for it to have the right date
firstFrameDir = dir([recName,'\*_0001.tiff']);
startDateTime = firstFrameDir.date;
save([recSavePrefix 'Local' stdStr '_corr.mat'],'startDateTime','timeVec', 'corrSpeckleContrast' , 'rawSpeckleContrast', 'meanVec', 'info','nOfChannels', 'recName','windowSize');

%% Plot
infoFields = fieldnames(info.name);
if ~isfield(info.name,'Gain')
    info.name.Gain = '';
end
firtsParamValue = info.name.(infoFields{1});
if ~ischar(firtsParamValue)
    firtsParamValue = num2str(firtsParamValue);
end
if isfield(info.name,'SDS')
    titleStr =  [ infoFields{1} firtsParamValue ' SDS=' num2str(info.name.SDS)  '; exp=' num2str(info.name.expT)  'ms; Gain='  num2str(info.name.Gain) 'dB' ];
else
    titleStr =  [ infoFields{1} firtsParamValue '; exp=' num2str(info.name.expT)  'ms; Gain='  num2str(info.name.Gain) 'dB' ];
end
    
[raw_SNR,  raw_FFT , raw_freq, raw_pulseFreq, raw_pulseBPM] = CalcSNR_Pulse(rawSpeckleContrast{1},frameRate);
[corr_SNR,  corr_FFT , corr_freq, corr_pulseFreq, corr_pulseBPM] = CalcSNR_Pulse(corrSpeckleContrast{1},frameRate);

if  plotFlag
    fig = figure('name',['SCOS ' recName ' Mono' num2str(nOfBits)],'Units','Normalized','Position',[0.1 0.1 0.8 0.8]);
    subplot(3,2,1);
        plot(timeVec,corrSpeckleContrast{1})
        ylabel('Corrected Contrast (var/I^2)')
%         title({ rawName , ['p2p = ' ]})
        title({titleStr, 'Corrected Signal'})
    subplot(3,2,2);
        plot(corr_freq,corr_FFT)
        ylabel(' FFT '); xlabel('f [Hz]')
        title(sprintf('Corrected FFT: SNR=%.2g Pulse=%.0fbpm',corr_SNR,corr_pulseBPM));
    subplot(3,2,3);
        plot(timeVec,rawSpeckleContrast{1})
        ylabel('Raw Contrast (var/I^2)')
        title('Raw Signal')
    subplot(3,2,4);
        plot(raw_freq,raw_FFT)
        ylabel(' FFT '); xlabel('f[Hz]')
        title(sprintf('Raw FFT: SNR=%.2g Pulse=%.0fbpm',raw_SNR,raw_pulseBPM));
    subplot(3,2,5);
        plot(timeVec,meanVec{1});
        ylabel('I [DU]')
        

        for plot_i=1:2:5
            subplot(3,2,plot_i);
            xlabel('Time [s]')
        end

    savefig(fig,[recSavePrefix 'Local' stdStr '_plot.fig']);
end

%% Plot rBFI
if any(corrSpeckleContrast{1} < 0 )
    warning('Error: There are negative values in the contrast !!!');
    warndlg('Error: There are negative values in the contrast !!! BFI has no meaning')
end

if plotFlag
    BFi = 1./corrSpeckleContrast{1};
    if timeVec(end) > 120
        timeToPlot = timeVec / 60; % convert to min
        xLabelStr = 'time [min]';
        rBFi = BFi/mean(BFi(1:round(10*frameRate))); % normalize by first 10 seconds        
    else
        timeToPlot = timeVec ; % convert to min
        xLabelStr = 'time [sec]';
        rBFi = BFi/prctile(BFi(1:round(10*frameRate)),5); % normalize by 5% percentile in first 10 sec
    end   

    fig7 = figure('Name',['rBFi: '  recName ],'Units','Normalized','Position',[0.1,0.1,0.4,0.4]); 
    subplot(2,1,1);    
    plot(timeToPlot,rBFi);
    ylim([0 10])
    title(titleStr)
    xlabel(xLabelStr)
    ylabel('rBFi');
    grid on
    hold on;
    set(gca,'FontSize',10);
    subplot(2,1,2)
    plot(timeToPlot,meanVec{1}); 
    xlabel(xLabelStr)
    ylabel('<I> [DU]');
    set(gca,'FontSize',10);
    grid on
    % tLim=10; subplot(2,1,1); xlim([0 tLim]); subplot(2,1,2); xlim([0 tLim])
    savefig(fig7,[recSavePrefix '_rBFi.fig']);
end
%%
toc(start_scos)
end

