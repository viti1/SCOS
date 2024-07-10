%  ---------------------------------------------------------------------------------------------------------
%  [ timeVec,  , rawSpeckleContrast , rawSpeckleVar, corrSpeckleVar , corrSpeckleContrast, meanVec , info] = PlotSCOSvsTime(recordName,windowSize,plotFlag,maskInput)
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

function [ timeVec, rawSpeckleContrast , rawSpeckleVar, corrSpeckleVar , corrSpeckleContrast, meanVec , info] = ...
    SCOSvsTime_WithNoiseSubtraction(recordName,backgroundName,windowSize,plotFlag,maskInput)
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
        lastF.recordName = [ fileparts(pwd) '\Records' ];
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
    save('.\lastRec.mat','recordName')
    
    maxWindowSize = 50; minWindowSize = 3;
    answer = inputdlg('Window Size','',[1 25],{'9'});
    windowSize = str2double(answer{1});
    if isnan(windowSize) || windowSize > maxWindowSize || windowSize < minWindowSize
        errordlg(['Window Size must be a number between ' num2str(minWindowSize)  ' and ' num2str( num2str(maxWindowSize) ) ]);
        error(['Window Size must be a number between ' num2str(minWindowSize)  ' and ' num2str( num2str(maxWindowSize) ) ])
    end
    clear answer
    
    dir_Background = [ dir([fileparts(recordName) , '\DarkIm*']) dir([fileparts(recordName) , '\Background*']) ] ;
    if isempty(dir_Background) || numel(dir_Background) > 1
        backgroundName = uigetdir( fileparts(recordName) ,'Please Select the background');    
    else
        backgroundName = fullfile(fileparts(recordName),dir_Background(1).name);
    end
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
        % [channels, masks, totMask, figIm] = CreateMask(recordName);
        % TBD!
        [ mask , circ ] = GetROI(mean(ReadRecord(recordName,20),3));
        masks{1} = mask;
        channels.Centers = circ.Center;
        channels.Radii = circ.Radius;
        totMask = mask;
        save(maskFile,'masks','channels','totMask');
        % close(figIm)
    else 
        M = load(maskFile);
        if isfield(M,'mask')
            masks{1} = M.mask;
            totMask = M.mask;
            channels.Centers = M.circ.Center;
            channels.Radii  = M.circ.Radius;
        else
            load(maskFile); %#ok<LOAD>
        end
    end
else
    
end



%% Check info
upFolders = strsplit(recordName,filesep);
shortRecName = strjoin(upFolders(end-2:end));

nOfFrames = GetNumOfFrames(recordName);
info = GetRecordInfo(recordName);
im1 = ReadRecord(recordName,1);

if backgroundName~=0
    info_background = GetRecordInfo(backgroundName);
    fields = {'expT','Gain'};
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
%% Get Background
disp('Load Background');

if ~isequal(backgroundName,0)
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
        end
    elseif endsWith(backgroundName,'.mat')
        bgS = load( backgroundName );
        background = bgS.meanIm;    
    end
else 
    background = zeros(size(im1));
end

if ~isequal(size(background),size(im1))
    error('The background should be the same picture size as the record. Record size is [%d,%d], but background size is [%d,%d]',size(im1,1), size(im1,2), size(background,1), size(background,1));
end
%% Get ReadoutNoise (for current Gain)
nOfBits = info.nBits;
% camDataFolder = [fileparts(mfilename('fullpath')) '\camerasData'];
% camDataFile = dir([ camDataFolder '\SN' num2str(info.cameraSN) '*readNoiseVsGain_Mono' num2str(nOfBits) '.mat']);
% if numel(camDataFile) == 0
%     error('Camera Data file was not found');
% elseif numel(camDataFile) > 1
%     disp(camDataFile.name)
%     error('more that one camera file was found');
% end
% 
% dData = load([camDataFolder filesep camDataFile(1).name]); 
maxCapacity = 10.5e3;% [e]
actualGain = ConvertGain(info.name.Gain,nOfBits,maxCapacity);
% readoutN   = interp1(ConvertGain(dData.gainArr,nOfBits,maxCapacity),dData.tempNoise, actualGain ,'spline');
dir_ReadoutNoise = dir([fileparts(recordName) , '\ReadNoise*']);
if isempty(dir_ReadoutNoise) || numel(dir_ReadoutNoise) > 1
    [folderReadNoise] = uigetdir(fileparts(recordName),'Please Enter Path to Read Noise Recording');
else
    folderReadNoise = fullfile(fileparts(recordName) , dir_ReadoutNoise(1).name);
end
disp('Calc ReadNoise Mat')
recReadNoise = ReadRecord(folderReadNoise);
readNoiseMat = std(recReadNoise,0,3);
readNoise = imboxfilt(readNoiseMat,windowSize) ;    

%% Calc spatialNoise 
disp('Calc Spatial Noise');
numFramesForSPNoise = 400;
if nOfFrames > 500 ;  numFramesForSPNoise=500; end
spRec = ReadRecord(recordName,numFramesForSPNoise);
spIm = mean(spRec,3) - background;
spVar = stdfilt( spIm ,true(windowSize)).^2;
[fitI_A,fitI_B] = FitMeanIm(spRec,totMask,windowSize);

%% Calc Specle Contrast
disp(['Calculating SCOS on "' recordName '" ... ']);
nOfChannels = numel(masks);
% init loop vars
[ rawSpeckleContrast , corrSpeckleContrast , meanVec] =InitNaN([nOfFrames 1],nOfChannels);
for i=1:nOfFrames
    if mod(i,50) == 0; fprintf('%d\t',i); end
    im = ReadRecord(recordName,1,i);
    
    stdIm = stdfilt(im,true(windowSize));
%     meanIm = imfilter(im, true(windowSize)/windowSize^2,'conv','same'); 
    
    for ch = 1:nOfChannels
        meanFrame = mean(im(masks{ch}));
        fittedI = fitI_A*meanFrame + fitI_B ;
        fittedISquare = fittedI.^2;
%         rawSpeckleContrast(i) = mean((stdIm_raw(masks{ch}) ./ meanIm(masks{ch}))).^2;
        
        rawSpeckleContrast{ch}(i) = mean((stdIm(masks{ch}).^2 ./ fittedISquare(masks{ch})));
        corrSpeckleContrast{ch}(i) = mean( ( stdIm(masks{ch}).^2 - actualGain.*fittedI(masks{ch})  - spVar(masks{ch}) - 1/12 - readNoise(masks{ch}).^2)./fittedISquare(masks{ch}) ); % - ( readoutN^2 )./fittedISquare(masks{ch}) );
        meanVec{ch}(i) = meanFrame;
    end 
end
fprintf('\n');
%% Create Time vector
timeVec = (0:(nOfFrames-1))'*(1/frameRate) ;   % FR = FrameRate
p2p_time = timeVec<timePeriodForP2P;
%% Save
stdStr = sprintf('Std%dx%d',windowSize,windowSize);
delete([recSavePrefix 'Local' stdStr '.mat']); % just for it to have the right date
save([recSavePrefix 'Local' stdStr '_corr.mat'],'timeVec', 'corrSpeckleContrast' , 'rawSpeckleContrast', 'meanVec', 'info','nOfChannels', 'recordName','windowSize');

%% Plot
infoFields = fieldnames(info.name);
if isfield(info.name,'SDS')
    titleStr =  [ infoFields{1} ' SDS=' num2str(info.name.SDS)  '; exp=' num2str(info.name.expT)  'ms; Gain='  num2str(info.name.Gain) 'dB' ];
else
    titleStr =  [ infoFields{1} '; exp=' num2str(info.name.expT)  'ms; Gain='  num2str(info.name.Gain) 'dB' ];
end
    
[raw_SNR,  raw_FFT , raw_freq, raw_pulseFreq, raw_pulseBPM] = CalcSNR_Pulse(rawSpeckleContrast{1},frameRate);
[corr_SNR,  corr_FFT , corr_freq, corr_pulseFreq, corr_pulseBPM] = CalcSNR_Pulse(corrSpeckleContrast{1},frameRate);

if  plotFlag
    fig = figure('name',['SCOS ' recordName ' Mono' num2str(nOfBits)]);
    subplot(3,2,1);
        plot(timeVec,corrSpeckleContrast{1})
        ylabel('Corrected Contrast (var/I^2)')
%         title({ rawName , ['p2p = ' ]})
        title({titleStr, 'Corrected Signal'})
    subplot(3,2,2);
        plot(corr_freq,corr_FFT)
        ylabel(' FFT ')
        title(sprintf('Corrected FFT: SNR=%.2g Pulse=%.0fbpm',corr_SNR,corr_pulseBPM));
    subplot(3,2,3);
        plot(timeVec,rawSpeckleContrast{1})
        ylabel('Raw Contrast (var/I^2)')
        title('Raw Signal')
    subplot(3,2,4);
        plot(raw_freq,raw_FFT)
        ylabel(' FFT ')
        title(sprintf('Raw FFT: SNR=%.2g Pulse=%.0fbpm',raw_SNR,raw_pulseBPM));
    subplot(3,2,5);
        plot(timeVec,meanVec{1});
        ylabel('I [DU]')
        

        for plot_i=1:5
            subplot(3,2,plot_i);
            xlabel('Time [s]')
        end

    savefig(fig,[recSavePrefix 'Local' stdStr '_plot_with_corrected.fig']);

end

%% Correct Jumps
nOfChannels = numel(corrSpeckleContrast);
params = GetParamsFromFileName(recordName);
frameRate = params.FR;
for ch =1:nOfChannels
    rawSpeckleContrast_jumpsCorrected{ch}  = CorrectJumps(rawSpeckleContrast{ch});
    corrSpeckleContrast_jumpsCorrected{ch} = CorrectJumps(corrSpeckleContrast{ch});   
end

if exist('rawSpeckleContrast_jumpsCorrected','var') 
    if ~exist('recSavePrefix','var')
        if exist(recordName,'file') == 7 % it's a folder
            recSavePrefix = [ recordName filesep ];
        else % it's a file
            recSavePrefix = [ recordName(1:find(recordName=='.',1,'last')-1) '_' ];
        end
    end
    stdStr = sprintf('Std%dx%d',windowSize,windowSize);
    save([recSavePrefix 'Local' stdStr '_corr.mat'],'timeVec', 'frameRate','corrSpeckleContrast' , 'rawSpeckleContrast','meanVec', 'info', 'recordName','windowSize','rawSpeckleContrast_jumpsCorrected','corrSpeckleContrast_jumpsCorrected');    

    for ch = 1:nOfChannels
        [raw_SNR2,  raw_FFT2 , raw_freq2, raw_pulseFreq2, raw_pulseBPM2] = CalcSNR_Pulse(rawSpeckleContrast_jumpsCorrected{ch},frameRate,false);
        [corr_SNR2,  corr_FFT2 , corr_freq2, corr_pulseFreq2, corr_pulseBPM2] = CalcSNR_Pulse(corrSpeckleContrast_jumpsCorrected{ch},frameRate,false);

        fig_jump_corrected = figure('Units','Normalized','Position',[0.01 0.4 0.95 0.3]);
        Nx=3; Ny=2;
        
        subplot(Ny,Nx,1);
        plot(timeVec,rawSpeckleContrast{ch})
        ylabel('Kraw^2')
        title(sprintf('Channel %d - Raw (<I>=%.0fDU) ',ch,mean(meanVec{ch})));
        xlim([0 timeVec(end)]);
        xlabel('Time [s]');

        subplot(Ny,Nx,2);
        plot(timeVec,rawSpeckleContrast_jumpsCorrected{ch})
        ylabel('Kraw^2')
        title(sprintf('Channel %d - Raw (<I>=%.0fDU) - Jumps corrected',ch,mean(meanVec{ch})));
        xlim([0 timeVec(end)]);
        xlabel('Time [s]');

        subplot(Ny,Nx,3);
        plot(raw_freq2,raw_FFT2)
        ylabel(' FFT')
        title(sprintf('FFT - jumps corrected: SNR=%.2g Pulse=%.0fbpm',raw_SNR2,raw_pulseBPM2));
        xlim([0 raw_freq2(end)]);
        xlabel('Frequency [Hz]')

        subplot(Ny,Nx,4);
        plot(timeVec,corrSpeckleContrast{ch})
        ylabel('Kraw^2')
        title(sprintf('Channel %d - Corr (<I>=%.0fDU) ',ch,mean(meanVec{ch})));
        xlim([0 timeVec(end)]);
        xlabel('Time [s]');
        
        subplot(Ny,Nx,5);
        plot(timeVec,corrSpeckleContrast_jumpsCorrected{ch})
        ylabel('Kraw^2')
        title(sprintf('Channel %d - Corr (<I>=%.0fDU) - Jumps corrected',ch,mean(meanVec{ch})));
        xlim([0 timeVec(end)]);
        xlabel('Time [s]');

        subplot(Ny,Nx,6);
        plot(corr_freq2,corr_FFT2)
        ylabel(' FFT')
        title(sprintf('FFT Corr - jumps corrected: SNR=%.2g Pulse=%.0fbpm',corr_SNR2,corr_pulseBPM2));
        xlim([0 corr_freq2(end)]);
        xlabel('Frequency [Hz]')

        savefig(fig_jump_corrected,[recSavePrefix 'Local' stdStr '_plot_jumps_corrected.fig']);
    end      
end

%% Plot Jump Corrected
%% Plot
Nx = 3;   
if plotFlag
    fig5 = figure('name',['SCOS ' shortRecName ' Mono' num2str(nOfBits) ],'Units','Normalized','Position',[0.01,1-0.16-nOfChannels*0.15,0.9,0.05+nOfChannels*0.15]);
    for Ch = 1:nOfChannels
        
        try
            [raw_SNR,  raw_FFT , raw_freq, raw_pulseFreq, raw_pulseBPM] = CalcSNR_Pulse(rawSpeckleContrast_jumpsCorrected{Ch},frameRate,false);
        catch err
            warning(err.message);
        end
        subplot(nOfChannels,Nx,Nx*(Ch-1)+1);
            plot(timeVec,meanVec{Ch});
            title(sprintf('Channel %d - mean I (<I>=%.0fDU  %.1%%)',Ch,mean(meanVec{Ch}),mean(meanVec{Ch})/2^nOfBits*100));
            xlim([0 timeVec(end)]);
            xlabel('Time [s]')
        subplot(nOfChannels,Nx,Nx*(Ch-1)+2);
            plot(timeVec,rawSpeckleContrast_jumpsCorrected{Ch})
            ylabel('Kraw^2')
            title(sprintf('Channel %d - Raw (<I>=%.0fDU  %.1f%%)',Ch,mean(meanVec{Ch})));
            xlim([0 timeVec(end)]);
            xlabel('Time [s]')
        subplot(nOfChannels,Nx,Nx*(Ch-1)+3);
            plot(raw_freq,raw_FFT)
            ylabel(' FFT ')
            title(sprintf('FFT: SNR=%.2g Pulse=%.0fbpm',raw_SNR,raw_pulseBPM));
            xlim([0 raw_freq(end)]);
            xlabel('Frequency [Hz]')
    end
    savefig(fig5,[recSavePrefix 'Contrast' stdStr '_jumpCorrected_plot_Raw.fig']);

    fig6 = figure('name',['SCOS ' shortRecName ' Mono' num2str(nOfBits)],'Units','Normalized','Position',[0.01,1-0.16-nOfChannels*0.15,0.9,0.05+nOfChannels*0.15]);
    for Ch = 1:nOfChannels        
        try
            [corr_SNR,  corr_FFT , corr_freq, corr_pulseFreq, corr_pulseBPM] = CalcSNR_Pulse(corrSpeckleContrast_jumpsCorrected{Ch},frameRate,false);
        catch err
            warning(err.message);
        end
        subplot(nOfChannels,Nx,Nx*(Ch-1)+1);
            plot(timeVec,meanVec{Ch});
            title(sprintf('Channel %d - mean I (<I>=%.0fDU)',Ch,mean(meanVec{Ch})));
            xlim([0 timeVec(end)]);
            xlabel('Time [s]')
        subplot(nOfChannels,Nx,Nx*(Ch-1)+2);
            plot(timeVec,corrSpeckleContrast_jumpsCorrected{Ch})
            ylabel('Kraw^2')
            title(sprintf('Channel %d - Corr (<I>=%.0fDU)',Ch,mean(meanVec{Ch})));
            xlim([0 timeVec(end)]);
            xlabel('Time [s]')
        subplot(nOfChannels,Nx,Nx*(Ch-1)+3);
            plot(corr_freq,corr_FFT)
            ylabel(' FFT ')
            title(sprintf('FFT: SNR=%.2g Pulse=%.0fbpm',corr_SNR,corr_pulseBPM));
            xlim([0 corr_freq(end)]);
            xlabel('Frequency [Hz]')
    end
    savefig(fig6,[recSavePrefix 'Contrast' stdStr 'jumpCorrected_plot_noiseSabtracted6.fig']);
            %     plot(timeVec,corrSpeckleContrast{k})
        %     ylabel('Kf^2')
        %     title(sprintf('Channel %d - Corrected partly ',k));
end