%  ---------------------------------------------------------------------------------------------------------
%%  [ timeVec,  , rawSpeckleContrast , rawSpeckleVar, corrSpeckleVar , corrSpeckleContrast, imMeanVec , info] = PlotSCOSvsTime(recordName,windowSize,plotFlag,maskInput)
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
    SCOSvsTimeMultiChannel(recordName,windowSize,plotFlag,resetMask)

addpath('.\baseFunc');
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
    
%     maxWindowSize = 50; minWindowSize = 3;
%     answer = inputdlg('Window Size','',[1 25],{'9'});
%     windowSize = str2double(answer{1});
%     if isnan(windowSize) || windowSize > maxWindowSize || windowSize < minWindowSize
%         errordlg(['Window Size must be a number between ' num2str(minWindowSize)  ' and ' num2str( num2str(maxWindowSize) ) ]);
%         error(['Window Size must be a number between ' num2str(minWindowSize)  ' and ' num2str( num2str(maxWindowSize) ) ])
%     end
    windowSize = 9;
    save('.\lastRec.mat','recordName')
    clear answer
end

if ~exist('resetMask','var') || isempty(resetMask)
    resetMask = false;
end
%% Create Mask
upFolders = strsplit(recordName,filesep);
rawName = strrep( strjoin(upFolders(end-2:end),'; '), '_',' ');
shortRecName = strjoin(upFolders(end-2:end));

if exist(recordName,'file') == 7 % it's a folder
    recSavePrefix = [ recordName filesep ];
else % it's a file
    recSavePrefix = [ recordName(1:find(recordName=='.',1,'last')-1) '_' ];
end
maskFile = [recSavePrefix 'Mask.mat'];

% delete(maskFile)
if  isstruct(resetMask) 
    masks{1} = resetMask.mask;
    totMask = resetMask.mask;
    channels.Centers = resetMask.circ.Center;
    channels.Radii   = resetMask.circ.Radius;
elseif  ~exist(maskFile,'file') || isequal(resetMask,1) || isequal(resetMask,true)

    savefig(figIm,[recordName '\ImWithChannels.fig']);
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

if ~isfield(info.name,'FR') || isnan(info.name.FR)
    error('Frame Rate must be part of the recording name as "FR"');
end

if ~isfield(info,'cameraSN') || isempty(info.cameraSN)
    if ~isequal(info.fileType,'.tiff')
        error('for this function file type should be .tiff');
    end
    error('Camera serial number was not found from the file name');
end
nOfBits = info.nBits;

%% Set params
% find read noise for current gain
nOfBits = info.nBits;
camDataFolder = [fileparts(mfilename('fullpath')) '\camerasData'];
camDataFile = dir([ camDataFolder '\SN' num2str(info.cameraSN) '*readNoiseVsGain_Mono' num2str(nOfBits) '.mat']);
if numel(camDataFile) == 0
    error('Camera Data file was not found');
elseif numel(camDataFile) > 1
    disp(camDataFile.name)
    error('more that one camera file was found');
end

dData = load([camDataFolder filesep camDataFile]); 
maxCapacity = 10.5e3;% [e]
actualGain = ConvertGain(info.name.Gain,nOfBits,maxCapacity);
readoutN   = interp1(ConvertGain(dData.gainArr,nOfBits,maxCapacity),dData.tempNoise, actualGain ,'spline');


%%  Calc Specle Contrast
disp(['Calculation SCOS on "' recordName '" ... '])
nOfFrames = GetNumOfFrames(recordName);
nOfChannels = numel(masks);

[ rawSpeckleContrast , corrSpeckleContrast , meanVec ] = InitNaN([nOfFrames 1],nOfChannels);
bt=tic;
at=tic;
for i=1:nOfFrames
%     fprintf('%d ',i);
    if mod(i,50)==0; fprintf('%d\t\t',i); disp(toc(at)); at=tic; end
    im = ReadRecord(recordName,1,i);
    stdIm = stdfilt(im,true(windowSize));
    meanIm = imfilter(im, true(windowSize)/windowSize^2,'conv','same'); % TBD!!! add Xiaojun algorithm
    meanImSquare = meanIm.^2;
    Kraw = (stdIm.^2)./meanImSquare ;
    for Ch = 1:nOfChannels
        rawSpeckleContrast{Ch}(i) = mean(Kraw(masks{Ch}));
        corrSpeckleContrast{Ch}(i) = mean(Kraw(masks{Ch}) - actualGain./meanIm(masks{Ch}) - (1/12 + readoutN^2)./meanImSquare(masks{Ch}) ); % readoutN^2
        meanVec{Ch}(i) = mean(im(masks{Ch}));
    end
end
fprintf('\n');
toc(bt)
%% Create Time vector

frameRate = info.name.FR; 

timeVec = (0:(nOfFrames-1))'*(1/frameRate) ;   % FR = FrameRate
p2p_time = timeVec<timePeriodForP2P;
stdStr = sprintf('Std%dx%d',windowSize,windowSize);

% infoFields = fieldnames(info.name);
% titleStr =  [ infoFields{1} '; exp=' num2str(info.name.expT)  'ms; Gain='  num2str(info.name.Gain) 'dB' ];

%% Save Data
save([recSavePrefix 'Local' stdStr '.mat'],'timeVec', 'shortRecName','corrSpeckleContrast' , 'rawSpeckleContrast','meanVec', 'info', 'recordName','windowSize');

%% Plot
Nx = 3;   
if plotFlag
    fig1 = figure('name',['SCOS ' shortRecName],'Units','Normalized','Position',[0.01,1-0.16-nOfChannels*0.15,0.9,0.05+nOfChannels*0.15]);
    for Ch = 1:nOfChannels
        
        try
            [raw_SNR,  raw_FFT , raw_freq, raw_pulseFreq, raw_pulseBPM] = CalcSNR_Pulse(rawSpeckleContrast{Ch},frameRate,false);
        catch err
            warning(err.message);
        end
        subplot(nOfChannels,Nx,Nx*(Ch-1)+1);
            plot(timeVec,meanVec{Ch});
            title(sprintf('Channel %d - mean I (<I>=%.0fDU)',Ch,mean(meanVec{Ch})));
            xlim([0 timeVec(end)]);
            xlabel('Time [s]')
        subplot(nOfChannels,Nx,Nx*(Ch-1)+2);
            plot(timeVec,rawSpeckleContrast{Ch})
            ylabel('Kraw^2')
            title(sprintf('Channel %d - Raw (<I>=%.0fDU)',Ch,mean(meanVec{Ch})));
            xlim([0 timeVec(end)]);
            xlabel('Time [s]')
        subplot(nOfChannels,Nx,Nx*(Ch-1)+3);
            plot(raw_freq,raw_FFT)
            ylabel(' FFT ')
            title(sprintf('FFT: SNR=%.2g Pulse=%.0fbpm',raw_SNR,raw_pulseBPM));
            xlim([0 raw_freq(end)]);
            xlabel('Frequency [Hz]')
    end
    savefig(fig1,[recSavePrefix 'Local' stdStr '_plot.fig']);

    fig1_corrected = figure('name',['SCOS ' shortRecName],'Units','Normalized','Position',[0.01,1-0.16-nOfChannels*0.15,0.9,0.05+nOfChannels*0.15]);
    for Ch = 1:nOfChannels        
        try
            [corr_SNR,  corr_FFT , corr_freq, corr_pulseFreq, corr_pulseBPM] = CalcSNR_Pulse(corrSpeckleContrast{Ch},frameRate,false);
        catch err
            warning(err.message);
        end
        subplot(nOfChannels,Nx,Nx*(Ch-1)+1);
            plot(timeVec,meanVec{Ch});
            title(sprintf('Channel %d - mean I (<I>=%.0fDU)',Ch,mean(meanVec{Ch})));
            xlim([0 timeVec(end)]);
            xlabel('Time [s]')
        subplot(nOfChannels,Nx,Nx*(Ch-1)+2);
            plot(timeVec,corrSpeckleContrast{Ch})
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
    savefig(fig1,[recSavePrefix 'Local' stdStr '_plot.fig']);
    savefig(fig1_corrected,[recSavePrefix 'Local' stdStr '_plot_corrected.fig']);

            %     plot(timeVec,corrSpeckleContrast{k})
        %     ylabel('Kf^2')
        %     title(sprintf('Channel %d - Corrected partly ',k));
end

% for plot_i=1:nOfChannels*Nx
%     subplot(nOfChannels,Nx,plot_i);    
% end
%% Correct Jumps
nOfChannels = numel(corrSpeckleContrast);
for channel_i =1:nOfChannels
    
%     jump_th = 0.065; %std(rawSpeckleContrast{channel_i}(1:min(20,end))) * 5;
    diff_rawSpeckleContrast =  diff(rawSpeckleContrast{channel_i});
    jump_th = median(abs(diff_rawSpeckleContrast)) * 10;
    jump_up_idx = find(diff(rawSpeckleContrast{channel_i}) > jump_th);
    jump_down_idx = find(diff(rawSpeckleContrast{channel_i}) < -jump_th);
    if ~isempty(jump_up_idx) && ~isempty(jump_down_idx)
        jumps_up   = diff_rawSpeckleContrast(jump_up_idx);
        jumps_down = diff_rawSpeckleContrast(jump_down_idx);

        rawSpeckleContrast_jumpsCorrected{channel_i} = rawSpeckleContrast{channel_i};
        corrSpeckleContrast_jumpsCorrected{channel_i} = rawSpeckleContrast{channel_i};
        if jump_up_idx(1) > jump_down_idx(1)
            if numel(jump_up_idx) < numel(jump_down_idx)
                jump_up_idx( end + 1 ) = numel(rawSpeckleContrast{channel_i}); %#ok<AGROW>
                jump_down_idx(numel(jump_up_idx)+1:end) = [];
            end
            for n=1:numel(jump_down_idx)
                rawSpeckleContrast_jumpsCorrected{channel_i}(jump_down_idx(n)+1:jump_up_idx(n)) = rawSpeckleContrast{channel_i}(jump_down_idx(n)+1:jump_up_idx(n) ) - jumps_down(n);
                corrSpeckleContrast_jumpsCorrected{channel_i}(jump_down_idx(n)+1:jump_up_idx(n)) = corrSpeckleContrast{channel_i}(jump_down_idx(n)+1:jump_up_idx(n) ) - jumps_down(n);                
            end
        else
            if numel(jump_down_idx) < numel(jump_up_idx)
                jump_down_idx( end + 1 ) = numel(rawSpeckleContrast{channel_i}); %#ok<AGROW>
                jump_up_idx(numel(jump_down_idx)+1:end) = [];
            end
            for n=1:numel(jump_up_idx)
                rawSpeckleContrast_jumpsCorrected{channel_i}(jump_up_idx(n)+1:jump_down_idx(n))  = rawSpeckleContrast{channel_i}(jump_up_idx(n)+1:jump_down_idx(n) ) - jumps_up(n);
                corrSpeckleContrast_jumpsCorrected{channel_i}(jump_up_idx(n)+1:jump_down_idx(n)) = corrSpeckleContrast{channel_i}(jump_up_idx(n)+1:jump_down_idx(n) ) - jumps_up(n);
            end    
        end
    end
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
%     save([recSavePrefix 'Local' stdStr '.mat'],'timeVec', 'frameRate','corrSpeckleContrast' , 'rawSpeckleContrast','meanVec', 'info', 'recordName','windowSize','rawSpeckleContrast_jumpsCorrected','corrSpeckleContrast_jumpsCorrected');    

    for channel_i = 1:nOfChannels
        [raw_SNR2,  raw_FFT2 , raw_freq2, raw_pulseFreq2, raw_pulseBPM2] = CalcSNR_Pulse(rawSpeckleContrast_jumpsCorrected{channel_i},frameRate,false);
        [corr_SNR2,  corr_FFT2 , corr_freq2, corr_pulseFreq2, corr_pulseBPM2] = CalcSNR_Pulse(corrSpeckleContrast_jumpsCorrected{channel_i},frameRate,false);

        fig_jump_corrected = figure('Units','Normalized','Position',[0.01 0.4 0.95 0.3]);
        Nx=3; Ny=2;
        
        subplot(Ny,Nx,1);
        plot(timeVec,rawSpeckleContrast{channel_i})
        ylabel('Kraw^2')
        title(sprintf('Channel %d - Raw (<I>=%.0fDU) ',channel_i,mean(meanVec{channel_i})));
        xlim([0 timeVec(end)]);
        xlabel('Time [s]');

        subplot(Ny,Nx,2);
        plot(timeVec,rawSpeckleContrast_jumpsCorrected{channel_i})
        ylabel('Kraw^2')
        title(sprintf('Channel %d - Raw (<I>=%.0fDU) - Jumps corrected',channel_i,mean(meanVec{channel_i})));
        xlim([0 timeVec(end)]);
        xlabel('Time [s]');

        subplot(Ny,Nx,3);
        plot(raw_freq2,raw_FFT2)
        ylabel(' FFT')
        title(sprintf('FFT - jumps corrected: SNR=%.2g Pulse=%.0fbpm',raw_SNR2,raw_pulseBPM2));
        xlim([0 raw_freq2(end)]);
        xlabel('Frequency [Hz]')

        subplot(Ny,Nx,4);
        plot(timeVec,corrSpeckleContrast{channel_i})
        ylabel('Kraw^2')
        title(sprintf('Channel %d - Corr (<I>=%.0fDU) ',channel_i,mean(meanVec{channel_i})));
        xlim([0 timeVec(end)]);
        xlabel('Time [s]');
        
        subplot(Ny,Nx,5);
        plot(timeVec,corrSpeckleContrast_jumpsCorrected{channel_i})
        ylabel('Kraw^2')
        title(sprintf('Channel %d - Corr (<I>=%.0fDU) - Jumps corrected',channel_i,mean(meanVec{channel_i})));
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
end

 