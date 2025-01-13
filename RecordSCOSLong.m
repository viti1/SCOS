clear 
clc
addpath('.\baseFunc');
%% Get User Input
skip_dark_frames = 0;
if exist('.\lastRec.mat','file')
    lastF = load('.\lastRec.mat');
else
    lastF.recName = [ fileparts(pwd) '\Records' ];
end
% folder = 'C:\SCOS\Records\Tests\TestFromMatlab_onHead_Hospital';
folder = uigetdir('C:\Users\viti_\OneDrive - Bar Ilan University\SCOS_Records\ShaareiZedek\12.01.2025','Where to Save?');
if folder==0 ; return; end
if ~exist(folder,'dir'); mkdir(folder); end
lastF.recordName = folder;
save('.\lastRec.mat','-struct','lastF')
nOfFrames = Inf;
nOfDarkFrames = 600;
nForSP = 600;
windowSize = 9;
frameRate = 1/100e-3; % Hz 
camParams.ExposureTime = 8000;
camParams.Gain = 16;
camParams.videoFormat = 'Mono10';
camParams.BlackLevel = 100;
camParams.TriggerSource = 'Line2';  % Hirose - Line3 or Line1, M8 - Line2
camParams.addToFilename.TriggerSource = false;
camParams.addToFilename.TriggerMode = false;
camParams.addToFilename.videoFormat = false;

setupParams.Laser = 'iBeam';
setupParams.LaserPower = 120; %mW
setupParams.Fiber = '90 deg 400um';
setupParams.SDS = 3;

saveTiff_flag = true;
showEveryNframes = 50;
%% Record Dark Recording
if ~skip_dark_frames
    uiwait(msgbox("Turn off Laser"));
    camParams.TriggerMode = 'Off'; % for dark recording
    warning('off','imaq:gentl:hardwareTriggerTriggerModeOff');
    [ darkImFull, darkVarImFull , darkRecName, infoDark ] = RecordFromCameraVarAndMean( nOfDarkFrames, camParams, [], folder, '.tiff','DarkIm', '', 1, 0);
    if size(darkImFull,4)~=1
        [ darkImFull, darkVarImFull , darkRecName, infoDark ] = RecordFromCameraVarAndMean( nOfDarkFrames, camParams, [], folder, '.tiff','DarkIm', '', 1, 0);
    end

    if abs(mean2(darkImFull) - infoDark.cam.BlackLevel) > 5
        warning(['Dark Image with suspicious level ' num2str(mean2(darkImFull))]);
        answer = input(['Dark Image with suspicious level ' num2str(mean2(darkImFull)) 10 'Do you want to continue? (Y/N) ' ],'s');
        if ~strcmpi(answer,'y')
            return;
        end       
    end
else
    % TBD
end
%% Create vid & src
vid = videoinput("gentl", 1, camParams.videoFormat);
vid.FramesPerTrigger = Inf; 
src = getselectedsource(vid);
triggerconfig(vid, 'hardware');
src.TriggerMode = 'On';

start(vid)
while(~vid.FramesAvailable); ; end
tic
getdata(vid, 1);
while(~vid.FramesAvailable); ; end
oneFrameTime = toc;
frameRateApproximation = round(1/oneFrameTime,1);
stop(vid);
if abs(frameRateApproximation-frameRate) > 4    
    errStr=sprintf('Please update frame rate: measured=%.4g , expected=%.4g',frameRateApproximation,frameRate);
    errordlg(errStr);
    delete(vid)
    error(errStr); %#ok<SPERR>
else
    fprintf('Approximated frame rate = %.3g\n',frameRateApproximation);
end
%%
% Create filename
%from Parameters Structs
[recName, recRawName,] = GenerateFileName(folder,camParams,[],'Main',['_FR' num2str(frameRate) 'Hz'],'.tiff',0,src);
disp(recName)
mkdir(recName);
%% Get Mask
uiwait(msgbox("Turn Laser On"));
start(vid);
while(~vid.FramesAvailable); pause(0.001); end
imagesBuff = getdata(vid, vid.FramesAvailable); 
im = mean(squeeze(imagesBuff),3) - camParams.BlackLevel;
stop(vid);
%
[totMaskFull, circ, figMask] = GetROI(im,windowSize);
meanI = round(mean(im(totMaskFull)));
pcntl = prctile(im(totMaskFull),[5 95]);
h2 = msgbox([ '<I>=' num2str(meanI,4) 'DU;   5%I=' num2str(pcntl(1),4) 'DU;   95%I=' num2str(pcntl(2),4) 'DU  '   ])
disp([ '<I>=' num2str(meanI) 'DU;   5%I=' num2str(pcntl(1),4) 'DU;   95%I=' num2str(pcntl(2),4) 'DU  '   ])
channels.Centers = circ.Center;
channels.Radii = circ.Radius;
        
savefig(figMask,[recName '\maskIm.fig']);
close(figMask);

%% Create info struct
info.setup = setupParams;
info.cameraSN = src.DeviceSerialNumber;
info.nBits = str2double(camParams.videoFormat(5:end));
[~,recShortName] = fileparts(recName);
info.name = GetParamsFromFileName(recShortName);
actualGain = GetActualGain(info);
save([recName '\info.mat'],'-struct','info');

%% Tiff Struct
tagstruct.ImageLength = size(im,1);
tagstruct.ImageWidth  = size(im,2);

if info.nBits == 8
    tagstruct.BitsPerSample = 8;
elseif ismember(info.nBits,9:16)
    tagstruct.BitsPerSample = 16;
else
    error(['Unsupported video format "' videoFormat '" for writing .tiffs'])
end
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software = 'MATLAB';
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.Compression = Tiff.Compression.None;
    
%% Decrease Image Size
[y,x] = find(totMaskFull) ;
half_win = floor(windowSize/2) + 20 ;%+2;
yLimits = [ min(y)-half_win , max(y)+half_win ];
xLimits = [ min(x)-half_win , max(x)+half_win ];

imSize = size(im);
if yLimits(1) < 1;  yLimits(1) = 1; end
if xLimits(1) < 1;  xLimits(1) = 1; end
if yLimits(2) > imSize(1); yLimits(2) = imSize(1); end
if xLimits(2) > imSize(2); xLimits(2) = imSize(2); end
          
set(vid, 'ROIPosition', [xLimits(1) yLimits(1) (diff(xLimits)+1) (diff(yLimits)+1)])
newPos = vid.ROIPosition;
xLimits = newPos(1) + [0 (newPos(3)-1)];
yLimits = newPos(2) + [0 (newPos(4)-1)];

totMask = totMaskFull(yLimits(1):yLimits(2)  , xLimits(1):xLimits(2)); 
masks{1} = totMask;
channels.Centers = circ.Center - [xLimits(1) yLimits(1)] + 1 ;
channels.Radii = circ.Radius;
save([recName '\Mask.mat'],'masks','channels','totMask');

save([recName '\ROI.mat'],'xLimits','yLimits');
darkIm = darkImFull(yLimits(1):yLimits(2)  , xLimits(1):xLimits(2)); 
darkVarIm = darkVarImFull(yLimits(1):yLimits(2)  , xLimits(1):xLimits(2)); 
darkVar = imboxfilt(darkVarIm,windowSize);
% TBD add save dark update image with dark full

fitI_A = ones(size(totMask));
fitI_B = zeros(size(totMask));
fitI_A_cut = ones(size(totMask));
fitI_B_cut = zeros(size(totMask));
%% Get images Sequence from Camera
fprintf('Recording "%s" ... \n',recName);
if ~isvalid(vid) 
    vid = videoinput("gentl", 1, camParams.videoFormat);
    set(vid, 'ROIPosition', [xLimits(1) yLimits(1) (diff(xLimits)+1) (diff(yLimits)+1)])
    vid.FramesPerTrigger = Inf; 
    src = getselectedsource(vid);
    triggerconfig(vid, 'hardware');
    src.TriggerMode = 'On';
end
start(vid);
ch =1;

% preallocate
if isinf(nOfFrames)
    nAlloc = 60000;
else
    nAlloc = nOfFrames;
end

rawSpeckleContrast{ch}  = nan(1,nAlloc);
corrSpeckleContrast{ch} = nan(1,nAlloc);
meanVec{ch} = nan(1,nAlloc);
    
% Start aquisition
k=1;
% h_waitbar = waitbar(0,'Recording ...');
ch = 1; imFig = [];
if camParams.videoFormat(end) == '8'
    spRec = uint8(nan(size(totMask,1),size(totMask,2),nForSP));
else
    spRec = uint16(nan(size(totMask,1),size(totMask,2),nForSP));
end
spSum = zeros(size(totMask));
fig_scos = figure('Name','SCOS Graph','Units','Normalized','Position',[0.31,0.2, 0.7, 0.45]); plot(0,0)
xlabel('time [min]');
ylabel('K_{corr}^2')
ax_scos = gca; 
timeVec = (0:(nAlloc-1))'*(1/frameRate)/60 ;   % FR = FrameRate

while  k<=nOfFrames
    while ~vid.FramesAvailable; pause(0.005); end
    im_raw = squeeze(getdata(vid, 1));
    
    %calc SCOS
    if k > nForSP
        im = double(im_raw) - darkIm;
        im_cut = im; %(roi.y,roi.x);
        stdIm = stdfilt(im_cut,true(windowSize));

        meanFrame = mean(im_cut(masks{ch}));
        fittedI = imboxfilt(im_cut,windowSize);
        fittedISquare = fittedI.^2;
    
        if k > length(rawSpeckleContrast{ch}) 
            rawSpeckleContrast{ch}  = [ rawSpeckleContrast{ch} nan(1,nAlloc) ];
            corrSpeckleContrast{ch} = [ rawSpeckleContrast{ch} nan(1,nAlloc) ];
            meanVec{ch} = [ meanVec{ch} nan(1,nAlloc) ];
            timeVec = [ timeVec (timeVec(end)+timeVec)]; 
        end

        rawSpeckleContrast{ch}(k) = mean((stdIm(masks{ch}).^2 ./ fittedISquare(masks{ch})));
        corrSpeckleContrast{ch}(k) = mean( ( stdIm(masks{ch}).^2 - actualGain.*fittedI(masks{ch})  - spVar(masks{ch}) - 1/12 - darkVar(masks{ch}))./fittedISquare(masks{ch}) ); % - ( readoutN^2 )./fittedISquare(masks{ch}) );
        meanVec{ch}(k) = meanFrame;
    else
        spRec(:,:,k) = im_raw;
        spSum = double(im_raw) + spSum;
    end
     
    % calc sp noise and fitting coefficients
    if k == nForSP
        stop(vid);
        pauseStart = clock();
        
        disp('Calc SP');        
        spIm = spSum/nForSP  - darkIm;
%         fig_spIm = my_imagesc(spIm); title(['Image average ' num2str(nForSP) ' frames'] );
%         imAx = gca;
%         imgH = findobj(imAx, 'Type', 'image');
%         savefig(fig_spIm, [recName '\spIm.fig']);
        spVar = stdfilt( spIm ,true(windowSize)).^2;
        save([recName  '\smoothingCoefficients.mat'],'spVar','fitI_A','fitI_B','spIm','totMask');
        clear fitI_A  fitI_B fitI_A_cut fitI_B_cut
        pauseEnd = clock();
        start(vid);
    end
    
    if k > nForSP && (k== nForSP+1 || mod(k,showEveryNframes) == 0)        
        fprintf('<I>=%.3gDU , K_raw = %.5g , Ks=%.5g , Kr=%.5g, Ksp=%.5g, Kq=%.5g, Kf=%.5g\n',meanFrame,rawSpeckleContrast{ch}(k), ...
            mean(actualGain.*fittedI(masks{ch})./fittedISquare(masks{ch})),mean(darkVar(masks{ch})./fittedISquare(masks{ch})),...
            mean(spVar(masks{ch})./fittedISquare(masks{ch})),mean(1./(12*fittedISquare(masks{ch}))),corrSpeckleContrast{ch}(k));
        figure(fig_scos);  plot(ax_scos,timeVec(1:k),corrSpeckleContrast{ch}(1:k)); 
    end
    
    % write tiff
    if saveTiff_flag
        t = Tiff([recName,sprintf('\\frame%0*d.tiff',3,k)],'w');
        setTag(t,tagstruct);
        write(t,im_raw);        
        close(t);
    end
    
    % print & increment k
    if mod(k,showEveryNframes)==0 
        if isvalid(h2); close(h2); end
        fprintf('frame %d : %d frames in buffer \n',k,vid.FramesAvailable);
        im = double(im_raw) - darkIm;

        if k==showEveryNframes
            [imFig , imgH ] = my_imagesc(im);
            imFig.Position = [0 0.2, 0.3 0.3 ]
        else            
            figure(imFig);
            imgH.CData = im;        
        end
        pcntl = prctile(im(masks{1}),[5 95]);
        title({['frame ' num2str(k) ] ,['<I>=' num2str(mean2(im),3) '  5%I=' num2str(pcntl(1),3) '  95%I=' num2str(pcntl(2),3) ''   ], ...
            [ 'K_{raw}^2=' num2str(rawSpeckleContrast{ch}(k),3)  '  K_{corr}^2=' num2str(corrSpeckleContrast{ch}(k),3 ) ]})         
    end
    if mod(k,1000) == 0; fprintf('\n'); end
    k = k + 1;
end
% close(h_waitbar)
fprintf('Finished\n');
stop(vid);
delete(vid);
rawSpeckleContrast{ch}(k:end) = [];
corrSpeckleContrast{ch}(k:end) = [];
meanVec{ch}(k:end) = [];

%% Calc First nForSP frames   
% h1 = msgbox(['Calculating first ' num2str(nForSP) ' frames...']);
disp( ['Calculating first ' num2str(nForSP) ' frames...'] );
for k = 1:nForSP
    im_cut = double(spRec(:,:,k)) - darkIm ;
    stdIm = stdfilt(im_cut,true(windowSize));
    
    meanFrame = mean(im_cut(masks{ch}));
    fittedI = imboxfilt(im_cut,windowSize);
    fittedISquare = fittedI.^2;
    
    rawSpeckleContrast{ch}(k) = mean((stdIm(masks{ch}).^2 ./ fittedISquare(masks{ch})));
    corrSpeckleContrast{ch}(k) = mean( ( stdIm(masks{ch}).^2 - actualGain.*fittedI(masks{ch})  - spVar(masks{ch}) - 1/12 - darkVar(masks{ch}))./fittedISquare(masks{ch}) ); % - ( readoutN^2 )./fittedISquare(masks{ch}) );
    meanVec{ch}(k) = meanFrame;
    if mod(k,100)==0; fprintf('%d  ',k); end
end
fprintf('\n');
% if isvalid(h1); close(h1); end
stdStr = sprintf('Std%dx%d',windowSize,windowSize);

nOfFrames = numel(rawSpeckleContrast{1});
timeVec = (0:(nOfFrames-1))'*(1/frameRate) ;   % FR = FrameRate

%% Save
saveName = [recName '\Local' stdStr '_corr.mat'];
if exist(saveName,'file'); delete(saveName); end % just for it to have the right date
nOfChannels = 1;
save(saveName,'timeVec', 'corrSpeckleContrast' , 'rawSpeckleContrast', 'meanVec', 'info','nOfChannels', 'recName','windowSize');

%% Plot
if isvalid(vid); delete(vid); end

infoFields = fieldnames(info.name);
if ~isfield(info.name,'Gain')
    info.name.Gain = '';
end
firtsParamValue = info.name.(infoFields{1});
if ~ischar(firtsParamValue)
    firtsParamValue = num2str(firtsParamValue);
end
if isfield(info.name,'SDS')
    SDSstr = ['SDS=' num2str(info.name.SDS) 'cm'];
elseif isfield(setupParams,'SDS')
    SDSstr = ['SDS=' num2str(setupParams.SDS) 'cm'];
else
    SDSstr = '';
end

titleStr =  [ infoFields{1} firtsParamValue SDSstr '; exp=' num2str(info.name.expT)  'ms; Gain='  num2str(info.name.Gain) 'dB' ];
  
% plot Contrast
[raw_SNR,  raw_FFT , raw_freq, raw_pulseFreq, raw_pulseBPM] = CalcSNR_Pulse(rawSpeckleContrast{1},frameRate);
[corr_SNR,  corr_FFT , corr_freq, corr_pulseFreq, corr_pulseBPM] = CalcSNR_Pulse(corrSpeckleContrast{1},frameRate);

    fig = figure('name',['SCOS ' recName ' Mono' num2str(info.nBits)],'Units','Normalized','Position',[0.1 0.1 0.8 0.8]);
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

    savefig(fig,[recName '\Local' stdStr '_plot.fig']);
% Plot rBFI
%% Calc rBFi
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
    
% Plot
fig8 = figure('Name',['rBFi: '  recName ],'Units','Normalized','Position',[0.1,0.1,0.4,0.4]);
subplot(2,1,1);
plot(timeToPlot,rBFi);
title(titleStr,'interpreter','none')
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
savefig(fig8,[recName '\_rBFi.fig']);
savefig(fig8,[fileparts(recName) '\_rBFi.fig']);

timingFile = fullfile(fileparts(recName),'timing.txt'); 
if exist(timingFile,'file')    
    markTiming(timingFile)
end