clear 
addpath('.\baseFunc');
%% Get User Input
if exist('.\lastRec.mat','file')
    lastF = load('.\lastRec.mat');
else
    lastF.recName = [ fileparts(pwd) '\Records' ];
end
folder = uigetdir(fileparts(lastF.recordName),'Where to Save?');
%  folder = '..\Records\VikaHead\Basler_1920\T3_short_40Hz';
% folder = 'C:\SCOS\Records\Tests\T1_FR5Hz';
if folder==0 ; return; end
lastF.recordName = folder;
save('.\lastRec.mat','-struct','lastF')
nOfFrames = 1000;
nForSP = 600;
windowSize = 9;
frameRate = 10; % Hz
camParams.ExposureTime = 10000;
camParams.Gain = 20;
camParams.BlackLevel = 100;
nOfDarkFrames = 1000;
camParams.videoFormat = 'Mono10';
camParams.TriggerSource = 'Line2';
camParams.addToFilename.TriggerSource = false;
camParams.addToFilename.TriggerMode = false;
camParams.addToFilename.videoFormat = false;

setupParams.Laser = 'iBeam';
setupParams.LaserPower = 96; %mW
setupParams.Fiber = 'Multicore_200um_x_19';
setupParams.SDS = 3;

%% Record Dark Recording
uiwait(msgbox("Turn off Laser"));
camParams.TriggerMode = 'Off'; % for dark recording
[ darkMean, darkVarIm , darkRecName, infoDark ] = RecordFromCameraVarAndMean( nOfDarkFrames, camParams, [], folder, '.tiff','DarkIm', '', 1, 0);
darkVar = imboxfilt(darkVarIm,windowSize);
mkdir(fileparts(darkRecName));
% save( darkRecName , 'meanIm','darkVar','darkVarIm'); % fit the other standart
% clear meanIm
if mean2(darkMean) - infoDark.cam.BlackLevel > 5
    warning(['Dark Image with suspicious level ' num2str(darkMean)]);
    answer = input(['Dark Image with suspicious level ' num2str(mean2(darkMean)) 10 'Do you want to continue? (Y/N) ' ],'s');
    if ~strcmpi(answer,'y')
        return;
    end       
end
%% Create vid & src
vid = videoinput("gentl", 1, camParams.videoFormat);
vid.FramesPerTrigger = Inf; 
src = getselectedsource(vid);
triggerconfig(vid, 'hardware');
src.TriggerMode = 'On';

% Create filename from Parameters Structs
[recName, recRawName,] = GenerateFileName(folder,camParams,[],'',['_FR' num2str(frameRate) 'Hz'],'.tiff',0,src);
mkdir(recName);
%% Get Mask
uiwait(msgbox("Turn Laser On"));
start(vid);
while(~vid.FramesAvailable); pause(0.001); end
imagesBuff = getdata(vid, vid.FramesAvailable); 
im = mean(squeeze(imagesBuff),3) - camParams.BlackLevel;
stop(vid);
%%
[totMask, circ, figMask] = GetROI(im,windowSize);
masks{1} = totMask;
meanI = round(mean(im(totMask)));
pcntl = prctile(im(totMask),[5 95]);
h2 = msgbox([ '<I>=' num2str(meanI) 'DU;   5%I=' num2str(pcntl(1)) 'DU;   95%I=' num2str(pcntl(2)) 'DU  '   ])
disp([ '<I>=' num2str(meanI) 'DU;   5%I=' num2str(pcntl(1)) 'DU;   95%I=' num2str(pcntl(2)) 'DU  '   ])
channels.Centers = circ.Center;
channels.Radii = circ.Radius;
save([recName '\Mask.mat'],'masks','channels','totMask');
        
savefig(figMask,[recName '\maskIm.fig']);
close(figMask);

%% Create info struct
info.setup = setupParams;
info.cameraSN = src.DeviceSerialNumber;
info.nBits = str2double(camParams.videoFormat(5:end));
[~,recShortName] = fileparts(recName);
info.name = GetParamsFromFileName(recShortName);
actualGain = GetActualGain(info)
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
[y,x] = find(totMask) ;
roi_lims  = [ min(y)-windowSize  , max(y)+windowSize
              min(x)-windowSize  , max(x)+windowSize ];

imageSize = size(totMask);
roi_lims( roi_lims < 1 ) = 1;        
if roi_lims(1,2) > imageSize(1)
    roi_lims(1,2) = imageSize(1);
end
if roi_lims(2,2) > imageSize(2)
    roi_lims(2,2) = imageSize(2);
end
roi.y = roi_lims(1,1):roi_lims(1,2);
roi.x = roi_lims(2,1):roi_lims(2,2);

masks_cut = cell(size(masks));
for ch = 1:numel(masks)    
    masks_cut{ch} = masks{ch}(roi.y  , roi.x); 
end
   
%% Get images Sequence from Camera
fprintf('Recording "%s" ... \n',recName);
if ~isvalid(vid) 
    vid = videoinput("gentl", 1, camParams.videoFormat);
    vid.FramesPerTrigger = Inf; 
    src = getselectedsource(vid);
    triggerconfig(vid, 'hardware');
    src.TriggerMode = 'On';
end
start(vid);
ch =1;

% preallocate
if isinf(nOfFrames)
    nAlloc = 10000;
else
    nAlloc = nOfFrames;
end

rawSpeckleContrast{ch}  = nan(1,nAlloc);
corrSpeckleContrast{ch} = nan(1,nAlloc);
meanVec{ch} = nan(1,nAlloc);
    
% Start aquisition
k=1;
h_waitbar = waitbar(0,'Recording ...');
ch = 1; imFig = [];
spRec = uint16(nan(size(im,1),size(im,2),nForSP));
while  k<=nOfFrames
    while ~vid.FramesAvailable; pause(0.005); end
    im_raw = squeeze(getdata(vid, 1));
    
    %calc SCOS
    if k > nForSP
        im = double(im_raw) - darkMean;
        im_cut = im(roi.y,roi.x);
        stdIm = stdfilt(im_cut,true(windowSize));

        meanFrame = mean(im_cut(masks_cut{ch}));
        fittedI = imboxfilt(im_cut,windowSize);
        fittedISquare = fittedI.^2;
    
        if k > length(rawSpeckleContrast{ch}) 
            rawSpeckleContrast{ch}  = [ rawSpeckleContrast{ch} nan(1,nAlloc) ];
            corrSpeckleContrast{ch} = [ rawSpeckleContrast{ch} nan(1,nAlloc) ];
            meanVec{ch} = [ meanVec{ch} nan(1,nAlloc) ];
        end

        rawSpeckleContrast{ch}(k) = mean((stdIm(masks_cut{ch}).^2 ./ fittedISquare(masks_cut{ch})));
        corrSpeckleContrast{ch}(k) = mean( ( stdIm(masks_cut{ch}).^2 - actualGain.*fittedI(masks_cut{ch})  - spVar(masks_cut{ch}) - 1/12 - darkVar(masks_cut{ch}))./fittedISquare(masks_cut{ch}) ); % - ( readoutN^2 )./fittedISquare(masks{ch}) );
        meanVec{ch}(k) = meanFrame;
    else
        spRec(:,:,k) = im_raw;
    end
     
    % calc sp noise and fitting coefficients
    if k == nForSP
        stop(vid);
        pauseStart = clock();
        
        disp('Calc SP');        
        spIm = mean(spRec,3) - darkMean;
        fig_spIm = my_imagesc(spIm); title(['Image average ' num2str(size(spRec,3)) ' frames'] );
        imAx = gca;
        imgH = findobj(imAx, 'Type', 'image');
        savefig(fig_spIm, [recName '\spIm.fig']);
        spVar = stdfilt( spIm ,true(windowSize)).^2;
        fitI_A = ones(size(totMask));
        fitI_B = zeros(size(totMask));
        fitI_A_cut = ones(size(masks_cut{ch}));
        fitI_B_cut = zeros(size(masks_cut{ch}));
        save([recName  '\smoothingCoefficients.mat'],'spVar','fitI_A','fitI_B','spIm','totMask');
    
        pauseEnd = clock();
        start(vid);
    end
    
    if k > nForSP && mod(k,100) == 0
        fprintf('<I>=%.3gDU , K_raw = %.5g , Ks=%.5g , Kr=%.5g, Ksp=%.5g, Kq=%.5g, Kf=%.5g\n',meanFrame,rawSpeckleContrast{ch}(k), ...
            mean(actualGain.*fittedI(masks_cut{ch})./fittedISquare(masks_cut{ch})),mean(darkVar(masks_cut{ch})./fittedISquare(masks_cut{ch})),...
            mean(spVar(masks_cut{ch})./fittedISquare(masks_cut{ch})),mean(1./(12*fittedISquare(masks_cut{ch}))),corrSpeckleContrast{ch}(k));
    end
    
    % write tiff
    t = Tiff([recName,sprintf('\\frame%0*d.tiff',3,k)],'w');
    setTag(t,tagstruct);
    write(t,uint16(im_raw));
    close(t);
    
    % print & increment k
    if mod(k,100)==0 
        if isvalid(h2); close(h2); end
        fprintf('frame %d : %d frames in buffer \n',k,vid.FramesAvailable);
        im = double(im_raw) - darkMean;

        if k==100
            [imFig , imgH ] = my_imagesc(im); 
        else            
            figure(imFig);
            imgH.CData = im;        
        end
        pcntl = prctile(im(totMask),[5 95]);
        title(['frame ' num2str(k) ': <I>=' num2str(mean2(im)) 'DU;   5%I=' num2str(pcntl(1)) 'DU;   95%I=' num2str(pcntl(2)) 'DU  '   ])         
    end
    if mod(k,1000) == 0; fprintf('\n'); end
    k = k + 1;
end
fprintf('Finished\n');

stop(vid);
delete(vid);
rawSpeckleContrast{ch}(k:end) = [];
corrSpeckleContrast{ch}(k:end) = [];
meanVec{ch}(k:end) = [];

%% Calc First nForSP frames   
h1 = msgbox(['Calculating first ' num2str(nForSP) ' frames...']);
for k = 1:nForSP
    im = double(spRec(:,:,k)) - darkMean ;
    im_cut = im(roi.y,roi.x);
    stdIm = stdfilt(im_cut,true(windowSize));
    
    meanFrame = mean(im_cut(masks_cut{ch}));
    fittedI = imboxfilt(im_cut,windowSize);
    fittedISquare = fittedI.^2;
    
    rawSpeckleContrast{ch}(k) = mean((stdIm(masks_cut{ch}).^2 ./ fittedISquare(masks_cut{ch})));
    corrSpeckleContrast{ch}(k) = mean( ( stdIm(masks_cut{ch}).^2 - actualGain.*fittedI(masks_cut{ch})  - spVar(masks_cut{ch}) - 1/12 - darkVar(masks_cut{ch}))./fittedISquare(masks_cut{ch}) ); % - ( readoutN^2 )./fittedISquare(masks{ch}) );
    meanVec{ch}(k) = meanFrame;
    if mod(k,100)==0; fprintf('%d  ',k); end
end
fprintf('\n');
if isvalid(h1); close(h1); end
stdStr = sprintf('Std%dx%d',windowSize,windowSize);


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

nOfFrames = numel(rawSpeckleContrast{1});
timeVec = (0:(nOfFrames-1))'*(1/frameRate) ;   % FR = FrameRate

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
    if timeVec(end) > 120
        timeToPlot = timeVec / 60; % convert to min
        xLabelStr = 'time [min]';
    else
        timeToPlot = timeVec ; % convert to min
        xLabelStr = 'time [sec]';
    end

    fig7 = figure('Name',['rBFi: '  recName ],'Units','Normalized','Position',[0.1,0.1,0.4,0.4]); 
    subplot(2,1,1);
    BFi = 1./corrSpeckleContrast{1};
    % rBFi = BFi/prctile(BFi(1:round(10*frameRate)),5); % normalize by 5% percentile in first 10 sec
    rBFi = BFi/mean(BFi(1:round(1*frameRate))); % normalize by first second
    plot(timeToPlot,rBFi); 
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
    savefig(fig7,[recName '\_rBFi.fig']);


