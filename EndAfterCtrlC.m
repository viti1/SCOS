if isvalid(vid);
    stop(vid);
    delete(vid);
    rawSpeckleContrast{ch}(k:end) = [];
    corrSpeckleContrast{ch}(k:end) = [];
    meanVec{ch}(k:end) = [];
end
msgbox("Turn off Laser");

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
saveName = [recName '\Local' stdStr '_corr_realtime.mat'];
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

    savefig(fig,[recName '\Local' stdStr '_plot_realtime.fig']);


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
ylim([0 10])
grid on
hold on;
set(gca,'FontSize',10);
subplot(2,1,2)
plot(timeToPlot,meanVec{1});
xlabel(xLabelStr)
ylabel('<I> [DU]');
set(gca,'FontSize',10);
grid on

timingFileDir = dir([ fileparts(recName),'\*timing*.txt' ]);
if numel(timingFileDir)==1  
    markTiming([fileparts(recName) '\' timingFileDir.name])
elseif numel(timingFileDir)>1
    warning('More that one timing file was found: ')
    disp({timingFileDir.name}')
end

% tLim=10; subplot(2,1,1); xlim([0 tLim]); subplot(2,1,2); xlim([0 tLim])
savefig(fig8,[recName '\_rBFi_realtime.fig']);
savefig(fig8,[fileparts(recName) '\rBFi_realtime.fig']);
