%% Load Data
% file = 'C:\Users\tarlevi\OneDrive - Bar-Ilan University - Students\PhD Research\SCOS\Records\Nisan Head\2024_07_22\Main_Mono12_expT15ms_Gain24_FR20Hz\LocalStd7x7_corr.mat';
% file ='C:\Users\tarlevi\OneDrive - Bar-Ilan University - Students\PhD Research\SCOS\Records\fromTomoya\FR20Hz_Gain16dB_expT10ms_BL400_right-20241009T063140Z-001\LocalStd7x7_corr.mat';
if exist('.\lastRec.mat','file')
    lastF = load('.\lastRec.mat');
else
    lastF.recName = [ fileparts(pwd) '\Records' ];
end

[ filename, foldername ] = uigetfile([lastF.recName '\Local*x*.mat']);%'..\Records\VikaHead\Local*x*.mat');
file = fullfile(foldername,filename);
[folderName,rawName]=fileparts(fileparts(file));
[~,subfolderName, ext]=fileparts(folderName);
subfolderName = [subfolderName ext];
titleStr = [ subfolderName ]; % ' ' strrep(rawName,'_',' ') ];
frameRate = ExtractParametersFromString(rawName,'FR'); 
load(file)


%% Calc rBFi
BFi = 1./corrSpeckleContrast{1};
if timeVec(end) > 120
    timeToPlot = timeVec / 60; % convert to min
    xLabelStr = 'time [min]';
        rBFi = BFi/mean(BFi([1:round(10*frameRate)])); % normalize by first 10 seconds
%     rBFi = BFi/mean(BFi(60*frameRate*7+ [1:round(10*frameRate)])); % normalize by first 10 seconds
else
    timeToPlot = timeVec ; % convert to min
    xLabelStr = 'time [sec]';
    rBFi = BFi/prctile(BFi(1:round(10*frameRate)),5); % normalize by 5% percentile in first 10 sec
end
   
%% Timing
if foldername(end)=='\'; foldername(end)=[]; end
timingDir = dir([ fileparts(foldername) '\*timing*.txt']);
if ~isempty(timingDir)

    timingFile = fullfile(fileparts(foldername),timingDir(1).name);
    T = readtable(timingFile);
    baselineTime = T.Time{1};
    if nnz(baselineTime==':') == 1
       baselineTime = [ baselineTime , ':00'] ;
    end

    if exist(timingFile,'file')
        if exist([foldername '\StartTime.mat'],'file')
            s = load([foldername '\StartTime.mat']);
            startTimeScos = s.startTime;
        else
            startTimeScos = getSCOSstartTime(foldername);
        end
    end
    baselineToScosStart = minutes(duration(startTimeScos) - duration(baselineTime));
else
    baselineToScosStart = 0;
end

%% Plot
fig8 = figure('Name',['rBFi: '  rawName ],'Units','Normalized','Position',[0.1,0.1,0.72 ,0.55]);
subplot(2,1,1);
plot(timeToPlot+baselineToScosStart,rBFi);
ylim([0 min(10,max(rBFi))])
title(titleStr,'interpreter','none')
xlabel(xLabelStr)
ylabel('rBFi');
grid on
hold on;
set(gca,'FontSize',10);

subplot(2,1,2)
plot(timeToPlot+baselineToScosStart,meanVec{1});
xlabel(xLabelStr)
ylabel('<I> [DU]');
set(gca,'FontSize',10);
grid on
% tLim=10; subplot(2,1,1); xlim([0 tLim]); subplot(2,1,2); xlim([0 tLim])

% Timing
for plt_i = 1:2
    subplot(2,1,plt_i);
    markTiming(timingFile);
    xlim([0 timeToPlot(end)+3])
end
%% Save
savefig(fig8,[fileparts(file) '\_rBFi.fig']);
save([fileparts(file) '\_rBFi.mat'],'timeVec','rBFi','timeToPlot','BFi');

save('.\lastRec.mat','recName');