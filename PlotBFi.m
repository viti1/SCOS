%% Load Data
% file = 'C:\Users\tarlevi\OneDrive - Bar-Ilan University - Students\PhD Research\SCOS\Records\Nisan Head\2024_07_22\Main_Mono12_expT15ms_Gain24_FR20Hz\LocalStd7x7_corr.mat';
% file ='C:\Users\tarlevi\OneDrive - Bar-Ilan University - Students\PhD Research\SCOS\Records\fromTomoya\FR20Hz_Gain16dB_expT10ms_BL400_right-20241009T063140Z-001\LocalStd7x7_corr.mat';
[ filename, foldername ] = uigetfile('..\Records\VikaHead\Local*x*.mat');
file = fullfile(foldername,filename);
[~,rawName]=fileparts(fileparts(file));
titleStr = strrep(rawName,'_',' ');
frameRate = ExtractParametersFromString(rawName,'FR'); 
load(file)

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
fig8 = figure('Name',['rBFi: '  recordName ],'Units','Normalized','Position',[0.1,0.1,0.4,0.4]);
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
savefig(fig8,[fileparts(file) '\_rBFi.fig']);
save([fileparts(file) '\_rBFi.mat'],'timeVec','rBFi','BFi');