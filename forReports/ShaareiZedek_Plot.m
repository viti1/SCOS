clear

recName = ['C:\Users\' getenv('username') '\OneDrive - Bar-Ilan University - Students\PhD Research\SCOS\Records\ShaareiZedek__\C4'];
% recName = 'C:\Users\tarlevi\OneDrive - Bar Ilan University\SCOS_Records\ShaareiZedek\B3 21.07.2024';
[recordsFolder , recRawName ]= fileparts(recName);
participantID = recRawName(1:2);

%% Timing
timingFileDir = dir([recName '\' participantID '*timing info.txt']); 
if isempty(timingFileDir)
    timingFileDir = dir([recordsFolder '\NIRS Data\' participantID '*timing info.txt']); 
    if isempty(timingFileDir)
        error('Could not find Timing file');
    end
    timingFile = fullfile(recordsFolder,'NIRS Data',timingFileDir.name);
else
    timingFile = fullfile(recName,timingFileDir.name);
end
Timing = readtable(timingFile);
baselineTime = Timing.Time(1);
%% NIRS
NIRSFolder = [ recordsFolder '\NIRS Data'];
nirsFile   = dir([NIRSFolder '\' participantID '*.csv']);
if ~isempty(nirsFile)    
    nirsFile = fullfile(NIRSFolder,nirsFile.name);
    NIRS = readtable(nirsFile);
    NIRS = renamevars(NIRS,2:7,{'Date','Time','Left','Right','SedlineEvents','Remarks'});
    if ischar(NIRS.Left)
        NIRS.Left = str2double(NIRS.Left);
        NIRS.Right = str2double(NIRS.Right);
    end
end
%% SCOS
SCOSdir = dir([recName '\Main*']);
if isempty(SCOSdir)
    SCOS_mat_dir = dir([recName , '\LocalStd*_corr.mat']);
    if isempty(SCOS_mat_dir)
        error(['Could not find ' recName , '\LocalStd*_corr.mat file']);
    end
    SCOSfile = fullfile(recName, SCOS_mat_dir(1).name );
    SCOS = load(SCOSfile);
    if ~exist('startTimeSCOS','var')
        tmp = strsplit(SCOS.startDateTime);
        startTimeSCOS = tmp{2};
    end
else
    SCOS_mat_dir = dir([recName filesep SCOSdir(1).name , '\LocalStd*_corr.mat']);
    SCOSfile = fullfile(recName,SCOSdir(1).name, SCOS_mat_dir(1).name );
    firstFrameDir = dir([recName,filesep,SCOSdir(1).name '\*_0001.tiff']);
    tmp = strsplit(firstFrameDir.date);
    startTimeSCOS = tmp{2};
    SCOS = load(SCOSfile);
end

%% Plot 
if isempty(nirsFile); nPlots = 2; else; nPlots =3; end
fig1 = figure('Units','Normalized','Position',[0.05 0.1 0.9 0.8]);
BFi = 1./SCOS.corrSpeckleContrast{1};
BFi(SCOS.corrSpeckleContrast{1} < 0 ) = NaN;
BFi(SCOS.meanVec{1} < 4 ) = NaN;
frameRate = SCOS.info.name.FR;
secForNormalization = 10;
% rBFi = BFi/prctile(BFi(1:round(10*frameRate)),5); % normalize by 5% percentile in first 10 sec
rBFi = BFi/mean(BFi(5:round(secForNormalization*frameRate)),'omitnan'); % normalize by secForNormalization

subplot(nPlots,1,1);
SCOStimeVec =  SCOS.timeVec/60 - minutes(duration(baselineTime,'InputFormat','hh:mm') - duration(startTimeSCOS));
plot(SCOStimeVec,rBFi); hold on;
xlabel('[min]');
ylabel('rBFI');
ylim([0 4]);
grid on; grid minor
markTiming(timingFile)

subplot(nPlots,1,2);
plot(SCOStimeVec,SCOS.meanVec{1}); hold on;
xlabel('[min]');
ylabel('Intensity');
grid on; grid minor
markTiming(timingFile);

if ~isempty(nirsFile)
    subplot(3,1,3);
    NIRStimeVec = minutes(NIRS.Time -  duration(baselineTime,'InputFormat','hh:mm'));
    % colors = {[0 0.4 0],[0.5 0.8 0.2]}; % green
    colors = {[0.8 0 0.2],[0.5 0 0]}; % red
    h(1) = plot(NIRStimeVec,NIRS.Left,'color',colors{1},'DisplayName','Left'); hold on;
    h(2) = plot(NIRStimeVec,NIRS.Right,'color',colors{2},'DisplayName','Right'); hold on;
    grid on; grid minor
    markTiming(timingFile);
    legend(h); % TBD check which is left and which is right
    xlabel('[min]');
    ylabel('Hb Saturation');
end
savefig(fig1, [ recName '\SignalsPlot.fig']);
