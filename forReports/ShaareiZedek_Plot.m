recName = ['C:\Users\' getenv('username') '\OneDrive - Bar Ilan University\SCOS_Records\ShaareiZedek\B4 21.07.2024'];
SCOSdir = dir([recName '\Main*']);
SCOSdir2 = dir([recName filesep SCOSdir.name , '\LocalStd*_corr.mat']);
% SCOSfile = fullfile(recName,SCOSdir(1).name, SCOSdir2(1).name );

[recFolder , recRawName ]= fileparts(recName);
participantID = recRawName(1:2);
NIRSFolder = [ recFolder '\NIRS Data'];
nirsFile   = dir([NIRSFolder '\' participantID '*.csv']); nirsFile = fullfile(NIRSFolder,nirsFile.name);
timingFile = dir([NIRSFolder '\' participantID '*timing info.txt']); timingFile = fullfile(NIRSFolder,timingFile.name);
firstFrameDir = dir([recName,filesep,SCOSdir(1).name '\*_0001.tiff']);
tmp = strsplit(firstFrameDir(1).date);
startTimeSCOS = tmp{2};

NIRS = readtable(nirsFile);
NIRS = renamevars(NIRS,2:7,{'Date','Time','Left','Right','SedlineEvents','Remarks'});
if ischar(NIRS.Left)
    NIRS.Left = str2double(NIRS.Left);
    NIRS.Right = str2double(NIRS.Right);
end
% firstIdx = find(~isnan(NIRS.Right) | ~isnan(NIRS.Left),1);
% startTimeNIRS = NIRS.Time(1);
% NIRStoSCOStimeDiff = etime(datevec(char(startTimeNIRS)),datevec(startTimeSCOS))/60; % Nirs - scos in min

% startTimeNIRS = NIRS.Time(firstIdx);
Timing = readtable(timingFile);
startTime = Timing.Time(1);

%% Plot
f = figure('Units','Normalized','Position',[0.05 0.1 0.9 0.8]);
% SCOS = load();
% figure('Units','Normalized','Position',[0.05 0.05 0.8 0.8]);
% subplot(3,1,1);
% plot(timeVec/60,1./corrSpeckleContrast_jumpsCorrected{1}); hold on;
% xlabel('[min]');
% ylabel('1/K^2');
% markTiming(timingFile,startTimeSCOS)

% subplot(3,1,2);
% plot(timeVec/60,meanVec{1}); hold on;
% xlabel('[min]');
% ylabel('Intensity');
% markTiming(timingFile);

subplot(3,1,3);
NIRStimeVec = minutes(NIRS.Time -  duration(startTime,'InputFormat','hh:mm'));
h(1) = plot(NIRStimeVec,NIRS.Left,'color',[0 0.4 0],'DisplayName','Left'); hold on;
h(2) = plot(NIRStimeVec,NIRS.Right,'color',[0.5 0.8 0.2],'DisplayName','Right'); hold on;
markTiming(timingFile);
legend(h); % TBD check which is left and which is right

