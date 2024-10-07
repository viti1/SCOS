recName = 'C:\Users\tarlevi\OneDrive - Bar Ilan University\SCOS_Records\ShaareiZedek\B4 21.07.2024';
SCOSdir = dir([recName '\Main*']);
SCOSdir2 = dir([recName filesep SCOSdir.name , '\LocalStd*_corr.mat']);
SCOSfile = fullfile(recName,SCOSdir(1).name, SCOSdir2(1).name );

[recFolder , recRawName ]= fileparts(recName);
participantID = recRawName(1:2);
NIRSFolder = [ recFolder '\NIRS Data'];
nirsFile   = dir([NIRSFolder '\' participantID '*.csv']); nirsFile = fullfile(NIRSFolder,nirsFile.name);
timingFile = dir([NIRSFolder '\' participantID '*timing info.txt']); timingFile = fullfile(NIRSFolder,timingFile.name);
firstFrameDir = dir([recName,filesep,SCOSdir(1).name '\*_0001.tiff']);
tmp = strsplit(firstFrameDir.date);
startTimeSCOS = tmp{2};

NIRS = readtable(nirsFile);
NIRS = renamevars(NIRS,2:5,{'Date','Time','Left','Right'});
NIRS.Left = str2double(NIRS.Left);
NIRS.Right = str2double(NIRS.Right);
firstIdx = find(~isnan(NIRS.Right) | ~isnan(NIRS.Left),1);
NIRStoSCOStimeDiff = etime(datevec(char(startTimeNIRS)),datevec(startTimeSCOS))/60; % Nirs - scos in min

startTimeNIRS = NIRS.Time(firstIdx);
Timing = readtable(timingFile);
SCOS = load(SCOSfile);
figure('Units','Normalized','Position',[0.05 0.05 0.8 0.8]);
subplot(2,1,1);
plot(timeVec/60,1./corrSpeckleContrast_jumpsCorrected{1}); hold on;
xlabel('[min]');
ylabel('1/K^2');
markTiming(timingFile,startTimeSCOS)

subplot(2,1,2);
plot(timeVec/60,meanVec{1}); hold on;
xlabel('[min]');
ylabel('Intensity');
markTiming(timingFile);

subplot(2,1,3);
NIRStimeVec = minutes(NIRS.Time -  duration(startTimeSCOS));
plot(NIRStimeVec,NIRS.Left,'color',[0 0.7 0]); hold on;
plot(timeVec/60,NIRS.Right,'color',[0.5 0.8 0.2]); hold on;
legend('Left','Right'); % TBD check which is left and which is right
