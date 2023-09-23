windowSize = 9; % for spatial local noise calculation
numOfFrames = 70;
noiseFolder = './Records/NoiseAndBackground';
detectorName = 'CameraData_Basler_acA1440-220um_SN40335401';

saveAnalysisFolder = [noiseFolder filesep detectorName];
if ~exist(saveAnalysisFolder,'dir'); mkdir(saveAnalysisFolder); end
set(groot,'defaultAxesFontSize',12)
addpath('./Code/func');

%% 0. How many frames to Average?
% use recording with uniform signal approximately at the middle of the dynamic range, with more than 80 frames
nFrames = 300;
% [ rec, recFile ] = RecordFromCamera(nFrames,5,20,100,30,'./Records/NoiseAndBackground/NumOfFrames/','.tiff','WhitePaper',['_' num2str(nFrames) 'frames'],[],"Mono8");
recFile = [noiseFolder '/NumOfFrames/WhitePaperTint5ms_Gain20dB_FR100Hz_BL30DU_300frames'];
rec = ReadRecord(recFile);

totRecSize = size(rec,3);
nOfFramesArr = [ 3:3:20 , 30, 40, 50, 60, 70, 80, 90:30:totRecSize ];

localSpatialNoise   = nan(size(nOfFramesArr));
globalSpatialNoise  = nan(size(nOfFramesArr));
temporalNoise       = nan(size(nOfFramesArr));

for i=1:numel(nOfFramesArr)
    nOfFramesArr(i)
    [localSpatialNoise(i), globalSpatialNoise(i), temporalNoise(i) , ~,~,~,meanImPerFrame] = CalcNoise(rec(:,:,1:nOfFramesArr(i)),windowSize);
end
%%
fig = figure('name','How many frames to Average?'); 
subplot(2,1,1)
plot(nOfFramesArr,localSpatialNoise,'.-'); hold on
plot(nOfFramesArr,globalSpatialNoise,'.-');
plot(nOfFramesArr,temporalNoise,'.-g'); 
legend('local Spatial', 'globalSpatial', 'temporal');
grid on
xlabel('number of frames');
ylabel('Noise [DU]')
[~,recName] = fileparts(recFile);
title([ recName ' ; Mean = ' num2str(mean2(rec(:,:,1)),3) 'DU' ],'interpreter','none');
grid minor

subplot(2,1,2);
plot(meanImPerFrame); xlabel('frame number'); ylabel('mean(I) [DU')

savefig(fig,[saveAnalysisFolder '\numOfFrames.fig'])
%% 1. Read Noise (minimum exposure time, with Cover), vs Gain
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
readNoiseFolder = [ noiseFolder '\ReadNoise\BlackLevel2'];
BL = 5;
[ ~ , gainArr ] = ChooseRecords( readNoiseFolder  ,['Cover*Tint*Gain*BL' num2str(BL) '*'] , 'Gain' );
gainArr = unique(gainArr);
% blackLevelArr = unique(blackLevelArr);

[ locSpatNoise , globSpatNoise, tempNoise , meanI]  = InitNaN([size(gainArr)]);

for gain_i = 1:numel(gainArr)
   
        filename =  ChooseRecords( readNoiseFolder  , sprintf('Cover*Gain%gdB*BL%gDU*.avi',gainArr(gain_i),BL ) );
        disp(filename{1});
        if numel(filename) > 1; error(' more that one option for specific gain (%g) and black level (%g)',gainArr(gain_i), num2str(BL)); end
        rec = ReadRecord(filename{1});
        [locSpatNoise(gain_i), globSpatNoise(gain_i), tempNoise(gain_i), ~,~,~, meanImPerFrame ] = CalcNoise(rec,windowSize,[],0);
        meanI(gain_i) = mean(meanImPerFrame);
end
%%
figure; 
plot(gainArr,[locSpatNoise,tempNoise],'*-');
legend({'local spatial','temporal'});
ylabel('Noise [DU]');
xlabel('Gain [dB]');
grid on ;
title([' Read Noise [e] vs Gain , BlackLevel=' , num2str(BL)])

%%
totNoise = sqrt(locSpatNoise(:,bl_i).^2 + tempNoise(:,bl_i).^2 )';
wellCapacity = 10.5e3;
nBits = 8;
actualGainArr = convertGain(gainArr,nBits,wellCapacity);
totNoiseE = totNoise./actualGainArr;
subplot(1,2,2);
plot(gainArr,tempNoise(:,bl_i)'./actualGainArr,'*-'); hold on
plot(gainArr,locSpatNoise(:,bl_i)'./actualGainArr,'*-');
plot(gainArr, totNoiseE ,'*-');
ylabel('Noise [e]');
xlabel('Gain [dB]');
grid on ;
title([' Read Noise [e] vs Gain , BlackLevel=' , num2str(blackLevel)])
% save([readNoiseBLFolder '\ReadNoise.mat'],'gainArr', 'readNoiseBLFolder',  'prefix' , 'blackLevelArr' ,'nOfFrames', 'frameRate' ,'tint', 'locSpatNoise' ,'globSpatNoise' ,'tempNoise','meanI');

%% 2. Dark Current Noise (slope of different exposure times, with Cover)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
currGain = 20;
tintUnits = 'ms';

[ records ] = ChooseRecords([ noiseFolder '\DarkCurrent'] ,['Cover*Gain' num2str(currGain) '*Tint*'] , 'Tint' );
[tintVec , offsetVec, localSpatialNoise , globalSpatialNoise, temporalNoise , meanI]  = InitNaN(size(records));
for rec_i = 1:numel(records)
    [rec, ~, currParams] = ReadRecord(records{rec_i},numOfFrames,{'Tint','Offset'},tintUnits);
    tintVec(rec_i) = currParams(1);
    offsetVec(rec_i) = currParams(2);
    [localSpatialNoise(rec_i), globalSpatialNoise(rec_i), temporalNoise(rec_i), ~,~,~, meanImPerFrame ] =  CalcNoise(rec,windowSize,[],0);
    meanI(rec_i) = mean(meanImPerFrame);
end

% Plot Noise vs Tint

figure; 
subplot(2,1,1); 

plot(tintVec,globalSpatialNoise.^2,'*-');
hold on; 
plot(tintVec,localSpatialNoise.^2,'*-');
plot(tintVec,temporalNoise.^2,'*-');
legend({'globalSpatial','localSpatial','Temporal'},'location','best');
title(['White Paper, Gain' num2str(currGain) ' : Noise vs Integration Time ']);
ylabel('Noise^2 [DU]')
xlabel(['Tint [' tintUnits ']']);


offset = ExtractParametersFromString( records , 'Offset');
subplot(2,1,2); 
plot(tintVec,meanI - offsetVec,'o-');
ylabel('mean(I) - Offset [DU]')
xlabel(['Tint [' tintUnits ']']);
%% 3. Shot noise (White paper, different illumination levels, different gains) 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

