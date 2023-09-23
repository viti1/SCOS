windowSize = 9; % for spatial local noise calculation
numOfFrames = 70;
noiseFolder = './Records/NoiseAndBackground';
detectorName = 'Basler_acA1440-220um_SN40335401';
%%
saveAnalysisFolder = [noiseFolder filesep detectorName];

%% 0. How many frames to Average?
% use recording with uniform signal approximately at the middle of the dynamic range, with more than 80 frames
nFrames = 300;
[ rec, recFile ] = RecordFromCamera(nFrames,5,20,100,30,'./Records/NoiseAndBackground/NumOfFrames/','.tiff','WhitePaper',['_' num2str(nFrames) 'frames'],[],"Mono8");
% recFile = [noiseFolder '/NumOfFrames/WhitePaperTint5ms_Gain20dB_FR100Hz_BL30DU_500frames'];
% rec = ReadRecord(recFile);

totRecSize = size(rec,3);
nOfFramesArr = [ 3:3:20 , 30, 40, 50, 60:30:totRecSize ];

localSpatialNoise   = nan(size(nOfFramesArr));
globalSpatialNoise  = nan(size(nOfFramesArr));
temporalNoise       = nan(size(nOfFramesArr));

for i=1:numel(nOfFramesArr)
    nOfFramesArr(i)
    [localSpatialNoise(i), globalSpatialNoise(i), temporalNoise(i) , ~,~,~,meanImPerFrame] = CalcNoise(rec(:,:,1:nOfFramesArr(i)),windowSize);
end

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

subplot(2,1,1);
plot(meanImPerFrame); xlabel('frame number'); ylabel('mean(I) [DU')

savefig(fig,[saveAnalysisFolder '\numOfFrames.fig']);
% save(fig,[saveAnalysisFolder '\numOfFrames.mat'],meanImPerFrame);

%% 1. Read Noise (minimum exposure time, with Cover), vs Gain
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

