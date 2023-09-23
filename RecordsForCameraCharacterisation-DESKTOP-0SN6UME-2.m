%% 1. Read Noise (minimum exposure time, with Cover)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% 2. Dark Current Noise (slope of different exposure times, with Cover) , different Gains
% -- Record --

gainArr = 0:5:30;
darkCurrentFolder = '.\Records\NoiseAndBackground\DarkCurrent';
prefix = 'Cover';
blackLevel = 30;
nOfFrames = 130;
frameRate = 100;
tintArr = [0.021 , 0.2 , 2 , 10 , 30 , 50 , 75, 100 , 125 , 150 , 200 ];
folder 
vid = videoinput("gentl", 1, "Mono8");
files = cell(numel(gainArr), numel(tintArr));
for gain = gainArr
    for tint = tintArr
        RecordFromCamera(nOfFrames,tint,gain,frameRate,blackLevel,darkCurrentFolder,'.avi',prefix,'',vid);
    end
end
delete(vid);

%% -- Plot --
tintUnits = 'ms';

[ recordsTemp , gainArr] = ChooseRecords([ noiseFolder '\DarkCurrent'] , 'Cover*Gain*Tint*', 'Gain');
gainArr = unique(gainArr);
[ ~, tintValues ] = SortRecords(recordsTemp, 'Tint');
tintValues = unique(tintValues);

[ tintVec, offsetVec, locSpatNoise , globSpatNoise, tempNoise , meanI]  = InitNaN(numel(gainArr), numel(tintValues));

for gain_i = 1:numel(gainArr)
    [ records , tintVec(gain_i,:) ] = ChooseRecords([ noiseFolder '\DarkCurrent'] ,['Cover*Gain' num2str(gainArr(gain_i)) '*Tint*'] , [] , 'Tint' );
    for tint_i = 1:numel(records)
        [ rec ,~, offsetVec(gain_i,tint_i)] = ReadRecord(records{tint_i},numOfFrames,tintUnits);
        [locSpatNoise(gain_i,tint_i), globSpatNoise(gain_i,tint_i), tempNoise(gain_i,tint_i), ~,~,~, meanImPerFrame ] = ...
            CalcNoise(rec,windowSize,[],0);
        meanI(gain_i,tint_i) = mean(meanImPerFrame);
    end
end

% A. Plot Noise vs Tint
gain = 20;
gain_i = find(gainArr==gain,1);
figure; 
subplot(2,1,1); 

plot(tintVec(gain_i,:),globSpatNoise(gain_i,:).^2,'*-');
hold on; 
plot(tintVec(gain_i,:),locSpatNoise(gain_i,:).^2,'*-');
plot(tintVec(gain_i,:),tempNoise(gain_i,:).^2,'*-');
legend({'globalSpatial','localSpatial','Temporal'},'location','best');
title(['White Paper, Gain' num2str(currGain) ' : Noise vs Integration Time ']);
ylabel('Noise^2 [DU]')
xlabel(['Tint [' tintUnits ']']);


subplot(2,1,2); 
plot(tintVec,meanI - offsetVec(gain_i,:),'o-');
ylabel('mean(I) - Offset [DU]')
xlabel(['Tint [' tintUnits ']']);

% Plot noises for different gains
figure(5,'name','Dark Current - for Different Gains';
subplot(2,1,1);
for gain_i = 1:numel(gainArr)
    plot(tintVec(gain_i,:),tempNoise(gain_i,:).^2,'*-');
    hold on;
end
legend(strcat({'Gain '} , num2cellstr(gainArr), 'dB'))
title('Temporal Noise')
subplot(2,1,1);
for gain_i = 1:numel(gainArr)
    plot(tintVec(gain_i,:),locSpatNoise(gain_i,:).^2,'*-');
    hold on;
end
legend(strcat({'Gain '} , num2cellstr(gainArr), 'dB'))
title('Local Spatial Noise')

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% 3. Shot noise (White paper, different illumination levels, different gains) 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

