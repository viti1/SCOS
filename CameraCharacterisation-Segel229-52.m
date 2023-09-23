%% Set Parameters 
windowSize = 9; % for spatial local noise calculation
numOfFrames = 70;
noiseFolder = '.\Records\NoiseAndBackground';
detectorName = 'Basler_acA1440-220um_SN40335401';
videoFormat = 'Mono8'; 

analysisFolder = [noiseFolder filesep detectorName filesep videoFormat];
if ~exist(analysisFolder,'dir'); mkdir(analysisFolder); end
set(groot,'defaultAxesFontSize',12)
mFolder = fileparts(mfilename('fullpath'));
addpath([fileparts(mfilename('fullpath')) '\func']);
bitDepth = sscanf(videoFormat,'%*[a-zA-Z]%d');

recordFlag = 0; % if off it will try to load the data from the predefined filename
rewriteFlag = 0; % if recordFlag == 0, and recordings already exist, rewrite it.

% Note : each section has two/three parts : Record , Calc and Plot
%% 0. How many frames to Average?
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Num Of Frames - Record and Calc
% use recording with uniform signal approximately at the middle of the dynamic range, with more than 300 frames
numOfFramesFolder = [ noiseFolder filesep 'NumOfFrames'];
matFileNumOfFrames = [ numOfFramesFolder filesep 'data.mat' ];
if recordFlag || ~exist(numOfFramesFolder,'dir')
    nFrames = 300;
   
    camParams.ExposureTime = 1000; % 1ms
    camParams.Gain   = 0;
    forceWrite = 1;
    setupParams.Background = 'WhitePaper'; setupParams.addToFilename = 0;
    camParams.videoFormat = videoFormat;  camParams.addToFilename.videoFormat = false;
    
    % Record
    [ rec, numOfFramesFile, info ] = RecordFromCamera(nFrames,camParams,setupParams,numOfFramesFolder,'.tiff','WhitePaper','',forceWrite);
    
    %     recFile = [numOfFramesFolder '/WhitePaperTint5ms_Gain20dB_FR100Hz_BL30DU_300frames'];
    %     rec = ReadRecord(recFile);

    totRecSize = size(rec,3);
    nOfFramesArr = [ 3:3:20 , 30, 40, 50, 60, 70, 80, 90:30:totRecSize ];

    locSpatNoise   = nan(size(nOfFramesArr));
    globSpatNoise  = nan(size(nOfFramesArr));
    tempNoise      = nan(size(nOfFramesArr));
    
    % Calc
    for i=1:numel(nOfFramesArr)
        fprintf('num of frames = %d\n',nOfFramesArr(i));
        [locSpatNoise(i), globSpatNoise(i), tempNoise(i) ,~,~,~,meanImPerFrame] = CalcNoise(rec(:,:,1:nOfFramesArr(i)),windowSize);
    end
    % Save
    save(matFileNumOfFrames,'locSpatNoise', 'globSpatNoise','tempNoise','nOfFramesArr','numOfFramesFile');
else
    load(matFileNumOfFrames);
end
%% Num Of Frames - Plot
fig = figure('name','How many frames to Average?'); 
subplot(2,1,1)
plot(nOfFramesArr,locSpatNoise,'.-'); hold on
plot(nOfFramesArr,globSpatNoise,'.-');
plot(nOfFramesArr,tempNoise,'.-g'); 
legend('local Spatial', 'globalSpatial', 'temporal');
grid on
xlabel('number of frames');
ylabel('Noise [DU]')
[~,recName] = fileparts(recFile);
title([ recName ' ; Mean = ' num2str(mean2(rec(:,:,1)),3) 'DU' ],'interpreter','none');
grid minor

subplot(2,1,2);
plot(meanImPerFrame); xlabel('frame number'); ylabel('mean(I) [DU')

savefig(fig,[analysisFolder '\numOfFrames.fig'])

%% 1. Read Noise (minimum exposure time, with Cover)
%% Read Noise vs BlackLevel - Record

gainArr = 10;
readNoiseBLFolder = [ analysisFolder filesep 'ReadNoise\vsBlackLevel'];
prefix = 'Cover';
suffix = '';
blackLevelArr = [0:6,10:10:30];
nOfFrames = 90;
frameRate = 100;
tint = 0.021 ; % the minimum avalable

% init
[ locSpatNoise , globSpatNoise, tempNoise , meanI]  = InitNaN([numel(gainArr), numel(blackLevelArr)]);
if ~exist(readNoiseBLFolder,'file'); mkdir(readNoiseBLFolder); end

% record and calc
files = cell([numel(gainArr), numel(blackLevelArr)]);
figRec = figure('name','Record');

vid = videoinput("gentl", 1, videoFormat);% vid =[];
gain = gainArr(1);
for bl_i = 1:numel(blackLevelArr)
    bl   = blackLevelArr(bl_i);
    [rec ,files{gain_i,bl_i}] = RecordFromCamera(nOfFrames,tint,gain,frameRate,bl,readNoiseBLFolder,'.avi',prefix,suffix,vid,videoFormat,forceRewrite,figRec,vid);
    %[locSpatNoise(gain_i,bl_i), globSpatNoise(gain_i,bl_i), tempNoise(gain_i,bl_i), ~,~,~, meanImPerFrame ] = CalcNoise(rec,windowSize,[],0);
    %meanI(gain_i,bl_i) = mean(meanImPerFrame);
end
delete(vid);

%%  Read Noise vs BlackLevel - Calc
[ records , gainArr ] = ChooseRecords( readNoiseBLFolder  ,['Cover*Tint*Gain*'] , 'Gain' );
[ records, blackLevelArr] = SortRecords(records,'BL');
gainArr = unique(gainArr);
blackLevelArr = unique(blackLevelArr);

[ locSpatNoise , globSpatNoise, tempNoise , meanI]  = InitNaN([numel(gainArr), numel(blackLevelArr)]);

%%  Read Noise vs BlackLevel - Plot

%%  Read Noise vs Gain - Record
%%  Read Noise vs Gain - Calc
%%  Read Noise vs Gain - Plot

%%
for gain_i = 4 %1:numel(gainArr)
    for bl_i = 1:numel(blackLevelArr)
        filename =  ChooseRecords( readNoiseBLFolder  , sprintf('Cover*Gain%gdB*BL%gDU*.avi',gainArr(gain_i),blackLevelArr(bl_i) ) );
        disp(filename{1});
        if numel(filename) > 1; error(' more that one option for specific gain (%g) and black level (%g)',gainArr(gain_i),blackLevelArr(bl_i)); end
        rec = ReadRecord(filename{1});
        [locSpatNoise(gain_i,bl_i), globSpatNoise(gain_i,bl_i), tempNoise(gain_i,bl_i), ~,~,~, meanImPerFrame ] = ...
            CalcNoise(rec,windowSize,[],0);
        meanI(gain_i,bl_i) = mean(meanImPerFrame);
    end
end
%%
save([readNoiseBLFolder '\ReadNoise.mat'],'gainArr', 'readNoiseBLFolder',  'prefix' , 'blackLevelArr' ,'nOfFrames', 'frameRate' ,'tint', 'locSpatNoise' ,'globSpatNoise' ,'tempNoise','meanI');

%%
load([readNoiseBLFolder '\ReadNoise.mat']);
%%
% 1A. plot - Find black level
fig1a = figure('name','black level influence on read noise','Units','Normalized','Position',[0.1 0.1 0.9 0.6]);
subplot(1,2,1);
for gain_i = 1:numel(gainArr)
    plot(blackLevelArr,locSpatNoise(gain_i,:),'*-'); hold on;
end
title('Read Noise - Local Spatial')
ylabel('Noise [DU]');
xlabel('Black Level [DU]')
legend(strcat({'Gain '} , num2cellstr(gainArr), 'dB'),'location','best');

subplot(1,2,2);
for gain_i = 1:numel(gainArr)
    plot(blackLevelArr,tempNoise(gain_i,:),'*-'); hold on
end
title('Read Noise - Temporal')
legend(strcat({'Gain '} , num2cellstr(gainArr), 'dB'),'location','best');

for i=1:2
    subplot(1,2,i);
    ylabel('Noise [DU]');
    xlabel('Black Level [DU]')
    set(gca,'YLim', [-0.1 1.3]);
end

savefig(fig1a, [analysisFolder '\ReadNoise vs BlackLevel.fig'])
%%
% 1B. plot  Read Noise in [DU]
blackLevel = 20 ;

bl_i = find(blackLevelArr == blackLevel);
fig1b = figure('name',['ReadNoise BL' num2str(blackLevel)],'Units','Normalized','Position',[0.3 0.2 0.63 0.4]);
subplot(1,2,1)
plot(gainArr,tempNoise(:,bl_i),'*-'); hold on
plot(gainArr,locSpatNoise(:,bl_i),'*-');
plot(gainArr,sqrt(locSpatNoise(:,bl_i).^2 + tempNoise(:,bl_i).^2 ) ,'*-');
ylabel('Noise [DU]');
xlabel('Gain [dB]')
title([' Read Noise [DU] vs Gain , BlackLevel=' , num2str(blackLevel)])
legend({'temporal', 'local spatial', 'total'},'location','northwest');
grid on;
grid  minor 


% 1C plot Read Noise in [e]
wellCapacity = 10.5e3;
nBits = 8;
actualGainArr = reshape( ConvertGain(gainArr,nBits,wellCapacity), [numel(gainArr) 1] );
totNoise = sqrt(locSpatNoise(:,bl_i).^2 + tempNoise(:,bl_i).^2 );

subplot(1,2,2);
plot(gainArr,tempNoise(:,bl_i)./actualGainArr,'*-'); hold on
plot(gainArr,locSpatNoise(:,bl_i)./actualGainArr,'*-');
plot(gainArr, totNoise./actualGainArr ,'*-');
ylabel('Noise [e]');
xlabel('Gain [dB]');
grid on ;
grid  minor 
title([' Read Noise [e] vs Gain , BlackLevel=' , num2str(blackLevel)])
legend({'temporal', 'local spatial', 'total'},'location','northeast');

savefig(fig1b, [analysisFolder '\ReadNoise vs Gain BL' num2str(blackLevel) '.fig'])

% clear 'gainArr', 'readNoiseBLFolder',  'prefix' , 'blackLevelArr' ,'nOfFrames', 'frameRate' ,'tint', 'locSpatNoise' ,'globSpatNoise' ,'tempNoise'
% 'meanI'
%% 2. Dark Current Noise (slope of different exposure times, with Cover) , different Gains
%% -- Record --

gainArr = 0:5:30;
darkCurrentFolder = '.\Records\NoiseAndBackground\DarkCurrent';
prefix = 'Cover';
blackLevel = 30;
nOfFrames = 130;
frameRate = 100;
tintArr = [0.021 , 0.2 , 2 , 10 , 30 , 50 , 75, 100 , 125 , 150 , 200 ];
 
vid = videoinput("gentl", 1, "Mono8");
files = cell(numel(gainArr), numel(tintArr));
for gain = gainArr
    for tint = tintArr
        [~,filename] = RecordFromCamera(nOfFrames,tint,gain,frameRate,blackLevel,darkCurrentFolder,'.avi',prefix,'',vid);
        filename
    end
end
delete(vid);

%% -- Plot --
tintUnits = 'ms';
recFormat = '*Cover*Tint*Gain*.avi';
[ recordsTemp , gainArr] = ChooseRecords([ noiseFolder '\DarkCurrent'] , recFormat , 'Gain');
gainArr = unique(gainArr);
[ ~, tintValues ] = SortRecords(recordsTemp, 'Tint');
tintValues = unique(tintValues);

[ tintVec, offsetVec, locSpatNoise , globSpatNoise, tempNoise , meanI]  = InitNaN([numel(gainArr), numel(tintValues)]);

for gain_i = 1:numel(gainArr)
    [ records , tintVec(gain_i,:) ] = ChooseRecords([ noiseFolder '\DarkCurrent'] ,['*Cover*Tint*Gain' num2str(gainArr(gain_i)) '*.avi'] , 'Tint' );
    for tint_i = 1:numel(records)
        disp(records{tint_i})
        [ rec ,~, offsetVec(gain_i,tint_i)] = ReadRecord(records{tint_i},numOfFrames,tintUnits);
        [locSpatNoise(gain_i,tint_i), globSpatNoise(gain_i,tint_i), tempNoise(gain_i,tint_i), ~,~,~, meanImPerFrame ] = ...
            CalcNoise(rec,windowSize,[],0);
        meanI(gain_i,tint_i) = mean(meanImPerFrame);
    end
end

save([analysisFolder '.\DarkNoise.mat'],'tintVec', 'gainArr', 'offsetVec', 'locSpatNoise' , 'globSpatNoise', 'tempNoise' , 'meanI');
%% 2A. Plot Noise vs Tint
plotVariance = 0; % weather to plot Noise(std) or Var (std^2)
gain = 20; % plot of single gain
gain_i = find(gainArr==gain,1);

if plotVariance
   yStr = 'Noise^2 [DU^2]';
else
   yStr = 'Noise [DU]'; 
end

fig2 = figure('name','Dark Noise - Single Gain'); 
% subplot(2,1,1); 

% plot(tintVec(gain_i,:),globSpatNoise(gain_i,:).^2,'*-');
hold on; 
if plotVariance
    plot(tintVec(gain_i,:),locSpatNoise(gain_i,:).^2,'*-'); %#ok<UNRCH>
    plot(tintVec(gain_i,:),tempNoise(gain_i,:).^2,'*-');
else
    plot(tintVec(gain_i,:),locSpatNoise(gain_i,:),'*-');
    plot(tintVec(gain_i,:),tempNoise(gain_i,:),'*-');
end
legend({'localSpatial','Temporal'},'location','best');
title(['White Paper, Gain' num2str(gain) ' : Noise vs Integration Time ']);
ylabel(yStr)
xlabel(['Tint [' tintUnits ']']);
grid on 

% subplot(2,1,2); 
% plot(tintVec,meanI(gain_i,:) - offsetVec(gain_i,:),'o-');
% ylabel('mean(I) - Offset [DU]')
% xlabel(['Tint [' tintUnits ']']);

savefig(fig2,[analysisFolder '.\Dark Noise Single Gain.fig']) 

%% 2B. Plot noises for different gains
fig3 = figure('name','Dark Current - for Different Gains');
subplot(1,2,2);
for gain_i = 1:numel(gainArr)
    if plotVariance
        plot(tintVec(gain_i,:),tempNoise(gain_i,:).^2,'*-'); %#ok<UNRCH>
    else
        plot(tintVec(gain_i,:),tempNoise(gain_i,:),'*-');
    end
    hold on;
end
ylabel(yStr)
xlabel(['Tint [' tintUnits ']']);
set(gca,'YLim',[0 2])
legend(strcat({'Gain '} , num2cellstr(gainArr), 'dB'))
title('Temporal Noise');
grid on 


subplot(1,2,1);
for gain_i = 1:numel(gainArr)
    if plotVariance            
        plot(tintVec(gain_i,:),locSpatNoise(gain_i,:).^2,'*-'); %#ok<UNRCH>
    else
        plot(tintVec(gain_i,:),locSpatNoise(gain_i,:),'*-');
    end
    hold on;
end
ylabel(yStr)
xlabel(['Tint [' tintUnits ']']);
legend(strcat({'Gain '} , num2cellstr(gainArr), 'dB'))
title('Local Spatial Noise')
set(gca,'YLim',[0 2])
grid on 
savefig(fig3,[analysisFolder '.\Dark Noise All Gain.fig']) 

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% 3. Shot noise (White paper, different illumination levels, different gains) 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%%
%% 4. Temperature influence
% 1. Shot Noise
% 2. Readout Noise
% 3. Dark Current Noise

src.DeviceTemperature
% set(findall(gcf,'-property','FontSize'),'FontSize',18)