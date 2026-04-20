setupParams.Phantom = '0.68%';

mainFolder = [  '..\Records\VikaHead\SNRvsGainVsTint' ];
if ~exist(mainFolder,'file'); mkdir(mainFolder); end


setupParams=struct();
setupParams.Camera = 'Basler_acA1440-220um_SN40335401';
setupParams.Objective = 'PlanNx20';
setupParams.ObjectiveLocation = 'Close';
% setupParams.Lens = '25mm';
% setupParams.Defocus = 'Middle';
setupParams.SDS = "1.7cm";
setupParams.FiberCollect = '400um';
setupParams.FiberFromLaser = '62.5um';
setupParams.LaserPulsed = 'Yes';
setupParams.Head = 'Vika';

fields = fieldnames(setupParams);
for fi = 1:numel(fields)
    setupParams.addToFilename.(fields{fi}) = false;
end
setupParams.addToFilename.SDS = true;
setupParams.addToFilename.Objective = true;

% setupParams.addToFilename.Lens = true;
% setupParams.addToFilename.Defocus = true;

nOfFrames = 150;
prefix = ''; %ObjectivePlanNx20Close';
recFormat = '.tiff';

camParams = struct();
% camParams.AcquisitionFrameRate = 100;
% camParams.addToFilename.AcquisitionFrameRate = false;
camParams.videoFormat = 'Mono8';
camParams.addToFilename.videoFormat = false;


overwriteFlag = false;
gainArr = [ 5 15];
expTArr = [1 3 5 7 10 15 20 25]*1e3; %[5 10 15 20 25]*1e3; 
suffix = 'FR20Hz';
prefis = 'ObjectiveClose';
setupParams.Table = 'OpticsTable';

camParams.Gain = 30;
camParams.ExposureTime = 20000;
camParams.TriggerMode = 'On';
im = RecordFromCamera(1, camParams, [],[],[],[],[],[],1);
% my_imagesc(im)
[ full_mask , circ ] = GetROI(im);
roi_half_size = ceil(circ.Radius) + 10 ;
roi = [ round(circ.Center(2)) + (-roi_half_size:roi_half_size) ;
        round(circ.Center(1)) + (-roi_half_size:roi_half_size) ];
mask        = full_mask(roi(1,:)  , roi(2,:));
fig=uifigure;
uialert(fig,'Are you Ready?', '','Icon','info','CloseFcn','uiresume(fig)')
uiwait(fig)
%%
for gain = gainArr(:)'
    for expt = expTArr(:)'
        camParams.Gain = gain;
        camParams.ExposureTime = expt;
        % gain
        % expt
        im = RecordFromCamera(1, camParams, [],[],[],[],[],[],1);
        meanI = mean2(im(full_mask));
        N = histcounts(im(full_mask),0:255);  
        if N(end) > numel(im)*0.0001 || meanI > 240 % check saturation
            break;
        end
        fprintf('Mean = %g\n',meanI);
        if meanI < 3
            continue;
        end
        RecordFromCamera(nOfFrames, camParams, setupParams , mainFolder , recFormat, prefix, suffix, overwriteFlag , 0 );            
    end
end
clear gain expt rec

