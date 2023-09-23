addpath('.\baseFunc')

setupParams.Camera = 'Basler_acA1440-220um_SN40335401';
setupParams.Objective = 'Edmunds_PlanN_x20';
setupParams.ObjectiveLocation = 'Close';
setupParams.Lens = 'None';
setupParams.Phantom = '0.68%';
setupParams.SDS = "1.5cm";
setupParams.FiberCollect = '400um';
setupParams.FiberFromLaser = '62.5um';
setupParams.LaserPulsed = 'No';
setupParams.addToFilename = false;

recFolder = [ '..\Records\NoiseAndBackground\Basler_1440GS_Vika01\Mono8\FixedPattern\Phantom_IL' setupParams.Phantom];


nOfFrames = 500;
prefix = 'ObjectivePlanNx20Close';
suffix = '';
recFormat = '.tiff';

camParams = struct();
camParams.AcquisitionFrameRate = 100; 
camParams.videoFormat = 'Mono8';
camParams.addToFilename.videoFormat = false;

%%
overwriteFlag = false;
plotFlag = true;
gainArr = 0:5:35;
expTArr = [1 5 10 15 20 25]*1e3; 
for gain = gainArr(:)'
    for expT = expTArr(:)'
        camParams.ExposureTime = expT;
        camParams.Gain = gain;
        im = RecordFromCamera(1, camParams);
        N = histcounts(im(:),256);  
        if N(end) > numel(im)*0.0001 % check saturation
            break;
        end
        fprintf('Mean = %g\n',mean2(im(613+200*[-1 1],852+200*[-1 1])));
        RecordFromCamera(nOfFrames, camParams, setupParams , recFolder , recFormat, prefix, suffix, overwriteFlag , plotFlag ); 
    end
end