function [localSpatialNoise, globalSpatialNoise, temporalNoise , localSpatialNoiseRel, globalSpatialNoiseRel, temporalNoiseRel, meanImPerFrame ] = ...
                CalcNoise2(rec,windowSize,mask,plot_flag)
    nOfFrames = size(rec,3);
    %nX = size(rec,2); nY = size(rec,1);
    
    if ~exist('mask','var') || isempty(mask)
        mask = CreateCircleMask([size(rec,1),size(rec,2)]);
    end
    
    if ~exist('plot_flag','var') 
        plot_flag = false;
    end
    
    meanIm = mean(rec,3);
    lpIm = medfilt2(meanIm,[windowSize,  windowSize]);

    temporalNoiseMat = std(rec,0,3);
    temporalNoise = mean(temporalNoiseMat(mask));
    temporalNoiseRel = mean(temporalNoiseMat(mask)./lpIm(mask));

    if plot_flag
        figure; 
        subplot(4,1,1);
        hist(temporalNoiseMat(:),500);
        title('temporal Noise histogram');
        xlabel('Noise [DU]')
        ylabel('number of Pixels')
    end

    % globalSpatialNoisePerFrame  = nan(nOfFrames,1);
    meanImPerFrame              = nan(nOfFrames,1);
    localSpatialNoisePerFrame   = nan(nOfFrames,1);
    for k=1:size(rec,3)
        im = rec(:,:,k);
        % globalSpatialNoisePerFrame(k) = std(im(mask));
        meanImPerFrame(k) = mean(im(mask));

        % stdIm = stdfilt(im,true(windowSize));

        % localSpatialNoisePerFrame(k) = mean(stdIm(mask));
    end
    % localSpatialNoisePerFrameRel = localSpatialNoisePerFrame./ meanImPerFrame;
    % localSpatialNoise = mean(localSpatialNoisePerFrame);
    % localSpatialNoiseRel = mean(localSpatialNoisePerFrameRel);
    % globalSpatialNoise = mean(globalSpatialNoisePerFrame);
    % globalSpatialNoiseRel = mean(globalSpatialNoisePerFrame./meanImPerFrame);

    stdIm = stdfilt(meanIm,true(windowSize));
    localSpatialNoise = mean(stdIm(mask));
    localSpatialNoiseRel = mean2(stdIm(mask)./lpIm(mask));
    globalSpatialNoise = std(meanIm(mask));
    globalSpatialNoiseRel = globalSpatialNoise/mean(lpIm(mask));
    
    
    if plot_flag
        subplot(4,1,2);
        plot(globalSpatialNoisePerFrame);
        title('Spatial Noise Per Frame');
        xlabel('Frame Num');
        ylabel('noise [DU]');
        
        subplot(4,1,3)
        plot(localSpatialNoisePerFrame)
        title('local Spatial Noise Per Frame');
        xlabel('Frame Num');
        ylabel(sprintf('std(I(%dx%d)) [DU]',windowSize,windowSize));
        
        subplot(4,1,4)
        plot(localSpatialNoisePerFrameRel)
        title('Normalized local Spatial Noise Per Frame');
        xlabel('Frame Num');
        ylabel(sprintf('std(I(%dx%d))/<I(%dx%d)>',windowSize,windowSize,windowSize,windowSize));

    end

    



