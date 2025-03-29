% ---------------------------------------------------------------------------------------------------------
% Function: SCOSvsTime_WithNoiseSubtraction_Ver2
%
% Description:
%   This function calculates the speckle contrast (SCOS) over time from an optical recording.
%   It subtracts background (dark) noise and accounts for spatial noise to compute both the raw and 
%   corrected speckle contrasts. Additionally, it calculates derived parameters such as the relative 
%   Blood Flow Index (rBFi) and provides FFT-based pulse frequency analysis.
%
% Usage:
%   [timeVec, rawSpeckleContrast, corrSpeckleContrast, meanVec, info] = ...
%       SCOSvsTime_WithNoiseSubtraction_Ver2(recordName, backgroundName, windowSize, plotFlag, maskInput)
%
% Input Arguments:
%   - recordName:   Path to the main recording (either a file or a folder)
%   - backgroundName:  Path to the dark/background recording (file or folder)
%   - windowSize:   Window size for local standard deviation filtering (used in stdfilt) 
%                   which influences speckle contrast estimation.
%   - plotFlag:     (Optional) Boolean flag to enable (true) or disable (false) plotting; default = true.
%   - maskInput:    (Optional) Either a boolean mask (with the same size as the image) or the boolean
%                   value true to indicate that the entire image is used as the mask.
%
% GUI Mode:
%   When no input arguments are provided, the function runs in GUI mode:
%     - The user selects the recording folder.
%     - The user is prompted to enter a window size.
%     - The user selects the background folder (or file) if not automatically detected.
%     - The first frame is displayed for the user to draw a Region of Interest (ROI). 
%       This ROI is then saved and reused on subsequent runs.
%
% File Naming Requirements:
%   The main and dark recording folder names should follow this format:
%     <Name>_Gain<X>dB_expT<>ms_FR<X>Hz_BL<>DU
%   where expT is the exposure time (ms), FR is the frame rate, and BL is the black level.
%
% Command Mode:
%   In command mode, recordName and windowSize must be provided.
%   The plotFlag is optional (default true), and maskInput can be a boolean value or a mask array.
%
% Output Arguments:
%   - timeVec:              Time vector corresponding to each frame.
%   - rawSpeckleContrast:   Cell array containing raw speckle contrast values (variance/I^2) per channel.
%   - corrSpeckleContrast:  Cell array containing corrected speckle contrast values after noise subtraction.
%   - meanVec:              Cell array of mean intensity values (in digital units, DU) for each frame.
%   - info:                 Structure containing metadata (frame rate, exposure time, gain, etc.) from the recording.
% ---------------------------------------------------------------------------------------------------------

function [timeVec, rawSpeckleContrast, corrSpeckleContrast, meanVec, info] = ...
    SCOSvsTime_WithNoiseSubtraction_Ver2(recordName, backgroundName, windowSize, plotFlag, maskInput)

% Set default plotting flag if not provided
if nargin < 3
    plotFlag = true;
end

% Add path to helper functions stored in the 'baseFunc' folder
addpath('.\baseFunc')

%% Constants
timePeriodForP2P = 2; % Duration [s] for any peak-to-peak related calculations

%% Input Parameter Checking and GUI Mode Handling
if nargin == 0  % If no input arguments, run in GUI mode
    plotFlag = 1;
    % Try to load last used recording path
    if exist('.\lastRec.mat', 'file')
        lastF = load('.\lastRec.mat');        
    else
        lastF.recordName = [fileparts(pwd) '\Records'];
    end
    
    % Prompt user to select the recording folder
    [recordName] = uigetdir(fileparts(lastF.recordName));
    if recordName == 0
        return;  % User canceled the selection
    end
    
    % If multiple .avi files are found, prompt the user to choose one
    if numel(dir([recordName, '\*.avi'])) > 1 
        [recordRawName, recordDir] = uigetfile([recordName '\*.avi']);
        if recordRawName == 0
            return;  % User canceled the selection
        end
        recordName = fullfile(recordDir, recordRawName);
    elseif ( numel(dir([recordName, '\*.avi'])) + numel(dir([recordName, '\*.tiff'])) + ...
             numel(dir([recordName, '\*.tif'])) + numel(dir([recordName, '\*.mat'])) ) < 1 
        errordlg(['No .avi or .tiff/.tif or .mat files found in ' recordName]);
        error(['No .avi or .tiff/.tif or .mat files found in ' recordName]);
    elseif numel(dir([recordName, '\*.avi'])) == 1 && ( numel(dir([recordName, '\*.tiff'])) + numel(dir([recordName, '\*.tif'])) > 1 )
        % If both .avi and .tiff files exist, assume the .avi file is the recording
        d = dir([recordName, '\*.avi']);
        recordName = fullfile(recordName, d(1).name);
    end
    
    % Save the last used recording name for next time
    save('.\lastRec.mat', 'recordName')
    
    % Prompt the user to input the window size
    maxWindowSize = 50; 
    minWindowSize = 3;
    answer = inputdlg('Window Size', '', [1 25], {'7'});
    windowSize = str2double(answer{1});
    if isnan(windowSize) || windowSize > maxWindowSize || windowSize < minWindowSize
        errordlg(['Window Size must be a number between ' num2str(minWindowSize) ' and ' num2str(maxWindowSize)]);
        error(['Window Size must be a number between ' num2str(minWindowSize) ' and ' num2str(maxWindowSize)]);
    end
    clear answer
end

% Determine if recordName refers to a file (file=2) or folder
isRecordFile = exist(recordName, 'file') == 2;

%% Handle Background (Dark) Recording Input
if nargin == 0 || isempty(backgroundName)
    % If a folder with '_dark' appended exists, assume it is the background folder
    if exist([recordName, '_dark'], 'dir')
        backgroundName = [recordName, '_dark'];
    else
        % Search for typical dark recording folder names (e.g., 'DarkIm*', 'background*', 'BG_*')
        dir_Background = [ dir([fileparts(recordName), '\DarkIm*']), ...
                           dir([fileparts(recordName), '\background*']), ...
                           dir([fileparts(recordName), '\BG_*']) ];
    
        if isempty(dir_Background) || numel(dir_Background) > 1 
            if isRecordFile
                backgroundName = uigetfile(fileparts(recordName), 'Please Select the background');
            else
                backgroundName = uigetdir(fileparts(recordName), 'Please Select the background'); 
            end
        else
            backgroundName = fullfile(fileparts(recordName), dir_Background(1).name);
        end
        if isequal(backgroundName, 0)
            disp('Aborting...')
            return;
        end
    end
end

%% Create Region of Interest (ROI) Mask
% Build a descriptive name based on folder structure for later use
upFolders = strsplit(recordName, filesep);
rawName = strrep(strjoin(upFolders(end-2:end), '; '), '_', ' ');

% Create prefix for saving related files (depending on file or folder input)
if exist(recordName, 'file') == 7  % recordName is a folder
    recSavePrefix = [recordName filesep];
else  % recordName is a file
    recSavePrefix = [recordName(1:find(recordName=='.', 1, 'last')-1) '_'];
end

% Compute an average frame from the first 20 frames to help define the ROI
[mean_frame] = mean(ReadRecord(recordName, 20), 3); 

maskFile = [recSavePrefix 'Mask.mat'];
if exist('maskInput', 'var')
    % If maskInput is provided, process it accordingly
    if isequal(maskInput, true)
        masks = {true(size(mean_frame))};
        totMask = masks{1};
        disp('Mask is the whole image')
    elseif iscell(maskInput)
        masks = maskInput;
        for k = 1:numel(masks)
            if ~isequal(size(masks{k}), size(mean_frame))
                error('wrong maskInput size, should be the same as the recording');
            end
        end
        totMask = masks{1};
        for k = 2:numel(masks)
            totMask = totMask | masks{k};
        end
        disp('Input mask is a cell array')
    elseif islogical(maskInput)
        if ~isequal(size(maskInput), size(mean_frame))
            error('wrong maskInput');
        end
        masks = {maskInput};
        totMask = masks{1};
        disp('Input single mask')
    else
        error('wrong maskInput type');
    end
    loadExistingFile_flag = 0;
else    
    % If no mask is provided, check for an existing mask file
    if ~exist(maskFile, 'file')
        loadExistingFile_flag = 0;    
    else
        if nargin == 0  % In GUI mode, ask whether to reuse the existing mask
            answer = questdlg('Mask file already exist, do you want to define it again?', '', 'Yes', 'No', 'No');
            loadExistingFile_flag = strcmp(answer, 'No');
        else  % In command mode, automatically load the saved mask
            loadExistingFile_flag = 1;
        end
    end
end

% Exclude margins from the mask based on half the window size to avoid border artifacts
ws2 = ceil(windowSize / 2);
if ~exist('masks', 'var')
    if ~loadExistingFile_flag
        % Allow the user to define the ROI by drawing a circle on the mean image
        [totMask, circ, figMask] = GetROI(mean(ReadRecord(recordName, 20), 3), windowSize);
        masks{1} = totMask;       
        channels.Centers = circ.Center;
        channels.Radii = circ.Radius;
        % Save the mask for future use
        save(maskFile, 'masks', 'channels', 'totMask');
        savefig(figMask, [recordName '\maskIm.fig'])
    else 
        % Load previously saved mask
        M = load(maskFile);
        if isfield(M, 'mask')
            masks{1} = false(size(M.mask));
            masks{1}(ws2+1:end-ws2, ws2+1:end-ws2) = M.mask(ws2+1:end-ws2);
            totMask = M.mask > 0;
            channels.Centers = M.circ.Center;
            channels.Radii  = M.circ.Radius;
        else
            load(maskFile);
            for k = 1:numel(masks)
                masks{k}([1:ws2, (end-ws2+1):end], :) = false;
                masks{k}(:, [1:ws2, (end-ws2+1):end]) = false;
            end
            totMask([1:ws2, (end-ws2+1):end], :) = false;
            totMask(:, [1:ws2, (end-ws2+1):end]) = false;            
        end
    end
else
    % Adjust any provided masks to remove edge regions where filtering is not valid
    for k = 1:numel(masks)
        masks{k}([1:ws2, (end-ws2+1):end], :) = false;
        masks{k}(:, [1:ws2, (end-ws2+1):end]) = false;
    end
    totMask([1:ws2, (end-ws2+1):end], :) = false;
    totMask(:, [1:ws2, (end-ws2+1):end]) = false;
end

%% Retrieve Recording Metadata and Validate Background Parameters
upFolders = strsplit(recordName, filesep);
shortRecName = strjoin(upFolders(end-2:end));

nOfFrames = GetNumOfFrames(recordName);
info = GetRecordInfo(recordName);
im1 = ReadRecord(recordName, 1);

% Ensure that critical background parameters (exposure time, gain, black level) match between
% the main recording and the background recording.
if backgroundName ~= 0
    info_background = GetRecordInfo(backgroundName);
    fields = {'expT', 'Gain', 'BL'};
    for fi = 1:numel(fields)
        param = fields{fi};
        if info_background.name.(param) ~= info.name.(param)
            error('Background.%s=%g   Record.%s=%g', param, info_background.name.(param), param, info.name.(param));
        end
    end
end

% Determine the frame rate from the camera information or the recording name
if isfield(info, 'cam') && isfield(info.cam, 'AcquisitionFrameRate')
    frameRate = info.cam.AcquisitionFrameRate;
else    
    if ~isfield(info.name, 'FR') || isnan(info.name.FR)
        error('Frame Rate must be part of the recording name as "FR"');
    end
    frameRate = info.name.FR; 
end

%% Load Background and Compute Background Noise Statistics
start_calib_time = tic;
disp('Load Background');
if ~exist(backgroundName, 'file')
    error([backgroundName, ' does not exist!']);
end

requiredBG_nOfFrames = 400;
if exist(backgroundName, 'file') == 7  % Background is provided as a folder
    if exist([backgroundName, '\meanIm.mat'], 'file')
        bgS = load([backgroundName, '\meanIm.mat']);
        if ~isfield(bgS, 'nOfFrames') 
            nOfFramesBG = GetNumOfFrames(backgroundName); 
        else
            nOfFramesBG = bgS.nOfFrames;
        end
        if nOfFramesBG < requiredBG_nOfFrames
            error('Not enough frames in background file. Required: %d, Exist: %d', requiredBG_nOfFrames, nOfFramesBG);
        end
        if isfield(bgS, 'recVar')
            background = bgS.recMean - info_background.name.BL;
            darkVar = bgS.recVar;
        else
            % Recalculate background statistics if necessary
            delete([backgroundName, '\meanIm.mat']);
            [background, darkVar] = ReadRecordVarAndMean(backgroundName);
            background = background - info_background.name.BL;
        end
        clear bgS
    else
        nOfFramesBG = GetNumOfFrames(backgroundName);
        if nOfFramesBG < requiredBG_nOfFrames
            error('Not enough frames in background file. Required: %d, Exist: %d', requiredBG_nOfFrames, nOfFramesBG);
        end
        [background, darkVar] = ReadRecordVarAndMean(backgroundName);
        background = background - info_background.name.BL;

        if abs(mean2(background)) > 3
            my_imagesc(background);
            title('Background');
            warning('Suspicious level of the background %gDU!', round(mean2(background), 2));
        end        
    end    
elseif endsWith(backgroundName, '.mat')
    % Load background information from a .mat file
    bgS = load(backgroundName);
    fields = fieldnames(bgS);
    if ismember('recMean', fields)
        background = bgS.recMean - info_background.name.BL;
    elseif startsWith(fields{1}, 'Video')
        darkRec = bgS.(fields{1}); 
        bgS.recMean = mean(darkRec, 3);
        background = bgS.recMean - info_background.name.BL;
        darkVar = std(darkRec, 0, 3).^2;
        bgS.recVar = darkVar;
        save(backgroundName, '-struct', 'bgS')
    else
        error('wrong fields in background file')
    end
    clear bgS
end

% Compute the local dark variance over a window using a box filter
darkVarPerWindow = imboxfilt(darkVar, windowSize);

% Verify that the background image and the recording frame have the same dimensions
if ~isequal(size(background), size(im1))
    error('The background should be the same size as the record. Record size is [%d,%d], but background size is [%d,%d]', ...
        size(im1,1), size(im1,2), size(background,1), size(background,2));
end

%% Retrieve Camera Gain and Black Level Information
nOfBits = info.nBits;
actualGain = GetActualGain(info);

if ~isfield(info.name, 'BL')
    BlackLevel = 0;
else
    BlackLevel = info.name.BL;
end

% Determine the file path for saving/loading spatial noise and smoothing coefficients
if isRecordFile
    smoothCoeffFile = [fileparts(recordName) '\smoothingCoefficients.mat'];
else
    smoothCoeffFile = [recordName '\smoothingCoefficients.mat'];
end

% Calculate spatial noise and smoothing coefficients if not already computed
if ~exist(smoothCoeffFile, 'file')
    disp('Calculating Spatial Noise and Smoothing Coefficients');
    numFramesForSPNoise = 400;
    if nOfFrames > 500
        numFramesForSPNoise = 500;
    end
    spRec = ReadRecord(recordName, numFramesForSPNoise) - BlackLevel;
    spIm = mean(spRec, 3) - background;
    fig_spIm = my_imagesc(spIm);
    title(['Image average over ' num2str(numFramesForSPNoise) ' frames']);
    savefig(fig_spIm, [recordName '\spIm.fig']);
    spVar = stdfilt(spIm, true(windowSize)).^2;
    [fitI_A, fitI_B] = FitMeanIm(spRec, totMask, windowSize);
    clear spRec
    save(smoothCoeffFile, 'spVar', 'fitI_A', 'fitI_B', 'spIm', 'totMask');
else
    disp('Loading Spatial Noise and Smoothing Coefficients');
    load(smoothCoeffFile);
end

disp('Calibration Time');
toc(start_calib_time);

%% Define a Reduced ROI for Faster Processing
% Find the coordinates of the ROI from the mask (totMask)
[y, x] = find(totMask);
roi_lims = [ min(y) - windowSize, max(y) + windowSize;
             min(x) - windowSize, max(x) + windowSize ];

% Ensure the ROI limits are within the valid image range
roi_lims(roi_lims < 1) = 1;        
if roi_lims(1,2) > info.imageSize(1)
    roi_lims(1,2) = info.imageSize(1);
end
if roi_lims(2,2) > info.imageSize(2)
    roi_lims(2,2) = info.imageSize(2);
end
roi.y = roi_lims(1,1):roi_lims(1,2);
roi.x = roi_lims(2,1):roi_lims(2,2);

% Crop the masks to the ROI to speed up processing
masks_cut = cell(size(masks));
for ch = 1:numel(masks)
    masks_cut{ch} = masks{ch}(roi.y, roi.x); 
end
   
fitI_A_cut = fitI_A(roi.y, roi.x);
fitI_B_cut = fitI_B(roi.y, roi.x);

%% Calculate Speckle Contrast (SCOS) for Each Frame
disp(['Calculating SCOS on "' recordName '" ...']);
disp(['Mono' num2str(nOfBits)]);
nOfChannels = numel(masks);
frameNames = dir([recordName, '\*.tiff']);

% Initialize cell arrays (one per channel) to store speckle contrast and mean intensity values
[rawSpeckleContrast, corrSpeckleContrast, meanVec] = InitNaN([nOfFrames 1], nOfChannels);

% Read the first frame to determine scaling factors (depending on bit depth)
if isRecordFile
    rec = ReadRecord(recordName);
    im1 = double(rec(:,:,1));
else
    im1 = double(imread([recordName, filesep, frameNames(1).name]));
end

% Adjust the intensity scaling if the recorded bit depth implies that the data is shifted
devide_by = 1;
if nOfBits == 12 && all(mod(im1(:), 2^4) == 0)
    devide_by = 2^4;
elseif nOfBits == 10 && all(mod(im1(:), 2^6) == 0)
    devide_by = 2^6;
end

start_scos = tic;
               
for i = 1:nOfFrames
    if i == 50
        time50frames = toc(start_scos);
        fprintf('\n Estimated Time = %g min (%d frames)\n', round(time50frames/50*nOfFrames/60,2), nOfFrames)
    end
    if mod(i,200) == 0 
        fprintf('%d\t', i);
        if mod(i,2000) == 0
            fprintf('\n');
        end
    end
    
    % Read the current frame and adjust intensity values
    if isRecordFile 
        im_raw = double(rec(:,:,i));
    else
        im_raw = double(imread([recordName, filesep, frameNames(i).name])) / devide_by;               
        im_raw = im_raw - BlackLevel;
    end
    
    % Subtract the background image from the raw frame
    im = im_raw - background;
    im_cut = im(roi.y, roi.x);
    
    % Compute the local standard deviation image using a window (stdfilt)
    stdIm = stdfilt(im_cut, true(windowSize));
    
    % Process each channel/ROI separately
    for ch = 1:nOfChannels
        % Calculate the mean intensity within the ROI for this channel
        meanFrame = mean(im_cut(masks_cut{ch}));
        % Use pre-computed calibration coefficients to estimate the fitted intensity image
        fittedI = fitI_A_cut * meanFrame + fitI_B_cut;         
        fittedISquare = fittedI.^2;
        
        % Compute the raw speckle contrast (variance normalized by squared intensity)
        rawSpeckleContrast{ch}(i) = mean((stdIm(masks_cut{ch}).^2 ./ fittedISquare(masks_cut{ch})));
        
        % Compute the corrected speckle contrast by subtracting known noise sources:
        % - Shot noise scaled by the actual gain (actualGain * fittedI)
        % - Spatial noise (spVar)
        % - Quantization noise (1/12)
        % - Local dark variance (darkVarPerWindow)
        corrSpeckleContrast{ch}(i) = mean((stdIm(masks_cut{ch}).^2 - actualGain .* fittedI(masks_cut{ch}) ...
            - spVar(masks_cut{ch}) - 1/12 - darkVarPerWindow(masks_cut{ch})) ./ fittedISquare(masks_cut{ch}));
        
        % Save the mean intensity for later analysis or plotting
        meanVec{ch}(i) = meanFrame;
        
        % For the first frame, print detailed information for debugging/calibration
        if i == 1
            fprintf('<I>=%.3gDU , K_raw = %.5g , Ks=%.5g , Kr=%.5g, Ksp=%.5g, Kq=%.5g, Kf=%.5g\n', ...
                meanFrame, rawSpeckleContrast{ch}(i), ...
                mean(actualGain .* fittedI(masks_cut{ch}) ./ fittedISquare(masks_cut{ch})), ...
                mean(darkVar(masks_cut{ch}) ./ fittedISquare(masks_cut{ch})), ...
                mean(spVar(masks_cut{ch}) ./ fittedISquare(masks_cut{ch})), ...
                mean(1./(12 * fittedISquare(masks_cut{ch}))), ...
                corrSpeckleContrast{ch}(i));
        end
    end 
end
fprintf('\n');

%% Create Time Vector Based on the Frame Rate
timeVec = (0:(nOfFrames-1))' * (1/frameRate);   % Time in seconds for each frame
p2p_time = timeVec < timePeriodForP2P;  % Logical index for frames within the first 'timePeriodForP2P' seconds

%% Save Processed Data to File
stdStr = sprintf('Std%dx%d', windowSize, windowSize);
if exist([recSavePrefix 'Local' stdStr '.mat'], 'file')
    delete([recSavePrefix 'Local' stdStr '.mat']); % Delete old file to ensure correct date stamp
end

% Obtain the start date/time from the first frame's file information (assumes TIFF naming convention)
firstFrameDir = dir([recordName, '\*0001.tiff']);
startDateTime = firstFrameDir.date;
save([recSavePrefix 'Local' stdStr '_corr.mat'], 'startDateTime', 'timeVec', 'corrSpeckleContrast', ...
    'rawSpeckleContrast', 'meanVec', 'info', 'nOfChannels', 'recordName', 'windowSize');

%% Plot SCOS Signals and FFT Analysis
% Construct a title string using metadata from the recording
infoFields = fieldnames(info.name);
if ~isfield(info.name, 'Gain')
    info.name.Gain = '';
end
firtsParamValue = info.name.(infoFields{1});
if ~ischar(firtsParamValue)
    firtsParamValue = num2str(firtsParamValue);
end
if isfield(info.name, 'SDS')
    titleStr = [infoFields{1} firtsParamValue ' SDS=' num2str(info.name.SDS) '; exp=' num2str(info.name.expT) 'ms; Gain=' num2str(info.name.Gain) 'dB'];
else
    titleStr = [infoFields{1} firtsParamValue '; exp=' num2str(info.name.expT) 'ms; Gain=' num2str(info.name.Gain) 'dB'];
end
    
% Calculate signal-to-noise ratio (SNR) and FFT parameters for both raw and corrected contrasts
[raw_SNR, raw_FFT, raw_freq, raw_pulseFreq, raw_pulseBPM] = CalcSNR_Pulse(rawSpeckleContrast{1}, frameRate);
[corr_SNR, corr_FFT, corr_freq, corr_pulseFreq, corr_pulseBPM] = CalcSNR_Pulse(corrSpeckleContrast{1}, frameRate);

if plotFlag
    fig = figure('name', ['SCOS ' recordName ' Mono' num2str(nOfBits)], 'Units', 'Normalized', 'Position', [0.1 0.1 0.8 0.8]);
    subplot(3,2,1);
        plot(timeVec, corrSpeckleContrast{1});
        ylabel('Corrected Contrast (var/I^2)');
        title({titleStr, 'Corrected Signal'});
    subplot(3,2,2);
        plot(corr_freq, corr_FFT);
        ylabel('FFT');
        xlabel('Frequency [Hz]');
        title(sprintf('Corrected FFT: SNR=%.2g Pulse=%.0fbpm', corr_SNR, corr_pulseBPM));
    subplot(3,2,3);
        plot(timeVec, rawSpeckleContrast{1});
        ylabel('Raw Contrast (var/I^2)');
        title('Raw Signal');
    subplot(3,2,4);
        plot(raw_freq, raw_FFT);
        ylabel('FFT');
        xlabel('Frequency [Hz]');
        title(sprintf('Raw FFT: SNR=%.2g Pulse=%.0fbpm', raw_SNR, raw_pulseBPM));
    subplot(3,2,5);
        plot(timeVec, meanVec{1});
        ylabel('I [DU]');
        xlabel('Time [s]');
    
    savefig(fig, [recSavePrefix 'Local' stdStr '_plot.fig']);
end

%% Plot Relative Blood Flow Index (rBFi)
if any(corrSpeckleContrast{1} < 0)
    warning('Error: There are negative values in the contrast !!!');
    warndlg('Error: There are negative values in the contrast !!! BFI has no meaning');
end

if plotFlag
    % Use minutes as time unit if the recording duration exceeds 2 minutes
    if timeVec(end) > 120
        timeToPlot = timeVec / 60;
        xLabelStr = 'Time [min]';
    else
        timeToPlot = timeVec;
        xLabelStr = 'Time [sec]';
    end

    fig7 = figure('Name', ['rBFi: ' recordName], 'Units', 'Normalized', 'Position', [0.1, 0.1, 0.4, 0.4]); 
    subplot(2,1,1);
    % Invert the corrected contrast to derive BFi and normalize using the first second of data
    BFi = 1 ./ corrSpeckleContrast{1};
    rBFi = BFi / mean(BFi(1:round(1 * frameRate)));
    plot(timeToPlot, rBFi); 
    title(titleStr);
    xlabel(xLabelStr);
    ylabel('rBFi');
    grid on;
    set(gca, 'FontSize', 10);
    
    subplot(2,1,2);
    plot(timeToPlot, meanVec{1});
    xlabel(xLabelStr);
    ylabel('<I> [DU]');
    grid on;
    set(gca, 'FontSize', 10);
    
    savefig(fig7, [recSavePrefix '_rBFi.fig']);
end

toc(start_scos)
end
