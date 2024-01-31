%  ---------------------------------------------------------------------------------------------------------
%  [ timeVec,  , rawSpeckleContrast , rawSpeckleVar, corrSpeckleVar , corrSpeckleContrast, imMeanVec , info] = PlotSCOSvsTime(recordName,windowSize,plotFlag,maskInput)
%  GUI mode:   - Choose the recording *folder*
%              - Choose widow size on which to do the std ( it will be used in stdfilter() function in order to calc the local std)
%              - First frame of the recording will appear, and then the
%                user should create a circle with the ROI (Region of Interest)
%                New figure will be created with speckleNoiseVec and speckleNoiseVecRel,
%                which is speckleNoiseVec devided by the the mean intensity inside the mask)
%                In the second time this function will run for the same
%                recording , the ROI is already saved , so ne need to
%                choose it again.
%                If the recording is .tiff files, FR (FrameRate) parameter should appeare in the name of the folder.
%                Example: murecord_FR20Hz 
%               
%  Command Mode: Same as GUI mode , but recordName and windowSize variables must be specified. 
%                plotFlag - [optional] create figure with 5 graphs, defualt
%                = true
%                
%                maskInput - true (then all the image is taken as mask) or
%                boolaen map (same size as record two first dimentions
%  ---------------------------------------------------------------------------------------------------------

function [ timeVec, rawSpeckleContrast , rawSpeckleVar, corrSpeckleVar , corrSpeckleContrast, imMeanVec , info] = ...
    SCOSvsTime(recordName,windowSize,plotFlag,maskInput)
if nargin <3
    plotFlag = true;
end
%% Constants
timePeriodForP2P = 2; % [s]

%% Check input parameters
if nargin == 0 % GUI mode
    plotFlag = 1;
    [recordName] = uigetdir();
    if recordName == 0; return; end % if 'Cancel' was pressed
    if numel(dir([recordName, '\*.avi' ])) > 1 
        [recordRawName, recordDir] = uigetfile([recordName '\*.avi']);
        if recordRawName == 0; return; end % if 'Cancel' was pressed
        recordName = fullfile(recordDir, recordRawName);
    elseif ( numel(dir([recordName, '\*.avi' ])) + numel(dir([recordName, '\*.tiff' ])) + numel(dir([recordName, '\*.tif' ])) +  numel(dir([recordName, '\*.mat' ])) ) < 1 
        errordlg(['No .avi or .tiff/.tif or .mat files found in ' recordName ])
        error(['No .avi or .tiff/.tif or .mat files found in ' recordName ]);
    elseif numel(dir([recordName, '\*.avi' ])) == 1 && ( numel(dir([recordName, '\*.tiff' ])) + numel(dir([recordName, '\*.tif' ])) ) > 1 
        % if in folder apear both .avi and .tiff files -> assume that .avi is the recording
        d = dir([recordName, '\*.avi' ]);
        recordName = fullfile( recordName , d(1).name );
    end
    
    maxWindowSize = 50; minWindowSize = 3;
    answer = inputdlg('Window Size','',[1 25],{'9'});
    windowSize = str2double(answer{1});
    if isnan(windowSize) || windowSize > maxWindowSize || windowSize < minWindowSize
        errordlg(['Window Size must be a number between ' num2str(minWindowSize)  ' and ' num2str( num2str(maxWindowSize) ) ]);
        error(['Window Size must be a number between ' num2str(minWindowSize)  ' and ' num2str( num2str(maxWindowSize) ) ])
    end
    
    clear answer
end

%% Create Mask
upFolders = strsplit(recordName,filesep);
rawName = strrep( strjoin(upFolders(end-2:end),'; '), '_',' ');

if exist(recordName,'file') == 7 % it's a folder
    recSavePrefix = [ recordName filesep ];
else % it's a file
    recSavePrefix = [ recordName(1:find(recordName=='.',1,'last')-1) '_' ];
end

[ first_frame ] = ReadRecord(recordName,1); 

maskFile = [recSavePrefix 'Mask.mat'];
if exist('maskInput','var')
    if isequal(maskInput, 1)
        mask = true(size(first_frame));
    elseif islogical(maskInput)
        if ~isequal( size(maskInput), size(first_frame) ) 
            error('wrong maskInput');
        end
        mask = maskInput;
    else
        error('wrong maskInput type');
    end
    loadExistingFile_flag = 0;
else    
    if ~exist(maskFile,'file')
        % TBD : add automatic circle recognition
        loadExistingFile_flag = 0;    
    else
        if nargin == 0 % GUI mode
            answer = questdlg('Mask file already exist, do you want to define it again?','','Yes','No','No');
            loadExistingFile_flag = strcmp(answer,'No');
        else      % Command line mode -> load the saved mask 
            loadExistingFile_flag = 1;
        end
    end
end

if ~exist('mask','var')    
    if ~loadExistingFile_flag
        %% ask the user to mark a circle on the image

        f = figure('position',[50,50,1200,800]); imagesc(first_frame); colormap gray; colorbar
        title(rawName,'interpreter','none');
        totMask  = false(size(first_frame,1),size(first_frame,2));    
        drawOneMore = true;
        k = 1; circles = struct();
        while drawOneMore
            circ = drawcircle('Color','r','FaceAlpha',0.2);
            circles(k).Center = circ.Center; 
            circles(k).Radius = circ.Radius;
            k=k+1;
            [x,y] = meshgrid(1:size(first_frame,2),1:size(first_frame,1));
            
            masks{k} = false(size(first_frame,1),size(first_frame,2)); %#ok<AGROW>
            masks{k}((x-circ.Center(1)).^2 + (y-circ.Center(2)).^2 < circ.Radius^2 ) = true; %#ok<AGROW>
            totMask = masks{k} | totMask ;
            answer = questdlg('One more channel?', '','Yes','No','Yes');
            drawOneMore = strcmp(answer,'Yes');
        end
        
        nOfChannels = numel(circles);
        save(maskFile,'mask','circles');
        close(f)
    else 
        load(maskFile)
    %     figure; imshowpair(first_frame,mask);  
    end
else
    
end

%% Read Record
disp(['Reading Record "' recordName '" ... '])
[ head_rec , info] = ReadRecord(recordName);

%%  Calc Specle Contrast
disp(['Calculation SCOS on "' recordName '" ... '])
nOfFrames = size(head_rec,3);

%% Load rellevant parameters from camera characterizaion file 
% dData = load('./BeslerAce1440-200u_.mat');  % detector Data
% readoutN   = interp1(detectorData.gainArr,detectorData.readoutNoise, info.Gain ,'spline');
% actualGain = 10^(info.Gain+detectorData.g0); 
% H = fspecial('average', [1 1]*windowSize);
% backgroundMFile = [ recSavePrefix 'background.mat' ]; 
% if exist(backgroundMFile,'file')
%     bgS = load(backgroundMFile);
%     background = bgS.recMean;
% else
%     if exist([ recSavePrefix 'background' ],'file')
%         backgroundRecFile = ReadRecord([ recSavePrefix 'background' ]);
%         if exist(RecordName,'file') ~= 7 % it's not a folder
%             [~,~,ext] = fileparts(RecordName);
%             backgroundRecFile = [ backgroundRecFile ext ];
%         end
%         backgroundRec = ReadRecord(backgroundRecFile);
%     %     background = PixFix(mean(backgroundRec,3),dData.BPMask);
%         background = mean(backgroundRec,3);
%     else
%         background = zeros(size(first_frame));
%     end
% end
% background = zeros(size(first_frame));
[rawSpeckleVar , rawSpeckleContrast , corrSpeckleVar , corrSpeckleContrast , imMeanVec] = InitNaN([nOfFrames 1]);
for i=1:nOfFrames
    im = head_rec(:,:,i);
    stdIm = stdfilt(im,true(windowSize));
    for k = 1:nOfChannels
        meanPerChannel{k} = mean(im(masks{k}));
        rawSpeckleVar{k}(i) = mean(stdIm(masks{k}))^2;
        rawSpeckleContrast{k}(i) = rawSpeckleVar(i)/meanFrame^2;
    end
    imMeanVec(i) = meanFrame; % mean(meanIm);
end
IMean = mean(imMeanVec);

%% Create Time vector
if ~isfield(info.name,'FR') || isnan(info.name.FR)
    error('Frame Rate must be part of the recording name as "FR"');
end
frameRate = info.name.FR; % TBD

timeVec = (0:(nOfFrames-1))'*(1/frameRate) ;   % FR = FrameRate
p2p_time = timeVec<timePeriodForP2P;

%% Plot
stdStr = sprintf('Std%dx%d',windowSize,windowSize);

infoFields = fieldnames(info.name);
if isfield(info.name,'SDS')
    titleStr =  [ infoFields{1} ' SDS=' num2str(info.name.SDS)  '; exp=' num2str(info.name.expT)  'ms; Gain='  num2str(info.name.Gain) 'dB' ];
else
    titleStr =  [ infoFields{1} '; exp=' num2str(info.name.expT)  'ms; Gain='  num2str(info.name.Gain) 'dB' ];
end
    
if plotFlag
        fig = figure('name',['SCOS ' recordName]);
    subplot(3,1,1);
        plot(timeVec,rawSpeckleVar)
        ylabel('Variance [DU^2]')
%         title({ rawName , ['p2p = ' ]})
        title(titleStr)
        subplot(3,1,2);
        plot(timeVec,rawSpeckleContrast)
        ylabel('Contrast (var/I^2)')
    subplot(3,1,3);
        plot(timeVec,imMeanVec);
        ylabel('I [DU]')
end

for plot_i=1:3
    subplot(3,1,plot_i);
    xlabel('Time [s]')
end
%% Save
if exist('fig','var')
    savefig(fig,[recSavePrefix 'Local' stdStr '_plot.fig']);
    save([recSavePrefix 'Local' stdStr '.mat'],'timeVec', 'corrSpeckleContrast' , 'rawSpeckleContrast','rawSpeckleVar','corrSpeckleVar', 'imMeanVec', 'info', 'recordName','windowSize');
end
 
