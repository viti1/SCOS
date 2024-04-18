%  ---------------------------------------------------------------------------------------------------------
%  [ timeVec,  , rawSpeckleContrast , rawSpeckleVar, imMeanVec , info] = SCOSvsTime_WithoutNoiseSubtraction (recordName,windowSize,plotFlag,maskInput)
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
    SCOSvsTime_WithoutNoiseSubtraction(recordName,windowSize,plotFlag,maskInput)

%% add path to baseFunc folder 
addpath('.\baseFunc\')

%% Define Constants
timePeriodForP2P = 2; % [s]

%% Check input parameters
if nargin <3
    plotFlag = true;
end

if nargin == 0 % GUI mode
    plotFlag = 1;
    if exist('.\lastRec.mat','file')
        lastF = load('.\lastRec.mat');        
    else
        lastF.recordName = '';
    end
    
    [recordName] = uigetdir(fileparts(lastF.recordName));
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
    save('.\lastRec.mat','recordName')

    
    maxWindowSize = 50; minWindowSize = 3;
    answer = inputdlg('Window Size','',[1 25],{'9'});
    windowSize = str2double(answer{1});
    if isnan(windowSize) || windowSize > maxWindowSize || windowSize < minWindowSize
        errordlg(['Window Size must be a number between ' num2str(minWindowSize)  ' and ' num2str( num2str(maxWindowSize) ) ]);
        error(['Window Size must be a number between ' num2str(minWindowSize)  ' and ' num2str( num2str(maxWindowSize) ) ])
    end
    
    clear answer    
end

%% Get Mask
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
%     elseif islogical(maskInput)
%         if ~isequal( size(maskInput), size(first_frame) ) 
%             error('wrong maskInput');
%         end
%         mask = maskInput;
    elseif isstruct(maskInput)
        if isfield(maskInput,'mask') && isfield(maskInput,'circ')
            save(maskFile,'-struct','maskInput');
            mask = maskInput.mask;
            circ = maskInput.circ;
        else
           loadExistingFile_flag = false;
        end    
    else
        error('wrong maskInput type');
    end
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
        % ask the user to mark a circle on the image
        [mask,circ ] = GetROI(first_frame);
        save(maskFile,'mask','circ');
    else 
        load(maskFile) %#ok<LOAD>
    end
end

%%  Calc Specle Contrast
disp(['Calculation SCOS on "' recordName '" ... '])
nOfFrames = GetNumOfFrames(recordName);
info = GetRecordInfo(recordName);

roi_half_size = ceil((circ.Radius+windowSize)) ;
roi = [ round(circ.Center(2)) + (-roi_half_size:roi_half_size) ;
        round(circ.Center(1)) + (-roi_half_size:roi_half_size) ];

cut_mask        = mask(roi(1,:)  , roi(2,:));
[rawSpeckleVar , rawSpeckleContrast , imMeanVec] = InitNaN([nOfFrames 1]);

for i=1:nOfFrames
    im = ReadRecord(recordName,1,i);
    if mod(i,50)==0; fprintf('%d\t\t',i); end

    cut_im   = im(  roi(1,:)  , roi(2,:));
   
    stdIm = stdfilt(cut_im,true(windowSize));  
    meanIm = imfilter(cut_im, true(windowSize)/windowSize^2,'conv','same'); 

    rawSpeckleVar(i)       = mean(stdIm(cut_mask).^2);
    rawSpeckleContrast(i)  = mean((stdIm(cut_mask)./meanIm(cut_mask)).^2);
    
    imMeanVec(i) = mean(cut_im(cut_mask)); % mean(meanIm);
end
IMean = mean(imMeanVec);
fprintf('\n');
%% Create Time vector
if ~isfield(info.name,'FR') || isnan(info.name.FR)
    error('Frame Rate must be part of the recording name as "FR"');
end
frameRate = info.name.FR; 
timeVec = (0:(nOfFrames-1))'*(1/frameRate) ;   % FR = FrameRate
p2p_time = timeVec<timePeriodForP2P;

%% Plot
stdStr = sprintf('%dx%d',windowSize,windowSize);
infoFields = fieldnames(info.name);
if isfield(info.name,'SDS')
    titleStr =  [ infoFields{1} ' SDS=' num2str(info.name.SDS)  '; exp=' num2str(info.name.expT)  'ms; Gain='  num2str(info.name.Gain) 'dB' ];
else
    titleStr =  [ infoFields{1} '; exp=' num2str(info.name.expT)  'ms; Gain='  num2str(info.name.Gain) 'dB' ];
end
      
if plotFlag
    fig1 = figure('name',['SCOS ' recordName]);
    subplot(3,1,1);
        plot(timeVec,rawSpeckleVar)
        ylabel('Raw Variance [DU^2]')
%         title({ rawName , ['p2p = ' ]})
        title({titleStr,'Raw Variance'});
    subplot(3,1,2);
        plot(timeVec,rawSpeckleContrast)
        ylabel('var(I)/I^2)')
        title('Raw Contrast');

    subplot(3,1,3);
        plot(timeVec,imMeanVec);
        ylabel('I [DU]')
        title('Intensity ');

    for plot_i=1:3
        subplot(3,1,plot_i);
        xlabel('Time [s]')
    end
end

%% Save
if exist('fig1','var')
    savefig(fig1,[recSavePrefix 'RawContrast' stdStr '_plot.fig']);
    save([recSavePrefix 'Local' stdStr '_RawOnly.mat'],'timeVec' , 'rawSpeckleContrast','rawSpeckleVar', 'imMeanVec', 'info', 'recordName','windowSize');
end
 
