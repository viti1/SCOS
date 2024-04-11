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
    SCOSvsTimeUpdated(recordName,windowSize,plotFlag,maskInput)
if nargin <3
    plotFlag = true;
end

addpath('.\baseFunc\')
%% Constants
timePeriodForP2P = 2; % [s]

%% Check input parameters
if nargin == 0 % GUI mode
    plotFlag = 1;
    if exist('.\lastRec.mat','file')
        lastF = load('.\lastRec.mat');        
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

%% 
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

%% Read Record
disp(['Reading Record "' recordName '" ... '])
[ head_rec , info] = ReadRecord(recordName);

%%  Calc Specle Contrast
disp(['Calculation SCOS on "' recordName '" ... '])

roi_half_size = ceil((circ.Radius+windowSize)) ;
roi = [ round(circ.Center(2)) + (-roi_half_size:roi_half_size) ;
        round(circ.Center(1)) + (-roi_half_size:roi_half_size) ];

nOfFrames = size(head_rec,3);

%% Load rellevant parameters from camera characterizaion file 
detectorFolder = [fileparts(fileparts(mfilename('fullpath'))) '\Records\NoiseAndBackground\Basler_1440GS_Vika01\Mono8'];
detectorName = 'Basler_1440GS_Vika01';
dData = load(['.\camerasData\' detectorName '\readNoiseVsGain.mat']);  % detector Data
readoutN   = interp1(dData.gainArr,dData.totNoise, info.name.Gain ,'spline');


cut_mask        = mask(roi(1,:)  , roi(2,:));
[rawSpeckleVar , rawSpeckleContrast , corrSpeckleVar , corrSpeckleContrast , imMeanVec] = InitNaN([nOfFrames 1]);

actualGain = ConvertGain(info.name.Gain,8,10.5e3);

for i=1:nOfFrames
    im = head_rec(:,:,i);
    cut_im   = im(  roi(1,:)  , roi(2,:));
   
    stdIm = stdfilt(cut_im,true(windowSize));  
    meanIm = imfilter(cut_im, true(windowSize)/windowSize^2,'conv','same'); 

    rawSpeckleVar(i)       = mean(stdIm(cut_mask))^2;
    rawSpeckleContrast(i)  = mean(stdIm(cut_mask)./meanIm(cut_mask))^2;
    corrSpeckleVar(i)      = mean((stdIm(cut_mask)- actualGain*meanIm(cut_mask) - readoutN -1/12))^2  ; 
    corrSpeckleContrast(i) = mean((stdIm(cut_mask)- actualGain*meanIm(cut_mask) - readoutN -1/12)./meanIm(cut_mask))^2  ;
    
    imMeanVec(i) = mean(cut_im(cut_mask)); % mean(meanIm);
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
stdStr = sprintf('%dx%d',windowSize,windowSize);
% ax1_exist_flag = ax1==1; % since ax1 and ax2 is rewritten to be an axes pointers, the satate beforehand need to be saved
% ax2_exist_flag = ax2==1;
% Ny = (ax1==1) + (ax2==1) + (ax3==1);
% 
% if ax1~=0 
%     if ax1==1 
%         fig = figure('position',[80 80    1124    250*Ny],'name',[rawName , ':  Std vs time' ]); %#ok<NASGU>
%         if ax2==1
%             ax1 = subplot(Ny,1,1);
%         else            
%             ax1 = axes;
%         end
%     end
%     p2pStd = max(speckleNoiseVec(p2p_time)) - min( speckleNoiseVec(p2p_time)) ; % measure peak-to-peak in the first seconds 
%     
%     plot(ax1,timeVec,speckleNoiseVec);
%     title({rawName, [stdStr ' vs Time (p2p= ' num2str(p2pStd,2) '[DU] )']})
%     xlabel('time [s]'); ylabel('<Std>[DU]'); grid on;
% end
% 
% if ax2~=0 
%     if ax2==1 
%         if ~exist('fig','var')
%             fig = figure('position',[80 80    1124    250*Ny],'name',[rawName , ':  Std vs time' ]);
%             ax2 = axes;
%         else        
%             ax2 = subplot(Ny,1,2);        
%         end
%     end
%     
%     p2pStdRel = max(speckleNoiseRelVec(p2p_time)) - min( speckleNoiseRelVec(p2p_time)) ; % measure peak-to-peak in the first 2 seconds
%     plot(ax2,timeVec,speckleNoiseRelVec);
%     if ~ax1_exist_flag
%         title({rawName, [ stdStr ' relative vs Time ( p2p = ' num2str(p2pStdRel,2) ' )' ]});
%     else
%         title([ stdStr ' relative vs Time ( p2p = ' num2str(p2pStdRel,2) ' )'] );
%     end
%     xlabel('time [s]'); ylabel('<Std>/<I>'); grid on;
% %     ylim(round(mean(speckleNoiseRelVec),2) + 0.015*[-1 1] )     
% end
% 
% if ax3~=0
%     if ax3==1 
%         if ~exist('fig','var')
%             fig = figure('position',[80 80    1124    250*Ny],'name',[rawName , ':  Std vs time' ]);
%             ax3 = axes;
%         else        
%             ax3 = subplot(Ny,1,Ny);        
%         end
%     end
% 
%     plot(ax3,timeVec,imMeanVec);     
%     p2pIntensity = max(imMeanVec(p2p_time)) - min( imMeanVec(p2p_time)) ; % measure peak-to-peak in the first 2 seconds 
%     if ~ax2_exist_flag && ~ax2_exist_flag
%         title({rawName, ' Intensity vs Time ( p2p = ' num2str(p2pIntensity,2) '[DU] )' });
%     else
%         title([' Intensity vs Time ( p2p = ' num2str(p2pIntensity,2) '[DU] , p2p/<<I>>=' num2str(p2pIntensity/IMean,2) ')'  ] );
%     end
%     xlabel('time [s]'); ylabel('<I> [DU]'); grid on;
%     grid on
% end
%

%% plot5
% if plotFlag
%     fig = figure('name',['SCOS ' recordName]);
%     subplot(5,1,1);
%         plot(timeVec,rawSpeckleVar)
%         ylabel('Raw Variance [DU^2]')
%         title({ rawName , ['p2p = ' ]})
%         subplot(5,1,2);
%         plot(timeVec,rawSpeckleContrast)
%         ylabel('Raw Contrast (var/I^2)')
%     subplot(5,1,3);
%         plot(timeVec,corrSpeckleVar)
%         ylabel('Corrected Variance [DU^2]')
%     subplot(5,1,4);
%         plot(timeVec,corrSpeckleContrast)
%         ylabel('Corrected Contrast (var/I^2)')
%     subplot(5,1,5);
%         plot(timeVec,imMeanVec);
%         ylabel('I [DU]')
% end

% for plot_i=1:5
%     subplot(5,1,plot_i);
%     xlabel('Time [s]')
% end
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
    
    fig2 = figure('name',['SCOS ' recordName]);
    subplot(3,1,1);
        plot(timeVec,corrSpeckleVar)
        ylabel('var(I) [DU^2]')
%         title({ rawName , ['p2p = ' ]})
        title({titleStr,'Corrected Variance'});
    subplot(3,1,2);
        plot(timeVec,corrSpeckleContrast);
        title('Corrected Contrast');
        ylabel('var(I)/I^2)')
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
    savefig(fig2,[recSavePrefix 'CorrContrast' stdStr '_plot.fig']);
    save([recSavePrefix 'Local' stdStr '.mat'],'timeVec', 'corrSpeckleContrast' , 'rawSpeckleContrast','rawSpeckleVar','corrSpeckleVar', 'imMeanVec', 'info', 'recordName','windowSize');
end
 
