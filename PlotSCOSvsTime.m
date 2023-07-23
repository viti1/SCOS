%  ---------------------------------------------------------------------------------------------------------
%  [ timeVec, speckleNoiseVec , speckleNoiseRelVec, IMean , info] = PlotSCOSvsTime(recordName,windowSize,ax1,ax2)
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
%               if ax1==1 Figure with plot of local Std (in the selected region) vs time is created.
%               if ax1 is an axis object, then the plot is created on this axes.
%               if ax2==1 Figure with plot of local Std (in the selected region)/I vs time is created (if ax1==1 then
%               the plot is added to the same figure)
%               if ax2 is an axis object, then the plot is created on this axes.
%               if ax1 or ax2 is 0 (false) then correcponding plot is not presented (defalt is 0)
%  ---------------------------------------------------------------------------------------------------------

function [ timeVec, speckleNoiseVec , speckleNoiseRelVec, imMeanVec , info] = PlotSCOSvsTime(recordName,windowSize,ax1,ax2,ax3)

%% Constants
timePeriodForP2P = 2; % [s]

%% Check input parameters
if nargin == 0 % GUI mode
    [recordName] = uigetdir();
    if recordName == 0; return; end % if 'Cancel' was pressed
    if numel(dir([recordName, '\*.avi' ])) > 1 
        [recordRawName, recordDir] = uigetfile([recordName '\*.avi']);
        if recordRawName == 0; return; end % if 'Cancel' was pressed
        recordName = fullfile(recordDir, recordRawName);        
    elseif ( numel(dir([recordName, '\*.avi' ])) + numel(dir([recordName, '\*.tiff' ])) + numel(dir([recordName, '\*.tif' ])) ) < 1 
        errordlg(['No .avi or .tiff/.tif files found in ' recordName ])
        error(['No .avi or .tiff/.tif files found in ' recordName ]);
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
    
    ax1=1;
    ax2=1;
    ax3=1;
end

if ~exist('ax1','var') || isempty('ax1')
    ax1 = 0;
elseif ~islogical(ax1) && ~( isnumeric(ax1) &&  ismember(ax1,[0,1]))  &&  ~( isgraphics(ax1) && strcmp(ax1.Type,'axes') )
    error('wrong input for ax1 - must be boolean or 0/1 or handle to exes');
end
if ~exist('ax2','var') || isempty('ax2')
    ax2 = 0;
elseif ~islogical(ax2) && ~( isnumeric(ax2) &&  ismember(ax2,[0,1]))  &&  ~( isgraphics(ax2) && strcmp(ax2.Type,'axes') ) 
    error('wrong input for ax2 - must be boolean or 0/1 or handle to exes');
end

if ~exist('ax3','var') || isempty('ax3')
    ax3 = 0;
elseif ~islogical(ax3) && ~( isnumeric(ax3) &&  ismember(ax3,[0,1]))  &&  ~( isgraphics(ax3) && strcmp(ax3.Type,'axes') ) 
    error('wrong input for ax2 - must be boolean or 0/1 or handle to exes');
end

%% Create Mask
upFolders = strsplit(recordName,filesep);
rawName = strrep( strjoin(upFolders(end-2:end),'; '), '_',' ');

if exist(recordName,'file') == 7 % it's a folder
    recSavePrefix = [ recordName filesep ];
else
    recSavePrefix = [ recordName(1:find(recordName=='.',1,'last')-1) '_' ];
end
maskFile = [recSavePrefix 'Mask.mat'];

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
   
[ first_frame ] = ReadRecord(recordName,1);     
if ~loadExistingFile_flag
    %% ask the user to mark a circle on the image

    f = figure('position',[50,50,1200,800]); imagesc(first_frame); colormap gray; colorbar
    title(rawName,'interpreter','none')
    circ = drawcircle('Color','r','FaceAlpha',0.2);
    c.Center = circ.Center; c.Radius = circ.Radius;
    mask = false(size(first_frame,1),size(first_frame,2));
    [x,y] = meshgrid(1:size(first_frame,2),1:size(first_frame,1));
    mask((x-circ.Center(1)).^2 + (y-circ.Center(2)).^2 < circ.Radius^2 ) = true;
    save(maskFile,'mask','c');
    close(f)
else 
    load(maskFile)
%     figure; imshowpair(first_frame,mask);  
end

%% Read Record
disp(['Reading Record "' recordName '" ... '])
[ head_rec , ~ , ~ , ~ , info] = ReadRecord(recordName,Inf, {'Tint','FR','Gain','Fiber'});

%%  Calc Specle Contrast
disp(['Calculation SCOS on "' recordName '" ... '])

roi_half_size = ceil((c.Radius+windowSize));
roi = [ round(c.Center(2)) + (-roi_half_size:roi_half_size) ;
        round(c.Center(1)) + (-roi_half_size:roi_half_size) ];

windowSize = 9;
nOfFrames = size(head_rec,3);
speckleNoiseVec = nan(nOfFrames,1);
imMeanVec = nan(nOfFrames,1);
for i=1:nOfFrames
    im = head_rec(:,:,i);
    cut_im   = im(  roi(1,:)  , roi(2,:));
    cut_mask = mask(roi(1,:)  , roi(2,:));
    %figure; imshowpair(cut_im,cut_mask);
    stdIm = stdfilt(cut_im,true(windowSize));
    speckleNoiseVec(i) = mean(stdIm(cut_mask));
    imMeanVec(i) = mean(cut_im(cut_mask));
end
speckleNoiseRelVec = speckleNoiseVec./imMeanVec;
IMean = mean(imMeanVec);

%% Create Time vector
if ~isfield(info,'FR') || isnan(info.FR.val)
    error('Frame Rate must be part of the recording name as "FR"');
end
timeVec = (0:(nOfFrames-1))'*(1/info.FR.val) ;   % FR = FrameRate
p2p_time = timeVec<timePeriodForP2P;

%% Plot
stdStr = sprintf('Std%dx%d',windowSize,windowSize);
ax1_exist_flag = ax1==1; % since ax1 and ax2 is rewritten to be an axes pointers, the satate beforehand need to be saved
ax2_exist_flag = ax2==1;
Ny = (ax1==1) + (ax2==1) + (ax3==1);

if ax1~=0 
    if ax1==1 
        fig = figure('position',[80 80    1124    250*Ny],'name',[rawName , ':  Std vs time' ]); %#ok<NASGU>
        if ax2==1
            ax1 = subplot(Ny,1,1);
        else            
            ax1 = axes;
        end
    end
    p2pStd = max(speckleNoiseVec(p2p_time)) - min( speckleNoiseVec(p2p_time)) ; % measure peak-to-peak in the first seconds 
    
    plot(ax1,timeVec,speckleNoiseVec);
    title({rawName, [stdStr ' vs Time (p2p= ' num2str(p2pStd,2) '[DU] )']})
    xlabel('time [s]'); ylabel('<Std>[DU]'); grid on;
end

if ax2~=0 
    if ax2==1 
        if ~exist('fig','var')
            fig = figure('position',[80 80    1124    250*Ny],'name',[rawName , ':  Std vs time' ]);
            ax2 = axes;
        else        
            ax2 = subplot(Ny,1,2);        
        end
    end
    
    p2pStdRel = max(speckleNoiseRelVec(p2p_time)) - min( speckleNoiseRelVec(p2p_time)) ; % measure peak-to-peak in the first 2 seconds
    plot(ax2,timeVec,speckleNoiseRelVec);
    if ~ax1_exist_flag
        title({rawName, [ stdStr ' relative vs Time ( p2p = ' num2str(p2pStdRel,2) ' )' ]});
    else
        title([ stdStr ' relative vs Time ( p2p = ' num2str(p2pStdRel,2) ' )'] );
    end
    xlabel('time [s]'); ylabel('<Std>/<I>'); grid on;
%     ylim(round(mean(speckleNoiseRelVec),2) + 0.015*[-1 1] )     
end

if ax3~=0
    if ax3==1 
        if ~exist('fig','var')
            fig = figure('position',[80 80    1124    250*Ny],'name',[rawName , ':  Std vs time' ]);
            ax3 = axes;
        else        
            ax3 = subplot(Ny,1,Ny);        
        end
    end

    plot(ax3,timeVec,imMeanVec);     
    p2pIntensity = max(imMeanVec(p2p_time)) - min( imMeanVec(p2p_time)) ; % measure peak-to-peak in the first 2 seconds 
    if ~ax2_exist_flag && ~ax2_exist_flag
        title({rawName, ' Intensity vs Time ( p2p = ' num2str(p2pIntensity,2) '[DU] )' });
    else
        title([' Intensity vs Time ( p2p = ' num2str(p2pIntensity,2) '[DU] , p2p/<<I>>=' num2str(p2pIntensity/IMean,2) ')'  ] );
    end
    xlabel('time [s]'); ylabel('<I> [DU]'); grid on;
    grid on
end

if exist('fig','var')
    savefig(fig,[recSavePrefix 'Local' stdStr '_plot.fig']);
    save([recSavePrefix 'Local' stdStr '.mat'],'timeVec', 'speckleNoiseVec' , 'speckleNoiseRelVec', 'imMeanVec', 'info', 'recordName','windowSize');
end

