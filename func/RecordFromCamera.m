function [ rec, filename ] = RecordFromCamera(nOfFrames,Tint,gain,frameRate,blackLevel,folder,saveFormat,prefix,suffix,vid,videoFormat,forceWrite)
%% [ rec, filename ] = RecordFromCamera(nOfFrames,Tint,gain,frameRate,blackLevel,folder,saveFormat,prefix,suffix,vid,videoFormat,forceWrite)
% Input : 
%   nOfFrames - desired number of frames. default = 1
%   Tint - integration (=exposure) time in ms
%   gain - analog camera gain in dB
%   frameRate - frames per second [Hz]
%   blackLevel - Offset of the image ( in [DU] ). For bester should be in range of [0,30]
%   folder - folder to which file should be saved. The parent folder should exist.
%            if not passed of is empty , record will not be saved.
%            Name of the file will contain all the parameters.
%   saveFormat - '.avi' or '.tiff' or '.mat'.  default = '.mat'
%                in case of .mat format, the file will start with 'Rec_' prefix
%   prefix - prefix for the filname
%   suffix - suffix for the filname 
%   vid    - preopened video object
%   videoFormat - "Mono8" or "Mono12" ( for additional formats - need to update the code )
%   forceWrite  - save recording even if recording with the same name already exist
%
%   * metadata is saved automatically in .mat file in variable called 'src' 
%     for '.mat' format  - inside the record matfile
%     for '.tiff' format - inside 'folder', file named 'src.mat'
%     for '.avi' format  - inside 'folder', the same name as the recording, but adeed '_src.mat' suffix instead of '.avi'
%
%   * any of the input parameters can be empty. In case of Tint,gain,frameRate,blackLevel - the parameters previously set to the camera are useed.
% Output :
%   rec      - 3D uint8 or uint16 (depending of the videoFormat) matrix
%   filename - output full file name (including path)
% ``````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````

%% open video object and set record parameters into video object
createVid= ~exist('vid','var') || isempty(vid);
if ~exist('videoFormat','var') || isempty(videoFormat)
    if ~createVid
        videoFormat = vid.VideoFormat;
    else
        videoFormat = "Mono8";
    end
elseif ~ismember({'Mono8','Mono12'},videoFormat)
    error(['Video format "' videoFormat '" is not supported. Must be "Mono8" or "Mono10"'])
end

if createVid
    vid = videoinput("gentl", 1, videoFormat);
end

src = getselectedsource(vid);
vid.FramesPerTrigger = Inf; 

if exist('Tint','var') && ~isempty(Tint)
    src.ExposureTime = Tint*1000;
end

if exist('gain','var') && ~isempty(gain)
    src.Gain = gain;
end

if exist('blackLevel','var') && ~isempty(blackLevel)
    src.BlackLevel = blackLevel;
end

if exist('frameRate','var') && ~isempty(frameRate)
    src.AcquisitionFrameRate = frameRate;
end

if ~exist('nOfFrames','var') || isempty(nOfFrames)
    nOfFrames = 1;
end

if ~exist('saveFormat','var') || isempty(saveFormat)
    saveFormat = '.mat';
end

if ~exist('forceWrite','var') || isempty(forceWrite)
    forceWrite = false;
end

%% 

start(vid);
while(~vid.FramesAvailable)
    pause(0.001)
end
% imagesBuff = getdata(vid, vid.FramesAvailable);
imagesBuff = getdata(vid, 1); % try this

oneFrame = squeeze( (imagesBuff(:,:,:,1)) );
rec = uint8(zeros(size(oneFrame,1),size(oneFrame,2),nOfFrames));


fig = figure;
imageAx = axes(fig); k=1;
imageAx.XLim = [0 size(oneFrame,2)];
imageAx.YLim = [0 size(oneFrame,1)];
axis equal
fprintf('Recording %s_Tint%gms_Gain%gdB_FR%gHz_BL%gDU%s ...\n',prefix,Tint,gain,frameRate,blackLevel,suffix);

while k <= nOfFrames 
        %currImage = getsnapshot(v);
        if vid.FramesAvailable
            fprintf('%d\t',k);
            if mod(k,50) == 0; fprintf('\n'); end
                
            currImagesBuff = getdata(vid, vid.FramesAvailable);
            buffSize = size(currImage,3);
            if strcmp(videoFormat,'Mono8')
                currImage = uint8(squeeze(currImagesBuff));
            elseif strcmp(videoFormat,'Mono12')
                currImage = uint16(squeeze(currImagesBuff));
            end
            rec(:,:,k:k+buffSize-1) = currImage;

            imagesc(currImage(:,:,1), 'Parent', imageAx); colormap gray;
%             axis equal
            title({sprintf('FPS=%.3g, Exposure=%.3gms, Gain=%.2gdB',frameRate,Tint,gain),sprintf('frame %d',k)})
            k = k + buffSize;
        end
end
fprintf('\n');
stop(vid);

if createVid
    delete(vid)
end

%% Save recording
if exist('folder','var') && ~isempty(folder)
    % -- check folder
    if exist(fileparts(folder),'dir') ~= 7
        error('The parent folder of the destination folder must exist');
    end
    if ~exist(folder,'dir')
        mkdir(folder);
    end
    
    % -- set recording name 
    recordName = sprintf('%sTint%gms_Gain%gdB_FR%gHz_BL%gDU%s',prefix,Tint,gain,frameRate,blackLevel,suffix);
    switch saveFormat
        case '.mat'
            recordName = ['Rec_' recordName '.mat'];
        case '.tiff'
            % check that folder is empty
            recordName = fullfile(folder,recordName);
        case '.avi'
            recordName = [recordName '.avi'];            
        otherwise
            error('wrong saveFormat, must be .tiff or .mat. or .avi')
    end
    filename = fullfile(folder,recordName);    

    % -- Check if recording already exist
    if strcmp(saveFormat,'.tiff')
        if exist(recordName,'dir') ~= 7
            mkdir(recordName)
        elseif ~isempty(dir([recordName '\*.tiff']))
            if forceWrite
                delete([recordName '\*.tiff'])
            else
                answer = userintput(['"' recordName '" Record exist. Do you want to rewrite it [y/n] ? '],"s");
                if ~strcmpi(answer,'y') % the answer is no (or any other)
                    disp('Aborting ...')
                    return
                else
                    delete([recordName '\*.tiff'])
                    disp(['Saving ' recordName ' ...'])
                end
            end
        end
    elseif ~forceWrite && exist(recordName,'file')        
        answer = userintput(['"' recordName '" Record exist. Do you want to rewrite it [y/n] ? '],"s");
        if ~strcmpi(answer,'y')
            disp('Aborting ...')
            return
        else
            disp(['Saving ' recordName ' ...'])
        end
    end
       
    % -- Save
    switch saveFormat
        case '.mat'
            save(fullfile(folder,recordName),'rec','src');
        case '.tiff'
            % prepare target struct
            tagstruct.ImageLength = size(rec,1);
            tagstruct.ImageWidth  = size(rec,2); 
            
            if strcmp(videoFormat,"Mono8")
                tagstruct.BitsPerSample = 8;
                tagstruct.SamplesPerPixel = 1;
            elseif strcmp(videoFormat,"Mono12")
                tagstruct.BitsPerSample = 16;
                tagstruct.SamplesPerPixel = 1;
            else 
                error('Unsupported format for writing .tiffs')
            end
            tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
            tagstruct.Software = 'MATLAB';
            tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
            
            %tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
            for k=1:size(rec,3)
                t = Tiff([destFolder,sprintf('\\frame%0*d.tiff',3,k)],'w');
                setTag(t,tagstruct);
                write(t,uint8(rec(:,:,k)));
                close(t);
            end
            save(fullfile(destFolder,'src.mat'),'src');
        case '.avi'
            v = VideoWriter( recordName ,'Grayscale AVI');
            open(v)
            if startsWith(videoFormat,"Mono")
                recToWrite = reshape(rec,[size(rec,1) size(rec,2) 1 size(rec,3) ]);
            else
                recToWrite = rec;
            end
            writeVideo(v,recToWrite)
            close(v)
            save([recordName(1:end-4) '_src.mat'],'src');            
        otherwise
            error('wrong saveFormat, must be .tiff or .mat. or .avi')
    end
    
end
end