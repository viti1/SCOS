function varargout = SCOS_GUI(varargin)
% SCOS_GUI MATLAB code for SCOS_GUI.fig
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SCOS_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @SCOS_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before SCOS_GUI is made visible.
function SCOS_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SCOS_GUI (see VARARGIN)

% Choose default command line output for SCOS_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SCOS_GUI wait for user response (see UIRESUME)
% uiwait(handles.fig_SCOS_GUI);

% ============ Start Code ===================================
addpath([fileparts(mfilename('fullpath')) '\baseFunc']);
handles.fig_SCOS_GUI.UserData.calcSCOS_Flag = false;
handles.fig_SCOS_GUI.UserData.scosData = nan(1,10000);
handles.fig_SCOS_GUI.UserData.scosTime = nan(1,10000);
handles.fig_SCOS_GUI.UserData.resetSCOS_Flag = true;
handles.fig_SCOS_GUI.UserData.calcSCOS_lastClock = 0;
handles.fig_SCOS_GUI.UserData.calcSCOS_beforeVideoStop_Flag = false;
handles.fig_SCOS_GUI.UserData.isVideoRunning = false;
handles.fig_SCOS_GUI.UserData.recFolder = [ fileparts(fileparts(mfilename('fullpath'))) '\Records\timeData' ];
if ~exist(handles.fig_SCOS_GUI.UserData.recFolder,'dir'); mkdir(handles.fig_SCOS_GUI.UserData.recFolder); end
handles.fig_SCOS_GUI.UserData.recName = [ handles.fig_SCOS_GUI.UserData.recFolder '\Rec_' strrep(char(datetime()),':','_') '.mat']; 
handles.fig_SCOS_GUI.UserData.saveEachNframes = 200;

% --- Outputs from this function are returned to the command line.
function varargout = SCOS_GUI_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on Calc Contrast start button press 
function btn_start_Callback(hObject, eventdata, handles)
if isequal(hObject.String,'Start SCOS') % 
    if handles.fig_SCOS_GUI.UserData.calcSCOS_Flag
        return; % do nothing - > theoretically not suppose to get here
    end
    handles.fig_SCOS_GUI.UserData.calcSCOS_Flag = true;
    handles.fig_SCOS_GUI.UserData.resetSCOS_clock = tic;
    set(hObject,'String','Stop SCOS')
else
    handles.fig_SCOS_GUI.UserData.calcSCOS_Flag = false;
    handles.fig_SCOS_GUI.UserData.calcSCOS_lastClock = handles.fig_SCOS_GUI.UserData.scosTime(find(isnan(handles.fig_SCOS_GUI.UserData.scosTime),1)-1)  ; 
    set(hObject,'String','Start SCOS');
    
    if ~isempty(handles.fig_SCOS_GUI.UserData.scosTime) && any(~isnan(handles.fig_SCOS_GUI.UserData.scosTime))
        scosData = handles.fig_SCOS_GUI.UserData.scosData;
        scosTime = handles.fig_SCOS_GUI.UserData.scosTime;
        save(handles.fig_SCOS_GUI.UserData.recName,'scosData','scosTime');
    end
    % handles.fig_SCOS_GUI.UserData.resetSCOS_clock = tic;
end

% --- Executes on button reset Contrast graph .
function btn_reset_Callback(hObject, eventdata, handles)
handles.fig_SCOS_GUI.UserData.resetSCOS_Flag = true;
plot(handles.ax_scos,0,0); % reset axes
handles.fig_SCOS_GUI.UserData.resetSCOS_clock = tic;
handles.fig_SCOS_GUI.UserData.calcSCOS_lastClock = 0;

function edt_exposureTime_Callback(hObject, eventdata, handles)
% hObject    handle to edt_exposureTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_exposureTime as text
%        str2double(get(hObject,'String')) returns contents of edt_exposureTime as a double
handles.fig_SCOS_GUI.UserData.src.ExposureTime = str2double(hObject.String)*1000; % should be in microseconds
% stop(handles.fig_SCOS_GUI.UserData.vid);
% start(handles.fig_SCOS_GUI.UserData.vid);

% --- Executes during object creation, after setting all properties.
function edt_exposureTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_exposureTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edt_gain_Callback(hObject, eventdata, handles)
% hObject    handle to edt_gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_gain as text
%        str2double(get(hObject,'String')) returns contents of edt_gain as a double
handles.fig_SCOS_GUI.UserData.src.Gain = str2double(hObject.String); 



% --- Executes during object creation, after setting all properties.
function edt_gain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edt_triggerDelay_Callback(hObject, eventdata, handles)
% hObject    handle to edt_triggerDelay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_triggerDelay as text
%        str2double(get(hObject,'String')) returns contents of edt_triggerDelay as a double
triggerDelayNum = str2double(hObject.String);
if isnan(triggerDelayNum)
    errordlg('Trigger Delay should be a number!');
    return;
end

if ~isvalid(handles.fig_SCOS_GUI.UserData.src)
    %stop(handles.fig_SCOS_GUI.UserData.vid)
    %delete(handles.fig_SCOS_GUI.UserData.vid);
    vid = videoinput("gentl", 1, "Mono8");
    src = getselectedsource(vid);
    src.TriggerDelay = striggerDelayNum;
    delete(vid);
    clear vid
else
    handles.fig_SCOS_GUI.UserData.src.TriggerDelay = triggerDelayNum; % should be in microseconds
end


% --- Executes during object creation, after setting all properties.
function edt_triggerDelay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_triggerDelay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in chk_externalTrigger.
function chk_externalTrigger_Callback(hObject, eventdata, handles)
    % -- set the external trigger field        
    if ~isvalid(handles.fig_SCOS_GUI.UserData.src)
        vid_was_valid = isvalid(handles.fig_SCOS_GUI.UserData.vid);
        if ~vid_was_valid
            vid = videoinput("gentl", 1, "Mono8");
        else
            vid = handles.fig_SCOS_GUI.UserData.vid;
        end
        src = getselectedsource(vid);

        if hObject.Value
            set( src, 'TriggerMode', 'on');
            src.TriggerMode = 'on';
            set( src, 'AcquisitionFrameRateEnable' , 'False'); % not sure it's really needed
        else
            set( src, 'TriggerMode', 'off');
            set( src, 'AcquisitionFrameRateEnable' , 'True'); 
            frameRate = get(src,'AcquisitionFrameRate');
        end
        if ~vid_was_valid
            delete(vid);
        end
    else
        if hObject.Value
            set( handles.fig_SCOS_GUI.UserData.src,'TriggerMode' , 'on'); 
            set( handles.fig_SCOS_GUI.UserData.src, 'AcquisitionFrameRateEnable','False');
        else
            set( handles.fig_SCOS_GUI.UserData.src ,'TriggerMode' ,'off' );
            set( handles.fig_SCOS_GUI.UserData.src, 'AcquisitionFrameRateEnable','True');  
            frameRate = get( handles.fig_SCOS_GUI.UserData.src,'AcquisitionFrameRate');
        end
    end   
    
    % --- change the FrameRate field accordingly 
    if hObject.Value        
        if isfield(handles.fig_SCOS_GUI.UserData,'last_external_trigger_frequency')
            default_ans = handles.fig_SCOS_GUI.UserData.last_external_trigger_frequency;
        else
            default_ans = '';
        end
        
        answer = inputdlg('What is the frequency of external trigger? [Hz]' , '', [1 45], {default_ans});
        if isnan(str2double(answer{1})); errordlg('frequency must be a number');
            set(hObject,'Value',0)
            return;
        end
        
        if str2double(answer{1})<=0; errordlg('frequency must be a positive number')
            set(hObject,'Value',0)
            return;            
        end
        handles.fig_SCOS_GUI.UserData.last_external_trigger_frequency = answer{1};
        
        set(handles.txt_frameRate,'String','External Trigger Frame Rate [Hz]')
        set(handles.edt_frameRate,'Value',str2double(answer{1}));
        set(handles.edt_frameRate,'String',answer{1})
    else % returned back to be camera inner trigger
        set(handles.txt_frameRate,'String','Frame Rate [Hz]')
        set(handles.edt_frameRate,'Value', frameRate );
        set(handles.edt_frameRate,'String', num2str(frameRate ))        
    end
    

function edt_frameRate_Callback(hObject, eventdata, handles)

frameRateNum = str2double(hObject.String);
if isnan(frameRateNum)
    errordlg('Frame Rate must be a number');
    set(handles.edt_frameRate,'String',num2str(hObject.Value));
    return;
end

set(handles.edt_frameRate,'Value',frameRateNum);

if handles.chk_externalTrigger.Value % we are working with external trigger
    handles.fig_SCOS_GUI.UserData.last_external_trigger_frequency = frameRate;
else     % we are working with internal camera frame rate 
    if isvalid(handles.fig_SCOS_GUI.UserData.src)
        handles.fig_SCOS_GUI.UserData.src.AcquisitionFrameRateEnable = 'True';
        handles.fig_SCOS_GUI.UserData.src.AcquisitionFrameRate = frameRateNum;
    else
        vid = videoinput("gentl", 1, "Mono8");
        src = getselectedsource(vid);
        src.AcquisitionFrameRateEnable = 'True';
        src.AcquisitionFrameRate = frameRateNum;
        delete(vid);
        %btn_startVideo_Callback(handles.btn_startVideo, eventdata, handles);
    end
end
% --- Executes during object creation, after setting all properties.
function edt_frameRate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_frameRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_updateROI.
function btn_updateROI_Callback(hObject, eventdata, handles)
% hObject    handle to btn_updateROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% stop Video
% handles.fig_SCOS_GUI.UserData.stopVideoFlag = true;

answer = questdlg('How do you want to choose ROI','Update ROI','Auto','Manual','Auto');
hData = findobj(handles.ax_image.Children,'Type','image');
im = double(get(hData,'Cdata'));


if strcmp(answer,'Auto')   
    x = 1 : size(im, 2); % Columns.
    y = 1 : size(im, 1); % Rows.
    [X, Y] = meshgrid(x, y);
    meanIm = mean(im(:));
    centerOfMassX = mean(im(:) .* X(:)) / meanIm;
    centerOfMassY = mean(im(:) .* Y(:)) / meanIm;
    % figure; imshowpair(CreateCircleMask(size(im),5,centerOfMassY,centerOfMassX),im)
    
    % Find MaxI
    maxI = mean(im(CreateCircleMask(size(im),5,centerOfMassY,centerOfMassX)));
    
    % Find Expected Radius
    relativeIforCircle = 0.3;
    numOfPixelsInSpot = nnz(im(:) > maxI*relativeIforCircle);
    approxRadius = sqrt(numOfPixelsInSpot/pi());
    
    [centers, radii] = imfindcircles(im > maxI*relativeIforCircle,round(approxRadius*[0.9 1.1]),'ObjectPolarity','bright','Sensitivity',0.99);
    if numel(radii) > 1 
        [circ.Radius,max_ind] = max(radii);
        circ.Center = centers(:,max_ind);
    elseif numel(radii) == 1 
        circ.Radius = radii;
        circ.Center = centers;
    else
        errordlg('Could not detect circle :(')
        return;
    end
    if isfield(handles.fig_SCOS_GUI.UserData,'ROI')
        delete(handles.fig_SCOS_GUI.UserData.ROI.circ_handle);
    end

    mask = false(size(im,1),size(im,2));
    [x,y] = meshgrid(1:size(im,2),1:size(im,1));
    mask((x-circ.Center(1)).^2 + (y-circ.Center(2)).^2 < circ.Radius^2 ) = true;
    
else
    axis(handles.ax_image);
    if isfield(handles.fig_SCOS_GUI.UserData,'ROI')
        delete(handles.fig_SCOS_GUI.UserData.ROI.circ_handle);
    end

    c_tmp = drawcircle('Color','r','FaceAlpha',0.2);
    circ.Center = c_tmp.Center;
    circ.Radius = c_tmp.Radius;
    delete(c_tmp);
    mask = false(size(im,1),size(im,2));
    [x,y] = meshgrid(1:size(im,2),1:size(im,1));
    mask((x-circ.Center(1)).^2 + (y-circ.Center(2)).^2 < circ.Radius^2 ) = true;
end


h = viscircles(handles.ax_image, circ.Center , circ.Radius,'Color','r','EnhanceVisibility',false,'LineWidth',1);
handles.fig_SCOS_GUI.UserData.ROI.mask = mask;
handles.fig_SCOS_GUI.UserData.ROI.circ = circ;
handles.fig_SCOS_GUI.UserData.ROI.circ_handle = h;



btn_AutoClim_Callback(handles.btn_AutoClim, eventdata, handles);

% --- Executes on button press in btn_record.
function btn_record_Callback(hObject, eventdata, handles)
% hObject    handle to btn_record (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
              
currFolder = fileparts(mfilename('fullpath'));
last_folder_matfile = [currFolder '\SCOS_GUI_lastrec.mat'];
if exist(last_folder_matfile,'file')
     tmp = load(last_folder_matfile);
     last_record_folder = tmp.recFolder;
     if ~exist(last_record_folder,'dir')
         last_record_folder = 0;
     end
else
    last_record_folder = 0;
end

if  last_record_folder == 0  
    last_record_folder = [fileparts(currFolder) filesep 'Records']; 
end

recFolder = uigetdir(last_record_folder,'Please Select Folder in which to save recording');
if recFolder==0
    return;
end
save(last_folder_matfile,'recFolder');

% stop running video and delete all opened vid
wasVideoRunning = handles.fig_SCOS_GUI.UserData.isVideoRunning;
if wasVideoRunning
    % stop video
    handles.fig_SCOS_GUI.UserData.stopVideoFlag = true;
    pause(0.02);
    if isvalid(handles.fig_SCOS_GUI.UserData.src) || isvalid(handles.fig_SCOS_GUI.UserData.vid)
        delete(handles.fig_SCOS_GUI.UserData.vid);
        objects = imaqfind;
        if numel(objects)>0;  delete(objects(1)); end
    end
end

if handles.chk_externalTrigger.Value
    %answer = inputdlg('What is the pulse generator frequency?','Pulse Generato frequence',[1 35]);
    frameRate = str2double(handles.edt_frameRate.String);
    if isnan(frameRate)
        errordlg('Wrong Frame Rate');
        return;
    end
else
    vid = videoinput("gentl", 1, "Mono8");
    src = getselectedsource(vid);
    frameRate = round(get(src,'AcquisitionFrameRate'),5);
    delete(vid);
end

nOfFrames = str2double(handles.edt_nOfFrames.String);
camParams.Gain = str2double(handles.edt_gain.String);
camParams.ExposureTime= str2double(handles.edt_exposureTime.String)*1000;

if handles.chk_externalTrigger.Value %triggerconfig(vid, 'hardware')
    camParams.TriggerMode = 'on';
else
    camParams.TriggerMode = 'off';
end

windowSize = str2double(handles.edt_windowSize.String); 

[~,recName]= RecordFromCamera(nOfFrames, camParams, [], recFolder ,'.tiff','',['FR' num2str(frameRate) 'Hz'],0,0); 
if isfield(handles.fig_SCOS_GUI.UserData,'ROI')
    SCOSvsTimeUpdated( recName ,windowSize, 1, handles.fig_SCOS_GUI.UserData.ROI)
else
    SCOSvsTimeUpdated( recName ,windowSize, 1  );
end

if wasVideoRunning
    set(handles.btn_startVideo,'String','Start Video');
    btn_startVideo_Callback(handles.btn_startVideo, eventdata, handles)
end


function edt_nOfFrames_Callback(hObject, eventdata, handles)
% hObject    handle to edt_nOfFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_nOfFrames as text
%        str2double(get(hObject,'String')) returns contents of edt_nOfFrames as a double


% --- Executes during object creation, after setting all properties.
function edt_nOfFrames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_nOfFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function btn_UpdateFolder_Callback(hObject, eventdata, handles)

function pushbutton10_Callback(hObject, eventdata, handles)

function pushbutton11_Callback(hObject, eventdata, handles)

function pushbutton12_Callback(hObject, eventdata, handles)

function ax_image_CreateFcn(hObject, eventdata, handles)        
    
% --- Executes on button press in btn_AutoClim.
function btn_AutoClim_Callback(hObject, eventdata, handles)
% hObject    handle to btn_AutoClim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 hData = findobj(handles.ax_image.Children,'Type','image');  
 im = get(hData,'CData');
 
 if isfield(handles.fig_SCOS_GUI.UserData,'ROI')
     mask = handles.fig_SCOS_GUI.UserData.ROI.mask;
 else
     mask = true(size(im));
 end
    
 [N,edges] = histcounts(im(mask),1000,'Normalization','cdf');
 
 upperLim = ceil(edges(find(N > 0.95,1)));
 lowerLim = floor(edges(find(N > 0.05,1)));
 minStretch = 6;
 if upperLim-lowerLim < minStretch
     lowerLim = max(0,mean([upperLim,lowerLim]) - minStretch/2) ;
     upperLim = upperLim + minStretch;
 end
 set(handles.ax_image,'CLim',[ lowerLim upperLim]);
        


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over btn_startVideo.
function btn_startVideo_ButtonDownFcn(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function btn_startVideo_CreateFcn(hObject, eventdata, handles)


function edt_windowSize_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edt_windowSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_windowSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_stopVideo.
function btn_stopVideo_Callback(hObject, eventdata, handles)
% hObject    handle to btn_stopVideo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles.fig_SCOS_GUI.UserData.vid,'vid')
    vid = handles.fig_SCOS_GUI.UserData.vid;
    stop(vid)
    delete(vid);
    clear(vid)
else
   objects = imaqfind;
   if numel(objects)>0;  delete(objects(1)); end
end

if ~isempty(handles.fig_SCOS_GUI.UserData.scosTime) && any(~isnan(handles.fig_SCOS_GUI.UserData.scosTime))
    scosData = handles.fig_SCOS_GUI.UserData.scosData;
    scosTime = handles.fig_SCOS_GUI.UserData.scosTime;
    save(handles.fig_SCOS_GUI.UserData.recName,'scosData','scosTime');
end

function ax_scos_CreateFcn(hObject, eventdata, handles)


% --- Executes on button press in btn_startVideo.
function btn_startVideo_Callback(hObject, eventdata, handles)
% hObject    handle to btn_startVideo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(get(hObject,'String'),'Stop Video')
    % stop calculating and plotting scos
    handles.fig_SCOS_GUI.UserData.calcSCOS_beforeVideoStop_Flag = handles.fig_SCOS_GUI.UserData.calcSCOS_Flag;
    if handles.fig_SCOS_GUI.UserData.calcSCOS_Flag  
        btn_start_Callback(handles.btn_start, eventdata, handles)
    end
    
    % Stop Video
    handles.fig_SCOS_GUI.UserData.isVideoRunning = false; % -> will stop the loop in previous instarnce of this function    

    % delete vid
    if isfield(handles.fig_SCOS_GUI.UserData.vid,'vid')
        vid = handles.fig_SCOS_GUI.UserData.vid;
        stop(vid)
        delete(vid);
        clear(vid)
    else
        vidobj = imaqfind;
        if numel(vidobj)>0;  stop(vidobj(1)); delete(vidobj(1)); end
    end
    
    % set(hObject,'String','Start Video') - > this will happen after video is stopped in other instance of the same function
else %  get(hObject,'String') == 'Start Video'    
    try
        if isfield(handles.fig_SCOS_GUI.UserData,'vid') && ~isvalid(handles.fig_SCOS_GUI.UserData.vid)
            delete(handles.fig_SCOS_GUI.UserData.vid);
        end

        vidobj = imaqfind;
        if ~isempty(vidobj); delete(vidobj(1)); end

        vid = videoinput("gentl", 1, "Mono8");
        vid.FramesPerTrigger = Inf; 
        handles.fig_SCOS_GUI.UserData.vid = vid;
    catch err
        if isequal(err.identifier,'imaq:videoinput:noDevices')
            errordlg('Connect camera then try again');
            return
        else
            errordlg({'Error While Starting videoinput ',err.message});
            return;
        end
    end
    handles.fig_SCOS_GUI.UserData.isVideoRunning = true;
    
    src = getselectedsource(vid);    
    handles.fig_SCOS_GUI.UserData.src = src;

    if strcmpi(src.TriggerMode,'On') %handles.chk_externalTrigger.Value
        triggerconfig(vid, 'hardware');
        % set(src,'TriggerMode','on');
    end

    % init camera params
    handles.edt_exposureTime.String = src.ExposureTime * 1e-3;
    handles.edt_gain.String = src.Gain;
    handles.edt_triggerDelay.String = src.TriggerDelay;
    handles.edt_frameRate.String = src.AcquisitionFrameRate;
    %handles.chk_externalTrigger.Value = strcmp(src.TriggerMode,'on');


    if handles.chk_externalTrigger.Value == 0
        handles.edt_frameRate.Visible = 1;
        handles.txt_frameRate.Visible = 1;
    end

    detectorFolder = [fileparts(fileparts(mfilename('fullpath'))) '\Records\NoiseAndBackground\Basler_1440GS_Vika01\Mono8'];
    dData = load([detectorFolder '\ReadNoise\vsGain\readNoiseVsGain.mat']);  % detector Data
   
    k = 1;    
    if ~isrunning(vid); start(vid); end

    pauseTriggerMode = 0.05;
    pauseUsualMode = 0.01;

    if strcmpi(src.TriggerMode,'On')
        pauseInterval = pauseTriggerMode;
    else
        pauseInterval = pauseUsualMode;
    end
    pause(pauseInterval);
    set(hObject,'String','Stop Video');
    scos_ind = find(isnan(handles.fig_SCOS_GUI.UserData.scosTime),1);
    if isempty(scos_ind)
        scos_ind = numel(handles.fig_SCOS_GUI.UserData.scosTime) + 1;        
    end
    currXlim = get(handles.ax_scos,'XLim');
    
    if handles.fig_SCOS_GUI.UserData.calcSCOS_beforeVideoStop_Flag
        btn_start_Callback(handles.btn_start, eventdata, handles)
    end
    while isvalid(handles.fig_SCOS_GUI) && handles.fig_SCOS_GUI.UserData.isVideoRunning  && isvalid(handles.fig_SCOS_GUI.UserData.src) && isvalid(vid) && vid.FramesAvailable  % VIKA TBD! 
        % ~handles.fig_SCOS_GUI.UserData.stopVideoFlag

        % --- Get one frame ------------------------
        currImagesBuff = getdata(vid, vid.FramesAvailable);
        %size(currImagesBuff,3)
        im = squeeze( currImagesBuff(:,:,end,end) );

        % --- Show in Axis -------------------------
        if k==1
            hImage = imagesc(handles.ax_image,im);
            [N,edges] = histcounts(im(:),1000,'Normalization','cdf');
            upperLim = ceil(edges(find(N > 0.99,1)));
            lowerLim = floor(edges(find(N > 0.01,1)));
            minStretch = 10;
            if upperLim-lowerLim < minStretch
                lowerLim  = mean([upperLim,lowerLim]) - minStretch/2 ;
                upperLim = upperLim + minStretch;
            end
            set(handles.ax_image,'CLim',[ lowerLim upperLim]);
            SetAxisEqual(handles.ax_image);
            colorbar;
            colormap gray;

            if isfield(handles.fig_SCOS_GUI.UserData,'ROI')
                circ = handles.fig_SCOS_GUI.UserData.ROI.circ;
                viscircles(handles.ax_image, circ.Center , circ.Radius,'Color','r','EnhanceVisibility',false,'LineWidth',1);
            end
        else
            set(hImage, 'CData', im);
        end
        k = k + 1;

        if isfield(handles.fig_SCOS_GUI.UserData,'ROI')
            mask = handles.fig_SCOS_GUI.UserData.ROI.mask;
        else
            mask = true(size(im)); %CreateCircleMask(size(im));
        end
        meanI = round(mean(im(mask)));
        pcntl = prctile(im(mask),[5 95]);
        set(handles.txt_avgI,'String',[ '<I>=' num2str(meanI) 'DU;   5%I  =' num2str(pcntl(1)) 'DU;   95%I =' num2str(pcntl(2)) 'DU  '   ]);        

        % --- Calc Scos ---------------------------
        if handles.fig_SCOS_GUI.UserData.calcSCOS_Flag
           if handles.fig_SCOS_GUI.UserData.resetSCOS_Flag
               handles.fig_SCOS_GUI.UserData.scosData = nan(1,10000);
               handles.fig_SCOS_GUI.UserData.scosTime = nan(1,10000);
               handles.fig_SCOS_GUI.UserData.resetSCOS_Flag = false;
               scos_ind=1;
%                timeJump = 0.5 ; %  min
               currXlim = [0 0.5]; 
           end
           if scos_ind > length(handles.fig_SCOS_GUI.UserData.scosData)
               handles.fig_SCOS_GUI.UserData.scosData = [ handles.fig_SCOS_GUI.scosData nan(1,10000)] ;
               handles.fig_SCOS_GUI.UserData.scosTime = [ handles.fig_SCOS_GUI.scosTime nan(1,10000)] ;
           end

           im = double(im);
           windowSize = str2double(handles.edt_windowSize.String);
           curr_gain = handles.fig_SCOS_GUI.UserData.src.Gain;
           readoutN   = interp1(dData.gainArr,dData.totNoise, curr_gain ,'spline');
           actualGain = ConvertGain(curr_gain,8,10.5e3);      

           stdIm = stdfilt(im,true(windowSize));
           meanIm = imfilter(im, true(windowSize)/windowSize^2,'conv','same');
           meanImSquare = meanIm.^2;
           Kraw = (stdIm.^2)./meanImSquare ;
           
           if mod(scos_ind,500) == 0
               disp(['Shot Noise^2 = ', num2str(actualGain*mean(meanIm(mask))) ])               
               disp(['Noise^2 = ' num2str(readoutN^2) ])               
               disp(['Var_raw = ' num2str(mean((stdIm(mask).^2))) ]);
               
%                stop(vid);
%                delete(vid);
%                %clear(vid)
%                 vid = videoinput("gentl", 1, "Mono8");
%                 vid.FramesPerTrigger = Inf;
%                 handles.fig_SCOS_GUI.UserData.vid = vid;
%                 src = getselectedsource(vid); 
%                 handles.fig_SCOS_GUI.UserData.src = src;
%                 
%                 if strcmpi(src.TriggerMode,'On') %handles.chk_externalTrigger.Value
%                     triggerconfig(vid, 'hardware');
%                 end
%                  
% pause(0.001);
%                 start(vid)
           end
           if mod(scos_ind, handles.fig_SCOS_GUI.UserData.saveEachNframes ) == 0 
               scosData = handles.fig_SCOS_GUI.UserData.scosData;
               scosTime = handles.fig_SCOS_GUI.UserData.scosTime;
               save(handles.fig_SCOS_GUI.UserData.recName,'scosData','scosTime');
           end

           handles.fig_SCOS_GUI.UserData.scosData(scos_ind) = mean(Kraw(mask) - actualGain./meanIm(mask) - ( readoutN^2 + 1/12)./meanImSquare(mask) );   % Kappa_corrected 
           %handles.fig_SCOS_GUI.UserData.scosData(scos_ind) =  mean( (stdIm(mask)- actualGain*meanIm(mask) - readoutN -1/12)./meanIm(mask) )^2  ;
           tm = toc(handles.fig_SCOS_GUI.UserData.resetSCOS_clock) ;

           handles.fig_SCOS_GUI.UserData.scosTime(scos_ind) =  tm/60 + handles.fig_SCOS_GUI.UserData.calcSCOS_lastClock  ; % /60 in order to convert to [min] from [sec]

           plot(handles.ax_scos,handles.fig_SCOS_GUI.UserData.scosTime(1:scos_ind),1./handles.fig_SCOS_GUI.UserData.scosData(1:scos_ind),'-');
           % set x and y axis labels
           ylabel(handles.ax_scos,'1/\kappa^2');
           xlabel(handles.ax_scos,'time [min]');
           set(handles.ax_scos,'XLim',currXlim  )
           
           % change XLim if it's not enough            
%            curr_xlim = get(handles.ax_scos,'XLim');
           if currXlim(2) <= handles.fig_SCOS_GUI.UserData.scosTime(scos_ind)
                currXlim(2) = 2*currXlim(2);
%                 timeJump = 2*curr_xlim(2);
%                 set(handles.ax_scos,'XLim',[0 timeJump*( mod(handles.fig_SCOS_GUI.UserData.scosTime(scos_ind),timeJump) + 1) ]); 
%                 set(handles.ax_scos,'XLim',[curr_xlim(1) (2*curr_xlim(2)) ] )
           end

           scos_ind = scos_ind + 1;
        end       

         %start(vid);
         pause(pauseInterval);

         % if still there is not frames -> while loop untill there are frames available
         temp_counter = 0;
         while  isvalid(handles.fig_SCOS_GUI) &&  handles.fig_SCOS_GUI.UserData.isVideoRunning && isvalid(vid) && temp_counter < Inf 
             if  ~isrunning(vid) || vid.FramesAvailable 
                 break;
             else
                 pause(pauseInterval/2); 
                 temp_counter = temp_counter + 1;
             end
         end
    %     trigger(vid);
    end

    if isvalid(handles.fig_SCOS_GUI);  handles.fig_SCOS_GUI.UserData.isVideoRunning = false; end
    if isvalid(vid) 
        stop(vid)
        delete(vid);
    end

    fprintf('\n')
    if isvalid(hObject); set(hObject,'String','Start Video'); end
end
             
% --- Executes when user attempts to close fig_SCOS_GUI.
function fig_SCOS_GUI_CloseRequestFcn(hObject, eventdata, handles)

handles.fig_SCOS_GUI.UserData.isVideoRunning = false;

if isfield(hObject.UserData,'vid') && isvalid(hObject.UserData.vid) && isrunning(hObject.UserData.vid)    
    stop(hObject.UserData.vid);
end
pause(0.065); % let it stop from the btn_startVideo_Callback
if isfield(hObject.UserData,'vid') && isvalid(hObject.UserData.vid)
    delete(hObject.UserData.vid);
end

%% Delete vid if there is one opened
objects = imaqfind;
if numel(objects)>0;  delete(objects(1)); end

%% Save
scosData = handles.fig_SCOS_GUI.UserData.scosData;
scosTime = handles.fig_SCOS_GUI.UserData.scosTime;
save(handles.fig_SCOS_GUI.UserData.recName,'scosData','scosTime');

%%
pause(0.1)
delete(hObject); %closes the figure


% --- Executes on button press in btn_open_recording.
function btn_open_recording_Callback(hObject, eventdata, handles)
    answer = questdlg('Open last saved recording?','','Open Last','Choose Recording','Open Last');
    if strcmp(answer,'Open Last')
        recFiles = dir([ handles.fig_SCOS_GUI.UserData.recFolder '\*.mat' ]);
        [~,last_rec_ind] = max([recFiles.datenum]);
        fileName = recFiles(last_rec_ind).name;
        
        D = load([handles.fig_SCOS_GUI.UserData.recFolder filesep fileName ]);
    else
        [fileName , filePath] =  uigetfile([ handles.fig_SCOS_GUI.UserData.recFolder '\*.mat' ]);
        if fileName == 0
            return;
        end
        D = load([filePath filesep fileName ]);
    end
    
    fig = figure('Name',[ 'Scos Result' fileName ]);
    plot(D.scosTime,1./D.scosData);
    ylabel('1/\kappa^2');
    xlabel('time [min]');
    title(['SCOS Results' fileName],'interpreter','none');
    

           
        
    