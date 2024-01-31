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
handles.fig_SCOS_GUI.UserData.resetSCOS_Flag = false;
handles.fig_SCOS_GUI.UserData.calcSCOS_lastClock = 0;

handles.fig_SCOS_GUI.UserData.isVideoRunning = false;


% --- Outputs from this function are returned to the command line.
function varargout = SCOS_GUI_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on Calc Contrast start button press 
function btn_start_Callback(hObject, eventdata, handles)
if handles.fig_SCOS_GUI.UserData.calcSCOS_Flag
    return; % do nothing
end
handles.fig_SCOS_GUI.UserData.calcSCOS_Flag = true;
handles.fig_SCOS_GUI.UserData.resetSCOS_clock = tic;

% --- Executes on button reset Contrast graph .
function btn_reset_Callback(hObject, eventdata, handles)
handles.fig_SCOS_GUI.UserData.resetSCOS_Flag = true;
% if handles.fig_SCOS_GUI.UserData.stopVideoFlag
    plot(handles.ax_scos,0,0); % reset
% end
handles.fig_SCOS_GUI.UserData.resetSCOS_clock = tic;
handles.fig_SCOS_GUI.UserData.calcSCOS_lastClock = 0;

% --- Executes on button stop Contrast calc
function btn_stop_Callback(hObject, eventdata, handles)
handles.fig_SCOS_GUI.UserData.calcSCOS_Flag = false;
handles.fig_SCOS_GUI.UserData.calcSCOS_lastClock = handles.fig_SCOS_GUI.UserData.scosTime(find(isnan(handles.fig_SCOS_GUI.UserData.scosTime),1)-1); 
% handles.fig_SCOS_GUI.UserData.resetSCOS_clock = tic;

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
            src.TriggerMode = 'on';
            set( src, 'AcquisitionFrameRateEnable' , 'False'); % not sure it's really needed
        else
            set( src,'TriggerMode', 'off');
            set( src, 'AcquisitionFrameRateEnable' , 'True'); 
            frameRate = get(src,'AcquisitionFrameRate');
        end
        if ~vid_was_valid
            delete(vid);
            clear vid
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
        
        answer = inputdlg('What is the frequency [Hz] of external trigger? ' , '', [1 45], {default_ans});
        if isnan(str2double(answer{1})); errordlg('frequency must be a number');
            set(hObject,'Value',0)
            return;
        end
        
        if str2double(answer{1})<=0; errordlg('frequency must be a positive number')
            set(hObject,'Value',0)
            return;            
        end
        handles.fig_SCOS_GUI.UserData.last_external_trigger_frequency = answer{1};
        
%         set(handles.edt_frameRate,'Enable','off'); % VIKA TBD : in edt_frameRate Callback add set(handles.edt_frameRate,'Value',___);
%                                                                 and set src only if we are not in external trigger mode
%                                                                 and change also handles.fig_SCOS_GUI.UserData.last_external_trigger_frequency
%                                                                 in calc SCOS we can take this value always
        set(handles.txt_frameRate,'String','External Trigger Frame Rate')
        set(handles.edt_frameRate,'Value',str2double(answer{1}));
        set(handles.edt_frameRate,'String',answer{1})
    else % returned back to be camera inner trigger
        set(handles.txt_frameRate,'String','Frame Rate')
        set(handles.edt_frameRate,'Value', frameRate );
        set(handles.edt_frameRate,'String', num2str(frameRate ))        
    end
    

function (hObject, eventdata, handles)
% hObject    handle to edt_frameRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_frameRate as text
%        str2double(get(hObject,'String')) returns contents of edt_frameRate as a double
if isvalid(handles.fig_SCOS_GUI.UserData.src)
    
    handles.fig_SCOS_GUI.UserData.src.AcquisitionFrameRateEnable = 'True';
    handles.fig_SCOS_GUI.UserData.src.AcquisitionFrameRate = str2double(hObject.String);
else
    vid = videoinput("gentl", 1, "Mono8");
    src = getselectedsource(vid);
    set( handles.fig_SCOS_GUI.UserData.src, 'AcquisitionFrameRateEnable','True');
    src.AcquisitionFrameRate = str2double(hObject.String);
    delete(vid);
    btn_startVideo_Callback(handles.btn_startVideo, eventdata, handles);
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


handles.fig_SCOS_GUI.UserData.stopVideoFlag = true;
pause(0.02);
if isvalid(handles.fig_SCOS_GUI.UserData.src)
    delete(handles.fig_SCOS_GUI.UserData.vid);
    objects = imaqfind;
    if numel(objects)>0;  delete(objects(1)); end
end
       
        
% handles.fig_SCOS_GUI.UserData.vid
currFolder = fileparts(mfilename('fullpath'));
last_folder_matfile = [currFolder '\SCOS_GUI_lastrec.mat'];
if exist(last_folder_matfile,'file')
     tmp = load(last_folder_matfile);
     last_record_folder = tmp.recFolder;
else
    last_record_folder = 0;
end

if  last_record_folder == 0  
    last_record_folder = [fileparts(currFolder) filesep 'Records']; 
end

recFolder = uigetdir(last_record_folder,'Please Select Folder in which to save recording');
save(last_folder_matfile,'recFolder');

if handles.chk_externalTrigger.Value
    answer = inputdlg('What is the pulse generator frequency?','Pulse Generato frequence',[1 35]);
    frameRate = str2double(answer);
else
    vid = videoinput("gentl", 1, "Mono8");
    src = getselectedsource(vid);
    frameRate = get(src,'AcquisitionFrameRate');
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
        
% --- Executes when user attempts to close fig_SCOS_GUI.
function fig_SCOS_GUI_CloseRequestFcn(hObject, eventdata, handles)

if isfield(hObject.UserData,'vid') && isvalid(hObject.UserData.vid)
    if isrunning(hObject.UserData.vid)
        stop(hObject.UserData.vid);
    end
    delete(hObject.UserData.vid);
end

objects = imaqfind;
if numel(objects)>0;  delete(objects(1)); end

delete(hObject); %closes the figure



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


% --- Executes on button press in btn_startVideo.
function btn_startVideo_Callback(hObject, eventdata, handles)
% hObject    handle to btn_startVideo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(get(hObject,'String'),'Stop Video')
    handles.fig_SCOS_GUI.UserData.isVideoRunning = false;
    if isfield(handles.fig_SCOS_GUI.UserData.vid,'vid')
        vid = handles.fig_SCOS_GUI.UserData.vid;
        stop(vid)
        delete(vid);
        clear(vid)
    else
        vidobj = imaqfind;
        if numel(vidobj)>0;  stop(vidobj(1)); delete(vidobj(1)); end
    end
    set(hObject,'String','Run Video')
else %  get(hObject,'String') == 'Run Video'    
    try
        if isfield(handles.fig_SCOS_GUI.UserData,'vid') && ~isvalid(handles.fig_SCOS_GUI.UserData.vid)
            delete(handles.fig_SCOS_GUI.UserData.vid);
        end

        vidobj = imaqfind;
        if ~isempty(vidobj); delete(vidobj(1)); end

        vid = videoinput("gentl", 1, "Mono8");
    %     vid.FramesPerTrigger =1;
        vid.FramesPerTrigger = Inf; 
        % triggerconfig(vid, 'immediate')

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

    % if succeded -> cannot start again
    set(hObject,'enable','off');

    src = getselectedsource(vid);
    % src.TriggerMode = 'On';
    src.TriggerMode
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

    k = 1; scos_ind = 1;
    handles.fig_SCOS_GUI.UserData.scosData = nan(1,10000);
    handles.fig_SCOS_GUI.UserData.scosTime = nan(1,10000);
    if ~isrunning(vid); start(vid); end

    pauseTriggerMode = 0.5;
    pauseUsualMode = 0.01;

    if strcmpi(src.TriggerMode,'On')
        pauseInterval = pauseTriggerMode;
    else
        pauseInterval = pauseUsualMode;
    end
    pause(pauseInterval);
    set(hObject,'String','Stop Video');
    while handles.fig_SCOS_GUI.UserData.isVideoRunning  && isvalid(handles.fig_SCOS_GUI.UserData.src) && isvalid(vid) && vid.FramesAvailable  % VIKA TBD! 
        % ~handles.fig_SCOS_GUI.UserData.stopVideoFlag

        % --- Get one frame ------------------------
        currImagesBuff = getdata(vid, vid.FramesAvailable);
        im = squeeze( currImagesBuff(:,:,1,end) );

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
        set(handles.txt_avgI,'String',[ '<I>=' num2str(meanI) 'DU'] );

        % --- Calc Scos ---------------------------
        if handles.fig_SCOS_GUI.UserData.calcSCOS_Flag
           if handles.fig_SCOS_GUI.UserData.resetSCOS_Flag
               handles.fig_SCOS_GUI.UserData.scosData = nan(1,10000);
               handles.fig_SCOS_GUI.UserData.scosTime = nan(1,10000);
               handles.fig_SCOS_GUI.UserData.resetSCOS_Flag = false;
               scos_ind=1;
           end
           if scos_ind > size(handles.fig_SCOS_GUI.UserData.scosData)
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
           handles.fig_SCOS_GUI.UserData.scosData(scos_ind) =  mean((stdIm(mask)- actualGain*meanIm(mask) - readoutN -1/12)./meanIm(mask))^2  ;
           tm = toc(handles.fig_SCOS_GUI.UserData.resetSCOS_clock) ;

           handles.fig_SCOS_GUI.UserData.scosTime(scos_ind) = tm + handles.fig_SCOS_GUI.UserData.calcSCOS_lastClock;

           plot(handles.ax_scos,handles.fig_SCOS_GUI.UserData.scosTime(1:scos_ind),handles.fig_SCOS_GUI.UserData.scosData(1:scos_ind),'-');
           scos_ind = scos_ind + 1;
           timeJump = 10; %set
           set(handles.ax_scos,'XLim',[0 (tm - mod(tm,timeJump) + timeJump)]) 
           if k==1
               ylabel(handles.ax_scos,'var(I)/<I>');
               xlabel(handles.ax_scos,'time [ms]')
           end       
        end       

         %start(vid);
         pause(pauseInterval);

         temp_counter = 0;
         while handles.fig_SCOS_GUI.UserData.isVideoRunning && isvalid(vid) && ~vid.FramesAvailable && temp_counter < 6 
             pause(pauseInterval/2); 
             temp_counter = temp_counter + 1;
         end
    %     trigger(vid);
    end

    handles.fig_SCOS_GUI.UserData.isVideoRunning = false;
    if isvalid(vid) 
        stop(vid)
        delete(vid);
        clear(vid)
    end

    fprintf('\n')
    set(hObject,'String','Run Video');
end