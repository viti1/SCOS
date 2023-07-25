clear 

% find testsFolderFunc
testsFolderFuncsList = jsondecode(fileread(fullfile(fileparts(mfilename('fullpath')),'testsFolder.json')));
testsFolderFunc  = [ testsFolderFuncsList.local filesep 'ReadRecord' ]; % create tests always in the local tets folder

% Create ReadRecord test Cases
recName1 = 'AviSingleFile_Tint3ms_Gain0_FR100Hz';
recName2 = 'Tiffs_f25mm_Tint0.1ms_Gain20_FrameRate20Hz';
recName3 = 'SeveralAviFiles\file1.avi';
recName4 = 'SeveralAviFiles';   % suppose to fail
recName5 = 'EmptyFolder';       % suppose to fail

n = 1;

% 1A. Avi File, all params
cases(n).input.file_or_folder = recName1;
cases(n).input.MaxNOfFrames = 10;
cases(n).input.parameters_names = {'Tint','Gain'};
cases(n).input.parameters_expected_units = {'ms',''};

cases(n).output = run_ReadRecord(cases(n).input, testsFolderFunc);
n = n + 1 ;

% 1B. Avi File, all params, Inf num of frames, wrong expected units
% cases(n).input.file_or_folder = recName1;
% cases(n).input.MaxNOfFrames = Inf;
% cases(n).input.parameters_names = {'Tint','Gain','X'};
% cases(n).input.parameters_expected_units = {'us','dB',''};
% 
% cases(n).output = run_ReadRecord(cases(n).input, testsFolderFunc);
% n= n + 1;
% 
% % 1C. Avi File, Inf num of frames, no parameters_expected_units
cases(n).input.file_or_folder = recName1;
cases(n).input.MaxNOfFrames = Inf;
cases(n).input.parameters_names = {'Tint','Gain','X'};

cases(n).output = run_ReadRecord(cases(n).input, testsFolderFunc);
n= n + 1;

% 1D. Avi File, 3 num of frames, no parameters_expected_units, no parameters_names
cases(n).input.file_or_folder = recName1;
cases(n).input.MaxNOfFrames = 3;

cases(n).output = run_ReadRecord(cases(n).input, testsFolderFunc);
n= n + 1;

% 1E. Avi File, 3 num of frames, no parameters_expected_units, no parameters_names
% cases(n).input.file_or_folder = recName1;
% cases(n).input.MaxNOfFrames = [];
% 
% cases(n).output = run_ReadRecord(cases(n).input,testsFolderFunc);
% n= n + 1;

% 1F. Avi File, no params
cases(n).input.file_or_folder = recName1;

cases(n).output = run_ReadRecord(cases(n).input,testsFolderFunc);
n= n + 1;

% 2A. Tiff Files, all Params
cases(n).input.file_or_folder = recName2;
cases(n).input.MaxNOfFrames = 3;
cases(n).input.parameters_names = {'Tint','Gain'};
cases(n).input.parameters_expected_units = {'ms',''};

cases(n).output = run_ReadRecord(cases(n).input,testsFolderFunc);
n = n + 1 ;

% 2B. Tiff Files, all Params
cases(n).input.file_or_folder = recName2;
cases(n).input.MaxNOfFrames = 20; % more than exist in folder
cases(n).input.parameters_names = {'Tint','Gain'};
cases(n).input.parameters_expected_units = {'ms',''};

cases(n).output = run_ReadRecord(cases(n).input,testsFolderFunc);
n = n + 1 ;

% 2C. Tiff Files, Inf frames
cases(n).input.file_or_folder = recName2;
cases(n).input.MaxNOfFrames = Inf;

cases(n).output = run_ReadRecord(cases(n).input,testsFolderFunc);
n = n + 1 ;

% 2C. Tiff Files, no frames
cases(n).input.file_or_folder = recName2;

% 3A. Single .Avi File
cases(n).input.file_or_folder = recName3;
cases(n).input.MaxNOfFrames = 5;

cases(n).output = run_ReadRecord(cases(n).input,testsFolderFunc);
n = n + 1 ;

% 3B. Single .Avi File Inf Frames
cases(n).input.file_or_folder = recName3;
cases(n).input.MaxNOfFrames = Inf;
cases(n).input.parameters_names = {'Tint','Gain'};

cases(n).output = run_ReadRecord(cases(n).input,testsFolderFunc);
n = n + 1 ;

% 3C. Single .Avi File no Frames
cases(n).input.file_or_folder = recName3;

cases(n).output = run_ReadRecord(cases(n).input,testsFolderFunc);
n = n + 1 ;

% 4. Several .avi files in folder
cases(n).input.file_or_folder = recName4;
cases(n).input.MaxNOfFrames = Inf;

try
    cases(n).output = run_ReadRecord(cases(n).input,testsFolderFunc);
catch err
    cases(n).output.error = err;
end
n = n + 1 ;

% 5. Empty Folder
cases(n).input.file_or_folder = recName5;
cases(n).input.MaxNOfFrames = Inf;

try
    cases(n).output = run_ReadRecord(cases(n).input,testsFolderFunc);
catch err
    cases(n).output.error = err;
end
n = n + 1 ;

%% Save .mat file with all scenarios
save(fullfile(testsFolderFunc,'testScenatios_ReadRecord.mat'),'cases','testsFolderFunc');

