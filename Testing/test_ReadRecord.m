clear
clc
runOnLocal = 1; % 1 or 0

%% Find tests Folder and test .mat file
[mFolder, mName ] = fileparts(mfilename('fullpath'));
functionName = mName(find(mName=='_',1)+1:end); %functionName ='ReadRecord';

testingFoldersList = jsondecode(fileread(fullfile(mFolder,'testsFolder.json'))); 
if runOnLocal
    testingFolder = testingFoldersList.local;   %#ok<UNRCH>
else
    testingFolder = testingFoldersList.global;  %#ok<UNRCH>
end

testingFolderFun = fullfile(testingFolder,functionName); % testing folder of the specific function 
if exist(testingFolderFun,'file')~=7
    error(['Testing Folder "' testingFolderFun '" does not exist!'])
end
casesStruct = load( fullfile( testingFolderFun ,['testScenatios_' functionName '.mat']) );

%% Run Over all Cases
passArr=true(1,numel(casesStruct.cases));
for test_i = 1:numel(casesStruct.cases)
    testCase = casesStruct.cases(test_i);
    disp([10 'Test Case Num ' num2str(test_i) ':'])
    disp(testCase.input);
    
    try
        currOutput = run_ReadRecord(testCase.input, testingFolderFun);
        if ~isequaln(currOutput,testCase.output)
            passArr(test_i) = false;
            warning(['Case ' num2str(test_i) ': OUTPUT is not as expected ! ']);
        else
            disp(' < PASSED >')
        end
    catch err
        if isfield(testCase.output,'error') 
            if isequaln(err.message,testCase.output.error.message) 
                disp('< PASSED >')
            else
                warning(['Case ' num2str(test_i) ': failed as expected but the error is wrong ! '])
                warning(['Expected Error : ' testCase.output.error.message ]);
                disp(testCase.output.error.stack)
                warning(['Actual   Error : ' err.message ])
                disp(err.stack)
            end
        else
            passArr(test_i) = false;
            warning(['Case ' num2str(test_i) ':  FAILED on running ! ']);
        end
    end
end

disp([10 '~~~~~~~~~~~~~~~~~~~~~~ Summary ~~~~~~~~~~~~~~~~~~~~~~~~~' ])
if all(passArr)
    disp([ functionName ' : Congradulations - All Tests Passed !' ])
else
    warning([ functionName '() function' 10 ' test passed : ' num2str(find(passArr)) 10 ' test failed : ' num2str(find(~passArr)) ])
end