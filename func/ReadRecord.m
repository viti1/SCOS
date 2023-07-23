%-------------------------------------------------------
% Input : full path of folder with .tiff or .tif or .avi files, or a file
% with  .tif/.tiff extention
% assuming gray scale image
%-------------------------------------------------------
function [Rec, parameters_names, parameters_values, parameters_units , info] = ReadRecord(file_or_folder,MaxNOfFrames,parameters_names,parameters_expected_units)
    if ~exist(file_or_folder,'file')
        error(['path ' file_or_folder ' don''t exist'])
    end

    if exist(file_or_folder,'file') == 7 % it's a folder
        folderpath = file_or_folder;
        % find all .tiff or .tif files
        tiff_files = [ dir([folderpath, '\*.tiff']) , dir([folderpath, '\*.tif']) ];
        avi_files  = dir([folderpath, '\*.avi']) ;
        if isempty(tiff_files) && isempty(avi_files)
            error(['There is no ''.tiff'' or ''.avi'' files in folder ''' folderpath '''']);
        elseif ~isempty(tiff_files) && ~isempty(avi_files)
            error([ '''' folderpath ''' contains both ''.tiff'' or ''.avi'' files. It must contain only one the the file types ' ]);        
        elseif ~isempty(avi_files)
            if numel(avi_files) > 1
                error([ folderpath ' must contain only one .avi file but contains ' numel(avi_files) ' files.']);
            end
            if exist('MaxNOfFrames','var')
                [Rec, avi_info] = Avi2Matrix( fullfile(folderpath,avi_files.name) ,MaxNOfFrames);
            else
                [Rec, avi_info] = Avi2Matrix( fullfile(folderpath,avi_files.name) );
            end
        elseif ~isempty(tiff_files)
            if exist('MaxNOfFrames','var')
                nOfFrames = min(MaxNOfFrames,numel(tiff_files));
            else
                nOfFrames = numel(tiff_files);
            end

            % get first image in order to find out the image size
            t = Tiff(fullfile(folderpath,tiff_files(1).name),'r');   
            im = read(t);
            close(t);

            % read all images
            Rec = nan(size(im,1),size(im,2),nOfFrames);
            for k = 1:nOfFrames
                t = Tiff(fullfile(folderpath,tiff_files(k).name),'r');
                Rec(:,:,k) = read(t);
                close(t);        
            end
        end
    else    
        [folderpath, ~ , ext] = fileparts(file_or_folder);
        if strcmp(ext,'.avi')
            if exist('MaxNOfFrames','var')
                [Rec, avi_info] = Avi2Matrix( file_or_folder ,MaxNOfFrames );
            else
                [Rec, avi_info] = Avi2Matrix( file_or_folder  );
            end
        elseif strcmp(ext,'.tiff') || strcmp(ext,'.tif')
            t = Tiff(file_or_folder,'r');
            Rec = read(t);
            close(t);
        else
            error(['Unsupported file type ' ext  ' . Supported types are .tif .tiff .avi '])
        end
    end

    % === Extract Parameters from record path =============================
    % set default parameters names, if not passed to the function 
    if ~exist('parameters_names','var')
        parameters_names = {'Tint','FrameRate','Gain','f'};
        parameters_expected_units = {'ms','Hz','','mm'};
    end

    sp = strsplit(folderpath,filesep); foldername = sp{end};
    [parameters_values, parameters_units] = ExtractParametersFromString(foldername,parameters_names);

    % check expected units
    if exist('expected_units','var')
        for k = 1:numel(parameters_names)
            if parameters_units{k} ~= parameters_expected_units{k} 
               error(['Parameter "' parameters_names '" units are expected to be "' parameters_expected_units{k} '" but appear as "' parameters_units{k} '"'])
            end
        end
    end

    info.FirstString = findFirstSubstring(foldername);
    for k = 1:numel(parameters_names)
        info.(parameters_names{k}).val =  parameters_values(k);
        info.(parameters_names{k}).units =  parameters_units{k};
    end

    % should add in case of avi file that FrameRate is correct
    if exist('avi_info','var')
        if ismember('FrameRate',parameters_names) && ~isnan(info.FrameRate.val)
            if avi_info.FrameRate ~= info.FrameRate.val 
                error('Wrong frame rate in the name of the file')
            end 
        elseif  ismember('FR',parameters_names) && ~isnan(info.FR.val)
            if avi_info.FrameRate ~= info.FR.val 
                error('Wrong frame rate in the name of the file')
            end 
        else
            info.FrameRate.val  = avi_info.FrameRate; % both names are appleacable
            info.FR.val         = avi_info.FrameRate;
            info.FR.units       = 'Hz';
            info.FrameRate.units= 'Hz';
        end
    end
end


% Convert parameters_values and parameters_units into struct
function substring = findFirstSubstring(str)
% Find the first underscore in the string
underscoreIndex = strfind(str, '_');

if isempty(underscoreIndex)
    % No underscore found
    substring = str;
else
    % Return substring from the beginning until the underscore
    substring = str(1:underscoreIndex(1)-1);
end
end