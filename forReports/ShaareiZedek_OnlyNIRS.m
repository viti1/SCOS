%% Init
mFolder = fileparts(mfilename('fullpath'));

recordsFolder = ['C:\Users\' getenv('username') '\OneDrive - Bar Ilan University\SCOS_Records\ShaareiZedek\NIRS Data']; 
recordsDir = [ dir([ recordsFolder '\A*.csv']); dir([ recordsFolder '\B*.csv']) ];
records = strcat(recordsFolder, filesep , {recordsDir.name})';

figureDir = fullfile(recordsFolder,'figs');

%% Create Figures
prefix = '\figs\';
figs = nan(size(numel(records)));
for ri = 1:numel(records)
    %%
    [~,recName ]  = fileparts(records{ri});
    participantID = recName(1:2);

    timingFile = [records{ri}(1:end-4) ' timing info.txt'];
    
    Timing = readtable(timingFile);
    startTime = Timing.Time{1};
    
    NIRS = readtable(records{ri});
    NIRS = renamevars(NIRS,2:7,{'Date','Time','Left','Right','SedlineEvents','Remarks'});
    if ~isnumeric(NIRS.Left)
        NIRS.Left  = str2double(NIRS.Left);        
    end
    if ~isnumeric(NIRS.Right)
        NIRS.Right = str2double(NIRS.Right);
    end

    figs(ri) = figure('Units','Normalized','Position',[0.05 0.1 0.7 0.5]);
    NIRStimeVec = minutes(NIRS.Time -  duration(startTime,'InputFormat','hh:mm'));
    % colors = {[0 0.4 0],[0.5 0.8 0.2]}; % green
    colors = {[1 0 0.8],[0.6 0 0.1]}; % red
    h(1) = plot(NIRStimeVec,NIRS.Left,'color',colors{1},'DisplayName','Left'); hold on;
    h(2) = plot(NIRStimeVec,NIRS.Right,'color',colors{2},'DisplayName','Right'); hold on;
    markTiming(timingFile);
    legend(h); % TBD check which is left and which is right
    xlabel('[min]');
    ylabel('rSO2 % ');
    xlims = xlim;
    endTime = minutes(duration(Timing.Time{end},'InputFormat','hh:mm') - duration(Timing.Time{1},'InputFormat','hh:mm'));
    xlim( [ min( NIRStimeVec(find(~isnan(NIRS.Left) & ~isnan(NIRS.Right),1)) , 0 ) , ...
            max( NIRStimeVec(find(~isnan(NIRS.Left) & ~isnan(NIRS.Right),1,'last')) , endTime*1.05 ) ] );
    title(recName,'interpreter','none');
    set(gca,'OuterPosition',[-0.1,0,1.2,1]);
    savefig(figs(ri),[figureDir '\NIRS ' recName '.fig']);
    saveas(figs(ri),[figureDir '\NIRS ' recName '.png']);
end

%% save to PPT
import mlreportgen.ppt.*

% Start PowerPoint
ppt = Presentation([ figureDir '/NIRS Summary.pptx' ]);
open(ppt);

% List of figure files
figureFiles = dir(fullfile(figureDir, '*.png'));

titleSlide = add(ppt,'Title Slide');
replace(titleSlide, 'Title', 'NIRS Results During Surgery');
replace(titleSlide, 'Subtitle', 'Saarei Zedek 2024');

% Loop through each figure file and add it to the presentation
for i = 1:length(figureFiles)
    % Add a new slide
    slide = add(ppt,'Title and Picture');
    % Insert the figure
    imagePath = fullfile(figureDir, figureFiles(i).name);
    replace(slide,'Picture',Picture(imagePath));
    
    participantID = figureFiles(i).name(6:8);
    replace(slide,'Title',participantID);
end

% Close PowerPoint
close(ppt);
rptview(ppt);