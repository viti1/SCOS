clear 
recordsFolder = '.\Records\Nisan Head\Nisan Head x20 Obj noLens 25mm SDS';
windowSize = 9;

%%
recordsNames = dir([ recordsFolder '\*']);
recordsNames(cellfun(@(x) ismember(x,{'.','..'}),{recordsNames.name})) = [];
recordsNames(cellfun(@(x) contains(x,'background'),{recordsNames.name})) = [];
records = fullfile(recordsFolder,{ recordsNames.name })';

tmp = strsplit(RecordsFolder,filesep);
% hirarchyForTitle = numel(tmp) - 2;
hirarchyForTitle = 1;
folderTitle = strjoin(tmp(end-hirarchyForTitle+1:end),' ');

%%
% f1 = figure('name','STD Relative vs Time','Position',[150 50 1400 1000]);
f2 = figure('name','STD vs Time','Position',[200 50 1400 1000]);
f3 = figure('name','STD Relative vs Time','Position',[250 50 1400 1000]);
%
Nx = 1; Ny = numel(records);
n = 1;
max_t = 0;
stdStr = sprintf('Std%dx%d',windowSize,windowSize);
for k=1:numel(records)
    recName =  recordsNames(k).name;
    
    % Calc
%     [ timeVec{k}, speckleNoiseVec{k} , speckleNoiseRelVec{k}, IVec{k} ,info{k}] = PlotSCOSvsTime(records{k},windowSize);

    % Plot Std
    figure(f2); subplot(Ny,Nx,n);
    plot(timeVec{k},speckleNoiseVec{k});
    title_curr = ['Fiber ' num2str(info{k}.Fiber.val) '\mum;  Gain=' num2str(info{k}.Gain.val) ' ' strrep(recName(strfind(recName,'defocus'):end),'_',' ') ';  Imean=' num2str(round(mean(IVec{k}))) '[DU]'];
    if n==1
        title({folderTitle,title_curr} );
    else
        title(title_curr);
    end
    xlabel('time [s]'); ylabel( [stdStr ' [DU]']); grid on;
    
    % Plot Std/<I>
    figure(f3); subplot(Ny,Nx,n);
    plot(timeVec{k},speckleNoiseRelVec{k});
    title(title_curr);
    xlabel('time [s]'); ylabel([stdStr ' / <I> ']); grid on;    
    ylim(mean(speckleNoiseRelVec{k}) + 0.01*[-1 1] )
    
    % increment
    n=n+1;
    if max_t< timeVec{k}(end)
        max_t = timeVec{k}(end);
    end
end
%%
max_t = round(max_t);
for n=1:Ny
    figure(f2); subplot(Ny,Nx,n);
    xlim([0 max_t]);
    figure(f3); subplot(Ny,Nx,n);
    xlim([0 max_t]); 
end
%%
savefig(f2,[recordsFolder sprintf('\\Std%dx%d vs Time ',windowSize,windowSize)]);
savefig(f3,[recordsFolder sprintf('\\Std%dx%dRelI vs Time ',windowSize,windowSize)]);