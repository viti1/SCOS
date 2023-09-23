clear 
recordsFolder = '.\Records\Vika Head\DirectCamera';
windowSize = 9;

%%
recordsNames = dir([ recordsFolder '\*']);
recordsNames(~[recordsNames.isdir]) = [];
recordsNames(cellfun(@(x) ismember(x,{'.','..'}),{recordsNames.name})) = [];
recordsNames(cellfun(@(x) contains(x,'background'),{recordsNames.name})) = [];
records = fullfile(recordsFolder,{ recordsNames.name })';

tmp = strsplit(recordsFolder,filesep);
% hirarchyForTitle = numel(tmp) - 2;
hirarchyForTitle = 1;
folderTitle = strjoin(tmp(end-hirarchyForTitle+1:end),' ');

%%
% f1 = figure('name','STD Relative vs Time','Position',[150 50 1400 1000]);
f2 = figure('name','Variance vs Time','Position',[200 50 1400 1000]);
f3 = figure('name','Contrast vs Time','Position',[250 50 1400 1000]);
%
Nx = 1; Ny = numel(records);
n = 1;
max_t = 0;
winStr = sprintf('%dx%d',windowSize,windowSize);
varStr = sprintf('Var%dx%d',windowSize,windowSize);

for k=1:numel(records)
    recName =  recordsNames(k).name;
    
    % Calc
    [ timeVec{k}, rawSpeckleContrast{k} , rawSpeckleVar{k}, ~,~,IVec{k} ,info{k}] = SCOSvsTime(records{k},windowSize,1);

    % Plot Std
    figure(f2); subplot(Ny,Nx,n);
    plot(timeVec{k},rawSpeckleVar{k});
%     title_curr = ['Fiber ' num2str(info{k}.name.Fiber.val) '\mum;  Gain=' num2str(info{k}.name.Gain.val) ' ' strrep(recName(strfind(recName,'defocus'):end),'_',' ') ';  Imean=' num2str(round(mean(IVec{k}))) '[DU]'];
    title_curr = [ 'SDS=' num2str(info{k}.name.SDS) '; exp=' num2str(info{k}.name.Exp) 'ms; Gain=' num2str(info{k}.name.Gain) 'dB ' strrep(recName(strfind(recName,'defocus'):end),'_',' ') ';  Imean=' num2str(round(mean(IVec{k}))) '[DU]'];
    if n==1
        title({folderTitle,title_curr} );
    else
        title(title_curr);
    end
    xlabel('time [s]'); ylabel( [varStr ' [DU]']); grid on;
    
    
    % Plot Std/<I>
    figure(f3); subplot(Ny,Nx,n);
    plot(timeVec{k},rawSpeckleContrast{k});
    title(title_curr);
    xlabel('time [s]'); ylabel([varStr ' / <I>^2 ']); grid on;    
%     ylim(mean(rawSpeckleContrast{k}) + 0.01*[-1 1] )
    
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
%% Save
savefig(f2,[recordsFolder sprintf('\\Variance%dx%d vs Time ',windowSize,windowSize)]);
savefig(f3,[recordsFolder sprintf('\\Contrast%dx%dRelI vs Time ',windowSize,windowSize)]);
matFile = [recordsFolder '\SCOS_' winStr '.mat'];
save(matFile,'timeVec', 'rawSpeckleContrast' , 'rawSpeckleVar' ,'IVec' ,'info','windowSize','folderTitle','recordsNames');
%%  Just plot
f2 = figure('name',[ 'Variance vs Time ' folderTitle ],'Position',[200 50 1400 1000]);
f3 = figure('name',[ 'Contrast vs Time ' folderTitle ],'Position',[250 50 1400 1000]);

Nx = 1; Ny = numel(recordsNames);
n = 1;
max_t = 0;
winStr = sprintf('%dx%d',windowSize,windowSize);

for k=1:numel(recordsNames)
    recName =  recordsNames(k).name;
    
    % Calc

    % Plot Std
    figure(f2); subplot(Ny,Nx,n);
    plot(timeVec{k},rawSpeckleVar{k});
%     title_curr = ['Fiber ' num2str(info{k}.name.Fiber.val) '\mum;  Gain=' num2str(info{k}.name.Gain.val) ' ' strrep(recName(strfind(recName,'defocus'):end),'_',' ') ';  Imean=' num2str(round(mean(IVec{k}))) '[DU]'];
    title_curr = [ 'exp=' num2str(info{k}.name.Exp) 'ms; Gain=' num2str(info{k}.name.Gain) 'dB ' strrep(recName(strfind(recName,'defocus'):end),'_',' ') ';  Imean=' num2str(round(mean(IVec{k}))) '[DU]'];
    if n==1
        title({folderTitle,title_curr} );
    else
        title(title_curr);
    end
    xlabel('time [s]'); ylabel( ['Var' winStr ' [DU]']); grid on;
    
    
    % Plot Std/<I>
    figure(f3); subplot(Ny,Nx,n);
    plot(timeVec{k},rawSpeckleContrast{k});
    if n==1
        title({folderTitle,title_curr} );
    else
        title(title_curr);
    end
    xlabel('time [s]'); ylabel(['Var' winStr ' / <I>^2 ']); grid on;    
%     ylim(mean(rawSpeckleContrast{k}) + 0.01*[-1 1] )
    
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