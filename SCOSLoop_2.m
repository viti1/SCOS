clear 
recordsFolder = [  'F:\nadav' ];
windowSize = 7;
plotFlag =1;
%%
recordsNames = dir([ recordsFolder '\*']);
recordsNames(~[recordsNames.isdir]) = [];
recordsNames(cellfun(@(x) ismember(x,{'.','..'}),{recordsNames.name})) = [];
recordsNames(cellfun(@(x) contains(x,'dark'),{recordsNames.name})) = [];
records = fullfile(recordsFolder,{ recordsNames.name })';

tmp = strsplit(recordsFolder,filesep);
% hirarchyForTitle = numel(tmp) - 2;
hirarchyForTitle = 1;
folderTitle = strjoin(tmp(end-hirarchyForTitle+1:end),' ');

timeVec = cell(1,numel(records));
corrSpeckleContrast = cell(1,numel(records));
meanVec = cell(1,numel(records));
rBFi = cell(1,numel(records));
for k=3:numel(records)
    [ timeVec{k}, ~ , corrSpeckleContrast{k}, meanVec{k} , rBFi{k}] = ...
     SCOSvsTime_WithNoiseSubtraction_Ver2(records{k},'',windowSize,plotFlag);    
end

%% Plot
figure('Units','Normalized','Position',[0 0.05 0.9 0.85])
N = numel(records);
for k=1:N
    subplot(N,1,k);
    plot(timeVec{k},rBFi{k});
    title(strrep(recordsNames(k).name,'_',' '));    
end



