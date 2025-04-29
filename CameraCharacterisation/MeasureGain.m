addpath('.\baseFunc')
shotNoiseFolder ='C:\SCOS\GainMeasurment\Lavi\Mono12_Gain8';
nBits = ExtractParametersFromString(shotNoiseFolder,'Mono');;
gainDB  = ExtractParametersFromString(shotNoiseFolder,'Gain');
SN = num2str(ExtractParametersFromString(shotNoiseFolder(20:end),'SN'));
records = dir([shotNoiseFolder '\*_I*']);
darkRecFile = dir([shotNoiseFolder '\*_Dark*']);
BlackLevel = ExtractParametersFromString(darkRecFile.name,'BL');
expT = ExtractParametersFromString(darkRecFile.name,'expT');

darkIm = mean( ReadRecord([shotNoiseFolder filesep darkRecFile.name ],68) , 3) - BlackLevel;

[IPerRec , spPerRec , spPerRec_std] = InitNaN([1,numel(records)]);
for ri = 1:numel(records)
    disp(records(ri).name)
    recName = fullfile(shotNoiseFolder,records(ri).name);
    rec = ReadRecord(recName,100);
    nFrames = size(rec,3);
    mean(rec(:));
        
    Ivec = nan(1,nFrames-1); spVec = nan(1,nFrames-1);  
    for i = 1:nFrames-1
       Ivec(i) = mean2(rec(:,:,i) - darkIm - BlackLevel ); 
       im_diff = rec(:,:,i+1)-rec(:,:,i); %rec_diff(:,:,i);
       spVec(i) = var(im_diff(:))/2;
    end
    spPerRec(ri) = mean(spVec);
    spPerRec_std(ri) = std(spVec);
    IPerRec(ri) = mean(Ivec);
    figure; 
    subplot(2,1,1); plot(Ivec); title([ records(ri).name ],'interpreter','none'); ylabel('<I> [DU]')
    subplot(2,1,2); plot(spVec); ylabel('var [DU]')    
end
% clear rec
%%
fig = figure; 
errorbar(IPerRec,spPerRec,spPerRec_std,'bo'); xlabel('<I> [DN]'); ylabel('\sigma^2 [DN^2]');
title(['AnalogGain = ' num2str(gainDB) 'dB, Exposure Time = ' num2str(expT) 'ms '])
p = polyfit(IPerRec,spPerRec,1);
hold on;
plot([0 IPerRec], polyval(p,[0 IPerRec]),'--r');
ylims = ylim;
xlims = xlim;
G = p(1) % 0.0238 , % 5.8750 for 24dB 20ms 12bits % 5.8617 for 24dB 15ms 12bits

% Gcalc = 2^8/10.5e3 % 0.0244    
Gcalc = 2^nBits/10.5e3 * 10^(gainDB/20)  

text(xlims(2)*0.1 , ylims(2)* 0.8 ,{['Slope = ' num2str(G,8)],['CalcG = 2^{' num2str(nBits) '}/capacity * 10^{g[dB]/20}= ' num2str(Gcalc,8)]});


save([ shotNoiseFolder filesep 'Results_Gain' num2str(gainDB) 'dB_expT' num2str(expT) 'ms.mat' ],...
    'IPerRec','spPerRec','spPerRec_std','expT','nBits','gainDB','Gcalc','G','BlackLevel','SN')
savefig(fig,[ shotNoiseFolder filesep 'Results_Gain' num2str(gainDB) 'dB_expT' num2str(expT) 'ms.fig' ]);


    