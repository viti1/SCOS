folder =  'C:\Users\tarlevi\OneDrive - Bar-Ilan University - Students\PhD Research\SCOS\Project';
figsNames = dir([folder '\*.fig']);

SNR = nan([numel(figsNames),1]);
SigFreq = nan([numel(figsNames),1]);
for k=1:numel(figsNames)
   [all_XData,all_YData,all_titles] = getDataFromFigure( [folder filesep figsNames(k).name] );
   freq = all_XData{1}{1};
   spectrum = all_YData{1}{1};
   figure; plot(freq,spectrum); title(all_titles{1})
   
   % cut the low frequency
   spectrum(freq < 0.4 ) = [];
   freq(freq < 0.4 ) = [];
   
   % find max power
   [ Sig , maxIdx ] = max(spectrum);
   SigFreq(k) = freq(maxIdx);
   sigIdx = freq( maxIdx );  
   Noise = mean( spectrum( freq > 2.5 ) );
   SNR(k) = Sig/Noise;
   %parameters = ExtractParametersFromString(figsNames(k).name,{'Tint','offset','Gain','Focus'});
end


Gain = ExtractParametersFromString({figsNames.name}','Gain');
Focus = ExtractParametersFromString({figsNames.name}','focus');

tbl = table(Focus,Gain,SNR)
