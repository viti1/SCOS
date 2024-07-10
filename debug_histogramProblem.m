analysisFolder = 'C:\Users\tarlevi\OneDrive - Bar-Ilan University - Students\PhD Research\SCOS\Records\NoiseAndBackground\Basler_1440GS_Vika01\Mono8\FixedPattern';
recName = [ analysisFolder '\ShotNoise\vsExpTvsGain_1\WhitePaper_BackgroundWhitePaper_Intensity15_FR100Hz_Gain0dB_expT3ms_000' ];

rec = ReadRecord(recName);
for k = 1:size(rec)
    im = rec(:,:,k);
    figure; h = histcounts(im(:),0:500); plot(h)
    title(sprintf('frame %d',k)); 
end
