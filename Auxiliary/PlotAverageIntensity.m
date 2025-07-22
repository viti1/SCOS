folder = 'E:\ShaareiZedek\12.01.2025\E2\Main_expT8ms_Gain18_FR10Hz_BL100';
tiffFiles = dir([ folder '\Basler_a2A1920-160umPRO__40513592__20250112_170130030_*']);
BL = ExtractParametersFromString( folder,'BL');

N = 100;
meanI = nan(1,N);
for i=1:N
    im = imread([folder , filesep , tiffFiles(i).name]); 
    meanI(i)= mean(double(im(:)))/64 - BL;
end

figure; plot(meanI); 
title('Average Intensity [DU]')