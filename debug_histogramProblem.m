% analysisFolder = 'C:\Users\tarlevi\OneDrive - Bar-Ilan University - Students\PhD Research\SCOS\Records\NoiseAndBackground\Basler_1440GS_Vika01\Mono8\FixedPattern';
% recName = [ analysisFolder '\ShotNoise\vsExpTvsGain_1\WhitePaper_BackgroundWhitePaper_Intensity15_FR100Hz_Gain0dB_expT3ms_000' ];

% recName = 'F:\SCOS\Records\NoiseAndBackground\Basler_1440GS_Menachem_SN40335410\Mono12\Gain24dB_expT15ms_I15';
% frame1=57; frame2=58;
recName = 'F:\SCOS\Records\Jumps\vika_skull_2.5cmSDS_fingertapping_Gain36dB_expT15ms_FR10Hz';


rec = ReadRecord(recName);
Ivec = nan(1,size(rec,3));
for k=1:size(rec,3)
    Ivec(k) = mean2(rec(:,:,k));
end
figure; plot(Ivec); 

%%
frame1 = 48; frame2=49;
rec = ReadRecord(recName,2,frame1);
for k = 1:size(rec)
    im = rec(:,:,k);
    figure; [ counts,bins] = histcounts(im(:),0:500); plot(bins(1:end-1),counts);
    title(sprintf('frame %d',k+frame1-1)); 
end
%%
im1 = rec(:,:,1);
im2 = rec(:,:,2);
[ counts1,bins1] = histcounts(im1(:),0:500); 
[ counts2,bins2] = histcounts(im2(:),0:500); 

figure; 
plot(bins1(1:end-1),counts1); hold on;
plot(bins2(1:end-1),counts2);
legend(['im' num2str(frame1) ],['im' num2str(frame2)] );

my_imagesc(im2-im1)