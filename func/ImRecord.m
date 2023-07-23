function [im, f,info] = ImRecord(recordName)

[ im , ~ , ~ , ~ , info] = ReadRecord(recordName,1); %, {'Tint','FR','Gain','Fiber'});

[~, rawName ]  = fileparts(recordName);
rawName = strrep(rawName,'_',' ');
f = figure('position',[50,50,1200,800]); imagesc(im); colormap gray; colorbar
title(rawName,'interpreter','none')
