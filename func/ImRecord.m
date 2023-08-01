function [im, f,info] = ImRecord(recordName)

[ im , ~ , ~ , ~ , info] = ReadRecord(recordName,1); %, {'Tint','FR','Gain','Fiber'});

[~, rawName ]  = fileparts(recordName);
rawName = strrep(rawName,'_',' ');
f = figure(); 
imagesc(im); colormap gray; colorbar
axis equal % axis equal should be before setting the limits
ax = gca();
ax.XLim = [ 0 size(im,2) ];
ax.YLim = [ 0 size(im,1) ];
title(rawName,'interpreter','none')
