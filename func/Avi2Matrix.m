function  [Rec, v] = Avi2Matrix(filename, MaxNOfFrames)
    v = VideoReader(filename);
    nOfFramesInRecord = v.Duration*v.FrameRate -1; % Not using the field 'NumFrames' because it is not consistent with different matlab versions
%     if isprop(v,'NumberOfFrames')
%         nOfFramesInRecord = v.NumberOfFrames;
%     else
%         nOfFramesInRecord = v.NumFrames;
%     end
    
    if exist('MaxNOfFrames','var')        
        nOfFramesToRead = min(nOfFramesInRecord,MaxNOfFrames);
    else
        nOfFramesToRead = nOfFramesInRecord;    
    end
    
    % Read 
    Rec  = zeros(v.Height,v.Width,nOfFramesToRead);
    for i = 1:nOfFramesToRead
        im3color = read(v,i);
        Rec(:,:,i) = im3color(:,:,1); % all three colors are the same
    end
