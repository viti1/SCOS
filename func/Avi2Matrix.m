function  [Rec, v] = Avi2Matrix(filename, MaxNOfFrames)
    v = VideoReader(filename);  
    if exist('MaxNOfFrames','var')
        nOfFramesToRead = min(v.NumFrames,MaxNOfFrames);
    else
        nOfFramesToRead = v.NumberOfFrames;    
    end
    
    % Read 
    Rec  = zeros(v.Height,v.Width,nOfFramesToRead);
    for i = 1:nOfFramesToRead
        im3color = read(v,i);
        Rec(:,:,i) = im3color(:,:,1); % all three colors are the same
    end
