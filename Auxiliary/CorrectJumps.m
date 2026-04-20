function ret = CorrectJumps(arr)
    ret = arr;
%     jump_th = 0.065; %std(rawSpeckleContrast{channel_i}(1:min(20,end))) * 5;
    diff_arr  =  diff(arr);

    jump_th = median(abs(diff_arr)) * 8;
    jump_up_idx = find(diff_arr > jump_th);
    jump_down_idx = find(diff_arr < -jump_th);
    if ~isempty(jump_up_idx) && ~isempty(jump_down_idx)
        jumps_up   = diff_arr(jump_up_idx);
        jumps_down = diff_arr(jump_down_idx);
       
        if jump_up_idx(1) > jump_down_idx(1)
            if numel(jump_up_idx) < numel(jump_down_idx)
                jump_up_idx( end + 1 ) = numel(arr); %#ok<AGROW>
                jump_down_idx(numel(jump_up_idx)+1:end) = [];
            end
            for n=1:numel(jump_down_idx)
                ret(jump_down_idx(n)+1:jump_up_idx(n))  =arr(jump_down_idx(n)+1:jump_up_idx(n)  ) - jumps_down(n);
            end
        else
            if numel(jump_down_idx) < numel(jump_up_idx)
                jump_down_idx( end + 1 ) = numel(arr); %#ok<AGROW>
                jump_up_idx(numel(jump_down_idx)+1:end) = [];
            end
            for n=1:numel(jump_up_idx)
                ret(jump_up_idx(n)+1:jump_down_idx(n))  = arr(jump_up_idx(n)+1:jump_down_idx(n) ) - jumps_up(n);
            end    
        end
    end
