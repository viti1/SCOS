function [out1, out2, out3, out4, out5, out6, out7, out8, out9, out10] = InitNaN(arrSize)
    nanArr = nan(arrSize);
    for i=1:nargout
       eval(['out' num2str(i) '=nanArr;']) 
    end