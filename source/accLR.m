function [AccLLR, LLR] = accLR(LLR, Go, RT)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2020 Shaoyu Qiao and Bijan Pesaran
% Pesaran Lab, New York University
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  [AccLLR, LLR] = accLR(LLR, Go, RT)


if nargin < 2; Go = 1; end
dt = 1;
Go;
if nargin<3;
    AccEnd = size(LLR,2);
    AccLLR = nan(size(LLR));
    for tr = 1:size(LLR,1)
        i = 0;
        for t = Go+1:dt:Go+AccEnd;
            i = i+1;
            AccLLR(tr,i) = sum(LLR(tr,Go:dt:t));
        end
    end
else
    AccLLR = nan(size(LLR,1),max(RT));
    for tr = 1:size(LLR,1)
        i = 0;
        for t = Go+1:dt:Go+RT(tr)
            i = i+1;
            AccLLR(tr,i) = sum(LLR(tr,Go:dt:t));
        end
    end
end
%sum(isnan(AccLLR))