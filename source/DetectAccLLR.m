function [p, ST] = DetectAccLLR(AccLLR, UpperLevel, LowerLevel)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2020 Shaoyu Qiao and Bijan Pesaran
% Pesaran Lab, New York University
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  [p, ST] = DetectAccLLR(AccLLR, UpperLevel, LowerLevel)
%
%   p(1) = Prob of upper level detections
%   p(2) = Prob of lower level detections
%   p(3) = Prob of don't knows
%
%   ST = Selection times for upper level detections.  nan indicates miss.
%
nTr = size(AccLLR,1);
Num1 = 0; Num2 = 0; Num3 = 0;
ST = zeros(1,nTr);

for iters = 1:nTr
    TSTemp1 = find(AccLLR(iters,:)>UpperLevel);
    TSTemp2 = find(AccLLR(iters,:)<LowerLevel);
    Flag1 = ~isempty(TSTemp1);
    Flag2 = ~isempty(TSTemp2);
    if Flag1 && ~Flag2    %  Hit upper Level1 only
        Num1 = Num1 + 1;
        ST(1,iters) = TSTemp1(1);
    elseif Flag2 && ~Flag1  %  Hit lower Level2 only
        Num2 = Num2 + 1;
        ST(1,iters) = nan;
    elseif Flag1 && Flag2   % Hit both levels
        if TSTemp2(1) > TSTemp1(1)
            Num1 = Num1 + 1;
            ST(1,iters) = TSTemp1(1);
        else
            Num2 = Num2 + 1;
            ST(1,iters) = nan;
        end
    else %  Hit no levels
        Num3 = Num3 + 1;
        ST(1,iters) = nan;
    end
end

p(1) = Num1./nTr;
p(2) = Num2./nTr;
p(3) = Num3./nTr;