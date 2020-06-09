function [LR, Res1, Res2] = likLR_Gaussian(Lfp, Evoked1, Evoked2, sigmaLfp1, sigmaLfp2)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2020 Shaoyu Qiao and Bijan Pesaran
% Pesaran Lab, New York University
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [LR, Res1, Res2] = likLR_Gaussian(Lfp, Evoked1, Evoked2, sigmaLfp1, sigmaLfp2) 

Evoked1Mat = repmat(Evoked1, size(Lfp,1),1);
Evoked2Mat = repmat(Evoked2, size(Lfp,1),1);

Res1 = (Lfp-Evoked1Mat);
Res2 = (Lfp-Evoked2Mat);

LR = log(sigmaLfp2)-log(sigmaLfp1)+Res2.^2/(2*sigmaLfp2.^2)-Res1.^2./(2*sigmaLfp1.^2);
