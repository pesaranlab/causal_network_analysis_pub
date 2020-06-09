function [Lfp, inds] = removeNoisyLfpTrials( Lfp, stdThresh )
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2020 Shaoyu Qiao and Bijan Pesaran
% Pesaran Lab, New York University
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [lfp, inds] = removeNoisyLfpTrials( lfp, stdThresh )
% removes noisy trials above a given standard deviation
% 
% in: Lfp trials x time
%     stdThresh number of standard deviations a signal must be within
% out: Lfp - Lfp with trials whose Lfp exceeds a stdThresh * std(Lfp(:)) 
%            removed trials x time
%      inds - indices of trials whose Lfp does not exceed the threshold

thresh = stdThresh*std(Lfp(:));
e = max(abs(Lfp'));
inds = find(e < thresh);
Lfp = Lfp(inds,:);

