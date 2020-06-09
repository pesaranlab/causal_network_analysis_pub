function [p, ST, Levels] = performance_levels(EventAccLLR, NullAccLLR, NumLevels,DetectType)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2020 Shaoyu Qiao and Bijan Pesaran
% Pesaran Lab, New York University
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [p, ST, Levels] = performance_levels(EventAccLLR, NullAccLLR, NumLevels, DetectType);
%
% This code is used to match the level for onset time detection between
% spikes and fields.
%   The level gives the same correct detect rates spikes and fields
% than rate of no detect in the context of baseline
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Input :       EventAccLLR: Accumulated log-likelihood ratios for the
%                           event distribution
%               NullAccLLR: Accumulated log-likelihood ratios for the null
%                       distribution
%               DetectType = String. 'Hit' or 'Correct'
%                      'Hit' only counts correctly classified Event trials
%                      'Correct' counts correctly classified Event and Null trials
%                           Defaults to 'Hit'
%   
%--------------------------------------------------------------------------
% Outputs:   Spike_Levels: The range of levels that were tested for spiking
%            Lfp_Levels: The range of levels that were tested for Lfps
%       p(:,1) = Prob of correct event detections for each level
%       p(:,2) = Prob of false miss detections for each level
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if nargin < 4 || isempty(DetectType)
    DetectType = 'Hit';
end

if nargin < 3 || isempty(DetectType)
    DetectType = 'Hit';
    NumLevels = 200;
end

[p, Levels] = characterize_level(EventAccLLR, NullAccLLR,NumLevels);

ST = zeros(1,length(Levels));

switch DetectType
    case 'Hit'
        for Levelind = 1:length(Levels)
            [Eventp, EventST] = DetectAccLLR(EventAccLLR, Levels(Levelind),-Levels(Levelind));
            if sum(~isnan(EventST))
                ST(Levelind) = mean(EventST(~isnan(EventST)));
            end
        end
    case {'Correct','Matched','CorrectMatched'}
        for Levelind = 1:length(Levels)
            [Eventp, EventST] = DetectAccLLR(EventAccLLR, Levels(Levelind),-Levels(Levelind));
            [Nullp, NullST] = DetectAccLLR(-NullAccLLR, Levels(Levelind),-Levels(Levelind));
            indEvent = find(~isnan(EventST));
            indNull = find(~isnan(NullST));
            if ~isempty(indEvent) & ~isempty(indNull)
                ST(Levelind) = (sum(EventST(indEvent))+sum(NullST(indNull)))./(length(indNull)+length(indEvent));
            elseif ~isempty(indEvent) & isempty(indNull)
                ST(Levelind) = sum(EventST(indEvent))./length(indEvent);
            elseif isempty(indEvent) & ~isempty(indNull)
                ST(Levelind) = sum(NullST(indNull))./length(indNull);
            end
        end
end
