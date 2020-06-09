function [LRRawEvent, LRRawNull, EventModel, NullModel,LPRawEvent,LPRawNull] = ...
    likRawNoHistModel(RawEvent, RawNull, Params, Filterflag)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2020 Shaoyu Qiao and Bijan Pesaran
% Pesaran Lab, New York University
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [LRRawEvent, LRRawNull, EventModel, NullModel] =
%                           likRawNoHistModel(RawEvent, RawNull, Params, Filterflag)
%
%   Inputs:
%       RawEvent = Raw events 
%       RawNull = Null events
%       Params = Parameters
%       Filterflag = flag for applying filters
%
%   Outputs:
%       LRRawEvent = Raw events filtered to Lfp events downsampled to 1 kHz
%       LRRawNull = = Raw null filtered to Lfp events downsampled to 1 kHz
%       EventModel
%       NullModel

if nargin < 4
    Filterflag = 1;
end

if(iscell(RawEvent))
    nSess = length(RawEvent);
else
    nSess = 1;
    myRawNull{1} = RawNull;
    myRawEvent{1} = RawEvent;
    RawNull = myRawNull; RawEvent = myRawEvent;
end

LPRawEvent = cell(1,nSess);
LPRawNull = cell(1,nSess);
sampling_rate = 1e3; % lfp

for iSess = 1:nSess
    if Filterflag
        tmpLPRawNull = mtfilter(RawNull{iSess},[0.03,100],3e4,0); % low-pass the raw signal at 100 Hz
        LPRawNull{iSess}  = tmpLPRawNull(:,1:3e4/sampling_rate:end); %downsample to 1 kHz
        
        tmpLPRawEvent = mtfilter(RawEvent{iSess},[0.03,100],3e4,0); % low-pass the raw signal at 100 Hz
        LPRawEvent{iSess} = tmpLPRawEvent(:,1:3e4/sampling_rate:end); %downsample to 1 kHz 
    else
        LPRawNull = RawNull;
        LPRawEvent = RawEvent;
    end
    
    LPRawMeanEvent{iSess} = mean(LPRawEvent{iSess},1);
    LPRawMeanNull{iSess} = mean(LPRawNull{iSess},1);
    
    nTrEvent = size(LPRawEvent{iSess},1);
    nTrNull = size(LPRawNull{iSess},1);
    nT = length(LPRawMeanEvent{iSess});
    
    ResNull = LPRawNull{iSess} - repmat(LPRawMeanNull{iSess}, nTrNull,1);
    sigmaNull(iSess) = std(ResNull(:));
    
    ResEvent = LPRawEvent{iSess} - repmat(LPRawMeanEvent{iSess},nTrEvent,1);
    sigmaEvent(iSess) = std(ResEvent(:));
    
    sigmaEventNull = (sigmaEvent(iSess)+sigmaNull(iSess))./2;
    LRRawEvent_tmp = zeros(nTrEvent,nT);
    ResidualRawEvent_tmp = zeros(nTrEvent,nT);
    for iTr = 1:nTrEvent
        LooTr = ~(ismember(1:nTrEvent,iTr));
        LPRawMeanEventTr = sum(LPRawEvent{iSess}(LooTr,:))./(nTrEvent-1);
        [LRRawEvent_tmp(iTr,:), ResidualRawEvent_tmp(iTr,:)] = ...
            likLR_Gaussian(LPRawEvent{iSess}(iTr,:), LPRawMeanEventTr, LPRawMeanNull{iSess}, sigmaEventNull, sigmaEventNull);
    end
    LRRawEvent{iSess} = LRRawEvent_tmp;
    ResidualRawEvent{iSess} = ResidualRawEvent_tmp;
    
    LRRawNull_tmp = zeros(nTrNull,nT);
    ResidualRawNull_tmp = zeros(nTrNull,nT);
    for iTr = 1:nTrNull
        LooTr = ~(ismember(1:nTrNull,iTr));
        LPRawMeanNullTr = sum(LPRawNull{iSess}(LooTr,:))./(nTrNull-1);
        [LRRawNull_tmp(iTr,:), ResidualRawNull_tmp(iTr,:)] = ...
            likLR_Gaussian(LPRawNull{iSess}(iTr,:), LPRawMeanEvent{iSess}, LPRawMeanNullTr, sigmaEventNull, sigmaEventNull);
    end
    LRRawNull{iSess} = LRRawNull_tmp;
    ResidualRawNull{iSess} = ResidualRawNull_tmp;
end

if nSess == 1
    EventModel.Residual = ResidualRawEvent{1};
    EventModel.Mean = LPRawMeanEvent{1};
    EventModel.sigma = sigmaEvent(1);
    
    NullModel.Mean = LPRawMeanNull{1};
    NullModel.sigma = sigmaNull(1);
    NullModel.Residual = ResidualRawNull{1};
    
    LRRawNull = LRRawNull{1}; LRRawEvent = LRRawEvent{1};
    LPRawNull = LPRawNull{1}; LPRawEvent = LPRawEvent{1};
else
    EventModel.Residual = ResidualRawEvent;
    EventModel.Mean = LPRawMeanEvent;
    EventModel.sigma = sigmaEvent;
    
    NullModel.Mean = LPRawMeanNull;
    NullModel.sigma = sigmaNull;
    NullModel.Residual = ResidualRawNull;
    
    LRRawEvent = averageLR(LRRawEvent);
    LRRawNull = averageLR(LRRawNull);
end
