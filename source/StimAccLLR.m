function [Results,Model] = StimAccLLR(data,AccLLRwin,bn,fs,TargCh,StimArtifactStartInd,StimArtifactBlankWin,Ind,LPRawNull,AccLLRwin_end_Null,NumLevels,Shuffle)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2020 Shaoyu Qiao and Bijan Pesaran
% Pesaran Lab, New York University
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This routine runs AccLLR for stim response detection
% [Results,Model] = StimAccLLR(data,AccLLRwin,bn,fs,TargCh,StimArtifactStartInd,StimArtifactBlankWin,Ind,RawNull,AccLLRwin_end_Null,Shuffle)
%
% Input:
%                     data: [trial, time] time series for specific channel
%                       bn: [t1,t2] time epoch for the trial/event  (ms)
%                AccLLRwin: time window for calculating AccLLR (ms)
%                       fs: sampling rate (Hz)
%                   TargCh: target channel
%     StimArtifactStartInd: index of the stimulus artifact start of specific pulse per trial
%     StimArtifactBlankWin: time window for blocking due to stimulation artifact (ms)
%                      Ind: target channel index
%                  RawNull: time series for the null event
%       AccLLRwin_end_Null: End of window of null event (ms)
%                NumLevels: Detection levels
%                  Shuffle: 1 = shuffle null and event events; 0 = no shuffle

% Output:
%                  Results: struct for results of AccLLR
%                    Model: struct for AccLLR models

% Written by Shaoyu Qiao, last updated Feb 27, 2020

if nargin < 12
    Shuffle = 0;
end

if nargin < 11
    NumLevels = AccLLRwin; % every 1 ms    
    Shuffle = 0;
end

if nargin < 10
    AccLLRwin_end_Null = -100;% in msec
    NumLevels = AccLLRwin; % every 1 ms    
    Shuffle = 0;
end

if nargin < 9
    LPRawNull = [];
    AccLLRwin_end_Null = -100;% in msec
    NumLevels = AccLLRwin; % every 1 ms
    Shuffle = 0;
end

if nargin < 8
    Ind = 1;
    LPRawNull = [];
    AccLLRwin_end_Null = -100;% in msec
    NumLevels = AccLLRwin; % every 1 ms
    Shuffle = 0;
end


%% AccLLR %%
AccLLRwin_Event = [StimArtifactBlankWin StimArtifactBlankWin+AccLLRwin]; % in msec
AccLLRwin_Null = [AccLLRwin_end_Null-AccLLRwin AccLLRwin_end_Null]; % in msec

RawEvent = data(:,StimArtifactStartInd+StimArtifactBlankWin*fs/1e3+1:end);
RawEvent = RawEvent(:,1:AccLLRwin*fs/1e3);
nTr_Event = size(RawEvent,1);

if isempty(LPRawNull)
    LPRawNull = data(:,StimArtifactStartInd-(AccLLRwin-AccLLRwin_end_Null)*fs/1e3+1:StimArtifactStartInd+AccLLRwin_end_Null*fs/1e3);
end
nTr_Null = size(LPRawNull,1);
min_num_trials = nTr_Null;

StartofAccumulationTime = 1;
StartTime = 0;

fs_lfp = 1e3; % sampling rate lfp, Hz
MaximumTimetoOnsetDetection = floor((AccLLRwin-1)*fs_lfp/1e3);
EndAcc = MaximumTimetoOnsetDetection*ones(1,min_num_trials);
StartofAccumulation = StartofAccumulationTime-StartTime;


[LRRawEvent,LRRawNull,EventModel,NullModel,LPRawEvent,LPRawNull] = likRawNoHistModel(RawEvent,LPRawNull,[]);
AccLLRRawNull = accLR(LRRawNull,StartofAccumulation, EndAcc);
AccLLRRawEvent = accLR(LRRawEvent, StartofAccumulation, EndAcc);


[p, ST, Levels] = performance_levels(AccLLRRawEvent, AccLLRRawNull,NumLevels);
ST = ST*1./(fs_lfp/1e3);

pCorrectDetect = (p(:,1,1)+p(:,2,2))./2;
pIncorrectDetect = (p(:,2,1)+p(:,1,2))./2;
pDiffCorrectVsIncorrectDetect = pCorrectDetect-pIncorrectDetect;
[PdiffmaxVal,PdiffmaxInd] = max(pDiffCorrectVsIncorrectDetect);

Level = Levels(PdiffmaxInd);
[Eventp, EventST] = DetectAccLLR(AccLLRRawEvent, Level,-Level);
[Nullp, NullST] = DetectAccLLR(-AccLLRRawNull, Level,-Level);

%% calculate ROC over AccLLR accmulation time
MaxTime = size(AccLLRRawNull,2);
ROC = zeros(1,MaxTime);
se = zeros(1,MaxTime);
Dt = (1/fs_lfp)*1e3; % msec
for iT = 1:MaxTime
    [ROC(iT),se(iT)] = myroc(AccLLRRawNull(:,iT)'+rand(1,size(AccLLRRawNull,1))*0.04, ...
        AccLLRRawEvent(:,iT)' + rand(1,size(AccLLRRawEvent,1))*0.04,0,0);
end
MaxTime_ms = MaxTime*(1/fs_lfp)*1e3;
t_ROC = Dt:Dt:MaxTime_ms;
t_ROC = t_ROC + StimArtifactBlankWin;

%% compile data
Results.pCorrectDetect = pCorrectDetect';
Results.pIncorrectDetect = pIncorrectDetect';
Results.pDiffCorrectVsIncorrectDetect = pDiffCorrectVsIncorrectDetect';
Results.ST = ST;
Results.Levels = Levels;
Results.Eventp = Eventp;
Results.EventST = EventST;
Results.Nullp = Nullp;
Results.NullST = NullST;
Results.Ch = TargCh;
Results.Ind = Ind;
Results.Fs = fs_lfp;
Results.bn.Trial = bn;
Results.bn.Null = AccLLRwin_Null;
Results.bn.Event = AccLLRwin_Event;
Results.AccLLRwin = AccLLRwin;
Results.NoHist.Null.LR = LRRawNull;
Results.NoHist.Event.LR = LRRawEvent;
Results.NoHist.Null.AccLLR = AccLLRRawNull;
Results.NoHist.Event.AccLLR = AccLLRRawEvent;
Results.StimArtifactBlankWin = StimArtifactBlankWin;
Results.StimArtifactStartInd = StimArtifactStartInd;
Results.StimArtifactStartInd = abs(bn(1)*fs_lfp/1e3);
Results.RawEvent = RawEvent;
Results.RawNull = LPRawNull;
Results.LPRawEvent = LPRawEvent;
Results.LPRawNull = LPRawNull;
Results.AccLLRwin_end_Null = AccLLRwin_end_Null;
Results.ROC.AUC = ROC;
Results.ROC.se = se;
Results.ROC.t = t_ROC;
Results.ROC.Dt = Dt;
Results.Perm = [];

Model.NoHist.Event = EventModel;
Model.NoHist.Null = NullModel;
Model.Type = 'lfp';
Model.Ind = Ind;
Model.shuffle = Shuffle;
