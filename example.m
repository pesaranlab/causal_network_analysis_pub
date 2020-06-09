clear variables
close all
clc

%% Add StimAccLLR to the path
addpath(genpath('./source'));


%% Load data
load('./example_data.mat');


%% Perform StimAccLLR for single-pulse stimulation response
% define input variables
AnalParams = data.Params.Anal;
ReceiverSig = data.Data.raw.receiver; % receiver recording
AccLLRwin = AnalParams.AccLLR.AccLLRwin; % AccLLR accumulation time
fs = data.Fs.raw; % sampling rate
bn = [-540 140]; % bout/trial epoch in ms
bnFs = bn*fs/1e3; % bout/trial epoch in samples
TargCh = data.receiverCh;
StimTimeInd = abs(bnFs(1));
StimArtifactBlankWin = data.Params.Stim.StimArtifactBlankWin;

% run StimAccLLR
[Results,Model] = StimAccLLR(ReceiverSig,AccLLRwin,bn,fs,TargCh,StimTimeInd,StimArtifactBlankWin);

% plot StimAccLLR results
StimTrials = data.StimTrials;
cmap = polarmap;
saveFigFlag = 0;
plotFigFlag = 1;
plotStimAccLLR(Results(1),StimTrials,bn,cmap,saveFigFlag,plotFigFlag)


%% Perform ROC analysis for modualtor decoding
% extract Hit and Miss events based on the receiver responses
EventST = Results.EventST;
HitIndx = find(~isnan(EventST));
MissIndx = find(isnan(EventST));

bn_baseline = [-510 -10];% baseline epoch in ms 
bnFs_spec = bn_baseline -  bn(1);
tindx_lfp_spec = linspace(bnFs_spec(1)+1,bnFs_spec(2),diff(abs(bnFs_spec)));

ModulatorSig = data.Data.lfp.modulator; % modulator baseline activity

X1 = ModulatorSig(HitIndx,tindx_lfp_spec); % Modulator baseline activity in Hit events
[X1,goodInd1] = removeNoisyLfpTrials(X1,5);
hitIndx = HitIndx(goodInd1);

X2 = ModulatorSig(MissIndx,tindx_lfp_spec); % Modulator baseline activity in Hit events
[X2,goodInd2] = removeNoisyLfpTrials(X2,5);
missIndx = MissIndx(goodInd2);


figure('NumberTitle', 'off','Name','Modulator Decoding')
% plot ROC curve of modulator activity (hit vs miss)
subplot(2,2,[1 3])
[auc,se,S1,S2,roc_Thresh] = calcRocSpecDiff(X1,X2,AnalParams);

S_all = [];
S_all(hitIndx) = S1;
S_all(missIndx) = S2;

% get different decoding rates across different ROC thresholds
[modDecodeHitRate,modDecodeMissRate,rocThresh,~] = runModulatorDecoder(S1,S2,S_all,EventST,roc_Thresh);

% find optimal decoding rate
ind = modDecodeHitRate > modDecodeMissRate;
indd = find(ind);
modDecodeHitPlusMissRate = (modDecodeHitRate(ind) + modDecodeMissRate(ind))/2;
[~,modDecodeHitPlusMissRateInd] = max(modDecodeHitPlusMissRate);
[~,iThresh] = min(abs(roc_Thresh-rocThresh(indd(modDecodeHitPlusMissRateInd))));
[optModDecodeHitRate,optModDecodeMissRate,optRocThresh,DecoderTrials] = runModulatorDecoder(S1,S2,S_all,EventST,roc_Thresh,iThresh);

% plot histogram 
subplot(2,2,2)
hg1=histogram(S1,'BinMethod','fd','Normalization','count','FaceColor','r');hold on
hg2=histogram(S2,'BinMethod','fd','Normalization','count','FaceColor','k');
plot([optRocThresh,optRocThresh],[0,max([hg1.Values hg2.Values])+1],'b--')
legend('Hit','Miss')
xlabel('Mean log \beta power');
ylabel('Count')
title('Modulator ROC decoding')
xlim([floor(rocThresh(1)) ceil(rocThresh(end))])

% plot decoding performance
subplot(2,2,4)
plot(rocThresh,(modDecodeHitRate+modDecodeMissRate)/2,'linewidth',2,'color','k');hold on
plot([optRocThresh,optRocThresh],[floor(min([modDecodeHitRate modDecodeMissRate])*10)/10,ceil(max([modDecodeHitRate modDecodeMissRate])*10)/10],'b--')
xlabel('Mean log \beta power');
ylabel('Decoding performance')
ylim([floor(min([modDecodeHitRate modDecodeMissRate])*10)/10,ceil(max([modDecodeHitRate modDecodeMissRate])*10)/10])
ylim([0.5 0.7])
