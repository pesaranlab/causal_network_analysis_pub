function plotStimAccLLR(Results,StimTrials,bn,cmap,saveFlag,plotFlag)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2020 Shaoyu Qiao and Bijan Pesaran
% Pesaran Lab, New York University
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This routine plots AccLLR results for electrical stimulation
% plotStimSyllableAccLLR(Results,StimTrials,bn,cmap,saveFlag,plotFlag,varargin)
%
%     Results: AccLLR Results structure
%  StimTrials: StimTrials structure
%          bn: time epoch for stim events in msec
%ERPplotRange: ERP plotting range
%       cmap : color map
%    saveFlag: 1 to save figures, 0 not
%    plotFlag: 1 to plot figures, 0 not

% Written by Shaoyu Qiao, July 20, 2016

if nargin < 6
    plotFlag = 0;
end

if nargin < 5
    saveFlag = 1;
    plotFlag = 0;
end

if nargin < 4
    saveFlag = 1;
    plotFlag = 0;
    cmap = hot;
end


StimPair = [StimTrials(1).Syllable.Pulse(1).Cathode StimTrials(1).Syllable.Pulse(1).Anode];
StimAmp = abs(StimTrials(1).Syllable.Pulse(1).AmplitudeCathode);
StimName = StimTrials(1).Name;
rec = StimTrials(1).Rec;
fs = Results(1).Fs;
day = StimTrials(1).Day;

bnFs = bn*fs/1e3;
t = linspace(bnFs(1)/(fs/1e3),bnFs(2)/(fs/1e3),sum(abs(bnFs)));

pCorrectDetect = Results.pCorrectDetect;
pIncorrectDetect = Results.pIncorrectDetect;
pDiffCorrectVsIncorrectDetect = Results.pDiffCorrectVsIncorrectDetect;
ST = Results.ST;
EventST = Results.EventST;
TargCh = Results.Ch;
LPRawEvent = Results.LPRawEvent;
LPRawNull = Results.LPRawNull;


AccLLRRawNull = Results.NoHist.Null.AccLLR;
AccLLRRawEvent = Results.NoHist.Event.AccLLR;
StimArtifactBlankWin = Results.StimArtifactBlankWin; % ms
StimArtifactStartInd = Results.StimArtifactStartInd;
AccLLRwin = Results.AccLLRwin;
nTr = length(Results.EventST);

[dum,ind] = sort(EventST);
nTr_CorrectDetect = sum(~isnan(dum));

mLPRawEvent = mean(LPRawEvent,1);
sdLPRawEvent = std(LPRawEvent,[],1);
semLPRawEvent = sdLPRawEvent/sqrt(nTr);

mLPRawNull = mean(LPRawNull,1);
sdLPRawNull = std(LPRawNull,[],1);
semLPRawNull = sdLPRawNull/sqrt(nTr);

mAccLLRRawEvent = mean(AccLLRRawEvent);
sdAccLLRRawEvent = std(AccLLRRawEvent);
semAccLLRRawEvent = sdAccLLRRawEvent/sqrt(nTr);
mAccLLRRawNull = mean(AccLLRRawNull);
sdAccLLRRawNull = std(AccLLRRawNull);
semAccLLRRawNull = sdAccLLRRawNull/sqrt(nTr);


%% calculate ROC choice probability
if isfield(Results,'ROC')
    AUC = Results.ROC.AUC;
    se = Results.ROC.se;
    t_ROC = Results.ROC.t;
    Dt = Results.ROC.Dt;
else
    MaxTime = size(AccLLRRawNull,2);
    AUC = zeros(1,MaxTime);
    se = zeros(1,MaxTime);
    Dt = 1; % msec
    for iT = 1:MaxTime
        [AUC(iT),se(iT)] = myroc(AccLLRRawNull(:,iT)'+rand(1,size(AccLLRRawNull,1))*0.04, ...
            AccLLRRawEvent(:,iT)' + rand(1,size(AccLLRRawEvent,1))*0.04);
    end
    t_ROC = 1:Dt:MaxTime;
    t_ROC = t_ROC + StimArtifactBlankWin;
end

%% %%%% plotting %%%%%%%%%%%%%%%%%
if plotFlag
    figure('visible','on','Position',[100 100 800 600], 'NumberTitle', 'off','Name','StimAccLLR')
else
    figure('visible','off','Position',[100 100 800 600], 'NumberTitle', 'off','Name','StimAccLLR')
end


if isequal(length(TargCh),2)
    suptitle([StimName, 'Stim-response @ e' num2str(TargCh(1)),'-e' num2str(TargCh(2)) ' [Stim: e' num2str(StimPair(1)) '(c)-e' num2str(StimPair(2)) '(a), Amp: ' num2str(StimAmp), '\muA]'])
else
    suptitle([StimName, 'Stim-response @ e' num2str(TargCh),' [Stim: e' num2str(StimPair(1)) '(c)-e' num2str(StimPair(2)) '(a), Amp: ' num2str(StimAmp), '\muA]'])
end

subplot(3,3,1)
plotSTA((1:length(mAccLLRRawEvent))/fs*1e3+StimArtifactBlankWin,mAccLLRRawEvent,semAccLLRRawEvent,[0 0 1]);
hold on
plotSTA((1:length(mAccLLRRawNull))/fs*1e3+StimArtifactBlankWin,mAccLLRRawNull,semAccLLRRawNull,[1 0 0]);
ylabel('Mean AccLLR')
xlim([0 StimArtifactBlankWin+AccLLRwin])
text(10,round(max(max([abs(mAccLLRRawEvent),abs(mAccLLRRawNull)])))/4*3,'LFP Event','color','b')
text(10,-round(max(max([abs(mAccLLRRawEvent),abs(mAccLLRRawNull)])))/4*3,'LFP Null','color','r')
box off;


h1=subplot(3,3,[2 5]);
plot(ST+StimArtifactBlankWin,pCorrectDetect,'x')
hold on;
plot(ST+StimArtifactBlankWin,pIncorrectDetect,'o')
xlabel('Selection time (ms)')
ylabel('Probability');
xlim([0 StimArtifactBlankWin+AccLLRwin])
legend('Hit','FA')
box off;
pos1 = get(h1,'Position');
hold off


h2 = subplot(3,3,8);
ERPplotRange = [-round(max(abs(LPRawEvent(:)))) round(max(abs(LPRawEvent(:))))];
imagesc((1:length(mAccLLRRawNull))/fs*1e3+StimArtifactBlankWin,1:nTr,LPRawEvent(ind,:),ERPplotRange);
box off;
xlabel('Time after stim onset (ms)');
ylabel('Sortetd trial')
xlim([0 StimArtifactBlankWin+AccLLRwin])
c = colorbar;
c.Label.String = 'Voltage (\muV)';
colormap(cmap)
for i = 1 : nTr_CorrectDetect
    text(dum(i)/fs*1e3+StimArtifactBlankWin,i,'x','color','k');
end
pos2 = get(h2,'Position');
set(h2,'Position',[pos2(1) pos2(2) pos1(3) pos2(4)]);


subplot(3,3,3)
plot([pCorrectDetect' pIncorrectDetect' pDiffCorrectVsIncorrectDetect'],'linewidth',2)
xlabel('Level of threshold for detection');ylabel('Probability')
ylim([0 1])
legend('Correct','InCorrect','Diff','Location','NorthEast')
box off;
legend boxoff


subplot(3,3,4)
tInd = StimArtifactStartInd+StimArtifactBlankWin*fs/1e3+1:StimArtifactStartInd+StimArtifactBlankWin*fs/1e3+size(LPRawEvent,2);
plotSTA(t(tInd),mLPRawEvent,semLPRawEvent,[0 0 1]); hold on
plotSTA(t(tInd),mLPRawNull,semLPRawNull,[1 0 0]);
xlim(round([0 t(tInd(end))]))
%xlabel('Time after stim onset (ms)');
ylabel('Voltage (\muV)');
title({'Evoked potential (\mu \pm sem)'})
box off;
hold off


subplot(3,3,7)
hold on;
h = patch([t_ROC,t_ROC(end:-Dt:1)],[AUC+2*se,AUC(end:-1:1)-2*se(end:-1:1)],0.7*[1,1,1]); % 95% CI
set(h,'FaceAlpha',.5,'EdgeAlpha',0,'Linestyle','none');
hold on;
plot(t_ROC,AUC,'k','Linewidth',2);
set(h,'Linewidth',2);
axis([0 StimArtifactBlankWin+AccLLRwin 0.4 1]);
h = line([0 StimArtifactBlankWin+AccLLRwin],[0.5,0.5]); set(h,'Linewidth',1,'Linestyle','--','Color','k');
xlabel('Time after stim onset (ms)');
ylabel('Choice probability')
set(gca,'ytick',[0.5 1])
hold off
