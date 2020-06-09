function [modDecodeHitRate,modDecodeMissRate,rocThresh,DecoderTrials] = runModulatorDecoder(S1,S2,S_all,EventST,roc_Thresh,iThresh)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2020 Shaoyu Qiao and Bijan Pesaran
% Pesaran Lab, New York University
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 6
    iThresh = [];
end

if isempty(iThresh)
    range = 1 : numel(roc_Thresh);
else
    range = iThresh;
end

[dum,ind] = sort(EventST);
nTr_CorrectDetect = sum(~isnan(dum));
receiverAccLLR_Hit_tr = sort(ind(1:nTr_CorrectDetect));

for i = range
    if isempty(iThresh)
        cutoff = roc_Thresh(i);
    else
        cutoff = roc_Thresh(iThresh);
    end
    if median(S1) > median(S2)
        modulatorDecoder_Hit_tr = [];
        modulatorDecoder_Miss_tr = [];
        modulatorDecoder_Miss_tr_bad = [];
        
        modulatorDecoder_Hit_tr = find(S_all >= cutoff);
        
        modulatorDecoder_Miss_tr = find(S_all < cutoff);
        modulatorDecoder_Miss_tr_bad = modulatorDecoder_Miss_tr(S_all(modulatorDecoder_Miss_tr)==0);
        modulatorDecoder_Miss_tr = modulatorDecoder_Miss_tr(S_all(modulatorDecoder_Miss_tr)~=0);
    else
        modulatorDecoder_Hit_tr = [];
        modulatorDecoder_Miss_tr = [];
        modulatorDecoder_Hit_tr_bad = [];
        
        modulatorDecoder_Hit_tr = find(S_all <= cutoff);
        modulatorDecoder_Hit_tr_bad = modulatorDecoder_Hit_tr(S_all(modulatorDecoder_Hit_tr)==0);
        modulatorDecoder_Hit_tr = modulatorDecoder_Hit_tr(S_all(modulatorDecoder_Hit_tr)~=0);
        
        modulatorDecoder_Miss_tr = find(S_all > cutoff);
    end
    % trials in modulatorDecoder_Hit_tr are in Receiver_AccLLR_Hit_tr
    ModHitInRecHit_tr = modulatorDecoder_Hit_tr(ismember(modulatorDecoder_Hit_tr,receiverAccLLR_Hit_tr));
    [dum1,temp1] = sort(EventST(ModHitInRecHit_tr));
    sorted_ModHitInRecHit_tr = ModHitInRecHit_tr(temp1);
    
    ModHitNotInRecHit_tr = setdiff(modulatorDecoder_Hit_tr,ModHitInRecHit_tr);
    
    % trials in modulatorDecoder_Miss_tr are in Receiver_AccLLR_Hit_tr
    ModMissInRecHit_tr = modulatorDecoder_Miss_tr(ismember(modulatorDecoder_Miss_tr,receiverAccLLR_Hit_tr));
    [dum2,temp2] = sort(EventST(ModMissInRecHit_tr));
    sorted_ModMissInRecHit_tr = ModMissInRecHit_tr(temp2);
    
    ModMissNotInRecHit_tr = setdiff(modulatorDecoder_Miss_tr,ModMissInRecHit_tr);
    
    modulatorDecodingHitRate(i) = numel(sorted_ModHitInRecHit_tr)/(numel(sorted_ModHitInRecHit_tr)+numel(ModHitNotInRecHit_tr));
    modulatorDecodingMissRate(i) = numel(ModMissNotInRecHit_tr)/(numel(sorted_ModMissInRecHit_tr)+numel(ModMissNotInRecHit_tr));
end

if isequal(numel(range),1)
    modDecodeHitRate = nonzeros(modulatorDecodingHitRate);
    modDecodeMissRate = nonzeros(modulatorDecodingMissRate);
    rocThresh = cutoff;
    
    DecoderTrials.sorted_ModHitInRecHit_tr = sorted_ModHitInRecHit_tr;
    DecoderTrials.ModHitNotInRecHit_tr = ModHitNotInRecHit_tr;
    DecoderTrials.sorted_ModMissInRecHit_tr = sorted_ModMissInRecHit_tr;
    DecoderTrials.ModMissNotInRecHit_tr = ModMissNotInRecHit_tr;
    DecoderTrials.modulatorDecoder_Hit_tr = modulatorDecoder_Hit_tr;
    DecoderTrials.modulatorDecoder_Miss_tr = modulatorDecoder_Miss_tr;
    DecoderTrials.dum1 = dum1;
    DecoderTrials.dum2 = dum2;
else
    CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
    CI_HitRate = CIFcn(modulatorDecodingHitRate,95);
    CI_MissRate = CIFcn(modulatorDecodingMissRate,95);
    
    indx = modulatorDecodingHitRate < CI_HitRate(2) & modulatorDecodingHitRate > CI_HitRate(1) & modulatorDecodingMissRate > CI_MissRate(1) & modulatorDecodingMissRate < CI_MissRate(2);
    
    modDecodeHitRate = modulatorDecodingHitRate(indx);
    modDecodeMissRate = modulatorDecodingMissRate(indx);
    rocThresh = roc_Thresh(indx);
    
    DecoderTrials = [];
end
