function [pval,se,aS_X1,aS_X2,roc_Thresh] = calcRocSpecDiff(X1,X2,AnalParams)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2020 Shaoyu Qiao and Bijan Pesaran
% Pesaran Lab, New York University
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [pval,se,aS_X1,aS_X2,roc_Thresh] = calcRocSpecDiff(X1,X2,AnalParams)
% This routine performs ROC analysis of LFP power in two conditions

% Written by Shaoyu Qiao, Aug 15, 2018


if(isfield(AnalParams,'Tapers'))
    N = AnalParams.Tapers(1);
    if length(AnalParams.Tapers)==3
        W = AnalParams.Tapers(2)./AnalParams.Tapers(1);
    else
        W = AnalParams.Tapers(2);
    end
end

if(isfield(AnalParams,'pad'))
    pad = AnalParams.pad; %
else
    pad = 2; % Hz
end

if(isfield(AnalParams.Spec.Test,'fk'))
    fk = AnalParams.Spec.Test.fk; %
else
    fk = [13 30]; %
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sampling = 1e3; % lfp sampling rate, Hz
p = N*W;
k = floor(2*p-1);
tapers = [N,p,k];
tapers(1) = floor(tapers(1).*sampling);
tapers = dpsschk(tapers);
K = length(tapers(1,:));
n = length(tapers(:,1));
nf = max(256,pad*2.^(nextpow2(n+1)));

if length(fk)==1
    fk = [0,fk]; 
end

nfk = floor(fk./sampling.*nf);
ntr1 = size(X1,1);
ntr2 = size(X2,1);

if ntr1 < 2 || ntr2 < 2
    pval = nan;
    S_X1 = nan; 
    S_X2 = nan;
    error(['Not enough trials ' num2str(ntr1) ':' num2str(ntr2)])
end


Xk1 = zeros(ntr1*K, diff(nfk));
Xk2 = zeros(ntr2*K, diff(nfk));
mX1 = sum(X1)./ntr1; mX2 = sum(X2)./ntr2;
X1 = X1 - mX1(ones(1,ntr1),:);
X2 = X2 - mX2(ones(1,ntr2),:);

for tr = 1:ntr1
    tmp1 = X1(tr,:)';
    xk = fft(tapers(:,1:K).*tmp1(:,ones(1,K)),nf)';
    Xk1((tr-1)*K+1:tr*K,:) = xk(:,nfk(1)+1:nfk(2));
end

for tr = 1:ntr2
    tmp2 = X2(tr,:)';
    xk = fft(tapers(:,1:K).*tmp2(:,ones(1,K)),nf)';
    Xk2((tr-1)*K+1:tr*K,:) = xk(:,nfk(1)+1:nfk(2));
end

S_X1 = Xk1.*conj(Xk1);
S_X2 = Xk2.*conj(Xk2);

% average across frequency bins
aS_X1 = sum(log(S_X1),2)./diff(nfk);
aS_X2 = sum(log(S_X2),2)./diff(nfk);

%run ROC analysis (pval = AUC; se = sd of residuals)
[pval,se,roc_Thresh] = myroc(aS_X1',aS_X2',0,1);
