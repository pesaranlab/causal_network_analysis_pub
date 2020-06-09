function [Y,ph] = mtfilter(X, tapers, sampling, f0, flag, complexflag)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2020 Shaoyu Qiao and Bijan Pesaran
% Pesaran Lab, New York University
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  MTFILTER Bandpass filter a time series using the multitaper method
%
%  [Y,PH] = MTFILTER(X, TAPERS, SAMPLING, F0, FLAG, COMPLEXFLAG)
%
%  Inputs:  X 		=  Time series array in [Space/Trials,Time] form
%	    TAPERS	=  Data tapers in [K,TIME], [N,P,K] or [N,W] form.
%				Defaults to [N,3,5] where N is duration of X
%	    SAMPLING 	=  Sampling rate of time series X in Hz.
%				Defaults to 1 Hz.
%	    F0		=  Center frequency of filter.
%				Defaults to 0 Hz.
%
%	    FLAG = 0:	Output data should be centered (Default)
%	    FLAG = 1:	Output data should not be centered.
%
%       COMPLEXFLAG = 0:  Output data should be real (Default)
%       COMPLEXFLAG = 1:  Output data should be complex.
%
%  Outputs: Y		=  Filtered version of X in [Space/Trials,Time] form.
%				If F0 is nonzero Y is complex.
%

%  Author:  Bijan Pesaran 04/15/98

if nargin < 3; sampling=1; end
if length(tapers) == 2
    n = tapers(1);
    w = tapers(2);
    p = n*w;
    k = floor(2*p-1);
    tapers = [n,p,k];
end
if length(tapers) == 3
    tapers(1) = tapers(1).*sampling;
    tapers = dpsschk(tapers);
end
if nargin < 4 f0 = 0; end
if nargin < 5 flag = 0; end
if nargin < 6 complexflag = 0; end

filt = mtfilt(tapers, sampling, f0);
N = length(filt);
%disp(['Filter is ' num2str(N) ' points long']);
szX = size(X);


if min(szX) > 1
    Y = zeros(szX, 'single')+i*zeros(szX, 'single');
    for ii = 1:szX(1)
        tmp = (conv(X(ii,:),filt));
        if ~flag
            Y(ii,:) = tmp(N./2:szX(2)+N./2-1);
        end
        
        if flag
            Y(ii,:) = tmp(N:szX(2)+N-1);
        end
        ph(ii,:) = phase(Y(ii,:));
        
        if complexflag == 0
            Y(ii,:) = real(Y(ii,:));
        end
        
        
    end
    
    
else
    Y = (conv(X,filt));
    if ~flag
        Y = Y(N./2:szX(2)+N./2-1);
    end
    
    if flag
        Y = Y(N:szX(2)+N-1);
    end
    ph = phase(Y);
    
    if complexflag == 0
        Y = real(Y);
    end
    
end


if f0 > 0
    Y = 2*Y;
end

