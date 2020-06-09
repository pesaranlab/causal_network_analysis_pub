function pr_op=dp_proj(tapers, sampling, f0)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2020 Shaoyu Qiao and Bijan Pesaran
% Pesaran Lab, New York University
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DP_PROJ Generate a prolate projection operator
%
%  PR_OP = DP_PROJ(TAPERS, SAMPLING, F0)
%
%  Inputs:  TAPERS 	=  Data tapers in [K,TIME], [N,P,K] or [N,W] form.
%	    SAMPLING 	=  Sampling rate for projection operator in Hz. 
%			   	Defaults to 1.
%	    F0		=  Center frequency for projection operator in Hz
%				Defaults to 0.
%
%  Outputs:  PR_OP	=  Projection operator in [TIME,K] form.
%

%   Author: Bijan Pesaran, version date 3/12/98.


if nargin < 2 	sampling = 1;  end
if length(tapers) == 2
   n = tapers(1);
   w = tapers(2);
   p = n*w;
   k = floor(2*p-1);
   tapers = [n,p,k];
end
if length(tapers) == 3 
   tapers(1) = round(tapers(1).*sampling); 
   tapers = dpsschk(tapers);
end
if nargin < 3   f0 = 0; end

% determine parameters and assign matrices
K = size(tapers,2);
N = size(tapers,1);
pr_op = zeros(N,K);

shifter = exp(-2.*pi.*j.*f0.*[1:N]./sampling);
for i = 1:K
	pr_op(:,i) = shifter'.*tapers(:,i);
end
