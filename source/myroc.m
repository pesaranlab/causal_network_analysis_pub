function [p,se,RESULT_S] = myroc(data1,data2,flag,figFlag)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2020 Shaoyu Qiao and Bijan Pesaran
% Pesaran Lab, New York University
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   [p,se,RESULT_S] = myroc(data1,data2,flag)
%
%	data1 = Vector of values for condition 1 for each sample
%	data2 = Vector of values for condition 2 for each sample
%   flag = 0/1 order data
%   figFlag = 1: figure; 0: nofigure
%	p = Choice probability

if nargin < 4
    figFlag = 0;
end

if nargin < 3
    figFlag = 0;
    flag = 0;
end

group = [ones(1,length(data1)) zeros(1,length(data2))];
data = [data1 data2];
c1 = 1; c2 = 0;
[myd,myg] = setup_roc(data',group',c1,c2,flag);

if figFlag
    [p,se,RESULT_S] = roc(myd,myg,'figure');
else
    [p,se,RESULT_S] = roc(myd,myg,'nofigure');
end
