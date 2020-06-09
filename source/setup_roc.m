function [rocdata,rocgroup] = setup_roc(data,group,c1,c2,flag)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2020 Shaoyu Qiao and Bijan Pesaran
% Pesaran Lab, New York University
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   [rocdata,rocgroup] = setup_roc(data,group,c1,c2,flag);

if nargin<5; flag = 0; end

ind1 = find(group==c1)';
ind2 = find(group==c2)';

rocdata = data([ind1,ind2],:);

if flag
    rocgroup(1:length(ind1)) = 1;
    rocgroup(length(ind1)+1:length(ind1)+length(ind2)) = -1;
else
    if mean(data(ind1))>mean(data(ind2))
        rocgroup(1:length(ind1)) = 1;
        rocgroup(length(ind1)+1:length(ind1)+length(ind2)) = -1;
    else
        rocgroup(1:length(ind1)) = -1;
        rocgroup(length(ind1)+1:length(ind1)+length(ind2)) = 1;
    end
end
