function plotSTA(t,sta,sem, plotColor)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2020 Shaoyu Qiao and Bijan Pesaran
% Pesaran Lab, New York University
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot stimulus/stimulation triggered averages (mean +/- sem)
usem = sta + sem;
lsem = sta - sem;

plot(t,sta,'linewidth',1, 'color', plotColor)
hold on
jbfill(t,usem,lsem,0.6*plotColor,0.6*plotColor);





