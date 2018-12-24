function eval_uneye(saveoption)
% evaluate 'uneye' s performance

if nargin < 1; saveoption = 1; end
foldername = 'Figure6_MIcrosaccade';
figpath = 'Z:\Katsuhisa\headfree_project\figures\';
datapath = 'Z:\Katsuhisa\headfree_project\dataset\uneye_pred';
addpath(genpath('Z:\Katsuhisa\code\integrated\'))

listings = dir(datapath);
listings(1:2) = [];
lenl = length(listings);
cols = cbrewer('qual', 'Set3', lenl/2);
close all;
amp = [];
peakv = [];
rate = zeros(1, lenl/2);
for i = 1:lenl/2
    pred = csvread([datapath '/' listings(2*i-1).name]);
    [ntr, nf] = size(pred);
    mat = csvread([datapath '/' listings(2*i).name]);
    nans = isnan(mat(1,:)) | isnan(mat(2,:));
    mat(:, nans) = [];
    mat(2, :) = mat(2, :)*500; % to velocity (deg/sec)
    rate(i) = 500*length(mat(1,:))/(ntr*nf);
%     figure(1);
%     subplot(4, 4, i)
%     plot(mat(1,:), mat(2,:), 'o', 'color', cols(i,:), 'markersize', 1)
%     set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
%     rr = corrcoef(mat(1,:), mat(2,:));
%     title({[num2str(i) ': r(pea) = ' num2str(rr(1,2))], [', n=' num2str(length(mat(1,:))) ', rate = ' num2str(rate(i))]})
    amp = [amp, mat(1,:)];
    peakv = [peakv, mat(2,:)];
end
figure(2);
[rr, pp] = corrcoef(amp, peakv);
% set(gca, 'XScale', 'log');set(gca, 'YScale', 'log');

% autosave figures
if saveoption==1
    savefig(figure(1), [figpath foldername '\raw_figs\ampVSpeakv_all_uneye.fig'])
end
    
% all collapsed
figure(2); 
[edgex, edgey, N] = ndhist(amp, peakv, 'bins', 3);
N(N==0) = nan;
imagesc(edgex, edgey, N)
set(gca,'YDir','normal')
cols = cbrewer('seq', 'YlGnBu', 99);
cols = [ones(1,3); cols];
xlabel('amplitude (deg)', 'fontsize', 6)
ylabel('peak velocity (deg/sec)', 'fontsize', 6)
title({['r(pea)=' num2str(rr(1,2)) ', p=' num2str(pp(1, 2))],...
    [', n= ' num2str(length(amp)) ', rate = ' num2str(mean(rate)) '+-' num2str(std(rate))]})
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
colormap(cols)
cb = colorbar();
cb.Ruler.Scale = 'log';
cb.Ruler.MinorTick = 'on';
cb.Location = 'eastoutside';
cb.Limits = [1 100];
set(gcf, 'Name', 'amplitude vs peak velocity', 'NumberTitle', 'off')

% autosave figures
if saveoption==1
    savefig(figure(2), [figpath foldername '\raw_figs\ampVSpeakv_uneye.fig'])
end