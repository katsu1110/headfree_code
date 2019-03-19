function uneye_goodex(saveoption)
% look for good example

if nargin < 1; saveoption = 0; end
foldername = 'Figure6_Microsaccade';
figpath = 'Z:/Katsuhisa/headfree_project/figures/';
datapath = 'Z:/Katsuhisa/headfree_project/dataset/uneye_pred';
eyepath = 'Z:/Katsuhisa/headfree_project/dataset/eyes_dis_csv';
addpath(genpath('Z:/Katsuhisa/code/integrated/cbrewer'))

% sessiondate = '2016.01.07';
sessiondate = '2016.01.06';

eyex = csvread([eyepath '/' sessiondate '_x.csv']);
eyey = csvread([eyepath '/' sessiondate '_y.csv']);
pred = csvread([datapath '/' sessiondate '_pred.csv']);

j = 605;
% j = 15;
saccc = 0;
while saccc < 2
    j = j + 1;
   [start_idx, end_idx] = consecutive_ones(pred(j,:));
    saccc = length(start_idx);
end
disp(j)
close all;
figure;
cols1 = cbrewer('qual', 'Pastel2', 8);
subplot(2,4, [1 2 5 6])
t = [1:size(eyex, 2)]/500;
plot(t, eyex(j,:), '-', 'color', cols1(7,:), 'linewidth', 0.5)
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
hold on;
% subplot(2,4, [5 6])
plot(t, eyey(j,:), '-', 'color', cols1(8,:), 'linewidth', 0.5)
xlabel('time after stimulus onset (sec)')
ylabel('eye position (deg)')
text(1.2, 0.5, 'x position', 'color', cols1(7,:), 'fontsize', 6)
text(1.2, 0.4, 'y position', 'color', cols1(8,:), 'fontsize', 6)
% set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
subplot(2,4, [3 7])
plot(eyex(j,:), eyey(j,:), '-', 'color', 0.5*[1 1 1], 'linewidth', 0.5)
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
subplot(2,4, [4 8])
vx = [0 diff(eyex(j,:))];
vy = [0 diff(eyey(j,:))];
plot(vx, vy, '-', 'color', 0.5*[1 1 1], 'linewidth', 0.5)
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

cols = cbrewer('qual', 'Dark2', saccc);
for i = 1:saccc
    idx = start_idx(i):end_idx(i);
    subplot(2,4, [1 2 5 6])
    hold on;
    plot(t(idx), eyex(j,idx), '-', 'color', cols(i,:), 'linewidth', 0.75)
%     subplot(2,4, [5 6])
    hold on;
    plot(t(idx), eyey(j,idx), '-', 'color', cols(i,:), 'linewidth', 0.75)
    subplot(2,4, [3 7])
    hold on;
    plot(eyex(j,idx), eyey(j,idx), '-', 'color', cols(i,:), 'linewidth', 0.75)
    subplot(2,4, [4 8])
    hold on;
    plot(vx(idx), vy(idx), '-', 'color', cols(i,:), 'linewidth', 0.75)
end
if saveoption == 1
    savefig(gcf, [figpath foldername '\raw_figs\' sessiondate '_' num2str(j) '.fig'])
end

function [start_idx, end_idx, count] = consecutive_ones(vector)
% https://de.mathworks.com/matlabcentral/fileexchange/34914-consecutive-ones
%[Function Description]
%Finds the number of consecutive ones in a binary signal. Returns the
%starting and ending points of the consecutive set and also the count. Since it
%does not use any for loop. It is pretty fast
%
%[Input]
%vector = A binary signal
%[Output]
%
%Author - Shreyes

temp = diff([0 vector 0]);
start_idx = find(temp == 1);
end_idx = find(temp == -1);
count = end_idx - start_idx;
end_idx = end_idx -1;