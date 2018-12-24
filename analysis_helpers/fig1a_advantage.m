function fig1a_advantage
% fig1a

addpath(genpath('Z:\Katsuhisa\code\integrated\cbrewer'))
xlabs = {'head-post surgery', 'waiting period (a few months) to ensure logevity of the implant', 'start training', 'with head-fixation'};
leg = {'training with head-fixation', 'training also without head-fixation'};
fz = 6;
ct = cbrewer('div', 'RdGy', 11);
figure;
p = plot([5 5.7], 0.2*[1 1], '-', 'color', ct(10,:), 'linewidth', 3);
p.Color(4) = 0.2;
hold on;
p = plot([0.3 3], 0.1*[1 1], '-', 'color', ct(2,:), 'linewidth', 3);
p.Color(4) = 0.2;
hold on;
p = plot([3.3 5.7], 0.1*[1 1], '-', 'color', ct(2,:), 'linewidth', 3);
p.Color(4) = 0.2;
text(0.5, 0.3, leg{1}, 'color', ct(10,:), 'fontsize', fz)
text(0.5, 0.25, leg{2}, 'color', ct(2,:), 'fontsize', fz)
xlim([0 6])
ylim([0 0.38])
xtickangle(45)
set(gca, 'XTick', [1 3 5], 'XTickLabel', xlabs, 'fontsize', fz)
set(gca, 'YTick', [])
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
set(gcf, 'position', [680   701   309   277])
savefig('Z:\Katsuhisa\headfree_project\figures\Figure1_HeadFreeSetup\raw_figs\advantage.fig')