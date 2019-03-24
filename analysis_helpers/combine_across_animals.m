function combine_across_animals

cols = [[197 25 125]; [77 146 33]; [0 0 0]];
mypath = 'Z:\Katsuhisa\headfree_project\figures\Figure2_ExperimentalParams\raw_figs';
names = {'fixmangoold', 'fixkiwiold', 'kaki_free_fix'};
types = {'ses_behav', 'expparams'};
lenn = length(names);
for t = 1:length(types)
    x = cell(lenn, 1);
    y = cell(lenn, 1);
    for n = 1:lenn
         fig = openfig([mypath '\' types{t} names{n} '.fig'], 'invisible');
         axesObjs = get(fig, 'Children');
         lena = length(axesObjs);         
         for p = 1:lena
            x{n}{p} = axesObjs(lena-p+1).Children.XData;
            y{n}{p} = axesObjs(lena-p+1).Children.YData;
            if contains(names{n}, 'kiwi')
                x{n}{p} = [x{n}{p}(1:4) nan(1, 11) x{n}{p}(5:end)+11];
                y{n}{p} = [y{n}{p}(1:4) nan(1, 11) y{n}{p}(5:end)];
                disp(x{n}{p})
            end
         end
         delete(fig);
    end
    figure;
    for n = 1:lenn
        for p = 1:lena
            subplot(1, lena, p)
            plot(x{n}{p}, y{n}{p}, '-o', 'color', cols(n, :)./255, 'linewidth', 0.1, ...
            'markersize', 1)
            hold on;
        end
    end
    savefig([mypath '/' types{t} '_all.fig'])
    close gcf
end