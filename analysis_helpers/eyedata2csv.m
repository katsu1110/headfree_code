function eyedata2csv
% convert eyedata to csv

% path =======================================
if mean(ismember('gpfs0', cd))==1
    mypath = '/gpfs01/nienborg/group';
else
    mypath = 'Z:';
end
eyepath = [mypath '/Katsuhisa/headfree_project/dataset/eyes_dis'];
savepath = [mypath '/Katsuhisa/headfree_project/dataset/eyes_dis_csv'];

listings = dir(eyepath);
listings(1:2) = [];
lenl = length(listings);
parfor i = 1:lenl
    % load data
    try
        e = load([eyepath '/' listings(i).name], 'eyemat');
        [x, y] = reshape_eyemat(e.eyemat);
        
        % date name
        idx = strfind(listings(i).name, '_');
        datename = listings(i).name(1:idx(end)-1);

        % to csv
        dlmwrite([savepath '/' datename '_x.csv'], ...
            x, 'delimiter', ',', 'precision', 9)
        dlmwrite([savepath '/' datename '_y.csv'], ...
            y, 'delimiter', ',', 'precision', 9)
        disp([listings(i).name ' saved as csv!'])
    catch
        disp([listings(i).name ' ERR'])
        continue
    end
end

function [x, y] = reshape_eyemat(eyemat)
unitr = unique(eyemat(end, :));
lentr = length(unitr);
trlens = zeros(1, lentr);
for i = 1:lentr
    trlens(i) = sum(eyemat(end,:)==unitr(i));
end
med = median(trlens);
x = nan(lentr, med); y = x; c = zeros(1, lentr);
for i = 1:lentr
    for k = 1:2
        v = eyemat(k, eyemat(end,:)==unitr(i));
        c(i) = length(v) - med;
        if length(v) > med
            v = v(1:med);
        else
            v = [v, v(end)*ones(1, med - length(v))];
        end
        if k==1
            x(i, :) = v;
        else
            y(i, :) = v;
        end
    end
end
close all;
figure;
histogram(c); hold on;
yy = get(gca, 'YLim');
plot(mean(c)*[1 1], yy, '-r')