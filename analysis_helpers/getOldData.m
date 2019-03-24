function getOldData
% extract data and perform basic analysis for the head-free project
% using the intial data from kiwi and mango (fixation task)
%

% path =======================================
if mean(ismember('gpfs0', cd))==1
    mypath = '/gpfs01/nienborg/group';
else
    mypath = 'Z:';
end

addpath(genpath([mypath '/Katsuhisa/headfree_project']))
addpath(genpath([mypath '/Katsuhisa/code/integrated/matlab_usefulfunc']))

% animals
animals = {'kiwi', 'mango'};
data = [];
for a = 1:length(animals)
    datapath = [mypath '/Katsuhisa/headfree_project/dataset/old_' animals{a} '_fix'];
    listings = dir(datapath);
    listings(1:2) = [];
    lenl = length(listings);
    params = cell(1, lenl);
    ids = cell(1, lenl);
    vals = zeros(1, lenl);
    eyeveclen = zeros(1, lenl);
    reward = cell(1, lenl);
    for l = 1:length(listings)
        % load
        try
            load([listings(l).folder '/' listings(l).name])
        catch
            disp([listings(l).name ' ..load error'])
        end
        % obtain experimental parameters
        try
           params{l} = getExparams(ex);
           sl = strfind(listings(l).name, '_');
           dots = strfind(listings(l).name, '.');
           ids{l} = listings(l).name(dots(1)-4:sl(1)-1);
           vals(l) = date2num(ids{l});
        catch
               disp([listings(l).name ' ...getExParams error'])
        end
        % obtain eyevec length
        for i = 1:length(ex.Trials)
            eyeveclen(l) = eyeveclen(l) + length(ex.Trials(i).Eye.v(1, ~isnan(ex.Trials(i).Eye.v(1,:))));
            reward{l} = [reward{l}, 1*(abs(ex.Trials(i).Reward) > 0)];
        end
    end
    % concatenation
    outs = isempty(params);
    params(outs) = [];
    ids(outs) = [];
    eyeveclen(outs) = [];
    reward(outs) = [];
    vals(outs) = [];
    univals = unique(vals);
    lenu = length(univals);
    for u = 1:lenu
        idx = find(vals == univals(u));
        p = params{idx(1)};
        e = eyeveclen(idx(1));
        r = reward{idx(1)};
        if length(idx) > 1
            for k = 2:length(idx)
                p = concatenate_params(p, params{idx(k)});
                e = e + eyeveclen(idx(k));
                r = [r, reward{idx(k)}];
            end
        end
        data.session(u).date = ids{idx(1)};
        data.session(u).params = p;
        data.session(u).eyedata.reward = r;
        data.session(u).eyedata.eyeveclen = e;
    end
    
%     c = 1;
%     d = ids{1};
%     idx = 1;
%     length(params)
%     for f = 2:length(params)
%         if strcmp(d, ids{f})
%             idx = [idx, f];
%         else
%             p = params{idx(1)};
%             e = eyeveclen(idx(1));
%             r = reward{idx(1)};
%             if length(idx) > 1
%                 for k = idx(2):idx(end)
%                     p = concatenate_params(p, params{k});
%                     e = e + eyeveclen(k);
%                     r = [r, reward{k}];
%                 end
%             end
%             data.session(c).date = d;
%             data.session(c).params = p;
%             data.session(c).eyedata.reward = r;
%             data.session(c).eyedata.eyeveclen = e;
%             c = c + 1;
% %             if f==length(params) && ~strcmp(d, ids{f})
% %                 data.session(c).date = idx{f};
% %                 data.session(c).params = params{idx(f)};
% %                 data.session(c).eyedata.reward = reward{idx(f)};
% %                 data.session(c).eyedata.eyeveclen = eyeveclen(idx(f));
% %             end
%             d = ids{f};
%             idx = f;
%         end
%     end
       
     % manually restore kiwi's session
    if strcmp(animals{a}, 'kiwi')
        dpp = 0.0166;
        dates = {'2014.4.28', '2014.4.30', '2014.5.1', '2014.5.2', '2014.5.6', '2014.5.7', '2014.5.8', ...
            '2014.5.9', '2014.5.12', '2014.5.13'};
        lenmanu = length(dates);
        reward = nan(lenmanu);
        eyeveclen = 500*60*[84, 24, 44, 63, 43, 53, 68, 66, 30, 132];
        len_alltr = [1732, 248, 878, 1107, 619, 100, 678, 700, 193, 576];
        len_tr = {[8 11 15 16 30 68 165 107 19 79 126 20], ...
            [9 7 15 2 8 37 18 11 40 9 17 71],  [6 22 20 16 3 54 17 20 32 15 64 144 27 15 236 56 160 2 89], ...
            [14 13 9 13 48 54 15 81 31 43 65 38 134 50 93 23 122 261], [9 25 5 27 17 30 76 331 44 5 28 22], [31 9 60], ...
            [58 153 175 68 39 192], [11 33 37 160 11 24 48 36 5 137 1 189], ...
            [1 1 50 17 124], [13 10 8 14 5 18 38 17 18 101 92 242]};
        fixdur = {[200 400 400 450 450 450 400 450 520 500 510 520]./1000, ...
            [200 300 400 400 450 450 500 520*ones(1, 5)]./1000, ...
            [150 400 400 400 500 500 500 540 540 560 580 600 620 620 640 660 680 700 700]./1000, ...
            [200 500 500 500 600 650 700 700 700 740 740 780 820 860 900 950 950 1000]./1000, ...
            [800 1000*ones(1, 5) 1100 1100 1100 800 1100 1100]./1000, ...
            [1 1 1], [1 1.1*ones(1, 5)], [0.9 1.1 1.2 1.3 1.4 1.4 1.5 1.6 1.7 1.7 1.5 1.5], ...
            [1 1 0.8*ones(1, 3)], [1 1.1 1.2 1.3 1.4 1.4 1.5*ones(1, 3) 1.6 1.7 1.7]};
        fixwin = {dpp*[100 100 60 60 60 60 60 60 60 60 60 70], ...
            dpp*[100 100 100 70 70 60 60 60 60 80 80 80], dpp*[100 100 85 70 70 70 75  75 80*ones(1, 11)], ...
            dpp*[100 100 85 75 70 70 70 70 80 80 70*ones(1, 8)], dpp*[100 100 50 60 70 70 70 80 80 100 100 80], ...
            dpp*[100 100 100], dpp*[100*ones(1, 4) 75 75], ...
            dpp*[100 80*ones(1,11)], dpp*[100 100 100 85 70], dpp*[100 100 100 100 100 85 85 80 100 100 100 100]};
        
        dates(5) = [];
        lenmanu = length(dates);
        reward(5) = [];
        eyeveclen(5) = [];
        len_alltr(5) = [];
        len_tr(5) = [];
        fixdur(5) = [];
        fixwin(5) = [];
        
        prestmdur = zeros(1, lenmanu);
        stmdur = fixdur;
        screenNum = ones(1, lenmanu);
        fixdotsz = nan(1, lenmanu);
        stmpos = repmat([3; 1], lenmanu);
        for i = 1:lenmanu
            datamanu.session(i).date = dates{i};
            datamanu.session(i).params.len_alltr = len_alltr(i);
            datamanu.session(i).params.len_tr = len_tr{i};
            datamanu.session(i).params.fixdur = fixdur{i};
            datamanu.session(i).params.fixwin = fixwin{i};
            datamanu.session(i).params.prestmdur = prestmdur(i);
            datamanu.session(i).params.stmdur = stmdur(i);
            datamanu.session(i).params.screenNum = screenNum(i);
            datamanu.session(i).params.fixdotsz = fixdotsz(i);
            datamanu.session(i).params.stmpos = stmpos(i);
            
            datamanu.session(i).eyedata.reward = reward(i);
            datamanu.session(i).eyedata.eyeveclen = eyeveclen(i);
        end
        ndata = length(data.session);
        data.session(ndata+1:ndata+lenmanu) = datamanu.session;
    end
    
    % sort chronologically
    data = sort_olddata(data);
    if strcmp(animals{a}, 'mango')
        data.session(1:3) = [];
    end
    
    % autosave
    fixdata = data;
    save([mypath '/Katsuhisa/headfree_project/dataset/fixdata_old_' animals{a} '.mat'], 'fixdata', '-v7.3')
    clearvars fixdata data
    disp([animals{a} ' data saved!'])
end

function num = date2num(dat)
dots = strfind(dat, '.');
num = 30*str2double(dat(dots(1)+1:dots(2)-1))...
        + str2double(dat(dots(2)+1:end));
    
function data = sort_olddata(data)
lend = length(data.session);
outs = zeros(1, lend);
for i = 1:lend
    if isempty(data.session(i).date)
        outs(i) = 1;
    end
end
data.session(outs==1) = [];
lend = length(data.session);
vals = zeros(1, lend);
for i = 1:lend
    vals(i) = date2num(data.session(i).date);
%     dots = strfind(data.session(i).date, '.');
%     vals(i) = 30*str2double(data.session(i).date(dots(1)+1:dots(2)-1))...
%         + str2double(data.session(i).date(dots(2)+1:end));
end
[~, sort_i] = sort(vals, 'ascend');
data.session = data.session(sort_i);
for i = 1:lend
    disp(data.session(i).date)
end
