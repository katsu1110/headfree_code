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
animals = {'mango', 'kiwi'};
data = [];
for a = 1:length(animals)
    datapath = [mypath '/Katsuhisa/headfree_project/dataset/old_' animals{a} '_fix'];
    listings = dir(datapath);
    listings(1:2) = [];
    lenl = length(listings);
    params = cell(1, lenl);
    ids = cell(1, lenl);
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
           ids{l} = listings(l).name(1:sl(1)-1);
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
    d = ids{1};
    c = 1;
    idx = 1;
    f = 2;
    while f < length(params)
        while strcmp(d, ids{f})
            idx = [idx, f];
            f = f + 1;
        end
        idx = unique(idx);
        data.session(c).date = d;
        p = params{idx(1)};
        e = eyeveclen(idx(1));
        r = reward{idx(1)};
        if length(idx) > 1
            for k = idx(2):idx(end)
                p = concatenate_params(p, params{k});
                e = e + eyeveclen(k);
                r = [r, reward{k}];
            end
        end
        data.session(c).params = p;
        data.session(c).eyedata.reward = r;
        data.session(c).eyedata.eyeveclen = e;
        idx = idx(end)+1;        
        d = ids{f+1};
    end
    
    % autosave
    fixdata = data;
    save([mypath '/Katsuhisa/headfree_project/dataset/fixdata_old_' animals{a} '.mat'], 'fixdata', '-v7.3')
    clearvars fixdata data
    disp([animals{a} ' data saved!'])
end