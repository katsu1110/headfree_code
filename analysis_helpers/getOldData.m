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
    e = eyeveclen(1);
    r = reward{1};
    p = params{1};
    d = ids{1};
    c = 1;
    for f = 2:length(params)
        if strcmp(d, ids{f})
            try
                p = concatenate_params(p, params{f});
            catch
                p = params{f};
            end
            e = e + eyeveclen(f);
            r = [r, reward{f}];
        else
            data.session(c).date = ids{f-1};
            data.session(c).params = p;
            data.session(c).eyedata.reward = r;
            data.session(c).eyedata.eyeveclen = e;
            c = c + 1;
            r = [];
            e = 0;
        end
    end
    
    % autosave
    fixdata = data;
    save([mypath '/Katsuhisa/headfree_project/dataset/fixdata_old_' animals{a} '.mat'], 'fixdata', '-v7.3')
    clearvars fixdata data
    disp([animals{a} ' data saved!'])
end