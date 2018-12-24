function ms_uncor = microsaccade_sessions(tasktype)
% detect microsaccades in the head_free setup (discrimination task)

if nargin < 1; tasktype = 'dis'; end

% path =======================================
if mean(ismember('gpfs0', cd))==1
    mypath = '/gpfs01/nienborg/group';
else
    mypath = 'Z:';
end
addpath(genpath([mypath '/Katsuhisa/code/integrated/matlab_usefulfunc']))

% list ==========================
datapath = [mypath '/Katsuhisa/headfree_project/dataset/eyes'];
switch tasktype
    case 'fix'
        startname = '2015.08.28';
        stopname = '2015.10.27';
    case 'dis'
        startname = '2015.12.16';
        stopname = '2016.01.13';
end
lists = meke_list(datapath, [startname '_eyemat.mat'], [stopname  '_eyemat.mat']);
lenl = length(lists);

% apply a detection algorithm ===========
out = cell(1, lenl);
% out2 = cell(1, lenl);
nodata = zeros(1, lenl);
name = 'engbert';
parfor i = 1:lenl
    % load eyedata
    e = load([datapath '/' lists(i).name]);
    
    % microsaccade detection
    try
        out{i} = microsaccade_detection(e.eyemat(1,:), e.eyemat(2,:), 500, name, 0);
%         out2{i} = microsaccade_detection(e.eyemat(6,:), e.eyemat(7,:), 500, name, 0);
    catch
        out{i} = nan;
%         out2{i} = nan;
        nodata(i) = 1;
    end
end
for i = 1:lenl
    ms_uncor.session(i).results = out{i};
%     ms_cor.session(i).results = out2{i};
end
ms_uncor.session(nodata==1) = [];
% ms_cor.session(nodata==1) = [];

% autosave
save([mypath  '/Katsuhisa/headfree_project/dataset/ms_uncor_' tasktype '.mat'], 'ms_uncor', '-v7.3');
% save([mypath  '/Katsuhisa/headfree_project/dataset/ms_cor_' tasktype '.mat'], 'ms_cor', '-v7.3');
disp('ms saved!')