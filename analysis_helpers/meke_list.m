function l = meke_list(datapath, startname, stopname)
% make a list of sessions
l = dir(datapath);
l(1:2) = [];
start_idx = 1;
stop_idx = length(l);
for i = 1:length(l)
    if strcmp(startname, l(i).name)
        start_idx = i;
    elseif strcmp(stopname, l(i).name)
        stop_idx = i;
    end
end
l = l(start_idx:stop_idx);