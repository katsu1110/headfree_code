function params = concatenate_params(p1, p2)
% concatenate params
fields = {'len_alltr', 'len_tr', 'fixdur', 'stmdur', 'prestmdur', ...
    'screenNum', 'fixwin', 'fixdotsz', 'stmpos'};
params = p1;
for f = 1:length(fields)
    params.(fields{f}) = [params.(fields{f}), p2.(fields{f})];
end

