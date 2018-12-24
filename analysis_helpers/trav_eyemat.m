function avgmat = trav_eyemat(eyemat)
% trial average of eye data

lentr = length(unique(eyemat(end, :)));
avgmat = nan(size(eyemat, 1), lentr);
for i = 1:lentr
    avgmat(:, i) = mean(eyemat(:, eyemat(end,:)==i), 2);
end