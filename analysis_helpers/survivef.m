function sv = survivef(x, nbin)
% make a survival function and perform a fitting
[N, sv.edges] = histcounts(x, nbin);
sumn = sum(N);
sv.survival = ones(1, nbin+1);
c = 0;
for n = 1:nbin
    c = c + N(n)/sumn;
    sv.survival(n+1) = 1 - c;
end
sv.survival(end) = sv.survival(end-1);

% fitting a Weibull distribution
% (extention of exponential distribution, no assumption of a constant hazard rate)
c = @(p) cost_wbl(p, sv.edges, sv.survival);
[~, idx] = min(abs(sv.survival - 0.5));
p0 = [sv.edges(idx), 1];
options = optimset('MaxFunEvals', 5000, 'MaxIter', 5000);
sv.fitparams = fminsearch(c, p0, options);

function y = survival_wbl(x, p)
y = 1 - cdf('wbl', x, p(1), p(2));

function c = cost_wbl(p, x, y)
c = sum((survival_wbl(x, p) - y).^2);