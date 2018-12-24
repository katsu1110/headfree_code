function fitted = fit_learning_curve(t, y, nrep)
% fit learning curve with a power low
% INPUT: t ... time
%        y ... performance
%        nrep ... the number of resampling for significance test
%

if nargin < 3; nrep = 0; end

% fitting a power low
c = @(p) cost(p, t, y);
p0 = [y(end) y(1) 0.5];
options = optimset('MaxFunEvals', 5000, 'MaxIter', 5000);
fitted.true_data.params = fminsearch(c, p0, options);

% resampling
if nrep > 0
    lent = length(t);
    fitted.surr_data.params = nan(nrep, 3);
    for r = 1:nrep
        % shuffle
        y_surr = y(randi(lent, 1, lent));
        % fit
        p0 = [y_surr(end) y_surr(1) 0.5];
        fitted.surr_data.params(r,:) = fminsearch(c, p0, options);
    end
    for i = 1:3
        fitted.permtest(i).p(1) = sum(fitted.true_data.params(i) < ...
            fitted.surr_data.params(:,i))/nrep;
        fitted.permtest(i).p(2) = sum(fitted.true_data.params(i) >= ...
            fitted.surr_data.params(:,i))/nrep;
    end
end

% subfunctions
function y = powerlow(x, p)
y = p(1) + p(2)*x.^(-p(3));

function c = cost(p, x, y)
c = sum((powerlow(x, p) - y).^2);