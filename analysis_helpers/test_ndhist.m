function test_ndhist
%%
% test visualization of outputs from ndhist
%

close all;

% high variance data
dh = normrnd(0, 5, [2, 1000]);

% low variance data
dl = normrnd(0, 2, [2, 1000]);

subplot(2,2,1)
[exh, eyh, ph] = ndhist(dh(1,:), dh(2,:));
title('high variance')
subplot(2,2,3)
imagesc(exh, eyh, ph./max(ph(:)))
subplot(2,2,2)
[exl, eyl, pl] = ndhist(dl(1,:), dl(2,:));
title('low variance')
subplot(2,2,4)
imagesc(exl, eyl, pl./max(pl(:)))