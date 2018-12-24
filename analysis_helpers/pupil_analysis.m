function ps = pupil_analysis(eyedata)
% pupil size analysis for a discrimination task
ns = 10000;
for i = 1:length(eyedata.stmstartidx)
    if abs(eyedata.reward(i)) > 0
        ns_cand = eyedata.stmstopidx(i) - eyedata.stmstartidx(i) + 1;
        if ns > ns_cand
            ns = ns_cand;
        end
    end
end

% available reward size
avrew = get_avrew(eyedata);

% matrix (sig, ch, avrew, rew, ps)
pmat = nan(sum(abs(eyedata.reward) > 0), 4+ns);
r = 1;
for i = 1:length(eyedata.stmstartidx)
    if abs(eyedata.reward(i)) > 0
        pmat(r, 1) = eyedata.dc(i);
        pmat(r, 2) = eyedata.ch(i);
        pmat(r, 3) = avrew(i);
        pmat(r, 4) = eyedata.reward(i);
        pstc = eyedata.filtered_p(eyedata.stmstartidx(i):eyedata.stmstopidx(i));
        pmat(r, 5:end) = pstc(1:ns);
        r = r + 1;
    end
end

% time-course with available reward size
modedc = 0;
ps.avrew(1).tc.ntr = sum(pmat(:,3)==0 & pmat(:,1)==modedc);
ps.avrew(1).tc.mean = nanmean(pmat(pmat(:,3)==0 & pmat(:,1)==modedc, 5:end), 1);
ps.avrew(1).tc.sem = nanstd(pmat(pmat(:,3)==0 & pmat(:,1)==modedc, 5:end), [], 1)...
    /sqrt(ps.avrew(1).tc.ntr);
ps.avrew(2).tc.ntr = sum(pmat(:,3)>0 & pmat(:,1)==modedc);
ps.avrew(2).tc.mean = nanmean(pmat(pmat(:,3)>0 & pmat(:,1)==modedc, 5:end), 1);
ps.avrew(2).tc.sem = nanstd(pmat(pmat(:,3)>0 & pmat(:,1)==modedc, 5:end), [], 1)...
    /sqrt(ps.avrew(2).tc.ntr);

% time-course with the time within session
idx = find(pmat(:,3)==0 & pmat(:,1)==modedc);
prc = prctile(idx, [0 20 40 60 80 100]);
for t = 1:5
    ps.timebin(t).tc.ntr = sum(idx >= prc(t) & idx <= prc(t+1));
    ps.timebin(t).tc.mean = nanmean(pmat(idx >= prc(t) & idx <= prc(t+1), 5:end), 1);
    ps.timebin(t).tc.sem = nanstd(pmat(idx >= prc(t) & idx <= prc(t+1), 5:end), [], 1)...
        /sqrt(ps.timebin(t).tc.ntr);
end

% ROC analysis
mindc = min(pmat(:,1));
maxdc = max(pmat(:,1));
for t = 1:ns
    ps.auc(t) = rocN(pmat(pmat(:,3)==0 & pmat(:,1)==mindc, 4+t), ...
        pmat(pmat(:,3)==0 & pmat(:,1)==maxdc, 4+t));
end
