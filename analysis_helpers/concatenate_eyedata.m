function eyedata = concatenate_eyedata(e1, e2)
% concatenate eyedata1 and 2
eyes = {'x','y','p'};
eyedata = e1;
for d = 1:length(eyes)
    eyedata.(eyes{d}) = [eyedata.(eyes{d}), e2.(eyes{d})];
end
trs = {'trstartidx', 'trstopidx', 'startfixidx', 'stmstartidx', 'stmstopidx'};
for d = 1:length(trs)
    eyedata.(trs{d}) = [eyedata.(trs{d}), e2.(trs{d})];
end
stm = {'dc','hdx','or', 'ch','reward','trlabel','rt'};
for d = 1:length(stm)
    if isfield(eyedata, stm{d})
        eyedata.(stm{d}) = [eyedata.(stm{d}), e2.(stm{d})];
    end
end
if isfield(eyedata, 'choicevec')
    j = length(eyedata.choicevec);
    for i = 1:length(e2.choicevec)
       eyedata.choicevec(j+i).x = e2.choicevec(i).x;
       eyedata.choicevec(j+i).y = e2.choicevec(i).y;
    end
end