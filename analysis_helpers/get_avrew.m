function avrew = get_avrew(eyedata)
compidx = abs(eyedata.reward) > 0;
rewards = eyedata.reward(compidx);
rewexp = zeros(1, sum(compidx));
for i = 4:sum(compidx)
    % moderate reward expectation trials (normally 0.30)
    if rewards(i-3)==1 && rewards(i-2)==1 && rewards(i-1)==1
        rewexp(i) = 1;
    end
    % high reward expectation trials (normally 0.60)
    if i > 4
        if rewards(i-4)==1 && rewards(i-3)==1 && rewards(i-2)==1 && rewards(i-1)==1
            rewexp(i) = 2;
        end  
    end
end
avrew = zeros(1, length(eyedata.reward));
avrew(compidx) = rewexp;