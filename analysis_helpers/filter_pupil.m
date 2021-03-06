function eyedata = filter_pupil(eyedata, method)
if nargin < 2; method = 'urai'; end

% preprocess pupil size (and eye)
samprate = 500;
trange = 100; % ms 
vlen = length(eyedata.x);
eyes = {'x','y','p'};

switch method
    case 'urai'
        % blink detection and interpolation
        blinkrange = [];
        blinks = find(abs(eyedata.p) > 4.8);
        for b = 1:length(blinks)
            start = blinks(b) - samprate*trange*0.001;
            if start <= 0 
                start = 1;
            end
            stop = blinks(b) + samprate*trange*0.001;
            if stop > vlen
                stop = vlen;
            end
            blinkrange = [blinkrange, start:stop];
        end    
        for d = 1:length(eyes)
            eyedata.(eyes{d})(blinkrange) = nan;
            eyedata.(eyes{d}) = nan_interp(eyedata.(eyes{d}));
        end
        % filter pupil size
        % dim = 3; cutoff = [0.05 4]; % de Gee et al., 2014
        dim = 2; cutoff = [0.01 10]; % Urai et al., 2017
        [B, A] = butter(dim, 2*cutoff/samprate);
        eyedata.filtered_p = filter(B, A, eyedata.p);

        % z-scoring
        stmps = [];
        for i = 1:length(eyedata.stmstartidx)
            if abs(eyedata.reward(i)) > 0
                stmps = [stmps, eyedata.filtered_p(eyedata.stmstartidx(i):eyedata.stmstopidx(i))];
            end
        end
        eyedata.filtered_p = (eyedata.filtered_p - nanmean(stmps))...
            /nanstd(stmps);
    case 'kawaguchi'
        eyedata.filtered_p = eyedata.p;
        ntr = sum(abs(eyedata.reward) > 0);
        % mean 
        me = nan(1, ntr);
        c = 1;
        for i = 1:length(eyedata.stmstartidx)
            if abs(eyedata.reward(i)) > 0
               me(c) = nanmean(eyedata.p(eyedata.stmstartidx(i):eyedata.stmstopidx(i)));
               c = c + 1;
            end
        end
        me = HiPaFi(me);
        c = 1;
        psvec = [];
        for i = 1:length(eyedata.stmstartidx)
            if abs(eyedata.reward(i)) > 0
                % high pass filter
               ps = eyedata.p(eyedata.stmstartidx(i):eyedata.stmstopidx(i)) - me(c);
               
               % low pass filter
               ps = LoPaFi(ps, 1);
               eyedata.filtered_p(eyedata.stmstartidx(i):eyedata.stmstopidx(i)) = ps;
               psvec = [psvec, ps];
               c = c + 1;
            end
        end
        psme = nanmean(psvec); pssd = nanstd(psvec);
        % z-scoring
        for i = 1:length(eyedata.stmstartidx)
            eyedata.filtered_p(eyedata.stmstartidx(i):eyedata.stmstopidx(i)) ...
                = (eyedata.filtered_p(eyedata.stmstartidx(i):eyedata.stmstopidx(i)) - psme)/pssd;
        end
end