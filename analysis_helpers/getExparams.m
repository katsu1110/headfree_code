function params = getExparams(ex)
% get experimental parameters from ex-file
len_alltr = length(ex.Trials);
tr = abs([ex.Trials.Reward]) > 0;
len_tr = sum(tr);

% session info
if isfield(ex,'fileID')
    filename = ex.fileName;
elseif isfield(ex,'Header')
    if isfield(ex.Header,'onlineFileName')
        filename = ex.Header.onlineFileName;
    elseif isfield(ex.Header,'Headers')
        filename = ex.Header.fileName;
    else
        try
            filename = ex.Header.fileName;
        catch
            filename = ex.fileName;
        end
    end
else
    filename = ex.fileName;
end

% fixation duration
if isfield(ex,'fix')
    fixDur = ex.fix.duration;
else
    fixDur = ex.fixduration;
end
% correct fixation duration for ORxTF files
if mean(ismember('grating.ORxTF',filename))==1
    fixDur = 2;
end
if fixDur==0.45
    fixDur = 2;
end
    
% degree per pixels
if isfield(ex,'setup')
    if ex.setup.monitorWidth==56
        screenNum = 1;
    else
        screenNum = 2;                
    end
else
    screenNum = ex.screen_number;
end

% degree per pixel
try
    dpp = compute_dpp(ex.setup.viewingDistance, ...
        ex.setup.monitorWidth, ex.setup.screenRect(3:4));
catch % old file?
    try
        dpp = compute_dpp(ex.viewingDistance, ...
            ex.monitorWidth, ex.screenRect(3:4));
    catch
        dpp = compute_dpp(101, 61.5, [1920 1080]); % initial mango sessions
    end
end
     
% fixation window
if isfield(ex,'fix')
   fixWin = dpp*[ex.fix.WinW; ex.fix.WinH];
   % size of fixation point
    fixpoint = dpp*ex.fix.PSz;
else
    fixWin = dpp*[ex.fixWinW; ex.fixWinH];
    % size of fixation point
    fixpoint = dpp*ex.fixPSz;
end

% stimulus position
if isfield(ex,'stim')
    if isfield(ex.stim,'vals')
        on = ex.stim.vals.co;
        if isfield(ex.stim.vals,'x0')
              if isempty(ex.stim.vals.x0) || on==0
                  ex.stim.vals.x0 = 0;
              end
        else
              ex.stim.vals.x0 = 0;
        end
        if isfield(ex.stim.vals,'y0')
              if isempty(ex.stim.vals.y0) || on==0
                  ex.stim.vals.y0 = 0;
              end
        else
              ex.stim.vals.y0 = 0;
        end
        stmpos = [ex.stim.vals.x0; ex.stim.vals.y0];
        stmpos(isempty(stmpos)) = 0;
    else
        stmpos = [3; 0];
    end
else
    stmpos = [0; 0];
end

% output structure
params.len_alltr = len_alltr;
params.len_tr = len_tr;
params.fixdur = fixDur;
if isfield(ex, 'fix')
    if isfield(ex.fix, 'preStimDuration')
        params.prestmdur = ex.fix.preStimDuration;
    else
        params.prestmdur = 0;
    end
    if isfield(ex.fix, 'stimDuration')
        params.stmdur = ex.fix.stimDuration;
    else
        params.stmdur = params.fixdur - params.prestmdur;
    end
else
    params.prestmdur = 0;
    params.stmdur = fixDur;
end
params.screenNum = screenNum;
params.fixwin = fixWin;
params.fixdotsz = fixpoint;
params.stmpos = stmpos;