function head_free_getData(animal)
% extract data and perform basic analysis for the head-free project
% INPUT: animal ... 'kaki', 'kaki2'
%
% session-by-session analysis:
% (fixation task)
% - fixation duration per trial
% - fixation window (area)
% - total fixation duration (working hours)
% - fixation breaks (%)
% - fixation precision
% - pupil size
%
% (orientation discrimination task)
% ... in addition ...
% - psychometric function
% - pupil size with rewards
% - saccadic choices
% - microsaccade

% input ==================
if nargin < 1; animal = 'kaki'; end

% path =======================================
if mean(ismember('gpfs0', cd))==1
    mypath = '/gpfs01/nienborg/group';
else
    mypath = 'Z:';
end

addpath(genpath([mypath '/Katsuhisa/headfree_project']))
addpath(genpath([mypath '/Katsuhisa/code/integrated/matlab_usefulfunc']))

% data =======================================
switch lower(animal)
    case 'kaki_free' % head-free
        datapath = cell(1,2);
        lists = cell(1,2);

        % fixation task
        datapath{1} = [mypath '/Katsuhisa/headfree_project/fixation task/'];
        startname = '2015.08.28';
        stopname = '2015.10.27';
        lists{1} = meke_list(datapath{1}, startname, stopname);

        % discrimination task
        datapath{2} = [mypath '/Katsuhisa/headfree_project/orientation discrimination task/'];
        % startname = '2015.11.04';
        startname = '2015.12.16';
        % startname = '2016.01.04';
        stopname = '2016.01.13';
        lists{2} = meke_list(datapath{2}, startname, stopname);
    case 'kaki_fixed' % head-fixed
        datapath = cell(1,1);
        lists = cell(1,1);

        % fixation task
        datapath{1} = [mypath '/Katsuhisa/headfree_project/fixation task/'];
        startname = '2016.02.01';
        stopname = '2016.06.04';
        lists{1} = meke_list(datapath{1}, startname, stopname); 
        
        % discrimination task
        datapath{2} = [mypath '/data/kaki/disparity discrimination task/'];
        startname = '2018.10.03';
        stopname = '2018.10.12';
        lists{2} = meke_list(datapath{2}, startname, stopname);
    case 'mango_fixed' % head-fixed
        datapath = cell(1,1);
        lists = cell(1,1);

        % fixation task
        datapath{1} = [mypath '/data/mango/FixTraining_DxTraining/'];
        startname = '2014.11.17';
        stopname = '2014.11.24';
        lists{1} = meke_list(datapath{1}, startname, stopname);  
        
        % discrimination task
        datapath{2} = [mypath '/data/mango/dx training/'];
        startname = '2017.11.01';
        stopname = '2017.11.10';
        lists{2} = meke_list(datapath{2}, startname, stopname);
   case 'kiwi_fixed' % head-fixed
        datapath = cell(1,1);
        lists = cell(1,1);

        % fixation task
        datapath{1} = [mypath '/Katsuhisa/headfree_project/kiwi_fixation/'];
        startname = '2016.05.10';
        stopname = '2016.05.20';
        lists{1} = meke_list(datapath{1}, startname, stopname);  
        
        % discrimination task
        datapath{2} = [mypath '/data/kiwi/dx_training/'];
        startname = '2015.10.12';
        stopname = '2015.10.23';
        lists{2} = meke_list(datapath{2}, startname, stopname);
end
        

% data extraction from sessions =======================
outfnames = {'dailyLog', 'edf', '_all'};
eyevecnames = {'x','y','p','filtered_p'};
% fixation task
for t = 1:length(lists)
    data = [];
%     if t==1
%         continue
%     end
    for l = 1:length(lists{t})
        disp(['working on... ' lists{t}(l).name ' ++++++++++++++++++'])
        sespath = [datapath{t} lists{t}(l).name];
        files = dir(sespath);
        files(1:2) = [];
        lenf = length(files);
        c = 1;
        params = cell(1, lenf);
        eyedata = cell(1, lenf);
        for f = 1:lenf
            filepath = [sespath '/' files(f).name];

            % check the filename
            ok = 1;
            for k = 1:3
                if ~isempty(strfind(files(f).name, outfnames{k}))
                    ok = 0;
                    break
                end
            end
            % load
            if ok==1
                try
                    load(filepath, 'ex')
                    if sum(abs([ex.Trials.Reward]) > 0) < 5
                        continue
                    end
                catch
%                     disp([filepath ' ...load error'])
                    continue
                end
            else
                continue
            end
            
            % obtain experimental parameters
            try
               params{c} = getExparams(ex);
            catch
%                disp([files(f).name ' ...getExParams error'])
               params{c} = nan;
            end
            
            % extract eye data
            try
               eyedata{c} = getEye(ex, t);
               c = c + 1;
            catch
                disp([files(f).name ' ...getEye error'])
                eyedata{c} = nan;
            end
        end

        % concatenation
        try
            nans = zeros(1, length(params));
            for f = 1:length(params)
                if ~isstruct(params{f}) || ~isstruct(eyedata{f})
                    nans(f) = 1;
                end
            end
            % find is needed somehow...
            params(find(nans)) = [];
            eyedata(find(nans)) = [];
            p = params{1};
            e = eyedata{1};
            for f = 2:length(eyedata)
                p = concatenate_params(p, params{f});
                e = concatenate_eyedata(e, eyedata{f});
            end
            data.session(l).date = lists{t}(l).name;
            data.session(l).params = p;
        catch
            data.session(l).params = nan;
            disp(['session ' num2str(l) ': concatenation error'])
        end
        
        % blink detection and filtering
        try
            data.session(l).eyedata = filter_pupil(e, 'kawaguchi');
        catch
            data.session(l).eyedata = e;
            disp(['session ' num2str(l) ': blink detection / filtering error'])
        end
        
        % fixation duration per session
        try
            data.session(l).eyedata.eyeveclen = length(data.session(l).eyedata.x); % raw eye vector length
        catch
            data.session(l).eyedata.eyeveclen = nan;
        end
        
        % survival analysis for fixation breaks
        try
            x = (data.session(l).eyedata.stmstopidx - ...
                data.session(l).eyedata.startfixidx)/500;
            x(x > max(data.session(l).params.fixdur)) = ...
                max(data.session(l).params.fixdur);
            data.session(l).survival = survivef(x, 10);
            data.session(l).raw_median = nanmedian(x);
        catch
            data.session(l).survival = nan;
            disp(['session ' num2str(l) ': survival analysis error'])
        end
        
        % eye data analysis
        try
            if isstruct(data.session(l).eyedata)
                try
                    fixwin = [max(data.session(l).params.fixwin(1,:)), ...
                        max(data.session(l).params.fixwin(2,:))];
                    fixframes = 0.75*min(data.session(l).params.stmdur)*500;
                catch
                    fixwin = [1 1];
                end
                nodata = zeros(1, length(data.session(l).eyedata.stmstartidx));
                eyex = [];
                eyey = [];
                ps = [];
                trs = [];
                meeye = [];
                c = 1;
                for k = 1:length(data.session(l).eyedata.stmstartidx)
                    if abs(data.session(l).eyedata.reward(k)) > 0
                        % eye positions during stimulus presentation period
                        idx1 = data.session(l).eyedata.stmstartidx(k);
                        idx2 = data.session(l).eyedata.stmstopidx(k);                        
                        treyex = data.session(l).eyedata.x(idx1:idx2);
                        treyey = data.session(l).eyedata.y(idx1:idx2);     
                        trps = data.session(l).eyedata.p(idx1:idx2);
                        
                        % check if this is really a completed trial
                        if idx2 - idx1 + 1 < fixframes || std(trps)==0
                            nodata(k) = 1;
                            continue
                        end                                                
                        
                        % smoothing
                        treyex = smooth_eye(treyex, 500);
                        treyey = smooth_eye(treyey, 500);
                                              
                        % adjust spurious eye positions      
                        treyex = treyex - nanmedian(treyex);
                        treyey = treyey - nanmedian(treyey);
                        trps = trps - nanmedian(trps);
                        treyex(abs(treyex) > fixwin(1)) = nan;
                        treyey(abs(treyey) > fixwin(2)) = nan;
                        if all(isnan(treyex)) || all(isnan(treyey))
                            nodata(k) = 1;
                            continue
                        end
                        treyex = nan_interp(treyex);
                        treyey = nan_interp(treyey);
                        trps = nan_interp(trps);
                        
                        % mean eye position
                        meeye = [meeye, [nanmean(treyex); nanmean(treyey)]];
                        
                        % eliminating calibration offset (e.g. Cherici et
                        % al., 2012)
                        treyex = treyex - mean(treyex);
                        treyey = treyey - mean(treyey);
                        
                        % concatenation
                        eyex = [eyex, treyex];
                        eyey = [eyey, treyey];
                        ps = [ps, trps];
                        trs = [trs, c*ones(1, length(treyex))];
                        c = c + 1;
                    end
                end                        
                % reassign completed trials
                fields = {'reward', 'trstartidx', 'startfixidx', 'stmstartidx', ...
                    'trstopidx', 'stmstopidx', 'ch', 'dc', 'hdx', 'or', 'rt', ...
                    'choicevec'};
                for f = 1:length(fields)
                    if isfield(data.session(l).eyedata, fields{f})
                        data.session(l).eyedata.(fields{f})(nodata==1) = [];
                    end
                end
                                
                % pupil artifact on eye positions
                [data.session(l).psarti, ceyes] = eyepos_vs_ps(eyex, eyey, ps);

                % store eye-data for further analysis such as 'U'n'Eye'
                % (https://github.com/berenslab/uneye#content)
                eyemat = [eyex; eyey; ps; ceyes; trs];
                save([mypath '/Katsuhisa/headfree_project/dataset/eyes/' animal '_' lists{t}(l).name '_eyemat.mat'], 'eyemat', '-v7.3')
                clearvars eyemat
            else
                continue
            end
        catch
            data.session(l).microsaccade = nan;
            disp(['session ' num2str(l) ': eye data analysis error'])
        end   
        
        % fixation precision
        try
            if isstruct(data.session(l).eyedata)              
                data.session(l).fixPrecision{1} = fixation_precision(eyex, eyey, 0.75);
                data.session(l).fixPrecision{2} = fixation_precision(ceyes(1,:), ceyes(2,:), 0.75);
                data.session(l).fixPrecision{3} = fixation_precision(ceyes(3,:), ceyes(4,:), 0.75);                
                data.session(l).fixPrecision{4} = fixation_precision(meeye(1,:), meeye(2, :), 0.75);
            else
                continue
            end
        catch
            data.session(l).fixPrecision = nan;
            disp(['session ' num2str(l) ': fixation precision error'])
        end   
        
        % discrimination task
        if t==2
            % PM fitting
            try
                if isfield(data.session(l).eyedata, 'or')
                    signs = ones(1, length(data.session(l).eyedata.or));
                    signs(data.session(l).eyedata.or<45) = -1;
                elseif isfield(data.session(l).eyedata, 'hdx')
                    signs = sign(data.session(l).eyedata.hdx);
                end
                dcs = signs.*data.session(l).eyedata.dc;
                unistm = unique(dcs);
                lend = length(unistm);
                py = nan(1, lend);
                ny = nan(1, lend);
                for d = 1:lend
                    ny(d) = sum(dcs==unistm(d) & abs(data.session(l).eyedata.reward) > 0);
                    py(d) = sum(data.session(l).eyedata.ch==1 & dcs==unistm(d) & ...
                        abs(data.session(l).eyedata.reward) > 0)/ny(d);
                end
                data.session(l).pm = fitPM(unistm,py,ny,'Gaussian','MLE',0);
            catch
                data.session(l).pm = nan;
                disp(['session ' num2str(l) ': PM error'])
            end
            
            % ps analysis
            try
                data.session(l).ps = pupil_analysis(data.session(l).eyedata);
            catch
                data.session(l).ps = nan;
                disp(['session ' num2str(l) ': PS error'])
            end  
            
        end    
        
        % remove eye vectors (as they are heavy...)
        for k = 1:4
            if isfield(data.session(l).eyedata, eyevecnames{k})
                data.session(l).eyedata = rmfield(data.session(l).eyedata,  eyevecnames{k});
            end
        end
        
    end
    
    % autosave 'data'
    if t==1
        fixdata = data;
        save([mypath '/Katsuhisa/headfree_project/dataset/fixdata_indvgain_' animal '.mat'], 'fixdata', '-v7.3')
        clearvars fixdata
    elseif t==2
        disdata = data;
        save([mypath '/Katsuhisa/headfree_project/dataset/disdata_indvgain_' animal '.mat'], 'disdata', '-v7.3')
    end    
    disp('data saved!')
end

%%
% % microsaccade detection for discrimination task
% tlabs = {'fix', 'dis'};
% for t = 1:2
% %     if t==1
% %         continue
% %     end
%     microsaccade_sessions(tlabs{t});
% end