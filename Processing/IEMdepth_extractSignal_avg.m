function [bEst, scans, stimLocs] = IEMdepth_extractSignal_avg(mr,whichTRs)
% extract the signal in each voxel on each trial, using specific timepoint
% in the BOLD signal after stimulus onset
% also return the conditions (attn task, side, contrast level, target?
% extracts signal from mr, at times of interest after trial start
% (returns one set of vox responses per trial)
% no need to compute design matrix here
%
% adapted from wmDrop_extractSignal_raw.m by TCS 2/6/2015
% - now takes in a set of TRs ([3 4]) and returns the average of those TRs
%   on each trial within each voxel


% bEst, etc is stacked from tpt 0:nTRs-1, use tpt to index into it,
% conds/scans are just repmat(conds/scans,nTRs,1)

root = mr.root;

tpts = whichTRs;
%avg_over_TRs = 1;  % can increase this if we want temporal smoothing


%     shift_TRs = mr.shift_TRs;
%     avg_over_TRs = mr.avg_over_TRs; %average over shiftTRs:shiftTRs+avg_over_TRs-1




% define some rtc creation stuff...  (if doing betas)
params.hpttp  = 5;      %     HRF time to peak for positive response (5)
params.hnttp = 15;      %     HRF5 time to peak for negative response (15)
params.hpnr = 6;        %     HRF positive/negative ratio (6)
params.hons = 0;        %     HRF onset (0)
params.hpdsp = 1;       %     HRF positive response dispersion (1)
params.hndsp = 1;       %     HRF negative response dispersion (1)
params.nvol = mr.NofTRs;%     number of volumes (data points)
params.prtr = mr.TR;    %     TR (in ms, 2000)
params.rnorm = 1;       %     Normalization (default: normalize HRF plateau to 1)
params.rcond = 0;       %     Conditions to remove (rest, 1)

data_orig=double(mr.data);

%index nan voxels
% voxGood = ~isnan(nanmean(data_orig,1)) & std(data_orig,[],1)~=0; % indexes nonNaN voxels
voxGood = ~isnan(nanmean(data_orig,1));
% % edit MMH 8/19/16 to fix BD141
% % this line gets rid of any columns with all zero elements
% % i am not sure where these columns from (something in the mask from localizer, the threshold might be low and it's taking extra voxels?)
% % but this seems to fix the problem, and matches an old data set? didn't
% % come up for other subjects...
% if strcmp(mr.subName,'BD141')
%     voxGood = ~isnan(nanmean(data_orig,1)) & sum(data_orig,1)~=0;
% end

data = mr.data;

runs = mr.r;
%    sortBy = mr.sortBy;

% build a design matrix for stims on the left and stims on the right,
% only based on all runs except the present run



startidx = 1;
%for tt = 1:length(tpts)
cnt = 1;
for ii=runs
    
    pname = [root, mr.subName, '/', mr.rp, num2str(ii), '/', mr.subName, '_', mr.rp, num2str(ii), '.prt'];
    %fprintf('loading PRT %s...\n',pname);
    prt = BVQXfile(pname);
    curD = (data(cnt*mr.NofTRs-mr.NofTRs+1:cnt*mr.NofTRs, :));
    
    mname = pname;
    mname(end-3:end) = '.mat';
    load(mname);
    
    n_trials = size(stim.stimLocs,1);
    
    if ii == 1
        
        nblank = n_trials * length(runs);
        bEst = nan(nblank,size(data,2));
        stimLocs = nan(nblank,3);
        %targPresent = nan(nblank,1);
        %conds = nan(nblank,size(stim.conditions,2)); % check dim2...
        scans = nan(nblank,1);
        %tpt   = nan(nblank,1);
        %targ_idx = nan(nblank,2);
        %tr_coord = nan(nblank,2);
        %tn_coord = nan(nblank,2);
        %tr = nan(nblank,1);
        %tn = nan(nblank,1);
        clear nblank;
        
    end
    
    thisidx = startidx:(startidx+n_trials-1);
    
    try
        for t = 1:length(prt.Cond)
            %braw(t,:) = curD(ceil(prt.Cond(t).OnOffsets(1)/mr.TR),:);
            trial_start_idx = ceil(prt.Cond(t).OnOffsets(1)/mr.TR);
            %disp(sprintf('stim dur = %i, t_idx = :',stim(t).stimExposeDur));
            t_idx = whichTRs + trial_start_idx;%:(trial_start_idx+avg_over_TRs-1)];
            %braw(t,:) = curD(round(prt.Cond(t).OnOffsets(1)/mr.TR)+shift_TRs,:);
            braw(t,:) = mean(curD(t_idx,:),1);
            clear trial_start_idx;
            
        end
        
    catch
        disp('whoops');
    end
    
  
    % TCS - changed from CreateRTC to CreateSDM for BVQXfile v 0.8
    % below is for betas...
    
    %rtc = prt.CreateRTC(params);
    %rtc = prt.CreateSDM(params);
    
    % build a design matrix to apply to data from the held out scan
    %X = [rtc.RTCMatrix, ones(mr.NofTRs,1)];
    
    % compute the beta vals for the tst set...
    %b = X\curD;  %inv(X'*X)*(X'*curD);
    %b = b(1:end-1,:);   % remove constant term
    
    
    %g = [g;gnew'];
    %conds(thisidx,:) = condnew;
    
    stimLocs(thisidx,:) = stim.stimLocs;

    
    bEst(thisidx,:) = braw;
    
    %         targ_idx(thisidx,:) = [stim.TRlocIdx stim.TNlocIdx];
    %         tr_coord(thisidx,:) = stim.TRcoordDeg;
    %         tn_coord(thisidx,:) = stim.TNcoordDeg;
    %         ang_offset(thisidx) = stim.allAngOffset;
    %
    %         tr(thisidx) = trnew;
    %         tn(thisidx) = tnnew;
    
    
    
    %bEst = [bEst; b];
    scans(thisidx) = ones(size(braw,1), 1)*ii ;
    %targPresent(thisidx) = stim.targPresent;
    %tpt(thisidx) = ones(size(braw,1),1)*tpts(tt);
    cnt = cnt + 1;
    startidx = thisidx(end)+1;
    clear stim;
    
end

%end
%stimLocs = stimLocs .* repmat([1 -1],size(stimLocs,1),1);
bEst = bEst(:,voxGood);
end