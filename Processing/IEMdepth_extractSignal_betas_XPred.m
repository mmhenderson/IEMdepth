function [bEst, predlabs, scans] = IEMdepth_extractSignal_betas_XPred(mr)
% extract the signal in each voxel on each trial, using trial-by-trial GLM
% also return the predicted channel response matrix for use during training
% (X_ret)


% define some rtc creation stuff...

root = mr.root;

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

voxGood = ~isnan(nanmean(data_orig,1)); % indexes nonNaN voxels

% edit MMH 8/19/16 to fix BD141
% this line gets rid of any columns with all zero elements
% i am not sure where these columns from (something in the mask from localizer, the threshold might be low and it's taking extra voxels?)
% but this seems to fix the problem, and matches an old data set? didn't
% come up for other subjects...
% if strcmp(mr.subName,'BD141')
%     voxGood = ~isnan(nanmean(data_orig,1)) & std(data_orig,[],1)~=0 & sum(data_orig,1)~=0;
% end


data = mr.data;

runs = mr.r;
%    sortBy = mr.sortBy;

% build a design matrix for stims on the left and stims on the right,
% only based on all runs except the present run
cnt = 1; bEst=[];scans=[];predlabs = [];

for ii=runs
    
    prtExt='_XPred';
    
    pname = [root, mr.subName, '/', mr.rp, num2str(ii), '/', mr.subName, '_', mr.rp, num2str(ii), prtExt,'.prt'];
    %fprintf('loading PRT: %s...\n',pname);
    prt = BVQXfile(pname);
    curD = (data(cnt*mr.NofTRs-mr.NofTRs+1:cnt*mr.NofTRs, :));
    
    
    
    mname = pname;
    mname(end-3:end) = '.mat';
    load(mname);
    
    
    
    rtc = prt.CreateSDM(params);
    
    % build a design matrix to apply to data from the held out scan
    X = [rtc.RTCMatrix, ones(mr.NofTRs,1)];
    
    % compute the beta vals for the tst set...
    b = X\curD;  %inv(X'*X)*(X'*curD);
    b = b(1:end-1,:);   % remove constant term
    
    predlabs = [predlabs; stim.posLabs];
    
    
    bEst = [bEst; b];
    scans = [scans; ones(size(b,1), 1)*ii];
    
    cnt = cnt + 1;
    
    clear p stim;
end

% CONVERT TO CARTESIAN COORDINATES: y * -1
%stimLocs = stimLocs .* repmat([1 -1],size(stimLocs,1),1);

bEst = bEst(:,voxGood);
%X_ret = X_ret/max(X_ret(:));
return