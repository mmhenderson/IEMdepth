function IEMdepth_fitRecons1D(varargin)
% the newer version of the fit function, based on the output of my average
% recon script (IEMdepth_average1DRecon.m)
% VAV 2/29/2016

% this version uses an exponentiated cosine function, implemented in
% eval_cos_1D_new
% also fits the subject-averaged recons, loaded from a previously averaged
% file allSubj_avgRecons1D
% also centers the recons based on their actual stimulus locations in
% stimLocs - saves these as allFitsCent, and saves the difference between
% center and real center as a fifth fit parameter

if nargin<4
    sessi = 1;
    sepconds = 1;
    absStat = 1;
    usestereox = 1;
end
%%
subj = {'AI','AP','BB','BC','BD','BJ','BM','BN','BO'};
VOIs = {'V1','V2','V3','V4','V3A','V3B','IPS0','LO1','LO2'};

root = '/usr/local/serenceslab/maggie/IEMdepth/';

sessStrs = {'_allGoodRuns','_allSess'};
sessStr = sessStrs{sessi+1};

sepcondstrs = {'_trnAllCond','_trnSepCond','_trnFixTstStim'};
sepcondstr = sepcondstrs{sepconds+1};

locSignStrs = {'posVoxOnly','allVoxAbs'};
locSignStr = locSignStrs{absStat+1};

stereostrs = {'_screenx','_stereox'};
stereostr = stereostrs{usestereox+1};

%%
% grid fit params
centers = -2.75:0.1:2.75;
sizes = 1.5:0.1:15;
[params1, params2] = meshgrid(centers,sizes);
params = [params1(:),params2(:)];

% uses fit constraints through gridfit_finetune_constr
% - constr: a 2 x size(bf_grid,2) matrix of [lower bounds; upper bounds] on
%   each optimizatin parameter. use -inf/inf for non-constrained
%   upper/lower bounds

% - gridstep: a vector (1 x n_fit_params - 2) specifying how much the fit
%   can deviate from the best fit in bf_grid. does not include amplitude &
%   baseline, which have separate constraints (passed in as constr)
% - constr: a 2 x 2 matrix of [lower bounds; upper bounds] on
%   the amplitude & baseline of the fits.
% fit_bounds = [-3, 0.1, 0, -5; 3, 20, 5, 5];   % changed to exclude baseline (fixed at zero)
ab_bounds = [0, -5; 5, 5];
gsteps = [centers(3)-centers(1), sizes(3)-sizes(1)];

%add an extra "subject" entry here so we can load the averaged recons, will
%fit them in the same way as single-subject
% subj=[subj,'allSubj'];

% when we center the recons - this is the number of points we take to left
% and right of the center point
% ptsUseCent = 13; % this is the maximum number to not get nans

%% loop over subjects
for ss = 1:length(subj)
    
    allbadfits = nan(length(VOIs),1);
    % load the averaged recons for this subject
    % BM_allVoxAbs_allDors_channelResp_1D_best100000ZVox_trnFixTstStim_allSess_stereox.mat
%     fn = sprintf('%sIEMdepth_chanResp/%s_%s_avgRecons1D_%s%s%s.mat',...
%             root,subj{ss},locSignStr,sepcondstr,sessStr,stereostr);
        
    for vv = 1:length(VOIs)
        
        fn = sprintf('%sIEMdepth_chanResp/%s_%s_avgRecons1D_%s%s%s%s.mat',...
            root,subj{ss},VOIs{vv},locSignStr,sepcondstr,sessStr,stereostr);
        fprintf('loading %s...\n\n',fn);
        load(fn);
        
        recdat = avgRec_1D;

        if ~exist('fitgrid','var')
            % theoretically these should be the same across all like
            % subjects & VOIs...
            basis_grid = recdat.basis_grid;
            fitgrid = make_grid(@eval_cos_1D, basis_grid', params);
            ndim = numel(recdat);
            ncond = size(recdat(1).avgRecons,3);
            nrecons = size(recdat(1).avgRecons,1);
            stimLocs = recdat(1).stimLocs;
        end

        allFits = nan(ndim, ncond, nrecons, length(basis_grid));
%         allFitsCent = nan(ndim, ncond, nrecons, ptsUseCent*2+1);        
        allFitParams = nan(ndim, ncond, nrecons, 4);
        allErr = nan(ndim, ncond, nrecons, 1);
        r2_cos = nan(ndim, ncond, nrecons, 1);
        r2_line = nan(ndim, ncond, nrecons, 1);
        lineFitParams = nan(ndim, ncond, nrecons, 2);

        
%         [bp(:,dd),be(dd),bfn(dd,:)] = gridfit(dat, fitgrid, params,...
%             ab_b(:,1), 0, 0);
%         [bfpar(dd,:),bferr(dd),bffcn(:,dd)] = gridfit_finetune_con(dat(:,dd),...
%             @eval_cos_1D, bp(dd,:), basis_grid, gsteps, ab_b, 0);
        
        for dd = 1:ndim               
            for cc = 1:ncond
                dat = squeeze(recdat(dd).avgRecons(:,:,cc))';
                
                ab_b = ab_bounds;
                ampTmp = max(dat)-min(dat);

                for ii = 1:size(dat,2)
                    ab_b(2,1) = ampTmp(ii);
                    [bp(ii,:),be(ii),bfn(ii,:)] = gridfit(dat(:,ii), fitgrid, params,...
                        ab_b(:,1), 0, 0);
                    [bfpar(ii,:),bferr(ii),bffcn(:,ii)] = gridfit_finetune_con(dat(:,ii),...
                        @eval_cos_1D, bp(ii,:), basis_grid, gsteps, ab_b, 0);
                end

%                 [bp,be,bfn] = gridfit(dat,fitgrid,params,...
%                     ab_bounds(:,1),0,1);
%                 [bfpar,bferr,bffcn] = gridfit_finetune_con(dat,...
%                     @eval_cos_1D, bp, basis_grid, gsteps, ab_bounds, 1);
                
                if any(bferr > be')
                    % error exploded when fine tuning...go back to best
                    % grid fit
                    badfits = find(bferr > be');
                    for i = 1:length(badfits)
                        bfpar(badfits(i),:) = bp(badfits(i),:);
                        bffcn(:,badfits(i)) = bfn(badfits(i),:);
                        bferr(badfits(i)) = be(badfits(i));
                    end
                    allbadfits(vv) = allbadfits(vv) + length(badfits);
                end
                allFits(dd,cc,:,:) = bffcn';
                allFitParams(dd,cc,:,:) = bfpar;
                allErr(dd,cc,:) = bferr;
                
                % also calculate R^2: 1 - (ssres/sstot)
                for r = 1:nrecons
                    ssres1 = sum((dat(:,r) - bffcn(:,r)).^2);
                    sstot1 = sum((dat(:,r) - mean(dat(:,r))).^2);
                    sse_cos(dd,cc,r) = ssres1;
                    r2_cos(dd,cc,r) = 1-(ssres1./sstot1);
                    clear ssres1 sstot1
                end

                clear bp bfpar bferr bffcn

                % also fit with a straight line for comparison, calc R^2
                X = [basis_grid', ones(length(basis_grid),1)];
                for r = 1:nrecons
                    lineFitParams(dd,cc,r,:) = X\dat(:,r);
                    pred = X(:,1)*squeeze(lineFitParams(dd,cc,r,1) ...
                        + lineFitParams(dd,cc,r,2));
                    ssres2 = sum((dat(:,r) - pred).^2);
                    sstot2 = sum((dat(:,r) - mean(dat(:,r))).^2);
                    sse_line(dd,cc,r) = ssres2;
                    r2_line(dd,cc,r) = 1-(ssres2./sstot2);
                    clear pred ssres2 sstot2
                end                
            end
            
        end

        fitdat(vv).allFitParams = allFitParams;
        fitdat(vv).allErr = allErr;
        fitdat(vv).allFits = allFits;
%         fitdat(vv).allFitsCent=allFitsCent;
        fitdat(vv).r2_cos = r2_cos;
        fitdat(vv).r2_line = r2_line;
        fitdat(vv).sse_cos = sse_cos;
        fitdat(vv).sse_line = sse_line;
        fitdat(vv).lineFitParams = lineFitParams;
        
    end

    fns = sprintf('%sIEMdepth_reconFits_vy/%s_%s_fits_1D%s%s%s.mat',...
            root,subj{ss},locSignStr,sepcondstr,sessStr,stereostr);
    save(fns,'fitdat','allbadfits','X','VOIs');
    disp([num2str(nansum(allbadfits)) ' bad fits for sub ' subj{ss}]);
end

end