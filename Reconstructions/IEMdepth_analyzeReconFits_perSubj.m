% Analyze reconstruction fits for all subjects, using bootstrapping to
% estimate variability of each parameter

% VAV 11/21/2018

clear

% This is the directory where the folders "IEMdepth_trialData" and
% "IEMdepth_chanResp" are located.
root = '/usr/local/serenceslab/maggie/IEMdepth/';

% set this to wherever you've put all the code folders (one level up from
% where this script lives, probably). Might be the same as root, but
% doesn't have to be. 
code_folder = '/usr/local/serenceslab/maggie/mFiles/IEMdepth/';

%%
subj = {'AI','AP','BB','BC','BD','BJ','BM','BN','BO'};
VOIs = {'V1','V2','V3','V4','V3A','V3B','IPS0','LO1','LO2'};
ns = numel(subj); nv = numel(VOIs);

dimlab = {'X','Z'};
sessi = 2;
sepconds = 2;
absStat = 2;
usestereox = 2;
sessStrs = {'_allGoodRuns','_allSess'};
sessStr = sessStrs{sessi};
sepcondstrs = {'_trnAllCond','_trnSepCond','_trnFixTstStim'};
sepcondstr = sepcondstrs{sepconds};
locSignStrs = {'posVoxOnly','allVoxAbs'};
locSignStr = locSignStrs{absStat};
stereostrs = {'_screenx','_stereox'};
stereostr = stereostrs{usestereox};

fitCond = 2;

fn = sprintf('%sIEMdepth_chanResp/%s_%s_avgRecons1D_%s%s%s%s.mat',...
                root,subj{1},VOIs{1},locSignStr,...
                sepcondstr,sessStr,stereostr);
load(fn,'avgRec_1D');
stimLocs = cat(2,avgRec_1D(:).stimLocs)';
[ndim,nrec] = size(stimLocs);
npar = 4;
    
%%

allFits = nan(nv,ns,ndim,nrec,npar);
for ss = 1:ns
    fns = sprintf('%sIEMdepth_reconFits_vy/%s_%s_fits_1D%s%s%s.mat',...
            root,subj{ss},locSignStr,sepcondstr,sessStr,stereostr);
    load(fns, 'fitdat');
    tmp = cat(5, fitdat(:).allFitParams);
    % allFits is nROI x nSubs x nDim x nLocs x nParams
    allFits(:,ss,:,:,:) = permute(tmp(:,fitCond,:,:,:),[5 2 1 3 4]);
    
    clear tmp fitdat
end

% fit error
sl = permute(repmat(stimLocs,1,1,nv,ns),[3 4 1 2]);
fitErr = abs(allFits(:,:,:,:,1) - sl);

%%
allFitCenters = squeeze(allFits(:,:,2,:,1));
roiLab = repmat([1:nv]',1,ns,6);
subLab = repmat(1:ns,nv,1,6);
posLab = repmat(permute(stimLocs(2,:),[3 1 2]),ns,nv,1);
fc = [allFitCenters(:), posLab(:), roiLab(:), subLab(:)];
csvwrite('fitCenters.csv',fc,1);

%% 2 way ANOVA on fit err

% To make the ANOVA table, make a column for the subject
dxz = table([1:ns]','VariableNames',{'Subj'});
% Now make a column for each repeated measure (i.e. ROI)
ii = 1;
for dd = 1:ndim
    for vv = 1:nv
        % Take the mean across the recon positions
        tmp = squeeze(mean(fitErr(vv,:,dd,:),4))';
        dxz = [dxz, table(reshape(tmp,[],1),'VariableNames',{sprintf('Y%d',ii)})];
        ii = ii + 1;
        clear tmp
    end
end

within = table(reshape(repmat({'X','Z'},nv,1),[],1),...
    [VOIs';VOIs'],'VariableNames',{'SpaceDim', 'ROI'});

rmxz = fitrm(dxz,'Y1-Y18~1','WithinDesign',within);

mauchly_xz = mauchly(rmxz)

% run my repeated measures anova here
ranovaxz = ranova(rmxz, 'WithinModel','SpaceDim*ROI')

%% Export a table for loading by RStudio

% To make the ANOVA table that will be exported, make a column for the subject
dp = table([1:ns]','VariableNames',{'Subj'});
% Now make a column for each repeated measure (i.e. ROI)
ii = 1;
for pp = 1:nrec
    for vv = 1:nv
        tmp = squeeze(fitErr(vv,:,2,pp))';
        dp = [dp, table(reshape(tmp,[],1),'VariableNames',{sprintf('Y%d',ii)})];
        ii = ii + 1;
        clear tmp
    end
end

within = table(reshape(repmat({'1','2','3','4','5','6'},nv,1),[],1),...
    reshape(repmat(VOIs',1,nrec),[],1),'VariableNames',{'Position', 'ROI'});

writetable(within,[code_folder 'MixedModels/fitErr2way.csv']);

%% 1-way ANOVA on fit error
% VAV 12/20/2018

oneway = struct('model',[],'mauchly',[],'results',[]);

for dd = 1:ndim
    % To make the ANOVA table, make a column for the subject
    d = table([1:ns]','VariableNames',{'Subj'});
    % Now make a column for each repeated measure (i.e. ROI)
    for vv = 1:nv
        % Take the mean across the recon positions
        tmp = squeeze(mean(fitErr(vv,:,dd,:),4))';
        d = [d, table(reshape(tmp,[],1),'VariableNames',{sprintf('y%d',vv)})];

        clear tmp
    end

    % This was to check ANOVA results against R.
    % writetable(X,'fitErr.csv');

    within = table(VOIs','VariableNames',{'ROI'});

    % fit the repeated measures model
    oneway(dd).model = fitrm(d,'y1-y9~1','WithinDesign',within);

    oneway(dd).mauchly = mauchly(oneway(dd).model);

    % run my repeated measures anova here
    oneway(dd).results = ranova(oneway(dd).model, 'WithinModel','ROI');

    % d2 = table2array(d);
    % [p,tbl,st] = kruskalwallis(d2(:,2:end))
end

%% 2 way ANOVA of fit center x ROI

% To make the ANOVA table, make a column for the subject
dp = table([1:ns]','VariableNames',{'Subj'});
% Now make a column for each repeated measure (i.e. ROI)
ii = 1;
for pp = 1:nrec
    for vv = 1:nv        
        % allFits is nROI x nSubs x nDim x nLocs x nParams
        tmp = squeeze(allFits(vv,:,2,pp,1));
        dp = [dp, table(reshape(tmp,[],1),'VariableNames',{sprintf('Y%d',ii)})];
        ii = ii + 1;
        clear tmp
    end
end

within = table(reshape(repmat({'1','2','3','4','5','6'},nv,1),[],1),...
    reshape(repmat(VOIs',1,nrec),[],1),'VariableNames',{'Position', 'ROI'});

rmp = fitrm(dp,'Y1-Y54~1','WithinDesign',within);

mauchly_mp = mauchly(rmp)

% run my repeated measures anova here
ranova_mp = ranova(rmp, 'WithinModel','Position*ROI')

% in_data: must be a matrix of num_measurements (i.e. subjects) x factor 1
% % (roi) x factor 2 (condition).
% md = permute(squeeze(squeeze(allFits(:,:,2,:,1))),[2,3,1]);
% % md is subj x Position x ROI
% [p_vals, ranova_table, ~] = get_f_dist_parfor(md,1000,1,1)

% mctable = multcompare(rmp,'ROI','By','Position');
% 
% % % each pairing is reported in both directions - look at only one direction
% [~,pmask] = fdr(mctable.pValue, .05);
% posDiffs = mctable.Difference > 0 & pmask;
% sigDiffs = mctable(posDiffs,:)


% dp2 = table2array(dp);
% [p,tbl,st] = kruskalwallis(dp2(:,2:end),within.Position)
% [p,tbl,st] = kruskalwallis(dp2(:,2:end),within.ROI)

%% FIT SUBJ SLOPES
X1 = [stimLocs(1,:)',ones(nrec,1)];
X2 = [stimLocs(2,:)',ones(nrec,1)];

indivSubSlopes = nan(nv,ns,ndim,2);
for vv = 1:nv
    for ss = 1:ns
        tmp = squeeze(allFits(vv,ss,:,:,1));
        indivSubSlopes(vv,ss,1,:) = X1 \ tmp(1,:)';
        indivSubSlopes(vv,ss,2,:) = X2 \ tmp(2,:)';
    end
end

%% 1-way ANOVA on the slopes
sd = table([1:ns]','VariableNames',{'Subj'});
% Now make a column for each repeated measure (i.e. ROI)
for vv = 1:nv
    % Take the mean across the recon positions
    sd = [sd, table(squeeze(indivSubSlopes(vv,:,2,1))',...
        'VariableNames',{sprintf('y%d',vv)})];
end

% non-parametric shuffle ANOVA
% rng(121);
% tic
% [p_vals, ranova_table, ~] = get_f_dist_parfor(squeeze(indivSubSlopes(:,:,2,1))',...
%     1000, 1, 1)
% toc

% This was to check ANOVA results against R.
% writetable(X,'fitErr.csv');

within = table(VOIs','VariableNames',{'ROI'});

% fit the repeated measures model
sloperm = fitrm(sd,'y1-y9~1','WithinDesign',within);

mauchly_tbl = mauchly(sloperm);

% run my repeated measures anova here
slopeanovatbl = ranova(sloperm, 'WithinModel','ROI');

sd2 = table2array(sd)
[p2,tbl2,st2] = kruskalwallis(sd2(:,2:end))

%% BOOTSTRAP THE FIT CENTERS
rs = 3791499;
rng(rs);
niters = 1000;

sorder = randsample(1:length(subj), niters*length(subj), 1);
sorder = reshape(sorder, [niters, length(subj)]);
iterStimLocs = nan(niters,ndim,nrec,nv);

% allFits is nv x ns x ndim x nrec x npar
parfor ii = 1:niters
    % center
    rsi = sorder(ii,:);
    thisf = permute(squeeze(mean(allFits(:,rsi,:,:,1),2)), [2 3 1]);
    % thisf is ndim x nrec x nv
    % save out resampled center fits for each iteration
    iterStimLocs(ii,:,:,:) = thisf;
end

% This shows if the fit centers are close to the real centers.
meanStimLocs = squeeze(mean(iterStimLocs));
ciStimLocs = prctile(iterStimLocs,[2.5,97.5]);
disp('V1');
disp(squeeze(ciStimLocs(:,2,:,1))');
disp('V3AB');
disp(squeeze(ciStimLocs(:,2,:,5))');

% meanDiffLocs is ndim x npos x nv
iterDiffLocs = abs(permute(iterStimLocs,[2,3,4,1]) - repmat(stimLocs,1,1,nv,niters));
meanDiffLocs = squeeze(mean(iterDiffLocs,4));
ciDiffLocs = prctile(iterDiffLocs,[2.5,97.5],4);
posMeanDiffLocs = squeeze(mean(meanDiffLocs,2));
ciPosDiffLocs = squeeze(prctile(mean(iterDiffLocs,2),[2.5,97.5],4));

figure;
for dd = 1:ndim
    d = squeeze(posMeanDiffLocs(dd,:));
    errorbar(1:nv,d,squeeze(ciPosDiffLocs(dd,:,1))-d,...
        squeeze(ciPosDiffLocs(dd,:,2))-d); hold on;
end
ax = gca;
ax.XTickLabel = VOIs;

% X error
[squeeze(posMeanDiffLocs(1,:))',squeeze(ciPosDiffLocs(1,:,1))',squeeze(ciPosDiffLocs(1,:,2))']
% Z error
[squeeze(posMeanDiffLocs(2,:))',squeeze(ciPosDiffLocs(2,:,1))',squeeze(ciPosDiffLocs(2,:,2))']

%% REGRESS SLOPE TO FIT STIM CENTERS

X1 = [stimLocs(1,:)',ones(nrec,1)];
X2 = [stimLocs(2,:)',ones(nrec,1)];
coefStimLocs = nan(niters,ndim,2,nv);

for vv = 1:nv
    t1 = nan(niters,2); t2 = nan(niters,2);
    parfor ii = 1:niters
        tmp = squeeze(iterStimLocs(ii,:,:,vv));
        t1(ii,:) = X1 \ tmp(1,:)';
        t2(ii,:) = X2 \ tmp(2,:)';
    end
    coefStimLocs(:,1,:,vv) = t1;
    coefStimLocs(:,2,:,vv) = t2;
    
    pvalSlopeZero(1,vv) = 2*min([mean(t1(:,1) > 0), mean(t1(:,1) < 0)]);
    pvalSlopeZero(2,vv) = 2*min([mean(t1(:,1) > 0), mean(t1(:,1) < 0)]);
end

meanCoefStimLocs = squeeze(mean(coefStimLocs)); % ndim x nCoefs x nv
ciCoefStimLocs = prctile(coefStimLocs,[2.5,97.5]);
disp('V1');
disp(squeeze(ciCoefStimLocs(:,2,:,1))');
disp('V3AB');
disp(squeeze(ciCoefStimLocs(:,2,:,5))');

%% MAKE A VIOLIN PLOT OF THE ERROR DISTRIBUTIONS 

figure;
h(1) = subplot(2,2,1);
violinplot(squeeze(mean(iterDiffLocs(1,:,:,:),2))',VOIs,...
    'ShowData',false,'ShowMean',false,'ViolinColor',[0.8,0.8,0.8]);
title('Recon X error');
ylim([0 0.6]);

h(2) = subplot(2,2,2);
violinplot(squeeze(mean(iterDiffLocs(2,:,:,:),2))',VOIs,...
    'ShowData',false,'ShowMean',false,'ViolinColor',[0.8,0.8,0.8]);
title('Recon Z error');
match_ylim(h(1:2),[0 1.15]);

%% Make a violin plot of the slope distributions (Figure 6B)
figure;
h(3) = subplot(1,2,1);
hl = line([0 10], [0 0]);
hl.Color = [0 0 0]; hl.LineWidth = 1; hold on;
hl = line([0 10], [1 1]);
hl.Color = [0 0 0]; hl.LineWidth = 1; hold on;
hv = violinplot(squeeze(coefStimLocs(:,1,1,:)),VOIs,...
    'ShowData',false,'ShowMean',false,'ViolinColor',[0.8,0.8,0.8]);
title('Recon X slopes');

h(4) = subplot(1,2,2);
hl = line([0 10], [0 0]);
hl.Color = [0 0 0]; hl.LineWidth = 1; hold on;
hl = line([0 10], [1 1]);
hl.Color = [0 0 0]; hl.LineWidth = 1; hold on;
hv = violinplot(squeeze(coefStimLocs(:,2,1,:)),VOIs,...
    'ShowData',false,'ShowMean',false,'ViolinColor',[0.8,0.8,0.8]);
title('Recon Z slopes');
match_ylim(h(3:4),[-0.25,1.1]);

prepFigForExport;
%%
% outfn = sprintf('%sfigs%sRecon1D_allIEMSlope_singleSubFits.eps', root, filesep);
% fprintf('Saving file...');
% print(gcf, '-depsc','-painters',outfn);
% fprintf('done!\n');

%% PLOT STIM LOC ERROR BARS AND INDIV SUBS (Figure 6A)
cm = plasma(10);

figure;
for dd = 1:2 % ndim
    subplot(1,2,dd);
%     squeeze(coefStimLocs(:,dd,1,:))    
    % squeeze(meanCoefStimLocs(dd,:,vv))'
    
    tmp = squeeze(allFits(:,:,dd,:,1));
    % plot single subj data
    for ss = 1:ns
        dat = mean(abs(squeeze(tmp(:,ss,:)) - repmat(stimLocs(dd,:),nv,1)),2);
        ht = scatter(1:nv, dat, [], cm(ss,:), 'filled'); hold on;
        alpha(0.25);
    end
    
	itererr = squeeze(mean(iterDiffLocs(dd,:,:,:),2));
    m = squeeze(mean(itererr,2));
    ci = (prctile(itererr,[2.5,97.5],2) - m)';
    he = errorbar(1:nv, m, ci(1,:), ci(2,:),'-o');
    he.Color = [0.25 0.25 0.25]; he.LineWidth = 1.5;
    he.LineStyle = 'none';
    he.Parent.XTick = 1:nv;
    he.Parent.XTickLabel = VOIs;
end

prepFigForExport;
% outfn = sprintf('%sfigs%sRecon1D_IEMErr_singleSubFits.eps', root, filesep);
% fprintf('Saving file...');
% print(gcf, '-depsc','-painters',outfn);
% fprintf('done!\n');

%% PAIRWISE VOI SLOPE COMPARISONS
pairedVOIs = combnk(1:nv,2);
fdr_mask = nan(size(pairedVOIs))';

for dd = 1:ndim
    for vi = 1:size(pairedVOIs,1)
        % coefStimLocs is niters x ndim x nCoefs x nv
        s1 = squeeze(coefStimLocs(:,dd,1,pairedVOIs(vi,1)));
        s2 = squeeze(squeeze(coefStimLocs(:,dd,1,pairedVOIs(vi,2))));
        diff_slope = s1 - s2;
        boot_ttest_slope(dd,vi) = 2*(min([mean(diff_slope < 0), ...
            mean(diff_slope > 0)]));
        
        clear s1 s2 diff_slope
    end
end

[~,fdr_mask] = fdr(boot_ttest_slope, 0.05);

diffx = find(fdr_mask(1,:));
fprintf('X encoding is significantly different in %s & %s\n', ...
    VOIs{pairedVOIs(diffx,:)'});

diffz = find(fdr_mask(2,:));
fprintf('Z encoding is significantly different in %s & %s\n', ...
    VOIs{pairedVOIs(diffz,:)'});
disp(' ');
ptmp = find(boot_ttest_slope(2,:) < .05);
fprintf('Without FDR, Z encoding is different in %s & %s\n', ...
    VOIs{pairedVOIs(ptmp,:)'});
disp(boot_ttest_slope(2,ptmp)) % p values

mean_diff_slope = arrayfun(@(vi) mean(squeeze(coefStimLocs(:,2,1,pairedVOIs(vi,1))-...
    squeeze(coefStimLocs(:,2,1,pairedVOIs(vi,2))))), ptmp);
ci_diff_slope = arrayfun(@(vi) prctile(squeeze(coefStimLocs(:,2,1,pairedVOIs(vi,1))-...
    squeeze(coefStimLocs(:,2,1,pairedVOIs(vi,2)))), [2.5,97.5]), ptmp,...
    'UniformOutput',0); 

%% PAIRWISE ERROR COMPARISONS
pairedVOIs = combnk(1:nv,2);
fdr_maskErr = nan(size(pairedVOIs))';

for dd = 1:ndim
    for vi = 1:size(pairedVOIs,1)
        % iterDiffLocs is niters x ndim x nCoefs x nv
        s1 = mean(squeeze(iterDiffLocs(dd,:,pairedVOIs(vi,1),:)),1);
        s2 = mean(squeeze(iterDiffLocs(dd,:,pairedVOIs(vi,2),:)),1);
        diff_err = s1 - s2;
        boot_err(dd,vi) = 2*(min([mean(diff_err < 0), ...
            mean(diff_err > 0)]));
        
        clear s1 s2 diff_err
    end
end

[~,fdr_maskErr] = fdr(boot_err, 0.05);

diffx = find(fdr_maskErr(1,:));
fprintf('X error is significantly different in %s & %s\n', ...
    VOIs{pairedVOIs(diffx,:)'});

diffz = find(fdr_maskErr(2,:));
fprintf('Z error is significantly different in %s & %s\n', ...
    VOIs{pairedVOIs(diffz,:)'});
disp(' ');
ptmp = find(boot_err(2,:) < .05);
fprintf('Without FDR, Z error is different in %s & %s\n', ...
    VOIs{pairedVOIs(ptmp,:)'});
disp(boot_err(2,ptmp)) % p values

mean_diff_err = arrayfun(@(vi) mean(mean(iterDiffLocs(2,:,pairedVOIs(vi,1),:)-...
    iterDiffLocs(2,:,pairedVOIs(vi,2),:),2)), ptmp);
ci_diff_err = arrayfun(@(vi) prctile(mean(iterDiffLocs(2,:,pairedVOIs(vi,1),:)-...
    iterDiffLocs(2,:,pairedVOIs(vi,2),:),2),[2.5,97.5]), ptmp,...
    'UniformOutput',0); 

%% PLOT STIM LOC ERROR BARS AND INDIV SUBS (Figure 5)
cm = plasma(10);
allX = cat(3, X1, X2);

figure; pp = 1;
for dd = 1:2 % ndim
%     figure;
    for vv = 1:nv
        hs(pp) = subplot(2,9,pp); hold all;
        tmp = squeeze(allFits(vv,:,dd,:,1));
        for ss = 1:ns
            ht = scatter(stimLocs(dd,:), tmp(ss,:), [], cm(ss,:), 'filled');
            alpha(0.25);
        end
        plot(stimLocs(dd,:), stimLocs(dd,:), 'Color', [0.8 0.8 0.8], ...
            'LineWidth', 1.5);
        sdat = squeeze(mean(iterStimLocs(:,dd,:,vv)));
        cidat = abs(squeeze(ciStimLocs(:,dd,:,vv)) - sdat');
        he = errorbar(stimLocs(dd,:), sdat, cidat(1,:), cidat(2,:));
        he.Color = [0.25 0.25 0.25]; he.LineWidth = 1.5;
        % plot the mean slope
        plot(stimLocs(dd,:), allX(:,:,dd) * squeeze(meanCoefStimLocs(dd,:,vv))', ...
                ':', 'Color', [0 0 0], 'LineWidth', 2);
            
        title(VOIs{vv});
        axis equal;
        
        if pp == 1
            legend([subj,{'true location', 'IEM estimate', 'line fit to IEM'}]);
        end
        if vv == 1
            ylabel(dimlab{dd});
        end
        
        pp = pp + 1;
    end
%     suptitle(dimlab{dd});
end

match_xlim(hs, [-2, 2]);
match_ylim(hs, [-2,2.25]);
prepFigForExport;
%%
% outfn = sprintf('%sfigs%sRecon1D_XZSlope_singleSubFits.eps', root, filesep);
% fprintf('Saving file...');
% print(gcf, '-depsc','-painters',outfn);
% fprintf('done!\n');

clear hs he

%% BOOTSTRAP THE OTHER PARAMS

% allFits is nv x ns x ndim x nrec x npar
parfor ii = 1:niters
    % center
    rsi = sorder(ii,:);
    thisf = permute(squeeze(mean(allFits(:,rsi,:,:,2:4),2)), [2 3 1 4]);
    % thisf is ndim x nrec x nv x npar
    iterPar(ii,:,:,:,:) = thisf;
end

ciPar = prctile(iterPar,[2.5,97.5]);

%% PLOT PAR ERROR BARS AND INDIV SUBS
parlab = {'center','sz','amp','base'};
plim = [0 16; 0 2; -0.75 1];
cm = plasma(10);
hs = []
% figure;
% ploti = 1;
for dd = 1:2 % ndim
    figure; ploti = 1;
    for pp = 1:3
        for vv = 1:nv
            hs(vv) = subplot(3,9,ploti);
            tmp = squeeze(allFits(vv,:,dd,:,pp+1));
            sdat = squeeze(mean(iterPar(:,dd,:,vv,pp)));
            cidat = abs(squeeze(ciPar(:,dd,:,vv,pp)) - sdat');
            he = errorbar(stimLocs(dd,:), sdat, cidat(1,:), cidat(2,:));
            he.Color = 'k'; he.LineWidth = 2; hold all;
            for ss = 1:ns
                ht = scatter(stimLocs(dd,:), tmp(ss,:), [], cm(ss,:), 'filled');
                alpha(0.25);
            end
            if pp == 1
                title(VOIs{vv});
            end
            if vv == 1
                ylabel(parlab{pp+1});
            end
            ploti = ploti + 1;
        end
        match_ylim(hs, plim(pp,:));
        suptitle(sprintf('%s fit parameters', dimlab{dd}));
        
        clear hs;
    end
    prepFigForExport;
end
