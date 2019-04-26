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
    
parlab = {'Cent','Size','Ampl','Base'};
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

%% Export a table of each param for loading by RStudio

X_err = [];X_amp = []; X_size = []; X_base = [];
ix = 0;
    
for pp = 1:nrec
    for vv = 1:nv
        for ss = 1:ns
            
            ix = ix+1;
            
            % allFits is nROI x nSubs x nDim x nLocs x nParams
            thiserr = allFits(vv,ss,2,pp,1) - stimLocs(2,pp);
            X_err(ix,:) = [abs(thiserr), vv, pp,ss];                        
              
            thisamp= allFits(vv,ss,2,pp,3);
            X_amp(ix,:) = [thisamp, vv, pp,ss];
            
            thissize= allFits(vv,ss,2,pp,2);
            X_size(ix,:) = [thissize, vv, pp,ss];

            thissize= allFits(vv,ss,2,pp,4);          
            X_base(ix,:) = [thissize, vv, pp,ss];
            
        end
    end
end

table_to_save = array2table(X_err, 'VariableNames',{'err','ROI','position','subject'});
if ~isfolder([code_folder 'MixedModels'])
    mkdir([code_folder 'MixedModels']);
end
writetable(table_to_save,[code_folder 'MixedModels/recon_err_tbl.txt']);

table_to_save = array2table(X_amp, 'VariableNames',{'amp','ROI','position','subject'});
if ~isfolder([code_folder 'MixedModels'])
    mkdir([code_folder 'MixedModels']);
end
writetable(table_to_save,[code_folder 'MixedModels/recon_amp_tbl.txt']);

table_to_save = array2table(X_size, 'VariableNames',{'size','ROI','position','subject'});
if ~isfolder([code_folder 'MixedModels'])
    mkdir([code_folder 'MixedModels']);
end
writetable(table_to_save,[code_folder 'MixedModels/recon_size_tbl.txt']);

table_to_save = array2table(X_base, 'VariableNames',{'base','ROI','position','subject'});
if ~isfolder([code_folder 'MixedModels'])
    mkdir([code_folder 'MixedModels']);
end
writetable(table_to_save,[code_folder 'MixedModels/recon_base_tbl.txt']);

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

%% 1-way ANOVA on the slopes (not significant)
sd = table([1:ns]','VariableNames',{'Subj'});
% Now make a column for each repeated measure (i.e. ROI)
for vv = 1:nv
    % Take the mean across the recon positions
    sd = [sd, table(squeeze(indivSubSlopes(vv,:,2,1))',...
        'VariableNames',{sprintf('y%d',vv)})];
end

within = table(VOIs','VariableNames',{'ROI'});

% fit the repeated measures model
sloperm = fitrm(sd,'y1-y9~1','WithinDesign',within);

mauchly_tbl = mauchly(sloperm);

% run my repeated measures anova here
slopeanovatbl = ranova(sloperm, 'WithinModel','ROI');

sd2 = table2array(sd);
% [p2,tbl2,st2] = kruskalwallis(sd2(:,2:end));


%% BOOTSTRAP ALL THE PARAMETERS
rs = 3791499;
rng(rs);
niters = 1000;

sorder = randsample(1:length(subj), niters*length(subj), 1);
sorder = reshape(sorder, [niters, length(subj)]);
iterPars = nan(niters,ndim,nrec,nv,npar);
% allFits is nROI x nSubs x nDim x nLocs x nParams
% iterPars is nIters x nDim x nRec x nROIs x nParams

parfor ii = 1:niters
    for pp=1:npar
        % center
        rsi = sorder(ii,:);
        % allFits is nROI x nSubs x nDim x nLocs x nParams
        thisf = permute(squeeze(mean(allFits(:,rsi,:,:,pp),2)), [2 3 1]);
        % thisf is ndim x nrec x nv
        % save out resampled center fits for each iteration
        iterPars(ii,:,:,:,pp) = thisf;
    end
end

%% Calculate error of the centers
% iterCenters is niter x ndim x nrec x nv 
% iterDiffLocs is nDim x npos x nv x niter
% meanDiffLocs is ndim x npos x nv
% ciDiffLocs is ndim x npos x nv x 2
iterCenters = squeeze(iterPars(:,:,:,:,1));
iterDiffLocs = abs(permute(iterCenters,[2,3,4,1]) - repmat(stimLocs,1,1,nv,niters));
meanDiffLocs = squeeze(mean(iterDiffLocs,4));
ciDiffLocs = prctile(iterDiffLocs,[2.5,97.5],4);

% posMeanDiffLocs = squeeze(mean(meanDiffLocs,2));
% ciPosDiffLocs = squeeze(prctile(mean(iterDiffLocs,2),[2.5,97.5],4));

% mean_ci_1 = squeeze(mean(ciDiffLocs(1,:,:,:),2));
% 
% mean_ci_2 = prctile(squeeze(mean(iterDiffLocs(1,:,:,:),2)), [2.5, 97.5], 2);

%% print out a couple of summary stats for putting in text

% amp averaged over position
iterVals = squeeze(mean(iterPars(:,:,:,:,3),3));
meanVals = squeeze(mean(iterVals,1));
ciVals = prctile(iterVals, [2.5, 97.5]);

avg_low_high_x = [meanVals(1,:)', squeeze(ciVals(:,1,:))'];
avg_low_high_z = [meanVals(2,:)', squeeze(ciVals(:,2,:))'];
% ci_x = squeeze(mean(ciAmps(1,:,::),2));

%% REGRESS SLOPE TO FIT STIM CENTERS

X1 = [stimLocs(1,:)',ones(nrec,1)];
X2 = [stimLocs(2,:)',ones(nrec,1)];
coefStimLocs = nan(niters,ndim,2,nv);

for vv = 1:nv
    t1 = nan(niters,2); t2 = nan(niters,2);
    parfor ii = 1:niters
        tmp = squeeze(iterCenters(ii,:,:,vv));
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

%% PLOT STIM LOC ERROR BARS AND INDIV SUBS (Figure 5)

% get the confidence intervals - iterCenters is niters x ndim x nrec x nv
ciCenters = prctile(iterCenters,[2.5,97.5]);

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
        sdat = squeeze(mean(iterCenters(:,dd,:,vv)));
        cidat = abs(squeeze(ciCenters(:,dd,:,vv)) - sdat');
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

%% PLOT RECON CENTER ERROR (W/ INDIV SUBS), BY ROI

cm = plasma(10);

figure;
for dd = 1:2 % ndim
    subplot(1,2,dd);

    % allFits is nROI x nSubs x nDim x nLocs x nParams
    tmp = squeeze(allFits(:,:,dd,:,1));
    % plot single subj data
    for ss = 1:ns
        dat = mean(abs(squeeze(tmp(:,ss,:)) - repmat(stimLocs(dd,:),nv,1)),2);
        ht = scatter(1:nv, dat, [], cm(ss,:), 'filled'); hold on;
        alpha(0.25);
    end
    
    % iterDiffLocs is nDim x npos x nv x niter
    iterVals = squeeze(mean(iterDiffLocs(dd,:,:,:),2))';
    meanVals = squeeze(mean(iterVals,1));
    ciVals = prctile(iterVals, [2.5, 97.5]) - repmat(mean(iterVals,1), 2, 1);
    m = meanVals; ci=ciVals;

    he = errorbar(1:nv, m, ci(1,:), ci(2,:),'-o');
    he.Color = [0.25 0.25 0.25]; he.LineWidth = 1.5;
    he.LineStyle = 'none';
    he.Parent.XTick = 1:nv;
    he.Parent.XTickLabel = VOIs;
end

prepFigForExport;

outfn = sprintf('%sfigs%sRecon1D_AbsErrorByROI_singleSubFits.eps', root, filesep);
fprintf('Saving file...');
print(gcf, '-depsc','-painters',outfn);
fprintf('done!\n');

%% PLOT RECON CENTER ERROR (W/ INDIV SUBS), BY POSITION

cm = plasma(10);

figure;
dd=2;
% for dd = 1:2 % ndim
%     subplot(1,2,dd);
%     squeeze(coefStimLocs(:,dd,1,:))    
    % squeeze(meanCoefStimLocs(dd,:,vv))'
    
    tmp = squeeze(allFits(:,:,dd,:,1));
    % plot single subj data
    for ss = 1:ns
        dat = mean(abs(squeeze(tmp(:,ss,:)) - repmat(stimLocs(dd,:),nv,1)),1);
        ht = scatter(1:nrec, dat, [], cm(ss,:), 'filled'); hold on;
        alpha(0.25);
    end
    
     % iterDiffLocs is nDim x npos x nv x niter
    iterVals = squeeze(mean(iterDiffLocs(dd,:,:,:),3))';
    meanVals = squeeze(mean(iterVals,1));
    ciVals = prctile(iterVals, [2.5, 97.5]) - repmat(mean(iterVals,1), 2, 1);
    m = meanVals; ci=ciVals;
    
    he = errorbar(1:nrec, m, ci(1,:), ci(2,:),'-o');
    he.Color = [0.25 0.25 0.25]; he.LineWidth = 1.5;
    he.LineStyle = 'none';
    he.Parent.XTick = 1:nrec;
    he.Parent.XTickLabel = stimLocs(2,:);
    xlim([0,7]);
% end

prepFigForExport;

outfn = sprintf('%sfigs%sRecon1D_AbsErrorByPos_singleSubFits.eps', root, filesep);
fprintf('Saving file...');
print(gcf, '-depsc','-painters',outfn);
fprintf('done!\n');
%% PLOT AMPLITUDE, SIZE, BASELINE (by ROI)
parlims = [NaN, NaN; 2.5,8.5; 0.25, 1.1; -0.4, 0.2];

for pp = 2:4

    cm = plasma(10);

    figure;
    for dd = 1:2 % ndim
        subplot(1,2,dd);
        
        % allFits is nROI x nSubs x nDim x nLocs x nParams
        tmp = squeeze(allFits(:,:,dd,:,pp));
        % plot single subj data
        for ss = 1:ns
            dat = mean(squeeze(tmp(:,ss,:)),2);
            ht = scatter(1:nv, dat, [], cm(ss,:), 'filled'); hold on;
            alpha(0.25);
        end

%         iterPars is nIters x nDim x nRec x nROIs x nParams
%         meanAmps is [nDim x nRec x nVOIs];
%         ciAmps is [2 x nDim x nRec x nVOIs];
        iterVals = squeeze(mean(iterPars(:,dd,:,:,pp),3));
        meanVals = squeeze(mean(iterVals,1));
        ciVals = prctile(iterVals, [2.5, 97.5]) - repmat(mean(iterVals,1), 2, 1);
        m = meanVals; ci=ciVals;

        he = errorbar(1:nv, m, ci(1,:), ci(2,:),'-o');
        he.Color = [0.25 0.25 0.25]; he.LineWidth = 1.5;
        he.LineStyle = 'none';
        he.Parent.XTick = 1:nv;
        he.Parent.XTickLabel = VOIs;
        he.Parent.XTickLabelRotation = 90;
        title(dimlab{dd});
        ylim(parlims(pp,:));
    end
    suptitle([parlab{pp} ' , plotted by ROI']);

    prepFigForExport;

    outfn = sprintf('%sfigs%sRecon1D_%sByROI_singleSubFits.eps', root, filesep, parlab{pp});
    fprintf('Saving file...');
    print(gcf, '-depsc','-painters',outfn);
    fprintf('done!\n');
end

%% PLOT AMPLITUDE, SIZE, BASELINE (by Position)
parlims = [NaN, NaN; 2.5,8.5; 0.25, 1.1; -0.4, 0.2];

for pp = 2:4

    % allFits is nROI x nSubs x nDim x nLocs x nParams
    % iterPars is nIters x nDim x nRec x nROIs x nParams
    % meanAmps is [nDim x nRec x nVOIs];
    % ciAmps is [2 x nDim x nRec x nVOIs];
    iterVals = squeeze(iterPars(:,:,:,:,pp));
    meanVals = squeeze(mean(iterVals,1));
    ciVals = prctile(iterVals, [2.5, 97.5]) - repmat(mean(iterVals,1), 2, 1,1,1);

    % same thing, but plotted by position instead of ROI
    cm = plasma(10);

    figure;
    dd=2;
%     for dd = 1:2 % ndim
%         subplot(1,2,dd);

        tmp = squeeze(allFits(:,:,dd,:,pp));
        % plot single subj data
        for ss = 1:ns
             dat = mean(squeeze(tmp(:,ss,:)),1);
            ht = scatter(1:nrec, dat, [], cm(ss,:), 'filled'); hold on;
            alpha(0.25);
        end

%         iterPars is nIters x nDim x nRec x nROIs x nParams
%         meanAmps is [nDim x nRec x nVOIs];
%         ciAmps is [2 x nDim x nRec x nVOIs];
        iterVals = squeeze(mean(iterPars(:,dd,:,:,pp),4));
        meanVals = squeeze(mean(iterVals,1));
        ciVals = prctile(iterVals, [2.5, 97.5]) - repmat(mean(iterVals,1), 2, 1);
        m = meanVals; ci=ciVals;

        he = errorbar(1:nrec, m, ci(1,:), ci(2,:),'-o');
        he.Color = [0.25 0.25 0.25]; he.LineWidth = 1.5;
        he.LineStyle = 'none';
        he.Parent.XTick = 1:nv;
        he.Parent.XTickLabel = stimLocs(2,:);
        title(dimlab{dd});
        xlim([0,7]);
%         ylim(parlims(pp,:));
%     end
    suptitle([parlab{pp} ' , plotted by position']);

    prepFigForExport;

    outfn = sprintf('%sfigs%sRecon1D_%sByPos_singleSubFits.eps', root, filesep, parlab{pp});
    fprintf('Saving file...');
    print(gcf, '-depsc','-painters',outfn);
    fprintf('done!\n');

end
%% PAIRWISE VOI SLOPE COMPARISONS ( some significant before FDR)
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

%% PAIRWISE ERROR COMPARISONS (bootstrap method, not significant)
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

