% run through all analyses, make all figures and tables for pairwise
% (six-way) Z position decoding.
% load the classifier output saved by the script IEMdepth_classify_Z6.m

% MMH 3/29/19

%% define subjects and flags for what to do
% close all

subj = {'AI','AP','BB','BC','BD','BJ','BM','BN','BO'};
vuse=[1:7,11:12];
VOIs={'V1','V2','V3','V4','V3A','V3B','IPS0','IPS1','IPS2','IPS3','LO1','LO2'};
for ss=1:length(subj)
    subj{ss} = [subj{ss} '_allROIs'];
end
nSubj=length(subj);
nVOIs=length(vuse);

saveFigs = 1;
figFolder='IEMdepth_figs';
ext='epsc';

typestr = 'Z2_oneVsOne';

%parameters for the classifier
kernelStr='linear';
% functStr='svmtrain (-t 0 -q)';
classStr='svmtrain';
% subMeanStr = 'subMean2';
subMeanStr = 'noSubMean';
usingA=1;
predStrs={'predB','predA'};
predStr=predStrs{usingA+1};
voxelStr = 'allVox';
% voxelStr = 'take150ZVox';
condStrs = {'trainStim','trainFixat'};
conduse = 2;

sigLevels = [0.05,0.01];

%% set up file info, other params

root='/usr/local/serenceslab/maggie/IEMdepth/';

nIter=1000;

folder='IEMdepth_classif';

zDispList = [38.592,25.612,11.124,-5.152,-23.567,-44.573];

posPairList = combnk(1:6,2);
dispPairList = reshape(zDispList(posPairList(:)),15,2);
dispDiffList = dispPairList(:,1)-dispPairList(:,2);
[distListSort,distOrder] = sort(dispDiffList,'ascend');

undist = unique(dispDiffList);
nDist = length(undist);
 

% preallocate arrays to store acc and d'
accs_allsub = zeros(nVOIs,nSubj,nDist);
accsrand_allsub = zeros(nVOIs,nSubj,nDist,nIter);

d_allsub = zeros(nVOIs,nSubj,nDist);
drand_allsub = zeros(nVOIs,nSubj,nDist,nIter);
 
accs_mean_allsub=nan(nVOIs,nSubj);
accsmeanrand_allsub=nan(nVOIs,nSubj,nIter);

d_mean_allsub=nan(nVOIs,nSubj);
dmeanrand_allsub=nan(nVOIs,nSubj,nIter);

% store p vals for significance of decoding 

pValsAcc_mean_allsub=nan(nVOIs,1);
pValsD_mean_allsub = nan(nVOIs,1);

pValsAcc_allsub = nan(nVOIs, nDist);
pValsD_allsub=  nan(nVOIs, nDist);

% this table will be all data in a different format, preparing for export to R script
X = zeros(nSubj*length(vuse)*nDist,4);
ix = 0;

% this array will be for making the pairwise dissimilarity matrix (RSA
% plot)
bigMat = zeros(nSubj,nVOIs,6,6);

%% loop over subs, load all data.

for ss=1:nSubj   
    
    fns=sprintf('%s%s/%s_%s_%s_%s_%s_%s_%s_%s.mat',...
                    root,folder,subj{ss},typestr,condStrs{conduse},voxelStr,predStr,classStr,kernelStr,subMeanStr);

    load(fns);
            
    for vv=1:nVOIs
        
         
        accs_mean_allsub(vv,ss) = mean(classStruct(vuse(vv)).accReal,1);
        accsmeanrand_allsub(vv,ss,:) = mean(classStruct(vuse(vv)).accRand,1);
        
        d_mean_allsub(vv,ss) =  mean(classStruct(vuse(vv)).dReal,1);
        dmeanrand_allsub(vv,ss,:) = mean(classStruct(vuse(vv)).dRand,1);
        
        for dd=1:nDist
            
            theseaccs = classStruct(vuse(vv)).accReal(distOrder(dd),:);
            accs_allsub(vv,ss,dd) = mean(theseaccs);
            
            theseaccs = classStruct(vuse(vv)).accRand(distOrder(dd),:);
            accsrand_allsub(vv,ss,dd,:) = mean(theseaccs,1);
            
            thesed = classStruct(vuse(vv)).dReal(distOrder(dd),:);
            d_allsub(vv,ss,dd) = mean(thesed);
            
            ix = ix+1;
            X(ix,:) = [mean(thesed), vv, distListSort(dd),ss];
            
            thesed = classStruct(vuse(vv)).dRand(distOrder(dd),:);
            drand_allsub(vv,ss,dd,:) = mean(thesed,1);
            
            % also put the values into this big matrix to make a
            % dissimilarity matrix
            gg1 = posPairList(dd,1);
            gg2 = posPairList(dd,2);
            bigMat(ss,vv,gg1,gg2) = classStruct(vuse(vv)).dReal(dd);
            bigMat(ss,vv,gg2,gg1) = classStruct(vuse(vv)).dReal(dd);
            
        end
        
    end
   
end

% export a table of all decoding d' scores that will be loaded by RStudio to
% perform linear mixed model analysis.
table_to_save = array2table(X, 'VariableNames',{'dprime','ROI','disparity','subject'});
writetable(table_to_save,'dprime_classifier_tbl_CORRECTED.txt')
%% compute p vals for each ROI, each cond separately 

for vv=1:nVOIs
        
    % all positions
    realAccs = squeeze(accs_mean_allsub(vv,:));
    nullAccs = squeeze(accsmeanrand_allsub(vv,:,:));

    realD = squeeze(d_mean_allsub(vv,:));
    nullD = squeeze(dmeanrand_allsub(vv,:,:));
   
    pValsAcc_mean_allsub(vv) = 2*min([mean(mean(realAccs)<mean(nullAccs)),mean(mean(realAccs)>mean(nullAccs))]);
    pValsD_mean_allsub(vv) = 2*min([mean(mean(realD)<mean(nullD)),mean(mean(realD)>mean(nullD))]);

     
    clear realAccs nullAccs realD nullD

    % and each position separately
    for dd=1:nDist

        
        realAccs = squeeze(accs_allsub(vv,:,dd));
        nullAccs = squeeze(accsrand_allsub(vv,:,dd,:));

        realD = squeeze(d_allsub(vv,:,dd));
        nullD = squeeze(drand_allsub(vv,:,dd,:));

        pValsAcc_allsub(vv,dd) = 2*min([mean(mean(realAccs)<mean(nullAccs)),mean(mean(realAccs)>mean(nullAccs))]);
        pValsD_allsub(vv,dd) = 2*min([mean(mean(realD)<mean(nullD)),mean(mean(realD)>mean(nullD))]);

         
        clear realAccs nullAccs realD nullD

    end
end

%% FDR correction of all these p-values
isSigAcc = zeros(size(pValsAcc_allsub,1),size(pValsAcc_allsub,2),length(sigLevels));
isSigD = zeros(size(pValsAcc_allsub,1),size(pValsAcc_allsub,2),length(sigLevels));

isSigAcc_mean = zeros(size(pValsAcc_mean_allsub,1),length(sigLevels));
isSigD_mean = zeros(size(pValsAcc_mean_allsub,1),length(sigLevels));

for aa=1:length(sigLevels)

    alpha=sigLevels(aa);
    
    % fdr correct d' within each position sep
    [p_fdr, p_masked] = fdr( pValsAcc_allsub(:,:), alpha);
    isSigAcc(:,:,aa)=p_masked; 
    [p_fdr, p_masked] = fdr( pValsD_allsub, alpha);
    isSigD(:,:,aa)=p_masked;
        
    % fdr correct the mean decoding d'
    [p_fdr, p_masked] = fdr( pValsAcc_mean_allsub, alpha);
    isSigAcc_mean(:,aa)=p_masked; 
    [p_fdr, p_masked] = fdr( pValsD_mean_allsub, alpha);
    isSigD_mean(:,aa)=p_masked; 
    
end

%% make plots of mean pairwise Z decoding (Figure 2B)
     
%NOTE - brackets for pairwise comparisons between these bars were added
%manually in illustrator, calculations for those comparisons are below.
figure;hold all;
ylabel('d-prime');
barMeans=nanmean(squeeze(d_allsub(:,:)),2);

if nSubj>1
    barErrs=nanstd(squeeze(d_allsub(:,:)),[],2)./sqrt(nSubj);
else
    barErrs=nan(size(barMeans));
end

colormap('jet')
bar(barMeans,'EdgeColor','none');
errorbar((1:size(barMeans,1)),barMeans,barErrs,'Marker','none',...
            'LineStyle','none','LineWidth',1,'Color',[0,0,0]);
set(gca, 'XTick', 1:nVOIs, ...
                'XTickLabel', VOIs(vuse),'XTickLabelRotation',90);

ylim([-0.1,0.5]);

horspacer=0.147;
verspacerbig = 0.03;
verspacersmall = 0.01;
markersize = 3;

astLocsEach=nan(nVOIs,2);
for vv=1:nVOIs
    if isSigD_mean(vv,1)
        if barMeans(vv)>0
            astLocsEach(vv)=barMeans(vv)+barErrs(vv)+verspacerbig;
        else
            astLocsEach(vv,1)=barMeans(vv)-barErrs(vv)-verspacerbig;
        end
    end
    if isSigD_mean(vv,2)
        astLocsEach(vv,2) = astLocsEach(vv,1);                    
    end
end

plot((1:nVOIs),astLocsEach(:,1),'o','Color','k','MarkerSize',markersize)
plot((1:nVOIs),astLocsEach(:,2),'o','Color','k','MarkerFaceColor','k','MarkerSize',markersize)

if saveFigs
    fnFig = [root figFolder filesep 'Dprime_allsubs_' typestr '_mean_' voxelStr];        
    fprintf('saving figure to %s...\n',fnFig);
    saveas(gcf,fnFig,ext);
end
           
%% Performing a one-way RM anova on the table of mean d' scores (not significant)

% switch from [nVOIs x nSubj] to [nSubj x nVOIs]
adat = d_mean_allsub';

% Create a table storing the respones
varNames = {'Y1','Y2','Y3','Y4','Y5','Y6','Y7','Y8','Y9'};
% create a table of all levels
t = array2table(data,'VariableNames',varNames);
factorNames = {'ROI'};
within = table(VOIs(vuse)','VariableNames',factorNames);

% fit the repeated measures model
rm = fitrm(t,'Y1-Y9~1','WithinDesign',within);

% this tells us that the assumption of sphericity is violated, and we must
% use the lower bound adjustment (last column of the rm output)
mauchly_tbl = mauchly(rm);

% run my repeated measures anova here
[ranovatbl] = ranova(rm, 'WithinModel','ROI');

% multiple comparisons
mctable = multcompare(rm,'ROI');

% each pairing is reported in both directions - look at only one direction
posDiffs = mctable{:,3}>0 & mctable{:,5}<alpha;
sigDiffs = mctable{posDiffs,[1,2]};

%% get bootstrapped dist of d' for each position and ROI

rs = 8754499;
rng(rs);

nPos = numel(unique(X(:,4)));
nIter = 1000;
bootPerf = nan(length(vuse), nPos, nIter);
parfor ii=1:nIter
    thisorder = datasample(1:nSubj, nSubj, 'Replace',true);
    smallmat = nan(length(vuse), nPos);
    for vv=1:length(vuse)
       for pp=1:nPos 

           dat = X(ismember(X(:,[2,4]), [vv,pp],'rows'), 1);
           % dat is [9 x 1] = this is the individual subj decoding acc for this
           % ROI and position.
           smallmat(vv,pp) = mean(dat(thisorder));

       end
    end
    bootPerf(:,:,ii) = smallmat;

end

%% summary stats on these bootstrapped dist

bootPerf_mean = zeros(length(vuse),1);
bootPerf_ci = zeros(length(vuse), 2);

for vv=1:length(vuse)

    mean_over_position = mean(squeeze(bootPerf(vv,:,:)),1);
    bootPerf_mean(vv) = mean(mean_over_position);
    bootPerf_ci(vv,:) = prctile(mean_over_position,[2.5,97.5]);    
end

%% use the bootstrapped distributions to perform a non-parametric t-test.

pairedVOIs = combnk(1:nv,2);
boot_ttest_perf = nan(size(pairedVOIs,1), 1);
boot_means = nan(size(pairedVOIs,1), 2);
boot_diffs = nan(size(pairedVOIs,1), 1);
for vi = 1:size(pairedVOIs,1)

    s1 = mean(squeeze(bootPerf(pairedVOIs(vi,1), :,:)),1);
    s2 = mean(squeeze(bootPerf(pairedVOIs(vi,2), :,:)),1);
    all_diffs = s1 - s2;
    boot_means(vi,:) = [mean(s1), mean(s2)];
    boot_diffs(vi,:) = mean(all_diffs);
    boot_ttest_perf(vi) = 2*(min([mean(all_diffs < 0), ...
        mean(all_diffs > 0)]));

end

[~,fdr_mask] = fdr(boot_ttest_perf, 0.05);

diffz = find(fdr_mask);
fprintf('average Z decoding is significantly different in %s & %s\n', ...
    VOIs{vuse(pairedVOIs(diffz,:)')});


%% fit slopes/intercepts to the d'/disparity diff relationship

rndseed = [879642];
rng(rndseed,'twister');

nIter = 1000;
allslopes = zeros(nVOIs,nIter);
allint = zeros(nVOIs,nIter);
allslopes_zeroint = zeros(nVOIs,nIter);

for ii=1:nIter    
    randorder = datasample(1:nSubj,nSubj,'Replace',true);
    for vv=1:nVOIs
        
        dprime_bydist = mean(squeeze(d_allsub(vv,randorder,:)),1)';        
        p = polyfit(distListSort,dprime_bydist,1);        
        %p= [slope, intercept]
        allslopes(vv,ii) = p(1);
        allint(vv,ii) = p(2);

        %fit slope with zero interept
        slope = distListSort\dprime_bydist;
        allslopes_zeroint(vv,ii) = slope;
    end        
end

medlines = zeros(nVOIs,2);

for vv=1:nVOIs
    diff=abs(allslopes(vv,:)-prctile(allslopes(vv,:),50));
    thisind = find(diff==min(diff),1);
    medlines(vv,:) = [allslopes(vv,thisind),allint(vv,thisind)]; 
end
%% make a plot showing the pairwise decoding (Figure 3A)

figure;hold all;
drange = [-.5,1];

ii=0;    
for vv=1:nVOIs
    ii=ii+1;
    ax=subplot(5,2,ii);hold all;

    title(sprintf('%s',VOIs{vuse(vv)}))

    if nSubj>1
        barErrsAll=nanstd(squeeze(d_allsub(vv,:,:)),[],1)./sqrt(nSubj);
        barMeansAll=nanmean(squeeze(d_allsub(vv,:,:)),1);
    else
        barMeansAll=squeeze(d_allsub(vv,:,:));
        barErrsAll=nan(size(barMeansAll));
    end

    h1=errorbar(ax,distListSort,barMeansAll(:),barErrsAll(:),'LineWidth',1,'Color',[0,0,0]);

    set(gca,'YLim', [drange], ...
        'XTick',0:10:90,'XTickLabel',{'0.0','10.0','20.0','30.0','40.0','50.0','60.0','70.0','80.0','90.0'},'XTickLabelRotation',90);

    xlabel('Disparity difference (in arcmin)')
    ylabel('d-prime');

    currxlim = get(gca,'XLim');
    line(currxlim,[0,0],'Color','k');

    xdat = (currxlim(1):currxlim(2))';
    ydat = medlines(vv,2)+medlines(vv,1)*xdat;

    line(xdat,ydat);

    astLocsEach=nan(nDist,2);
    for dd=1:nDist
        if isSigAcc(vv,dd,1)
            if barMeansAll(dd)>0
                if nSubj>1
                    astLocsEach(dd,1)=barMeansAll(dd)+barErrsAll(dd)+verspacerbig;
                else
                    astLocsEach(dd,1)=barMeansAll(dd)+verspacerbig;
                end
            else
                if nSubj>1
                    astLocsEach(dd,1)=barMeansAll(dd)-barErrsAll(dd)-verspacerbig;
                else
                    astLocsEach(dd,1)=barMeansAll(dd)-verspacerbig;
                end
            end
        end
        if isSigAcc(vv,dd,2)
            astLocsEach(dd,2) = astLocsEach(dd,1);                    
        end
    end

    plot(distListSort,astLocsEach(:,1),'o','Color','k','MarkerSize',markersize)
    plot(distListSort,astLocsEach(:,2),'o','Color','k','MarkerFaceColor','k','MarkerSize',markersize)

    uistack(h1,'top')

end

if saveFigs
    set(gcf,'Position',get(0,'ScreenSize'))
    set(gcf,'Color','w')
    fnFig = [root figFolder filesep 'Dprime_2Z_oneVsOne_byDisp_wLine_' voxelStr];
    fprintf('saving figure to %s...\n',fnFig);
    saveas(gcf,fnFig,ext);
end

%% Make a plot of slope distribution (Figure 3B)

%NOTE - some aspects of this plot are manually added in illustrator. These
%include making all circle markers filled except for V1, to indicate the V1
%distribution overlaps w/ zero. Also adding brackets for pairwise
%comparisons of slopes b/w ROIs. 
% the calculations corresponding to those elements are below.
figure;hold all;

title('Slope')

ylabel('Slope');

vp = violinplot(allslopes', VOIs,'ShowData',false, 'ShowMean',false, 'ViolinColor',[0.8,0.8,0.8]);

set(gca,'XTickLabel',VOIs(vuse))
line(get(gca,'XLim'), [0,0]);

%% Compare all slopes to zero, and FDR correct

fdrThresh = zeros(2,2);
fdrLabs = {'slope: q=0.05','int: q=0.05';'slope: q=0.01','int q=0.01'};

pValsSlope = min([mean(allslopes<0,2),mean(allslopes>0,2)],[],2);     

isSigSlope = zeros(nVOIs,2);

for aa=1:length(sigLevels)
    alpha=sigLevels(aa);

    [p_fdr, p_masked] = fdr(pValsSlope, alpha);
    isSigSlope(:,aa)=p_masked; 

    fdrThresh(aa,1) = p_fdr;
end

%% Perform all pairwise slope comparisons b/w ROIs
pList=[];
vCompList=[];
ii=0;
meanList=[];

for vv1=1:nVOIs
    for vv2 = vv1+1:nVOIs
            
        set1 = allslopes(vv1,:)';
        set2 = allslopes(vv2,:)';

        p = min([mean(set1<set2),mean(set2<set1)]);
        
        pList = cat(1,pList,p);
        meanList=cat(1,meanList,[mean(allslopes(vv1,:)),mean(allslopes(vv2,:))]);
        vCompList = cat(1,vCompList,[VOIs(vuse(vv1)),VOIs(vuse(vv2))]);              
    end
end

[fdrThreshBetween,pListCorr] = fdr(pList,0.01);
mctable = table(vCompList,meanList,pList,pListCorr);

%% Make RSA plot (extended data, Figure 3-1)
 
groupStrs = {'38.6','25.6','11.1','-5.2','-23.6','-44.6'};
confMat = zeros(nSubj,length(vuse),6,6);
nPoss = zeros(6,6);

% put the significance values into matrix format
sigMat = zeros(nVOIs, 6, 6);
for vv=1:nVOIs    
    for pp=1:nDist
        gg1 = posPairList(pp,1);
        gg2 = posPairList(pp,2);
        % use distOrder vector here to make sure we're grabbing the correct
        % point...ordering here is not straightforward.
        if isSigD(vv,distOrder(pp),1)
            sigMat(vv,gg1,gg2) = 1;
        end
        if isSigD(vv,distOrder(pp),2)
            sigMat(vv,gg1,gg2) = 2;
        end
    end
end

% make the actual plot
ax = [];
nGroups = 6;

figure('units','normalized','outerposition',[0,0,1,1]);
hold all;

for vv=1:nVOIs
    
    thisMat = squeeze(mean(bigMat(:,vv,:,:),1));

    for ii=1:6
        for jj=ii:6
            thisMat(ii,jj) =nan;
        end       
    end
    
    subplot(3,3,vv);hold all
   
    sanePColor([thisMat]);
    axis equal;

    xlim([.5, size(thisMat,1)+.5])
    ylim([.5, size(thisMat,1)+.5])


    set(gca,'XTick', 1:nGroups, ...
                'XTickLabel', groupStrs,'XTickLabelRotation',90);
    set(gca,'YTick', 1:nGroups, ...
        'YTickLabel', groupStrs);
    xlabel('Disparity (arcmin)');
    ylabel('Disparity (arcmin)');
    
    set(gca,'YDir','reverse')
    
    % add significance stars
    [x,y] = ind2sub([6,6], find(squeeze(sigMat(vv,:,:))==1));    
    plot(x, y, 'ko')
    
    [x,y] = ind2sub([6,6], find(squeeze(sigMat(vv,:,:))==2));   
    plot(x, y, 'k.', 'MarkerSize',10)

    title(sprintf('%s: Z discriminability (d'')',VOIs{vuse(vv)}));
    colormap(viridis)
    colorbar();
    
    ax = [ax, gca];
end

match_clim(ax);

fnFig = [root figFolder filesep 'RSA_6Z_oneVsOne'];      
fprintf('saving figure to %s...\n',fnFig);
print(fnFig,'-depsc','-r0')
