function IEMdepth_plotClassifier_Z2_mean()

%plot the results of decoding within each VOI: compute mean accuracy and
%significance at subject level, FDR correct the p-vals

%% define subjects and flags for what to do


subj = {'AI','AP','BB','BC','BD','BJ','BO','BN','BM'};

nSubj=length(subj);
    
VOIs = {'V1','V2','V3','V4','V3A','V3B','IPS0','IPS1','IPS2','IPS3','LO1','LO2'};

vuse = [1:8,11:12];
nVOIs=length(vuse);

plotSig=1;

plotVoxNum=0;

plotAcc=0;
plotD=1;
% 
% numVoxUse=50;
% voxelStr=sprintf('take%dZ1WayVox',numVoxUse);
voxelStr = 'allVox';
classStr = 'svmtrain';
kernelStr='linear';
predStr = 'predA';
subMeanStr = 'noSubMean';
% subMeanStr = 'subMean2';
tirange=[1];

typestr='Z2_oneVsOne';
titlestr='6-way - oneVsOne';

ymaxes=[1,1];
drange=[-.2,0.6];
acclims = [0,0.6];
chanceVals=[1/2];

dimStrs = {'X position','Z position'};
nDim= length(dimStrs);

condStrs = {'trainStim','trainFixat'};
conduse = 2;

nIter=1000;

% close all

%% set up file info, other params

root='/usr/local/serenceslab/maggie/IEMdepth/';
folder='IEMdepth_classif';

saveFig = 1;
figFolder='IEMdepth_figs';
ext='epsc';

sigLevels=[0.05,0.01];

horspacer=0.147;
verspacerbig = 0.03;
verspacersmall = 0.01;
markersize = 3;

for ti=tirange

    chanceVal=chanceVals(ti);
    ymax=ymaxes(ti);

    % arrays to store acc and d' 
    accs_allsub=nan(nVOIs,nSubj);
    accsrand_allsub=nan(nVOIs,nSubj,nIter);
    d_allsub=nan(nVOIs,nSubj);
    drand_allsub=nan(nVOIs,nSubj,nIter);
     
    % p vals for significance of decoding in each condition alone
    pValsAcc_allsub=nan(nVOIs,1);
    pValsD_allsub = nan(nVOIs,1);


%% loop over subs
for ss=1:nSubj   

    fns=sprintf('%s%s/%s_allROIs_%s_%s_%s_%s_%s_%s_%s.mat',...
                    root,folder,subj{ss},typestr,condStrs{conduse},voxelStr,predStr,classStr,kernelStr,subMeanStr);

    load(fns);

    for vv=1:length(vuse)
        
        accs_allsub(vv,ss) = mean(classStruct(vuse(vv)).accReal,1);
        accsrand_allsub(vv,ss,:) = mean(classStruct(vuse(vv)).accRand,1);
        
        d_allsub(vv,ss) =  mean(classStruct(vuse(vv)).dReal,1);
        drand_allsub(vv,ss,:) = mean(classStruct(vuse(vv)).dRand,1);
        
    end

end

%% compute p vals for each roi, each cond separately 

if plotSig
for vv=1:nVOIs

    %% all trials
    realAccs = squeeze(accs_allsub(vv,:));
    nullAccs = squeeze(accsrand_allsub(vv,:,:));

    realD = squeeze(d_allsub(vv,:));
    nullD = squeeze(drand_allsub(vv,:,:));

    pValsAcc_allsub(vv) = 2*min([mean(mean(realAccs)<mean(nullAccs)),mean(mean(realAccs)>mean(nullAccs))]);
    pValsD_allsub(vv) = 2*min([mean(mean(realD)<mean(nullD)),mean(mean(realD)>mean(nullD))]);

end

%% FDR correction
isSigAcc = zeros(size(pValsAcc_allsub,1),length(sigLevels));
isSigD = zeros(size(pValsAcc_allsub,1),length(sigLevels));

for aa=1:length(sigLevels)
    alpha=sigLevels(aa);

    [p_fdr, p_masked] = fdr( pValsAcc_allsub, alpha);
    isSigAcc(:,aa)=p_masked; 

    [p_fdr, p_masked] = fdr( pValsD_allsub, alpha);
    isSigD(:,aa)=p_masked; 
end

end
%% make plot for all subs

if plotAcc
       
    figure;hold all;
    
    title(sprintf('Accuracy, %s',titlestr))
    
    barMeans=nanmean(squeeze(accs_allsub(:,:)),2);

    if nSubj>1
        barErrs=nanstd(squeeze(accs_allsub(:,:)),[],2)./sqrt(nSubj);
    else
        barErrs=nan(size(barMeans));
    end  
        
    colormap('jet')
    bar(barMeans);
    errorbar((1:size(barMeans,1)),barMeans,barErrs,'Marker','none',...
                'LineStyle','none','LineWidth',1,'Color',[0,0,0]);
    set(gca,'YLim', [acclims], 'XTick', 1:nVOIs, ...
                    'XTickLabel', VOIs(vuse),'XTickLabelRotation',90);
        line([0,nVOIs+1],[chanceVal,chanceVal],'Color','k');

%         legend(dimStrs,'Location','EastOutside');

    if plotSig
        
        astLocsEach=nan(nVOIs,2);
        for vv=1:nVOIs;

            if isSigAcc(vv,1)
                if barMeans(vv)>0
                    astLocsEach(vv,1)=barMeans(vv)+barErrs(vv)+verspacerbig;
                else
                    astLocsEach(vv,1)=barMeans(vv)-barErrs(vv)-verspacerbig;
                end
            end
            if isSigAcc(vv,2)
                astLocsEach(vv,2) = astLocsEach(vv,1);                    
            end

        end
        
        
        plot((1:nVOIs),astLocsEach(:,1),'o','Color','k','MarkerSize',markersize)
        plot((1:nVOIs),astLocsEach(:,2),'o','Color','k','MarkerFaceColor','k','MarkerSize',markersize)

     end
              
end
%%

if plotD
       
    figure;hold all;
    
    title(titlestr)
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
% 
%     for vv=1:nVOIs
%         scatter(ones(nSubj,1)*vv,d_allsub(vv,:,2));
%     end

%     legend(dimStrs,'Location','EastOutside');

    ylim([-0.1,0.5]);

    if plotSig

        astLocsEach=nan(nVOIs,2);
        for vv=1:nVOIs;

            if isSigD(vv,1)
                if barMeans(vv)>0
                    astLocsEach(vv)=barMeans(vv)+barErrs(vv)+verspacerbig;
                else
                    astLocsEach(vv,1)=barMeans(vv)-barErrs(vv)-verspacerbig;
                end
            end
            if isSigD(vv,2)
                astLocsEach(vv,2) = astLocsEach(vv,1);                    
            end

        end
        
        
        plot((1:nVOIs),astLocsEach(:,1),'o','Color','k','MarkerSize',markersize)
        plot((1:nVOIs),astLocsEach(:,2),'o','Color','k','MarkerFaceColor','k','MarkerSize',markersize)

    end
     
     
     if saveFig
        fnFig = [root figFolder filesep 'Dprime_allsubs_' typestr '_mean_' voxelStr];
        
        fprintf('saving figure to %s...\n',fnFig);
        saveas(gcf,fnFig,ext);
    end
           
end

    %% plot nVox
    
if plotVoxNum
    
    figure();hold all
    title('Number of voxels per VOI per sub')
    plot(1:nVOIs,nVox_allsub,'*');
    set(gca, 'XTick', 1:nVOIs, ...
                'XTickLabel', VOIs(1:nVOIs),'XTickLabelRotation',90);
    legend(subj,'Location','EastOutside');
end
  

end

    %% tilefigs

if length(get(0,'Children')) > 1
    tilefigs();
end


end
    

