function IEMdepth_plotClassifier_XZ2()

% plot the result of 2-way classification in X and in Z - this is a coarse
% division of all trials into two groups (left-right or front-back)

% using all voxels, not using IPS2-3

%% define subjects and flags for what to do


subj = {'AI','AP','BB','BC','BD','BJ','BO','BN','BM'};
% subj = {'AP','BB','BC','BD','BJ','BM','BN','BO'};
nSubj=length(subj);

VOIs = {'V1','V2','V3','V4','V3A','V3B','IPS0','IPS1','IPS2','IPS3','LO1','LO2'};

vuse = [1:7,11:12];
nVOIs=length(vuse);

plotSig=1;
plotAcc=0;
plotD=1;
% 
% numVoxUse=150;
% voxelStr=sprintf('take%dZ1WayVox',numVoxUse);
voxelStr = 'allVox';
classStr = 'svmtrain';
kernelStr='linear';
predStr = 'predA';
subMeanStr = 'noSubMean';

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


chanceVal=1/2;
drange = [-.5,4];
ymax=1;

typestr = 'XZ2_singleTrialPreds';
titlestr='2-way classification';
    
% arrays to store acc and d' 
accs_allsub=nan(nVOIs,nSubj,nDim);
accsrand_allsub=nan(nVOIs,nSubj,nDim,nIter);
d_allsub=nan(nVOIs,nSubj,nDim);
drand_allsub=nan(nVOIs,nSubj,nDim,nIter);

% p vals for significance of decoding in each condition alone
pValsAcc_allsub=nan(nVOIs,nDim);
pValsD_allsub = nan(nVOIs,nDim);

nVox = zeros(nVOIs,nSubj);

%% loop over subs
for ss=1:nSubj   
        
    fns=sprintf('%s%s/%s_allROIs_%s_%s_%s_%s_%s_%s_%s.mat',...
                    root,folder,subj{ss},typestr,condStrs{conduse},voxelStr,predStr,classStr,kernelStr,subMeanStr);
    load(fns);

    for vv=1:length(vuse)
        
        nVox(vv,ss) = classStruct(vuse(vv)).numVoxTot;
        
        accs_allsub(vv,ss,1) = classStruct(vuse(vv)).accRealX;
        accsrand_allsub(vv,ss,1,:) = classStruct(vuse(vv)).accRandX;
        
        accs_allsub(vv,ss,2) = classStruct(vuse(vv)).accRealZ;
        accsrand_allsub(vv,ss,2,:) = classStruct(vuse(vv)).accRandZ;
        
        d_allsub(vv,ss,1) = classStruct(vuse(vv)).dRealX;
        drand_allsub(vv,ss,1,:) = classStruct(vuse(vv)).dRandX;
        
        d_allsub(vv,ss,2) = classStruct(vuse(vv)).dRealZ;
        drand_allsub(vv,ss,2,:) = classStruct(vuse(vv)).dRandZ;
    end

end

%% compute p vals for each roi, each cond separately 

if plotSig
for vv=1:nVOIs
    for cc=1:nDim

        %% all trials
        realAccs = squeeze(accs_allsub(vv,:,cc));
        nullAccs = squeeze(accsrand_allsub(vv,:,cc,:));
 
        realD = squeeze(d_allsub(vv,:,cc));
        nullD = squeeze(drand_allsub(vv,:,cc,:));

        pValsAcc_allsub(vv,cc) = 2*min([mean(mean(realAccs)<mean(nullAccs)),mean(mean(realAccs)>mean(nullAccs))]);
        pValsD_allsub(vv,cc) = 2*min([mean(mean(realD)<mean(nullD)),mean(mean(realD)>mean(nullD))]);


    end
end

%% FDR correction

isSigAcc = zeros(size(pValsAcc_allsub,1),size(pValsAcc_allsub,2),length(sigLevels));
isSigD = zeros(size(pValsAcc_allsub,1),size(pValsAcc_allsub,2),length(sigLevels));

for aa=1:length(sigLevels)
    alpha=sigLevels(aa);

    [p_fdr, p_masked] = fdr( pValsAcc_allsub, alpha);
    isSigAcc(:,:,aa)=p_masked; 

    [p_fdr, p_masked] = fdr( pValsD_allsub, alpha);
    isSigD(:,:,aa)=p_masked; 
end

end
%% make plot for all subs

if plotAcc
       
    figure;hold all;
    
    title(sprintf('Accuracy, %s',titlestr))
    
    barMeansAll = zeros(nVOIs,nDim);
    barErrsAll = zeros(nVOIs,nDim);
    
    for cc=1:nDim
    
        barMeans=nanmean(squeeze(accs_allsub(:,:,cc)),2);

        if nSubj>1
            barErrs=nanstd(squeeze(accs_allsub(:,:,cc)),[],2)./sqrt(nSubj);
        else
            barErrs=nan(size(barMeans));
        end

        barMeansAll(:,cc) = barMeans;
        barErrsAll(:,cc) = barErrs;
    end
        
    
        
    colormap('jet')
    bar(barMeansAll);
    errorbar((1:size(barMeansAll,1))-horspacer,barMeansAll(:,1),barErrsAll(:,1),'Marker','none',...
                'LineStyle','none','LineWidth',1,'Color',[0,0,0]);
    errorbar((1:size(barMeansAll,1))+horspacer,barMeansAll(:,2),barErrsAll(:,2),'Marker','none',...
        'LineStyle','none','LineWidth',1,'Color',[0,0,0]);
    set(gca,'YLim', [0,ymax], 'XTick', 1:nVOIs, ...
                    'XTickLabel', VOIs(vuse),'XTickLabelRotation',90);
        line([0,nVOIs+1],[chanceVal,chanceVal],'Color','k');

    legend(dimStrs,'Location','EastOutside');

    if plotSig
        astLocsEach=nan(nVOIs,nDim,2);
        for vv=1:nVOIs;
            for cc=1:nDim
                if isSigAcc(vv,cc,1)
                    if barMeansAll(vv,cc)>0
                        astLocsEach(vv,cc,1)=barMeansAll(vv,cc)+barErrsAll(vv,cc)+verspacerbig;
                    else
                        astLocsEach(vv,cc,1)=barMeansAll(vv,cc)-barErrsAll(vv,cc)-verspacerbig;
                    end
                end
                if isSigAcc(vv,cc,2)
                    astLocsEach(vv,cc,2) = astLocsEach(vv,cc,1);                    
                end
            end
        end
        
        
        plot((1:nVOIs)-horspacer,astLocsEach(:,1,1),'o','Color','k','MarkerSize',markersize)
        plot((1:nVOIs)-horspacer,astLocsEach(:,1,2),'o','Color','k','MarkerFaceColor','k','MarkerSize',markersize)
        plot((1:nVOIs)+horspacer,astLocsEach(:,2,1),'o','Color','k','MarkerSize',markersize)
        plot((1:nVOIs)+horspacer,astLocsEach(:,2,2),'o','Color','k','MarkerFaceColor','k','MarkerSize',markersize)
        
    end
     
    if saveFig
        fnFig = [root figFolder filesep 'Acc_allsubs_' typestr '_' voxelStr];
        
        fprintf('saving figure to %s...\n',fnFig);
        saveas(gcf,fnFig,ext);
    end
              
end
%%

if plotD
       
    figure;hold all;
    
    title([titlestr ' ' voxelStr])
    ylabel('d-prime');
    
    barMeansAll = zeros(nVOIs,nDim);
    barErrsAll = zeros(nVOIs,nDim);
    
    for cc=1:nDim
    
        barMeans=nanmean(squeeze(d_allsub(:,:,cc)),2);

        if nSubj>1
            barErrs=nanstd(squeeze(d_allsub(:,:,cc)),[],2)./sqrt(nSubj);
        else
            barErrs=nan(size(barMeans));
        end

        barMeansAll(:,cc) = barMeans;
        barErrsAll(:,cc) = barErrs;
    end
                
    colormap('jet')
    bar(barMeansAll,'EdgeColor','none');
    errorbar((1:size(barMeansAll,1))-horspacer,barMeansAll(:,1),barErrsAll(:,1),'Marker','none',...
                'LineStyle','none','LineWidth',1,'Color',[0,0,0]);
    errorbar((1:size(barMeansAll,1))+horspacer,barMeansAll(:,2),barErrsAll(:,2),'Marker','none',...
        'LineStyle','none','LineWidth',1,'Color',[0,0,0]);
    set(gca,'YLim', drange, 'XTick', 1:nVOIs, ...
                    'XTickLabel', VOIs(vuse),'XTickLabelRotation',90);
        line([0,nVOIs+1],[0,0],'Color','k');

    legend(dimStrs,'Location','EastOutside');

    if plotSig

        astLocsEach=nan(nVOIs,nDim,2);
        for vv=1:nVOIs;
            for cc=1:nDim
                if isSigD(vv,cc,1)
                    if barMeansAll(vv,cc)>0
                        astLocsEach(vv,cc,1)=barMeansAll(vv,cc)+barErrsAll(vv,cc)+verspacerbig;
                    else
                        astLocsEach(vv,cc,1)=barMeansAll(vv,cc)-barErrsAll(vv,cc)-verspacerbig;
                    end
                end
                if isSigD(vv,cc,2)
                    astLocsEach(vv,cc,2) = astLocsEach(vv,cc,1);                    
                end
            end
        end
        
        
        plot((1:nVOIs)-horspacer,astLocsEach(:,1,1),'o','Color','k','MarkerSize',markersize)
        plot((1:nVOIs)-horspacer,astLocsEach(:,1,2),'o','Color','k','MarkerFaceColor','k','MarkerSize',markersize)
        plot((1:nVOIs)+horspacer,astLocsEach(:,2,1),'o','Color','k','MarkerSize',markersize)
        plot((1:nVOIs)+horspacer,astLocsEach(:,2,2),'o','Color','k','MarkerFaceColor','k','MarkerSize',markersize)
        
    end
     
     
    if saveFig
        fnFig = [root figFolder filesep 'Dprime_allsubs_' typestr '_' voxelStr];
        
        fprintf('saving figure to %s...\n',fnFig);
        saveas(gcf,fnFig,ext);
        
        
    end
%          
              
end

    %% tilefigs

if length(get(0,'Children')) > 1
    tilefigs();
end

end
    

