% run through all analyses, make all figures and tables for two-way
% decoding of X and Z position
% load the classifier output saved by the script IEMdepth_classify_XZ2.m

% MMH 3/29/19

clear

%% set up some path stuff
% This is whatever directory contains the folder "IEMdepth_classif"
root = '/usr/local/serenceslab/maggie/IEMdepth/';
% This is the path where the classifier output files are located.
load_folder = [root 'IEMdepth_classif'];

%%
% set up a list of areas, subjects, etc.
subj = {'AI','AP','BB','BC','BD','BJ','BM','BN','BO'};
vuse=[1:7,11:12];
VOIs={'V1','V2','V3','V4','V3A','V3B','IPS0','IPS1','IPS2','IPS3','LO1','LO2'};
nSubj=length(subj);
nVOIs=length(vuse);

%parameters for the classifier, tell the script the data file to load.
kernelStr='linear';
typestr = 'XZ2_singleTrialPreds';
classStr='svmtrain';
subMeanStr = 'noSubMean';
predStr='predA';
voxelStr = 'allVox';
condStrs = {'trainStim','trainFixat'};
conduse = 2;

dimStrs = {'X position','Z position'};
nDim= length(dimStrs);

nIter=1000;

sigLevels=[0.05,0.01];

% arrays to store acc and d' 
accs_allsub=nan(nVOIs,nSubj,nDim);
accsrand_allsub=nan(nVOIs,nSubj,nDim,nIter);
d_allsub=nan(nVOIs,nSubj,nDim);
drand_allsub=nan(nVOIs,nSubj,nDim,nIter);

% p vals for significance of decoding in each condition alone
pValsAcc_allsub=nan(nVOIs,nDim);
pValsD_allsub = nan(nVOIs,nDim);

%% loop over subs
for ss=1:nSubj   
        
    fns=sprintf('%s/%s_allROIs_%s_%s_%s_%s_%s_%s_%s.mat',...
                    load_folder,subj{ss},typestr,condStrs{conduse},voxelStr,predStr,classStr,kernelStr,subMeanStr);
    load(fns);

    for vv=1:length(vuse)
        
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

for vv=1:nVOIs
    for cc=1:nDim

        %all trials
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

fdrThresh = zeros(2,2);
fdrLabs = {'Acc - 0.05','Dprime - 0.05'; 'Acc - 0.01','Dprime - 0.01'};

for aa=1:length(sigLevels)
    alpha=sigLevels(aa);

    [p_fdr, p_masked] = fdr( pValsAcc_allsub, alpha);
    isSigAcc(:,:,aa)=p_masked; 
    fdrThresh(aa,1) = p_fdr;

    [p_fdr, p_masked] = fdr( pValsD_allsub, alpha);
    isSigD(:,:,aa)=p_masked;
    fdrThresh(aa,2) = p_fdr;
end

%% Make the plot (Figure 2A)

figure;hold all;

horspacer=0.147;
verspacerbig = 0.03;
verspacersmall = 0.01;
markersize = 3;

drange = [-.5,4];


% title([titlestr ' ' voxelStr])
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
       