function IEMdepth_anaClassifier_Z6_mean()

%plot the results of decoding within each VOI: compute mean accuracy and
%significance at subject level, FDR correct the p-vals

%% define subjects and flags for what to do


subj = {'AI','AP','BB','BC','BD','BJ','BO','BN','BM'};

nSubj=length(subj);
    
VOIs = {'V1','V2','V3','V4','V3A','V3B','IPS0','IPS1','IPS2','IPS3','LO1','LO2'};

vuse = [1:7,11:12];
nVOIs=length(vuse);

plotSig=1;

plotVoxNum=0;

plotAcc=1;
plotD=0;

% numVoxUse=150;
% voxelStr=sprintf('take%dZVox',numVoxUse);
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


figFolder='IEMdepth_figs';

sigLevels=[0.05,0.01];


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

fdrThresh = zeros(2,2);
fdrLabs = {'Acc - 0.05','Dprime - 0.05'; 'Acc - 0.01','Dprime - 0.01'};

for aa=1:length(sigLevels)
    alpha=sigLevels(aa);

    [p_fdr, p_masked] = fdr( pValsAcc_allsub, alpha);
    isSigAcc(:,aa)=p_masked; 

    fdrThresh(aa,1) = p_fdr;
    
    [p_fdr, p_masked] = fdr( pValsD_allsub, alpha);
    isSigD(:,aa)=p_masked; 
    
    fdrThresh(aa,2) = p_fdr;
end

%% save output

rowLabs = VOIs;
colLabs = {'Mean accuracy','SE',...
    'Mean dprime','SE','p'}; 

outTable = zeros(nVOIs, 5);

outTable(:,1) = nanmean(squeeze(accs_allsub(:,:)),2);
outTable(:,2) = nanstd(squeeze(accs_allsub(:,:)),[],2)./sqrt(nSubj);

outTable(:,3) = nanmean(squeeze(d_allsub(:,:)),2);
outTable(:,4) = nanstd(squeeze(d_allsub(:,:)),[],2)./sqrt(nSubj);
outTable(:,5) = pValsD_allsub(:);

fnTab = [root figFolder filesep 'DecTable_' typestr '_mean_' voxelStr '.mat'];
save(fnTab, 'outTable','colLabs','rowLabs','fdrThresh','fdrLabs');
