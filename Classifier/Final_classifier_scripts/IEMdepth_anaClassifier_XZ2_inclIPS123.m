function IEMdepth_anaClassifier_XZ2_inclIPS123()

% make a table of the result of two-way classification of x and z position

% ignore IPS2-3

% MMH 9/20/17

%% define subjects and flags for what to do


subj = {'AI','AP','BB','BC','BD','BJ','BO','BN','BM'};

nSubj=length(subj);

VOIs = {'V1','V2','V3','V4','V3A','V3B','IPS0','IPS1','IPS2','IPS3','LO1','LO2'};

vuse = [1:12];
nVOIs=length(vuse);

figFolder='IEMdepth_figs';

% numVoxUse=150;
% voxelStr=sprintf('take%dZVox',numVoxUse);
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

%% set up file info, other params

root='/usr/local/serenceslab/maggie/IEMdepth/';
folder='IEMdepth_classif';

sigLevels=[0.05,0.01];

typestr = 'XZ2_singleTrialPreds';

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
        
    fns=sprintf('%s%s/%s_allROIs_%s_%s_%s_%s_%s_%s_%s.mat',...
                    root,folder,subj{ss},typestr,condStrs{conduse},voxelStr,predStr,classStr,kernelStr,subMeanStr);
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

%% save the table
% in same folder where the figs are

rowLabs = VOIs;
colLabs = {'Mean accuracy - X','SE','p',...
    'Mean accuracy - Z','SE','p',...
    'Mean dprime - X','SE','p',...
    'Mean dprime - Z','SE','p'}; 

outTable = zeros(nVOIs, length(colLabs));

outTable(:,1) = nanmean(squeeze(accs_allsub(:,:,1)),2);
outTable(:,2) = nanstd(squeeze(accs_allsub(:,:,1)),[],2)./sqrt(nSubj);
outTable(:,3) = pValsAcc_allsub(:,1);

outTable(:,4) = nanmean(squeeze(accs_allsub(:,:,2)),2);
outTable(:,5) = nanstd(squeeze(accs_allsub(:,:,2)),[],2)./sqrt(nSubj);
outTable(:,6) = pValsAcc_allsub(:,2);

outTable(:,7) = nanmean(squeeze(d_allsub(:,:,1)),2);
outTable(:,8) = nanstd(squeeze(d_allsub(:,:,1)),[],2)./sqrt(nSubj);
outTable(:,9) = pValsD_allsub(:,1);

outTable(:,10) = nanmean(squeeze(d_allsub(:,:,2)),2);
outTable(:,11) = nanstd(squeeze(d_allsub(:,:,2)),[],2)./sqrt(nSubj);
outTable(:,12) = pValsD_allsub(:,2);


fnTab = [root figFolder filesep 'DecTable_' typestr '_' voxelStr '_inclIPS123.mat'];
save(fnTab, 'outTable','rowLabs','colLabs','fdrThresh','fdrLabs');


