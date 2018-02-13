function IEMdepth_anaClassifier_Z6_oneVsOne()

% make table of classifier performance for 6Z one vs one

% MMH 9/14/17

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

saveFig = 1;
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

zDispList = [-35.8,-22.7,-8,8.5,27.2,48.4];

posPairList = combnk(1:6,2);
dispPairList = reshape(zDispList(posPairList(:)),15,2);
dispDiffList = dispPairList(:,2)-dispPairList(:,1);
[distListSort,distOrder] = sort(dispDiffList,'ascend');

undist = unique(dispDiffList);
nDist = length(undist);
    
accs_allsub = zeros(nVOIs,nSubj,nDist);
accsrand_allsub = zeros(nVOIs,nSubj,nDist,nIter);

d_allsub = zeros(nVOIs,nSubj,nDist);
drand_allsub = zeros(nVOIs,nSubj,nDist,nIter);


%% loop over subs
for ss=1:nSubj   
    
    fns=sprintf('%s%s/%s_%s_%s_%s_%s_%s_%s_%s.mat',...
                    root,folder,subj{ss},typestr,condStrs{conduse},voxelStr,predStr,classStr,kernelStr,subMeanStr);

    load(fns);
            
    for vv=1:nVOIs
        
        for dd=1:nDist
            
            theseaccs = classStruct(vuse(vv)).accReal(distOrder(dd),:);
            accs_allsub(vv,ss,dd) = mean(theseaccs);
            
            theseaccs = classStruct(vuse(vv)).accRand(distOrder(dd),:);
            accsrand_allsub(vv,ss,dd,:) = mean(theseaccs,1);
            
            thesed = classStruct(vuse(vv)).dReal(distOrder(dd),:);
            d_allsub(vv,ss,dd) = mean(thesed);
            
            thesed = classStruct(vuse(vv)).dRand(distOrder(dd),:);
            drand_allsub(vv,ss,dd,:) = mean(thesed,1);
        end
        
    end
   
end

%% compute p vals for each roi, each cond separately 
tic
% if plotSig
for vv=1:nVOIs
    for dd=1:nDist

        %% all trials
        realAccs = squeeze(accs_allsub(vv,:,dd));
        nullAccs = squeeze(accsrand_allsub(vv,:,dd,:));

        realD = squeeze(d_allsub(vv,:,dd));
        nullD = squeeze(drand_allsub(vv,:,dd,:));

        pValsAcc_allsub(vv,dd) = 2*min([mean(mean(realAccs)<mean(nullAccs)),mean(mean(realAccs)>mean(nullAccs))]);
        pValsD_allsub(vv,dd) = 2*min([mean(mean(realD)<mean(nullD)),mean(mean(realD)>mean(nullD))]);


    end
end
toc

%% FDR correction
isSigAcc = zeros(size(pValsAcc_allsub,1),size(pValsAcc_allsub,2),length(sigLevels));
isSigD = zeros(size(pValsAcc_allsub,1),size(pValsAcc_allsub,2),length(sigLevels));

fdrThresh = zeros(2,2);
fdrLabs = {'Acc - 0.05','Dprime - 0.05'; 'Acc - 0.01','Dprime - 0.01'};


for aa=1:length(sigLevels)

    alpha=sigLevels(aa);

    [p_fdr, p_masked] = fdr( pValsAcc_allsub(:,:), alpha);
    isSigAcc(:,:,aa)=p_masked; 

    fdrThresh(aa,1) = p_fdr;
    
    [p_fdr, p_masked] = fdr( pValsD_allsub, alpha);
    isSigD(:,:,aa)=p_masked;
    
    fdrThresh(aa,2) = p_fdr;
    
end

% end

meanAccs = squeeze(nanmean(accs_allsub,2));


%% save the output

% in same folder where the figs are

rowLabs = VOIs;

outTableAcc = zeros(nVOIs, 3*nDist);
outTableD = zeros(nVOIs, 3*nDist);

outTable = zeros(nVOIs,5*3,5);

ii=0;
jj=0;
pp=1;
colLabs = cell(5*3,5);
for dd=1:nDist
    
    ii=ii+1;
    colLabsAcc{ii} = sprintf('Mean accuracy: disp %.4f',distListSort(dd));
    outTableAcc(:,ii) = nanmean(squeeze(accs_allsub(:,:,dd)),2);
    colLabsD{ii} = sprintf('Mean dprime: disp %.4f',distListSort(dd));
    outTableD(:,ii) = nanmean(squeeze(d_allsub(:,:,dd)),2);
    ii=ii+1;
    colLabsAcc{ii} = sprintf('SE: disp %.4f',distListSort(dd));
    outTableAcc(:,ii) = nanstd(squeeze(accs_allsub(:,:,dd)),[],2)./sqrt(nSubj);
    colLabsD{ii} = sprintf('SE: disp %.4f',distListSort(dd));
    outTableD(:,ii) = nanstd(squeeze(d_allsub(:,:,dd)),[],2)./sqrt(nSubj);
    ii=ii+1;
    colLabsAcc{ii} = sprintf('p: disp %.4f',distListSort(dd));
    outTableAcc(:,ii) = pValsAcc_allsub(:,dd);
    colLabsD{ii} = sprintf('p: disp %.4f',distListSort(dd));
    outTableD(:,ii) = pValsD_allsub(:,dd);
   
    
    jj=jj+1;
    colLabs{jj,pp} = sprintf('Mean accuracy: disp %.4f',distListSort(dd));
    outTable(:,jj,pp) = nanmean(squeeze(accs_allsub(:,:,dd)),2);
    jj=jj+1;
    colLabs{jj,pp} = sprintf('SE: disp %.4f',distListSort(dd));
    outTable(:,jj,pp) = nanstd(squeeze(accs_allsub(:,:,dd)),[],2)./sqrt(nSubj);
    jj=jj+1;
    colLabs{jj,pp} = sprintf('Mean dprime: disp %.4f',distListSort(dd));
    outTable(:,jj,pp) = nanmean(squeeze(d_allsub(:,:,dd)),2);    
    jj=jj+1;
    colLabs{jj,pp} = sprintf('SE: disp %.4f',distListSort(dd));
    outTable(:,jj,pp) = nanstd(squeeze(d_allsub(:,:,dd)),[],2)./sqrt(nSubj);
    jj=jj+1;
    colLabs{jj,pp} = sprintf('p: disp %.4f',distListSort(dd));
    outTable(:,jj,pp) = pValsD_allsub(:,dd);
    
    if ~mod(dd,3);
        jj=0;
        pp=pp+1;
    end
        
end

fnTab = [root figFolder filesep 'DecTable_' typestr '_' voxelStr '.mat'];
save(fnTab, 'outTable','colLabs','outTableAcc','outTableD','colLabsAcc','colLabsD','rowLabs','fdrThresh','fdrLabs');



