function IEMdepth_anaSlopeDist_Z6()

% load resampled confidence intervals for slope and intercept of the
% d'/disparity relationship in each ROI, save in a table

% MMH 10/3/17
%% define subjects and flags for what to do

subj = {'AI','AP','BB','BC','BD','BJ','BO','BN','BM'};
VOIs = {'V1','V2','V3','V4','V3A','V3B','IPS0','IPS1','IPS2','IPS3','LO1','LO2'};
nSubj = length(subj);
% VOIs={'V1','V2','V3','V4','V3A','V3B','IPS0','IPS1','IPS2','IPS3','IPS1-3','allDors'};

% vuse = 1:nVOIs;

vuse = [1:7,11:12];
nVOIs=length(vuse);

plotSig=1;

plotVoxNum=0;

plotAcc=0;
plotD=1;

% numVoxUse=150;
% voxelStr=sprintf('take%dZVox',numVoxUse);
voxelStr = 'allVox';
classStr = 'svmtrain';
kernelStr='linear';
predStr = 'predA';
subMeanStr = 'noSubMean';
% subMeanStr = 'subMean2';
% tirange=[2];

typestr='Z2_oneVsOne';
titlestr='Slope of dprime versus disp';

ymaxes=[1,1];
drange=[-.2,0.6];

chanceVals=[1/2];

dimStrs = {'X position','Z position'};
nDim= length(dimStrs);

condStrs = {'trainStim','trainFixat'};
conduse = 2;

nIter=1000;

close all

%% set up file info, other params

root='/usr/local/serenceslab/maggie/IEMdepth/';
folder='IEMdepth_classif';

saveFig = 0;
figFolder='IEMdepth_figs';
ext='epsc';

sigLevels=[0.05,0.01];

horspacer=0.147;
verspacerbig = 0.0001;
verspacersmall = 0.01;
markersize = 3;

%% load the slope data

fnout=sprintf('%s%s/allSubj_allROIs_%s_%s_slopeByDisp.mat',...
                    root,folder,typestr,voxelStr);

load(fnout)
                
medslopes = prctile(allslopes,50,2);
lbslopes = prctile(allslopes,2.5,2);
ubslopes = prctile(allslopes,97.5,2);

medint = prctile(allint,50,2);
lbint = prctile(allint,2.5,2);
ubint = prctile(allint,97.5,2);

%% Compare each slope value to zero, and FDR correct

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

%% Compare each int value to zero, and FDR correct

pValsInt = min([mean(allint<0,2),mean(allint>0,2)],[],2);     

isSigInt = zeros(nVOIs,2);

for aa=1:length(sigLevels)
    alpha=sigLevels(aa);

    [p_fdr, p_masked] = fdr(pValsInt, alpha);
    isSigInt(:,aa)=p_masked; 

    fdrThresh(aa,2) = p_fdr;
end

outTable = table(medslopes,lbslopes,ubslopes,pValsSlope,isSigSlope,medint,lbint,ubint,pValsInt,isSigInt);

%% Compare slopes among ROIs
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

%% save output

fnTab = [root figFolder filesep 'SlopeInt_dprimeVsDisp_' typestr '_' voxelStr '.mat'];
save(fnTab,'mctable','outTable','fdrThresh','fdrLabs','fdrThreshBetween');

