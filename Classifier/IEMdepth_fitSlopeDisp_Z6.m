function IEMdepth_fitSlopeDisp_Z6()

% fit slopes/intercept to the decoding peformance/disparity difference relationship
% in each ROI, bootstrapping across subs

% This analysis does not include IPS2-3, but does include LO1-2

% MMH 9/14/17

%% define subjects and flags for what to do

rndseed = [879642];
rng(rndseed,'twister');

subj = {'AI','AP','BB','BC','BD','BJ','BM','BN','BO'};
VOIs={'V1','V2','V3','V4','V3A','V3B','IPS0','IPS1','IPS2','IPS3','LO1','LO2'};
vuse=[1:7,11:12];

nSubj=length(subj);
nVOIs=length(vuse);

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
alpha=0.01;

nVox = zeros(nSubj,nVOIs);
nIter=1000;

%% set up file info, other params

root='/usr/local/serenceslab/maggie/IEMdepth/';

folder='IEMdepth_classif';

zDispList = [38.592,25.612,11.124,-5.152,-23.567,-44.573];

posPairList = combnk(1:6,2);
dispPairList = reshape(zDispList(posPairList(:)),15,2);
dispDiffList = dispPairList(:,1)-dispPairList(:,2);
[distListSort,distOrder] = sort(dispDiffList,'ascend');

undist = unique(dispDiffList);
nDist = length(undist);
    
% accs_allsub = zeros(nVOIs,nSubj,nDist);
d_allsub = zeros(nVOIs,nSubj,nDist);

% thisind=0;
% X = zeros(nSubj*length(vuse)*nDist,4);

%% loop over subs
for ss=1:nSubj   
    
    fns=sprintf('%s%s/%s_allROIs_%s_%s_%s_%s_%s_%s_%s.mat',...
                    root,folder,subj{ss},typestr,condStrs{conduse},voxelStr,predStr,classStr,kernelStr,subMeanStr);

    load(fns);

    for vv=1:length(vuse)
        nVox(ss,vv) = classStruct(vuse(vv)).numVoxUseActual;
        
        for dd=1:nDist
            d_allsub(vv,ss,dd) = classStruct(vuse(vv)).dReal(distOrder(dd));
        end
        
    end
   
end

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

%% save the output

fnout=sprintf('%s%s/allSubj_allROIs_%s_%s_slopeByDisp.mat',...
                    root,folder,typestr,voxelStr);

save(fnout, 'allslopes','allint','allslopes_zeroint')

