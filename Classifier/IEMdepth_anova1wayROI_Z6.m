function IEMdepth_anova1wayROI_Z6()

% run a one-way anova for two-way Z decoding accuracies
% factor is ROIs, categorical
% uses ranova function in matlab for RM test, and multcompare for multiple
% comparisons 

% This analysis does not include IPS2-3, but does include LO1-2

% MMH 9/19/17

%% define subjects and flags for what to do
% close all

subj = {'AI','AP','BB','BC','BD','BJ','BM','BN','BO'};
% subj = {'AP','BB','BC','BD','BJ','BM','BN','BO'};
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

%% set up file info, other params

root='/usr/local/serenceslab/maggie/IEMdepth/';
figFolder='IEMdepth_figs';

folder='IEMdepth_classif';

thisind=0;
X = zeros(nSubj*length(vuse),3);

% dall= zeros(nSubj,nVOIs);

%% loop over subs
for ss=1:nSubj   
    
    fns=sprintf('%s%s/%s_allROIs_%s_%s_%s_%s_%s_%s_%s.mat',...
                    root,folder,subj{ss},typestr,condStrs{conduse},voxelStr,predStr,classStr,kernelStr,subMeanStr);

    load(fns);

    for vv=1:length(vuse)

        thisd = mean(classStruct(vuse(vv)).dReal);

        thisind = thisind+1;
        X(thisind,:) = [thisd, vv, ss];       

    end
   
end

%% code from John

ii=0;
data = zeros(nSubj,nVOIs);

for vv=1:nVOIs    
    ii=ii+1;
    data(:,ii) = X(X(:,2)==vv,1);
end

% Create a table storing the respones
varNames = {'Y1','Y2','Y3','Y4','Y5','Y6','Y7','Y8','Y9'};

t = array2table(data,'VariableNames',varNames);
% Create a table reflecting the within subject factors 'TestCond', 'Attention', and 'TMS' and their levels
factorNames = {'ROI'};

within = table(VOIs(vuse)','VariableNames',factorNames);

% fit the repeated measures model
rm = fitrm(t,'Y1-Y9~1','WithinDesign',within);

mauchly_tbl = mauchly(rm);

% run my repeated measures anova here
[ranovatbl] = ranova(rm, 'WithinModel','ROI')

%% multiple comparisons
% compare which ROIs are significantly different - averaging across X/Z
mctable = multcompare(rm,'ROI');

% % each pairing is reported in both directions - look at only one direction
posDiffs = mctable{:,3}>0 & mctable{:,5}<alpha;
sigDiffs = mctable{posDiffs,[1,2]}

%% save the output

fnSave = [root figFolder filesep 'Dprime_anova1wayROI_Z6_' voxelStr '.mat'];
save(fnSave,'ranovatbl','mctable','sigDiffs','mauchly_tbl');


