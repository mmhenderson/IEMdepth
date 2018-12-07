function IEMdepth_anovaROIxDim_XZ2()

% run an anova on d' values for 2-way pairwise X and Z classification
% 2-way RM anova with ROI and dimension as categorical factors
% using ranova for anova, and multcompare for pairwise comparisons

% This analysis does not include IPS2-3, but does include LO1-2

% MMH 9/19/17

%% define subjects and flags for what to do
% close all

subj = {'AI','AP','BB','BC','BD','BJ','BM','BN','BO'};
VOIs={'V1','V2','V3','V4','V3A','V3B','IPS0','IPS1','IPS2','IPS3','LO1','LO2'};
vuse=[1:7,11:12];

nSubj=length(subj)
nVOIs=length(vuse)

typestr = 'XZ2_singleTrialPreds';

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

figFolder='IEMdepth_figs';


%% set up file info, other params

root='/usr/local/serenceslab/maggie/IEMdepth/';

folder='IEMdepth_classif';

nDim=2;

thisind=0;
X = zeros(nSubj*length(vuse)*nDim,4);

%% loop over subs
for ss=1:nSubj   
    
    fns=sprintf('%s%s/%s_allROIs_%s_%s_%s_%s_%s_%s_%s.mat',...
                    root,folder,subj{ss},typestr,condStrs{conduse},voxelStr,predStr,classStr,kernelStr,subMeanStr);

    load(fns);

    for vv=1:length(vuse)

        thisd = classStruct(vuse(vv)).dRealX;

        thisind = thisind+1;
        X(thisind,:) = [thisd, vv, 1,ss];
        
        thisd = classStruct(vuse(vv)).dRealZ;
        
        thisind = thisind+1;
        X(thisind,:) = [thisd, vv, 2,ss];

    end
   
end

%% code from John

ii=0;
data = zeros(nSubj,18);
for dd=1:2
    for vv=1:nVOIs    
        ii=ii+1;
        data(:,ii) = X(X(:,2)==vv & X(:,3)==dd,1);
    end
end

% Create a table storing the respones
varNames = {'Y1','Y2','Y3','Y4','Y5','Y6','Y7','Y8','Y9','Y10',...
    'Y11','Y12','Y13','Y14','Y15','Y16','Y17','Y18',};
t = array2table(data,'VariableNames',varNames);
% Create a table reflecting the within subject factors 'TestCond', 'Attention', and 'TMS' and their levels
factorNames = {'SpaceDim', 'ROI'};

within = table({'X';'X';'X';'X';'X';'X';'X';'X';'X';...
    'Z';'Z';'Z';'Z';'Z';'Z';'Z';'Z';'Z'},...
    [VOIs(vuse)';VOIs(vuse)'],'VariableNames',factorNames);

% fit the repeated measures model
rm = fitrm(t,'Y1-Y18~1','WithinDesign',within);

mauchly_tbl = mauchly(rm)

% run my repeated measures anova here
[ranovatbl] = ranova(rm, 'WithinModel','SpaceDim*ROI')

%% multiple comparisons
% compare which ROIs are significantly different 
mctable = multcompare(rm,'ROI','By','SpaceDim');

%% save the output

fnSave = [root figFolder filesep 'Dprime_anovaROIxDim_XZ2_' voxelStr '.mat'];
save(fnSave,'ranovatbl','mctable','mauchly_tbl');

%% old

% % each pairing is reported in both directions - look at only one direction
% posDiffs = mctable{:,3}>0 & mctable{:,5}<alpha;
% sigDiffs = mctable{posDiffs,[1,2]}

%%
%subjects as random factor - treating disp as a continuous variable
% [~,ftab,stats] = anovan(X(:,1),{X(:,2),X(:,3),X(:,4)},'random',[3],'model','interaction','display','on');
% [~,ftab,stats] = anovan(X(:,1),{X(:,2),X(:,3)},'model','interaction','display','on');

% ftab
% multcompare(stats,'Display','on')
% [~,ftab] = anovan(X(:,1),{X(:,3),X(:,2),X(:,4)},'random',[3],'continuous',[1],'model','full','display','on');
