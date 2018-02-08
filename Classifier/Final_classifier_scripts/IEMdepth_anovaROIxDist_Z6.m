function IEMdepth_anovaROIxDist_Z6()

% run an anova on d' values for 2-way pairwise Z classification
% Fixed effects are ROI (categorical) and disparity difference (continuous)
% rm with ranova

% This analysis does not include IPS2-3, but does include LO1-2

% MMH 9/14/17

%% define subjects and flags for what to do
% close all

subj = {'AI','AP','BB','BC','BD','BJ','BM','BN','BO'};
VOIs={'V1','V2','V3','V4','V3A','V3B','IPS0','IPS1','IPS2','IPS3','LO1','LO2'};
vuse=[1:8,11:12];

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

%% set up file info, other params

root='/usr/local/serenceslab/maggie/IEMdepth/';

folder='IEMdepth_classif';
figFolder='IEMdepth_figs';

zDispList = [-35.8,-22.7,-8,8.5,27.2,48.4];

posPairList = combnk(1:6,2);
dispPairList = reshape(zDispList(posPairList(:)),15,2);
dispDiffList = dispPairList(:,2)-dispPairList(:,1);
[distListSort,distOrder] = sort(dispDiffList,'ascend');

undist = unique(dispDiffList);
nDist = length(undist);
    
% accs_allsub = zeros(nVOIs,nSubj,nDist);
d_allsub = zeros(nVOIs,nSubj,nDist);

thisind=0;
X = zeros(nSubj*length(vuse)*nDist,4);

%% loop over subs
for ss=1:nSubj   
    
    fns=sprintf('%s%s/%s_allROIs_%s_%s_%s_%s_%s_%s_%s.mat',...
                    root,folder,subj{ss},typestr,condStrs{conduse},voxelStr,predStr,classStr,kernelStr,subMeanStr);

    load(fns);

    for vv=1:length(vuse)
        nVox(ss,vv) = classStruct(vuse(vv)).numVoxUseActual;
        
        for dd=1:nDist
            
%             theseaccs = classStruct(vv).accReal(distOrder(dd),:);
%             accs_allsub(vv,ss,dd) = mean(theseaccs);
            
            thesed = classStruct(vuse(vv)).dReal(distOrder(dd),:);
            d_allsub(vv,ss,dd) = mean(thesed);
            
            thisind = thisind+1;
            X(thisind,:) = [mean(thesed), vv, distListSort(dd),ss];
        end
        
    end
   
end

%% code from John (ranova)

ii=0;
data = zeros(nSubj,nVOIs*nDist);
for dd=1:nDist
    for vv=1:nVOIs    
        ii=ii+1;
        distCells(ii) = distListSort(dd);
        vCells{ii} = VOIs{vuse(vv)};
        data(:,ii) = X(X(:,2)==vv & X(:,3)==distListSort(dd),1);
        varNames{ii} = sprintf('Y%d',ii);
    end
end

% Create a table storing the respones
% varNames = {'Y1','Y2','Y3','Y4','Y5','Y6','Y7','Y8','Y9','Y10',...
%     'Y11','Y12','Y13','Y14','Y15','Y16','Y17','Y18','Y19','Y20'};
t = array2table(data,'VariableNames',varNames);
% Create a table reflecting the within subject factors 'TestCond', 'Attention', and 'TMS' and their levels
factorNames = {'Disp', 'ROI'};

within = table(distCells',vCells','VariableNames',factorNames);

% fit the repeated measures model
rm = fitrm(t,'Y1-Y150~1','WithinDesign',within);

% run my repeated measures anova here
[ranovatbl] = ranova(rm, 'WithinModel','Disp*ROI')

%% multiple comparisons
% compare which ROIs are significantly different - averaging across X/Z
mctableROI = multcompare(rm,'ROI');

% each pairing is reported in both directions - look at only one direction
posDiffs = mctableROI{:,3}>0 & mctableROI{:,5}<0.05;
sigDiffs = mctableROI{posDiffs,[1,2]};

%% multiple comparisons
% compare which ROIs are significantly different - averaging across X/Z
mctableDisp = multcompare(rm,'Disp');

% each pairing is reported in both directions - look at only one direction
posDiffs = mctableDisp{:,3}>0 & mctableDisp{:,5}<0.05;
sigDiffs = mctableDisp{posDiffs,[1,2]};


%% save the output

fnSave = [root figFolder filesep 'Dprime_anovaROIxDist_Z6_' voxelStr '.mat'];
save(fnSave,'ranovatbl');

%% rm_anova2

% stats = rm_anova2(X(:,1),X(:,4),X(:,2),X(:,3),{'ROI','dim'});
% 
% [thisF,thisP,factornames] = my_RMAOV2(X,alpha);


%% anovan
%subjects as random factor - treating disp as a continuous variable
% [~,ftab,stats] = anovan(X(:,1),{X(:,2),X(:,3),X(:,4)},'random',[3],'continuous',[2],'model','full','display','on');

% multcompare(stats,'Display','on')
% [~,ftab] = anovan(X(:,1),{X(:,3),X(:,2),X(:,4)},'random',[3],'continuous',[1],'model','full','display','on');
