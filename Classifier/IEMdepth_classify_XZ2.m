function IEMdepth_classify_XZ2(subj,VOIs,nVox2Use,nIter,subMean2)

% Perform SVM classification analysis for X and Z position, using a two-way 
% scheme where all trials are grouped into either front-back, or left-right

% Test the significance of decoding performance by performing a permutation 
% test, where the labels of both training and testing sets are shuffled.

% Output gets saved as a .mat file for each subj.

% All main analyses were run using the default values for input arguments.

% MMH 9/14/17

%% set up some path stuff
% This is whatever directory contains the folder "IEMdepth_trialData"
root = '/usr/local/serenceslab/maggie/IEMdepth/';
% This is where trialdata structures should be.
load_folder = [root 'IEMdepth_trialData'];
% This is the path where the classifier output files are saved.
save_folder = [root 'IEMdepth_classif'];
if ~isfolder(save_folder)
    mkdir(save_folder)
end
%%
% Make sure we have the libSVM version of "svmtrain" on the path. There is
% a matlab version with the same name and we don't want to mistakenly use that. 
% libSVM 3.1 can be downloaded from:
% https://www.csie.ntu.edu.tw/~cjlin/libsvm/
svm_version = which('svmtrain');
if ~strcmp(svm_version(end-14:end), 'svmtrain.mexa64')
    error('Verify that you are using the libSVM (3.1) version of "svmtrain". If you think you are using the correct version, you can comment this out and proceed.')
end

%% define subjects, VOIs, etc

rndseed = [442445];
rng(rndseed,'twister');

if nargin <5
    % subtract mean over the voxel dimension? Usually set to no.
    subMean2=0;
end

if nargin<4 || isempty(nIter)
    nIter=1000;
end

if nargin<3 
    nVox2Use = [];
end

if isempty(nVox2Use)
    nVox2Use=10000;
    voxStr = 'allVox';
else
    voxStr = sprintf('take%dZ1WayVox',nVox2Use);
end

if nargin<2 || isempty(VOIs)
    VOIs={'V1','V2','V3','V4','V3A','V3B','IPS0','IPS1','IPS2','IPS3','LO1','LO2'};
end

if nargin==0 || isempty(subj)
    subj = {'AI','AP','BB','BC','BD','BJ','BM','BN','BO'};
end

%set conditions to do: 1:2 is stim:fix
cc=2;
condStrs={'trainStim','trainFixat'};
 
%parameters for the classifier
kernelStr='linear';
functStr='svmtrain (-t 0 -q)';
classStr='svmtrain';

classSaveStr = 'XZ2_singleTrialPreds';

% parameters for the data to load
sessStr='_allSess';
locSignStr='allVoxAbs';

subMeanStrs = {'noSubMean','subMean2'};
subMeanStr=subMeanStrs{subMean2+1};

% use average signal over 3-4 TRs? alternative is to use a beta weight.
usingA=1;
predStrs={'predB','predA'};
predStr=predStrs{usingA+1};

%% loop over subjects

for ss = 1:size(subj,2)     

    % load trialData file
    trialData_fn = sprintf('%s/%s_allROIs_%s%s.mat',load_folder,...
        subj{ss},locSignStr,sessStr);
    fprintf('Loading trial data for %s from: %s...\n',subj{ss},trialData_fn);
    load(trialData_fn);

    %% define trial labels that are same for all VOIs
    %define trial labels that are same for all VOIs
    condAll = trialDataAll(1).condAll;
    runsAll = trialDataAll(1).runsAll;
    nruns = length(unique(runsAll));
    stimLocsAll = round(trialDataAll(1).stimLocsAll,4);
    xAll=round(stimLocsAll(:,1),1);
    zAll=round(stimLocsAll(:,3),1);

    groupRealZ=zeros(size(zAll,1),1);
    groupRealX=zeros(size(xAll,1),1);

    [~,~,groupReal72]=unique(stimLocsAll,'rows');

    groupRealZ(zAll<0,1)=1;
    groupRealZ(zAll>0,1)=2;
    groupRealX(xAll<0,1)=1;
    groupRealX(xAll>0,1)=2;

    %make an array to easily convert between positions- will
    %shuffle across all 72, then simplify those indices to just
    %the x, z position groups we want.
    c72toZ=zeros(72,1);
    c72toX=zeros(72,1);

    for xz=1:length(unique(groupRealZ))
        tmpinds=unique(groupReal72(groupRealZ==xz));
        c72toZ(tmpinds)=xz;
        tmpinds2=unique(groupReal72(groupRealX==xz));
        c72toX(tmpinds2)=xz;
    end

    if sum(groupRealZ>0)~=numel(groupRealZ) || sum(groupRealX>0)~=numel(groupRealX) || sum(groupReal72>0)~=numel(groupReal72)...
            || length(unique(groupRealX))~=2 || length(unique(groupRealZ))~=2 || length(unique(groupReal72))~=72
        warning('error in group assignments');
        pause;
    end
    fprintf('Subject %s, labeling trials for two-way class: %d/%d trials assigned\n',subj{ss},sum(groupRealZ>0),numel(groupRealZ));

    %just take the trials for this condition - separate out the
    %group labels we want, and also the run labels we want
    groupRealZ=groupRealZ(condAll==cc,:);
    groupRealX=groupRealX(condAll==cc,:);
    groupReal72=groupReal72(condAll==cc,:);
    runsAll=runsAll(condAll==cc,:);  
    fprintf('separating by conditions before training: %d/%d trials in %s\n',size(groupRealZ,1),size(condAll,1),condStrs{cc});

    for vv = 1:length(VOIs)


        %% define the (nTrials x nVoxels) matrix of response patterns

        %define the predictors for each VOI here
        if usingA
            predAll = trialDataAll(vv).aAll;
        else
            predAll = trialDataAll(vv).bAll;
        end 

        %just take the trials for this condition
        predAll=predAll(condAll==cc,:); 

        %% use the p-values for each voxel to decide which set to use

        voxelindsuse_all=[];

        %% run the decoder with real data
        [a,d] = IEMdepth_classifier(predAll,groupRealX,runsAll,classStr,voxelindsuse_all,subMean2);
        accRealX = a;
        dRealX = d;
        [a,d] = IEMdepth_classifier(predAll,groupRealZ,runsAll,classStr,voxelindsuse_all,subMean2);
        accRealZ = a;
        dRealZ = d;

        fprintf('    ...accuracy on real data: %.2f in X, %.2f in Z. Starting shuffle over %d iterations...\n',accRealX, accRealZ ,nIter);

        %% run the decoder of nIter of shuffled data

        accRandX=nan(nIter,1);
        accRandZ=nan(nIter,1);
        dRandX = nan(nIter,1);
        dRandZ = nan(nIter,1);

        parfor ii=1:nIter

            %shuffle the labels (within runs, all positions at same
            %time)
            unruns=unique(runsAll);

            groupRand72=nan(length(groupReal72),1);

            for cv=1:length(unruns)
                theseinds=runsAll==unruns(cv);
                theselabs=groupReal72(theseinds);
                groupRand72(theseinds)=theselabs(randperm(length(theselabs)));        
            end

            groupRandX=c72toX(groupRand72);
            groupRandZ=c72toZ(groupRand72);


            %% run the decoder with rand data
            [a,d] = IEMdepth_classifier(predAll,groupRandX,runsAll,classStr,voxelindsuse_all,subMean2);
            accRandX(ii) = a;
            dRandX(ii) = d;
            [a,d] = IEMdepth_classifier(predAll,groupRandZ,runsAll,classStr,voxelindsuse_all,subMean2);
            accRandZ(ii) = a;
            dRandZ(ii) = d;

        end

        %% save the output
        classStruct(vv).accRealX=accRealX;
        classStruct(vv).accRandX=accRandX;
        classStruct(vv).accRealZ=accRealZ;
        classStruct(vv).accRandZ=accRandZ;
        classStruct(vv).dRealX=dRealX;
        classStruct(vv).dRandX=dRandX;
        classStruct(vv).dRealZ=dRealZ;
        classStruct(vv).dRandZ=dRandZ;
        classStruct(vv).predStr=predStr;
        classStruct(vv).kernelStr=kernelStr;
        classStruct(vv).functStr=functStr;
        classStruct(vv).classStr=classStr;
        classStruct(vv).numVoxTot=size(predAll,2);
        classStruct(vv).numCV=nruns;
        classStruct(vv).voxStr=voxStr;

    end

    %save the results - all VOIs, this subject, this cond
    accDist_fn = sprintf('%s/%s_allROIs_%s_%s_%s_%s_%s_%s_%s.mat',...
            save_folder,subj{ss},classSaveStr,condStrs{cc},voxStr,predStr,classStr,kernelStr,subMeanStr);
    save(accDist_fn, 'classStruct');
    fprintf('Saving to: %s...\n',accDist_fn);

end %end subj loop 

end

