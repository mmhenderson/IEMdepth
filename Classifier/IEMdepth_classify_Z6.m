function IEMdepth_classify_Z6(subj,VOIs,nVox2Use,nIter,subMean2)

% SVM classification for Z position - perform every possible 2-way
% comparison of Z positions, collapsing over the X dimension. 
% There are 6C2 = 15 comparisons made

% Test the significance of decoding performance by performing a permutation 
% test, where the labels of both training and testing sets are shuffled.

% Output gets saved as a .mat file for each subj.

% All analyses in the text were run using the default values for input arguments.

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

rndseed = [211572];
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

% list all the pairwise comparisons between the 6 Z positions.
posPairList = combnk(1:6,2);
for pp=1:length(posPairList)
    posPairStrs{pp} = sprintf('%dvs%d',posPairList(pp,1),posPairList(pp,2));
end

classSaveStr = sprintf('Z2_oneVsOne');

%set conditions to do: 1:2 is stim:fix
cc=2;
condStrs={'trainStim','trainFixat'};
 
%parameters for the classifier
kernelStr='linear';
functStr='svmtrain (-t 0 -q)';
classStr='svmtrain';

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

    condAll = trialDataAll(1).condAll;
    runsAll = trialDataAll(1).runsAll;
    nruns = length(unique(runsAll));
    stimLocsAll = round(trialDataAll(1).stimLocsAll,4);

    zAll=round(stimLocsAll(:,3),1);

    groupRealZ=zeros(size(zAll,1),1);

    zAll=round(zAll,1);

    groupRealZ(zAll==-1.5)=1;
    groupRealZ(zAll==-0.9)=2;
    groupRealZ(zAll==-0.3)=3;
    groupRealZ(zAll== 0.3)=4;
    groupRealZ(zAll== 0.9)=5;
    groupRealZ(zAll== 1.5)=6;

    if sum(groupRealZ>0)~=numel(groupRealZ) ||  length(unique(groupRealZ))~=6 
       error('error in group assignments');

    end

    for vv = 1:length(VOIs)

        %define the predictors for each VOI here
        if usingA
            predAll = trialDataAll(vv).aAll;
        else
            predAll = trialDataAll(vv).bAll;
        end 
       
        % using all voxels
        voxelindsuse_all=[];

        % initialize matrices to store results
        accReal = zeros(length(posPairList),1);
        dReal = zeros(length(posPairList),1);

        accRand = zeros(length(posPairList),nIter);
        dRand = zeros(length(posPairList),nIter);

        %% loop over all 15 pairwise comparisons
        for pp=1:length(posPairList)

            % select the trials to use - 2/6 z positions, and one condition
            indsuse = groupRealZ==posPairList(pp,1) | groupRealZ==posPairList(pp,2);
            indsuse = indsuse & condAll==cc;

            % define the (nTrials x nVoxels) matrix of response patterns
            predUse = predAll(indsuse,:);
            groupUse = groupRealZ(indsuse,:);
            runLabsUse = runsAll(indsuse,:);

            if any(unique(groupUse)'~=posPairList(pp,:)) || any(unique(runLabsUse)'~=1:nruns)
                error('error selecting runs to use')
            end

            fprintf('%s %s %s - %s - using %d/%d trials\n',subj{ss},VOIs{vv},condStrs{cc},posPairStrs{pp},length(groupUse),length(groupRealZ))

            %% run the decoder with real data
            [a,d] = IEMdepth_classifier(predUse,groupUse,runLabsUse,classStr,voxelindsuse_all,subMean2);
            accReal(pp) = a;
            dReal(pp) = d;

            fprintf('    ...accuracy on real data: %.2f. Starting shuffle over %d iterations...\n',accReal(pp) ,nIter);

            %% run the decoder of nIter of shuffled data

            parfor ii=1:nIter

                %shuffle the labels (within runs)
                unruns=unique(runsAll);

                groupRand = zeros(size(groupUse));

                for cv=1:length(unruns)
                    theseinds=runLabsUse==unruns(cv);
                    theselabs=groupUse(theseinds);
                    groupRand(theseinds)=theselabs(randperm(length(theselabs)));        
                end

                % run the decoder with shuffled data
                [a,d] = IEMdepth_classifier(predUse,groupRand,runLabsUse,classStr,voxelindsuse_all,subMean2);
                accRand(pp,ii) = a;
                dRand(pp,ii) = d;

            end
        end

        %% save the output
        classStruct(vv).accReal=accReal;
        classStruct(vv).accRand=accRand;

        classStruct(vv).dReal=dReal;
        classStruct(vv).dRand=dRand;

        classStruct(vv).predStr=predStr;
        classStruct(vv).kernelStr=kernelStr;
        classStruct(vv).functStr=functStr;
        classStruct(vv).classStr=classStr;               
        classStruct(vv).voxStr=voxStr;

        classStruct(vv).posPairList = posPairList;
        classStruct(vv).posPairStrs = posPairStrs;

        classStruct(vv).numVoxTot=size(predAll,2);

        classStruct(vv).numCV=nruns;

    end

    %save the results - all VOIs, this subject, this cond
    accDist_fn = sprintf('%s/%s_allROIs_%s_%s_%s_%s_%s_%s_%s.mat',...
            save_folder,subj{ss},classSaveStr,condStrs{cc},voxStr,predStr,classStr,kernelStr,subMeanStr);
    save(accDist_fn, 'classStruct');
    fprintf('Saving to: %s...\n',accDist_fn);


end %end subj loop  

end
