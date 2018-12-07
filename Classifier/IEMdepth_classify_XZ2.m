function IEMdepth_classify_XZ2(subj,VOIs,nVox2Use,nIter,subMean2)

% Perform SVM classification analysis for X and Z position, using a two-way 
% scheme where all trials are grouped into either front-back, or left-right

% Test the significance of decoding performance by performing a permutation 
% test, where the labels of both training and testing sets are shuffled.

% Output can be visualized using IEMdepth_plotClassifier_XZ2

% All main analyses were run using the default values for input arguments.

% MMH 9/14/17

%% define subjects, VOIs, etc

rndseed = [442445];
rng(rndseed,'twister');

if nargin <5
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

nSubj=length(subj);
nVOIs=length(VOIs);

%what are the classifier types we want to compare?
classTypeStrs={'Front-Back','Left-Right'};

classSaveStr = 'XZ2_singleTrialPreds';

%set conditions to do: 1:2 is stim:fix
cc=2;

%parameters for the classifier
kernelStr='linear';
functStr='svmtrain (-t 0 -q)';
classStr='svmtrain';
% subMean=1;
subMeanStrs = {'noSubMean','subMean2'};
subMeanStr=subMeanStrs{subMean2+1};
nBalanceIter=100;
shuffleStr='Shuffle all pred labels';
usingA=1;
predStrs={'predB','predA'};
predStr=predStrs{usingA+1};

useAll=1;
absStat=1;
%parameters for the subjects/VOIs
if useAll
    sessStr='_allSess';
else
    sessStr='_allGoodRuns';
end

% VOIs={'V1','V2','V3','V4','V3A','V3B','IPS0','IPS1','IPS2','IPS3','IPS1-3','allDors'};

locSignStrs={'posVoxOnly','allVoxAbs'};
locSignStr=locSignStrs{absStat+1};

% stereostrs={'_screenx','_stereox'};
% stereostr=stereostrs{usestereox+1};

%locations of files
root = '/usr/local/serenceslab/maggie/IEMdepth/';
folder='IEMdepth_classif';
% voxelfolder='IEMdepth_voxelFStats';
% termind=2;  % always using the second term; the selectivity for Z position

condStrs={'trainStim','trainFixat'};
cc = 2;

%% loop over the classifier types (x, z, etc), conditions, and subjects

        for ss = 1:size(subj,2)     

            % load trialData file
            trialData_fn = sprintf('%sIEMdepth_trialData/%s_allROIs_%s%s.mat',root,...
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

            if strcmp(classTypeStrs{1},'Front-Back') && strcmp(classTypeStrs{2},'Left-Right')
                                       
                termind=2;
                [unique72,~,groupReal72]=unique(stimLocsAll,'rows');
                                
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

            end
            
           if sum(groupRealZ>0)~=numel(groupRealZ) || sum(groupRealX>0)~=numel(groupRealX) || sum(groupReal72>0)~=numel(groupReal72)...
                    || length(unique(groupRealX))~=2 || length(unique(groupRealZ))~=2 || length(unique(groupReal72))~=72
                warning('error in group assignments');
                pause;
            end
            fprintf('Subject %s, grouping according to %s: %d/%d trials assigned\n',subj{ss},shuffleStr,sum(groupRealZ>0),numel(groupRealZ));

            %just take the trials for this condition - separate out the
            %group labels we want, and also the run labels we want
            groupRealZ=groupRealZ(condAll==cc,:);
            groupRealX=groupRealX(condAll==cc,:);
            groupReal72=groupReal72(condAll==cc,:);
            runsAll=runsAll(condAll==cc,:);  
            fprintf('separating by conditions before training: %d/%d trials in %s\n',size(groupRealZ,1),size(condAll,1),condStrs{cc});

%             %load voxel F scores for each CV
%             fn = sprintf('%s%s/IEMdepth_voxelAnovaTwoWay_%s.mat',...
%             root,voxelfolder,subj{ss});              
%             load(fn);

             %load voxel F scores for each CV
%             fn = sprintf('%s%s/IEMdepth_voxelAnovaOneWay6Z_%s.mat',...
%             root,voxelfolder,subj{ss});              
%             load(fn);
            
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
%                 for rr=1:nruns
% 
%                     %choose the voxels to use on this CV
%                     [~,pOrderAsc]=sort(squeeze(voxelAnova(vv).voxelP(rr,:))','ascend');
%                     
%                     if length(pOrderAsc)>=nVox2Use;
%                         voxelindsuse=pOrderAsc(1:nVox2Use);
%                     else
%                         voxelindsuse=pOrderAsc;
%                     end
%                     
%                     %store this sort order to save time later
%                     voxelindsuse_all=cat(1,voxelindsuse_all,voxelindsuse');
%                end
                
                %% run the decoder with real data
                [a,d] = IEMdepth_classifier(predAll,groupRealX,runsAll,classStr,nBalanceIter,voxelindsuse_all,subMean2);
                accRealX = a;
                dRealX = d;
                [a,d] = IEMdepth_classifier(predAll,groupRealZ,runsAll,classStr,nBalanceIter,voxelindsuse_all,subMean2);
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
                    [a,d] = IEMdepth_classifier(predAll,groupRandX,runsAll,classStr,nBalanceIter,voxelindsuse_all,subMean2);
                    accRandX(ii) = a;
                    dRandX(ii) = d;
                    [a,d] = IEMdepth_classifier(predAll,groupRandZ,runsAll,classStr,nBalanceIter,voxelindsuse_all,subMean2);
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
                classStruct(vv).shuffleStr=shuffleStr;
                classStruct(vv).classTypeStrs=classTypeStrs;
                classStruct(vv).numVoxTot=size(predAll,2);
                classStruct(vv).numVoxUseSet=nVox2Use;
                classStruct(vv).numVoxUseActual=length(voxelindsuse);
%                 classStruct(vv).numTrialsTotal=size(predAll,1);
%                 classStruct(vv).totalPredRank=rank(predAll);
                classStruct(vv).numCV=nruns;
                classStruct(vv).voxStr=voxStr;
                
            end
            
            %save the results - all VOIs, this subject, this cond
            accDist_fn = sprintf('%s%s/%s_allROIs_%s_%s_%s_%s_%s_%s_%s.mat',...
                    root,folder,subj{ss},classSaveStr,condStrs{cc},voxStr,predStr,classStr,kernelStr,subMeanStr);
            save(accDist_fn, 'classStruct');
            fprintf('Saving to: %s...\n',accDist_fn);
            
            
        end %end subj loop  

end

