function IEMdepth_classify_Z6(subj,VOIs,nVox2Use,nIter,subMean2)

% SVM classification for Z position - perform every possible 2-way
% comparison of Z positions, collapsing over the X dimension. 
% There are 6C2 = 15 comparisons made

% Test the significance of decoding performance by performing a permutation 
% test, where the labels of both training and testing sets are shuffled.

% Output can be visualized using IEMdepth_plotClassifier_XZ2

% This script requires the output of IEMdepth_findTopVox in order to run -
% if nVox2Use is set to a non-empty value, then nVox2Use is the maximum
% number of voxels per ROI to use for classification, sorted by F score in
% a two-way anova for X and Z position (performed on an independent data
% partition from the test set)
% If nVox2Use is empty, then all voxels will be used. 

% subMean2 is a flag for whether to subtract the mean pattern across voxels
% from each trial pattern before training classifier

% MMH 9/14/17


%% define subjects, VOIs, etc

rndseed = [211572];
rng(rndseed,'twister');

if nargin <5
    subMean2=0;
end

if nargin<4 || isempty(nIter)
    nIter=2;
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
% classTypeStrs={'Front-Back','Left-Right'};

posPairList = combnk(1:6,2);
for pp=1:length(posPairList)
    posPairStrs{pp} = sprintf('%dvs%d',posPairList(pp,1),posPairList(pp,2));
end

classSaveStr = sprintf('Z6_oneVsOne');

%set conditions to do: 1:2 is stim:fix
cc=2;
condStrs={'trainStim','trainFixat'};
 
%parameters for the classifier
kernelStr='linear';
functStr='svmtrain (-t 0 -q)';
classStr='svmtrain';
% subMean=1;
subMeanStrs = {'noSubMean','subMean2'};
subMeanStr=subMeanStrs{subMean2+1};
nBalanceIter=100;
% shuffleStr='Shuffle all pred labels';
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
voxelfolder='IEMdepth_voxelFStats';
% termind=2;  % always using the second term; the selectivity for Z position


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
%             xAll=round(stimLocsAll(:,1),1);
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
%                         
%             %load voxel F scores for each CV
%             fn = sprintf('%s%s/IEMdepth_voxelAnovaTwoWay_%s.mat',...
%             root,voxelfolder,subj{ss});              
%             load(fn);

             %load voxel F scores for each CV
            fn = sprintf('%s%s/IEMdepth_voxelAnovaOneWay6Z_%s.mat',...
            root,voxelfolder,subj{ss});              
            load(fn)
            
            for vv = 1:length(VOIs)
                
                %define the predictors for each VOI here
                if usingA
                    predAll = trialDataAll(vv).aAll;
                else
                    predAll = trialDataAll(vv).bAll;
                end 
                %% use the p-values for each voxel to decide which set to use

                voxelindsuse_all=[];
                for rr=1:nruns
                    %figure out which trials we're using for model estimation
                    %(training) and to compute channel responses (testing)
                    traininds = runsAll~=rr;
                    testinds = ~traininds;

                    %choose the voxels to use on this CV
                    [~,pOrderAsc]=sort(squeeze(voxelAnova(vv).voxelP(rr,:))','ascend');
                    if length(pOrderAsc)>=nVox2Use;
                        voxelindsuse=pOrderAsc(1:nVox2Use);
                    else
                        voxelindsuse=pOrderAsc;
                    end

                    %store this sort order to save time later
                    voxelindsuse_all=cat(1,voxelindsuse_all,voxelindsuse');
                    %fprintf('%s, crossval run %d/%d: training mat is %d by %d with rank=%d\n',VOIs{vv},rr,nruns,size(predAll(traininds,voxelindsuse),1),size(predAll(traininds,voxelindsuse),2),rank(predAll(traininds,voxelindsuse)));
                end
                
                % initialize matrices to store results
                accReal = zeros(length(posPairList),1);
                dReal = zeros(length(posPairList),1);
                
                accRand = zeros(length(posPairList),nIter);
                dRand = zeros(length(posPairList),nIter);
                
                
                for pp=1:length(posPairList)
                    
                    %% select the trials to use - 2/6 z positions, and one condition
                    
                    indsuse = groupRealZ==posPairList(pp,1) | groupRealZ==posPairList(pp,2);
                    indsuse = indsuse & condAll==cc;
                    
                    %% define the (nTrials x nVoxels) matrix of response patterns

                    predUse = predAll(indsuse,:);
                    groupUse = groupRealZ(indsuse,:);
                    runLabsUse = runsAll(indsuse,:);
                    
                    if any(unique(groupUse)'~=posPairList(pp,:)) || any(unique(runLabsUse)'~=1:nruns)
                        error('error selecting runs to use')
                    end
                    
                    fprintf('%s %s %s - %s - using %d/%d trials\n',subj{ss},VOIs{vv},condStrs{cc},posPairStrs{pp},length(groupUse),length(groupRealZ))
                    
%                     if pp==1
%                         % in each comparison, 1/3 of the trials get used
%                         realLabs = zeros(length(posPairList),length(groupUse));
%                         predLabs = zeros(length(posPairList),length(groupUse));
%                         predLabsRand = zeros(length(posPairList),length(groupUse),nIter);
%                     end
%                     
%                     realLabs(pp,:) = groupUse;
                    
                    %% run the decoder with real data
                    [a,d] = IEMdepth_classifier(predUse,groupUse,runLabsUse,classStr,nBalanceIter,voxelindsuse_all,subMean2);
                    accReal(pp) = a;
                    dReal(pp) = d;
%                     predLabs(pp,:) = pl;

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
                        
                        %% run the decoder with shuffled data
                        [a,d] = IEMdepth_classifier(predUse,groupRand,runLabsUse,classStr,nBalanceIter,voxelindsuse_all,subMean2);
                        accRand(pp,ii) = a;
                        dRand(pp,ii) = d;
%                         predLabsRand(pp,:,ii) = pl;

                    end

                end

                %% save the output
                classStruct(vv).accReal=accReal;
                classStruct(vv).accRand=accRand;

                classStruct(vv).dReal=dReal;
                classStruct(vv).dRand=dRand;
%                 
%                 classStruct(vv).realLabs = realLabs;
%                 classStruct(vv).predLabs = predLabs;
%                 classStruct(vv).predLabsRand = predLabsRand;

                classStruct(vv).predStr=predStr;
                classStruct(vv).kernelStr=kernelStr;
                classStruct(vv).functStr=functStr;
                classStruct(vv).classStr=classStr;               
%                 classStruct(vv).shuffleStr=shuffleStr;
                classStruct(vv).voxStr=voxStr;
                
                classStruct(vv).posPairList = posPairList;
                classStruct(vv).posPairStrs = posPairStrs;
                
                classStruct(vv).numVoxTot=size(predAll,2);
                classStruct(vv).numVoxUseSet=nVox2Use;
                classStruct(vv).numVoxUseActual=length(voxelindsuse);

                classStruct(vv).numCV=nruns;
                classStruct(vv).nBalanceIter=nBalanceIter;
                
            end
            
            %save the results - all VOIs, this subject, this cond
            accDist_fn = sprintf('%s%s/%s_allROIs_%s_%s_%s_%s_%s_%s_%s.mat',...
                    root,folder,subj{ss},classSaveStr,condStrs{cc},voxStr,predStr,classStr,kernelStr,subMeanStr);
            save(accDist_fn, 'classStruct');
            fprintf('Saving to: %s...\n',accDist_fn);
            
    
        end %end subj loop  
%     end %end cond loop
% end %end type loop


end
