%% do ANOVA test for each voxel's position selectivity
% these F scores can be used to select voxels for use in classifier script
% two way anova using anovan - use 6 x and 6 z positions
% save out a structure array voxelAnova for each subject

% MMH 9/14/17

close all;clear

%% set up subjects and data root, etc
sessStr='_allSess';

subj = {'AI','AP','BB','BC','BD','BJ','BM','BN','BO'};

VOIs={'V1','V2','V3','V4','V3A','V3B','IPS0','IPS1','IPS2','IPS3','LO1','LO2'};

dimStrs={'Z'};

absStat=1; %take the absolute value of the localizer vmp for each voxel? this tells which MRs and trialdata structs to load
locSignStrs={'posVoxOnly','allVoxAbs'};
locSignStr=locSignStrs{absStat+1};

nDim=length(dimStrs);
nSubj=length(subj);
nVOIs=length(VOIs);

root = '/usr/local/serenceslab/maggie/IEMdepth/';

savefolder='IEMdepth_voxelFStats';

sigLevel=0.01;

nTerms=1;

% numUse=100;

termNames={'six-way Z pos'};

%% load data for each subj, do anova

for ss = 1:nSubj  
    trialData_fn = sprintf('%sIEMdepth_trialData/%s_allROIs_%s_allSess.mat',root,...
        subj{ss},locSignStr);
    fprintf('Loading trial data for %s from: %s...\n',subj{ss},trialData_fn);
    load(trialData_fn);

    for vv = 1:nVOIs
        fprintf('Starting single-voxel ANOVA for subject %s, %s\n',subj{ss},VOIs{vv});           
                   
        %look at data in this VOI - do anova for x and z position       
        vdat = trialDataAll(vv);          
        voxResp = vdat.aAll;
        stimLocsAll = vdat.stimLocsAll;
        xAll=round(stimLocsAll(:,1),1);
        zAll=round(stimLocsAll(:,3),1);
        runsAll=vdat.runsAll;
        nruns=length(unique(runsAll));
        
        clear vdat
        
        nVox=size(voxResp,2);      

        %these cutoffs just divide the space into thirds for
        %left, right, middle (or back, front, middle)
%         groupX=ones(size(voxResp,1),1);
%         groupX(xAll>-1.3)=2;
%         groupX(xAll>-0.6)=3;
%         groupX(xAll>=0)=4;
%         groupX(xAll>=0.7)=5;
%         groupX(xAll>=1.3)=6;
        
        groupZ=ones(size(voxResp,1),1);
        groupZ(zAll>-1.2)=2;
        groupZ(zAll>-0.9)=3;
        groupZ(zAll>=0)=4;
        groupZ(zAll>=0.8)=5;
        groupZ(zAll>=1.2)=6;
              
        
        %do an F test for each voxel
        pValsEachVox=nan(nruns,nVox);
        fValsEachVox=nan(nruns,nVox);  
        
%         bestInds=zeros(nruns,nTerms,nVox);
        
        for rr=1:nruns
            fprintf('    CV %d/%d\n',rr,nruns);
            % figure out which trials we're using for model estimation
            % (training) and to compute channel responses (testing)
            traininds = runsAll~=rr;
            testinds = ~traininds;
            
            for vx=1:nVox
                [~,stats] = anovan(voxResp(traininds,vx), {groupZ(traininds,:)},'model','full','display','off');
                pValsEachVox(rr,vx)=cell2mat(stats(2,7));
                fValsEachVox(rr,vx)=cell2mat(stats(2,6));
            end

        end   
                       
        %store in a structure to load later
        voxelAnova(vv).voxelF=fValsEachVox;
        voxelAnova(vv).voxelP=pValsEachVox;   
%         voxelAnova(vv).bestInds=bestInds;
        voxelAnova(vv).termNames=termNames;
        
    end  %end VOI loop   

    %save the results
    fn = sprintf('%s%s/IEMdepth_voxelAnovaOneWay6Z_%s.mat',...
            root,savefolder,subj{ss});
    save(fn, 'voxelAnova');
    fprintf('Saving results for all ROIs to: %s\n\n',fn);
    
end %end subj loop

