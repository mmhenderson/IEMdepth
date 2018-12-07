function IEMdepth_mergeTrialData(subj,varargin)
%Merge all the trialData files from this subject into one file, which is
%indexed by visual area, condition, and session
%choose whichever sessions we want to combine, but they have to be from
%same localizer

% version 2 - uses different ROI list, including LO1-LO2

% 9/5/17, MMH

%this version 8/19/16 has the option of selecting only some of the runs -
%if we want to ignore runs, or whole sessions, where the goggles
%malfunctioned or where accuracy was low
%this runs all subjects, can index just one if you want (subjInds below)

%% set up the subjects, sessions etc.

subjInds=1:9;

subj = {'AI','AP','BB','BC','BD','BJ','BM','BN','BO'};
    
VOIs={'V1','V2','V3','V4','V3A','V3B','IPS0','IPS1','IPS2','IPS3','LO1','LO2'};

sessionsAll={{'141','143'},...
            {'141','143','144'},...                
            {'141','142','143'},...
            {'141'},...
            {'141'},...
            {'141'},...
            {'141','142'},...
            {'141','142'},...
            {'141','142'}}; 
  
sessStr='_allSess';
goodRunsAll={{[1:2,5],[1:8]},...
            {[1:4,6:8],[1:6],[1:7]},...
            {[1:7],[1:6],[1:8]},...
            {[1:7]},...
            {[1:8]},...
            {[1:8]},...
            {[1:6],[1:6]},...
            {[1:7],[1:7]},...
            {[1:7],[1:9]}};
        
%this includes every run, except the AI141 runs (3,4) with the
%reversed grid scaling because the positions won't match with
%our model, and one run where AP141 coughed

subj=subj(subjInds);
sessionsAll=sessionsAll(subjInds);
goodRunsAll=goodRunsAll(subjInds);

root='/usr/local/serenceslab/maggie/IEMdepth/';

absStat=1; %take the absolute value of the localizer vmp for each voxel? this tells which MRs and trialdata structs to load
locSignStrs={'posVoxOnly','allVoxAbs'};
locSignStr=locSignStrs{absStat+1};

trialDataStr = {'DepthStim','DepthFix'};

for ss = 1:length(subj)   % loop over subjects
    
    theseSessions=sessionsAll{ss};
    theseGoodRuns=goodRunsAll{ss};
    %initialize cell array for this subject, each cell is a VOI
%     trialDataAll=cell(length(VOIs),1);
     
    for vv=1:length(VOIs)
        
        %initialize arrays for this VOI in this subject
        aAll=[];
        bAll=[];
        stimLocsAll=[];
        condAll=[];
        runsAll=[];
        dirAll = [];
        respAll = [];
        fixChgAll = [];
        
        bXAll=[];
        predlabsXAll=[];
        condlabsXAll=[];
        runlabsXAll=[];

        bZAll=[];
        predlabsZAll=[];
        condlabsZAll=[];
        runlabsZAll=[];
        
        runsPrevSess=0;

        for se = 1:length(theseSessions)
            runsUse=theseGoodRuns{se};

%             runsUse=runsUse{1};
            for cc = 1:length(trialDataStr)

                %% get spinDir
            
                expFileName = 'IEMdepth';
                sn = [subj{ss},sessionsAll{ss}{se}];
                for rr=runsUse
                    
                    fn = sprintf( '%s%s/%s_Behav/%s_%s_sess01_run%02.f_task%i.mat' , root, sn, sn, sn,expFileName, rr, cc );
                    if ~exist(fn)
                        dbstop
                    end
                    
                    load(fn)
                    dirAll = [dirAll;p.hasRotation(~isnan(p.hasRotation))];
                    fixChgAll = [fixChgAll;p.hasFixChange(~isnan(p.hasRotation))];
                        
                    if ~strcmp(sn,'AI141') && ~strcmp(sn,'AP141')
                        respAll = [respAll;p.resp(~isnan(p.hasRotation))];
                    else
                        fn2 = sprintf( '%s%s/%s_Behav/%s_%s_sess01_run%02.f_task%i_output.mat' , root, sn, sn, sn,expFileName, rr, cc );
                        load(fn2);
                        respAll = [respAll;r.resp(~isnan(p.hasRotation))];
                    end
                end
                
                %% load trialData files
                
                trialData_fn = sprintf('%sIEMdepth_trialData/%s%s_%s_%s_%s.mat',root,subj{ss},theseSessions{se},locSignStr,VOIs{vv},trialDataStr{cc});
                fprintf('Loading trial data from: %s...\n',trialData_fn);
                load(trialData_fn);

                if length(unique(scans))>length(runsUse)
                    %need to remove a few runs and relabel the rest
                    thisRunNum=runsPrevSess;
                    for rr=runsUse
                        thisRunNum=thisRunNum+1;
%                         fprintf('se==%d,vv==%d,cc==%d,rr==%d,thisRunNum==%d\n',se,vv,cc,rr,thisRunNum);
   
                        aAll=cat(1,aAll,a_resp(scans==rr,:));
                        bAll=cat(1,bAll,b_resp(scans==rr,:));
                        stimLocsAll=cat(1,stimLocsAll,stimLocs(scans==rr,:));
                        condAll=cat(1,condAll,repmat(cc,sum(scans==rr),1));
                        runsAll=cat(1,runsAll,repmat(thisRunNum,sum(scans==rr),1));
                        
                        bXAll=cat(1,bXAll,betasX(runlabsX==rr,:));
                        predlabsXAll=cat(1,predlabsXAll,predlabsX(runlabsX==rr,:));
                        condlabsXAll=cat(1,condlabsXAll,repmat(cc,sum(runlabsX==rr),1));
                        runlabsXAll=cat(1,runlabsXAll,repmat(thisRunNum,sum(runlabsX==rr),1));

                        bZAll=cat(1,bZAll,betasZ(runlabsZ==rr,:));
                        predlabsZAll=cat(1,predlabsZAll,predlabsZ(runlabsZ==rr,:));
                        condlabsZAll=cat(1,condlabsZAll,repmat(cc,sum(runlabsZ==rr),1));
                        runlabsZAll=cat(1,runlabsZAll,repmat(thisRunNum,sum(runlabsZ==rr),1));

                        
                    end
                else

                    
                    if se>1 && size(aAll,2)~=size(a_resp,2)
                        error('size')
                    end
                    
                    %use all the runs
                    aAll=cat(1,aAll,a_resp);
                    
                    
                    bAll=cat(1,bAll,b_resp);
                    stimLocsAll=cat(1,stimLocsAll,stimLocs);
                    condAll=cat(1,condAll,repmat(cc,size(a_resp,1),1));
                    runsAll=cat(1,runsAll,scans+runsPrevSess);
                    
                    bXAll=cat(1,bXAll,betasX);
                    predlabsXAll=cat(1,predlabsXAll,predlabsX);
                    condlabsXAll=cat(1,condlabsXAll,repmat(cc,size(betasX,1),1));
                    runlabsXAll=cat(1,runlabsXAll,runlabsX+runsPrevSess);
                    
                    bZAll=cat(1,bZAll,betasZ);
                    predlabsZAll=cat(1,predlabsZAll,predlabsZ);
                    condlabsZAll=cat(1,condlabsZAll,repmat(cc,size(betasZ,1),1));
                    runlabsZAll=cat(1,runlabsZAll,runlabsZ+runsPrevSess);
                    
                end
            end 
            runsPrevSess=max(runsAll);
        end  % end session loop
        
        %put all arrays for this VOI into the correct cell
        trialDataAll(vv).aAll=aAll;
        trialDataAll(vv).bAll=bAll;
        trialDataAll(vv).stimLocsAll=stimLocsAll;
        trialDataAll(vv).condAll=condAll;
        trialDataAll(vv).runsAll=runsAll;
        trialDataAll(vv).spinDir = dirAll;
        trialDataAll(vv).resp = respAll;
        trialDataAll(vv).fixChgDir = fixChgAll;
        
        trialDataXPred(vv).betas=bXAll;
        trialDataXPred(vv).predlabs=predlabsXAll;
        trialDataXPred(vv).condlabs=condlabsXAll;
        trialDataXPred(vv).runlabs=runlabsXAll;
        
        trialDataZPred(vv).betas=bZAll;
        trialDataZPred(vv).predlabs=predlabsZAll;
        trialDataZPred(vv).condlabs=condlabsZAll;
        trialDataZPred(vv).runlabs=runlabsZAll;

        %make sure you have the right number of predictors...
        if size(bAll,1)~=3*size(bXAll,1) || size(aAll,1)~=6*size(bZAll,1)
            error('mismatch in trial numbers after merging')
        end
        
    end  % end VOI loop
      
    %save the new file for this subject 
    trialDataAll_fn = sprintf('%sIEMdepth_trialData/%s_allROIs_%s%s', root,subj{ss},locSignStr,sessStr);
    fprintf('Saving %d runs over %d sessions for %s to %s\n\n', runsPrevSess, length(theseSessions),subj{ss}, trialDataAll_fn);
    save(trialDataAll_fn','trialDataAll','trialDataXPred','trialDataZPred');
    
    clear trialDataAll trialDataXPred trialDataZPred
    
end
        


