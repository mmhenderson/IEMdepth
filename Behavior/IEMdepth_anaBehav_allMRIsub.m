function IEMdepth_plotBehav_allMRIsub(subName,makePlot)
% process behavior for scanner subjects 

%%
clear 
close all

root = '/usr/local/serenceslab/maggie/IEMdepth/';

subj = {'AI','AP','BB','BC','BD','BJ','BM','BN','BO'};
    
subsIgnore = {'AI131','AI142','AP142','BDIPSbig'};
% VOIs={'V1','V2','V3','V4','V3A','V3B','IPS0','IPS1','IPS2','IPS3','LO1','LO2'};

nSubj = length(subj);

expName = 'IEMdepth';
expPrefix = {'Stimulus','Fixation'};
nCond = length(expPrefix);

allaccs = [];
runlist = [];

allRT = [];
runlistRT = [];

dlims = [-1,5];
acclims = [0,100];

alld = [];

taskAccs = zeros(nSubj,2);

for ss=1:nSubj
    subName=subj{ss};
   
    sessfolders = dir([root filesep subj{ss} '*']);
    
    for ff=1:length(sessfolders)
        if ~any(strcmp(subsIgnore, sessfolders(ff).name))            
  
            for tt=1:length(expPrefix)

                subName = sessfolders(ff).name;
                
                taskStr = sprintf('task%01.f',tt);    

                behavRoot=[root subName filesep subName '_Behav' filesep];

                %count the runs, then loop over them
                fn1=dir([behavRoot,subName,'_' expName, '*_run*_', taskStr, '.mat']);
                [nRuns1,~]=size(fn1);
%                 
%                 fn2=dir([behavRoot,subName,'_' expName, '*_run*_', taskStr, '_output.mat']);
%                 [nRuns2,~]=size(fn2);
%                 
%                 if nRuns1~=nRuns2
%                     error('input and output file structures do not match')
%                 end

                fprintf('Subj %s: found %.f runs for %s\n',subName,nRuns1,expPrefix{tt});

                for rr=1:nRuns1

                    thisfn = [behavRoot fn1(rr).name];
%                     fprintf('loading %s...\n',thisfn);
                    load(thisfn);
                    
                    if strcmp(subName,'AI141') || strcmp(subName,'AP141')
                        thisfn2 = thisfn(1:end-4);
                        thisfn2 = [thisfn2 '_output.mat'];
                        load(thisfn2);
                        
                        allaccs = cat(1,allaccs,r.accuracy);

                        predlabs = r.resp;
                        thisRT = r.responseTime;
                    else
                        allaccs= cat(1,allaccs,p.accuracy);
                        predlabs = p.resp;
                        thisRT = p.responseTime;
                    end

                    
                    %get d' for this run - responses on all trials after the first,
                    %where a response was made
                    
                    if tt==1
                        reallabs = p.hasRotation;
                    else
                        reallabs = p.hasFixChange;
                    end
                    
                    if any(isnan(predlabs))
                        inds = ~isnan(predlabs);
                        predlabs = predlabs(inds);
                        reallabs = reallabs(inds);
                    end

                    thisd = get_dprime(predlabs,reallabs,unique(reallabs));
                    alld = cat(1,alld,thisd);

                    
                    thisRT = thisRT(~isnan(thisRT));
                    allRT = cat(1,allRT,thisRT);
                    runlistRT= cat(1,runlistRT,repmat([ss,tt,ff,rr],length(thisRT),1));
                    
                    runlist = cat(1,runlist,[ss,tt,ff,rr]);

                end

              
            end
        end
    end       
end

for ss=1:nSubj
    for tt=1:length(expPrefix)
        taskAccs(ss,tt) = mean(allaccs(runlist(:,1)==ss & runlist(:,2)==tt));
    end
end

fprintf('mean accuracy on target task: %.2f +/- %.2f\n',mean(taskAccs(:,1)),std(taskAccs(:,1)));

fprintf('mean accuracy on fixation task: %.2f +/- %.2f\n',mean(taskAccs(:,2)),std(taskAccs(:,2)));