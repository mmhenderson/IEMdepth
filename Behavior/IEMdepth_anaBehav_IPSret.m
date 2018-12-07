% Process behavioral performance for all IPS mapping runs
% some subjects have several sessions - this script will analyze
% performance over all sessions in all subjects.

%%
clear 
close all

root = '/usr/local/serenceslab/retBV/functional_ret/';

subj = {'AI','AP','BB','BC','BD','BJ','BM','BN','BO'};
    
subsIgnore={'APIPSbig_CC'};

nSubj = length(subj);

% possible prefixes that the runs might have in the file name
expPrefix = {'wedge','IPS3'};

allaccs = [];
runlist = [];

subAccs = zeros(nSubj,1);

for ss=1:nSubj
    subName=subj{ss};
   
    sessfolders = dir([root filesep subj{ss} filesep '*IPS*']);
    
    for ff=1:length(sessfolders)
        
       
        if ~any(strcmp(subsIgnore, sessfolders(ff).name))            
            
            behavRoot = [root subj{ss} filesep sessfolders(ff).name filesep sessfolders(ff).name '_Behav'];
            if ~isdir(behavRoot)
                behavRoot = [root subj{ss} filesep sessfolders(ff).name filesep sessfolders(ff).name '_behav'];
                if ~isdir(behavRoot)
                    error('no behav folder found')
                end
            end
            
            behavFileList = cat(1,dir([behavRoot filesep '*' expPrefix{1} '*']), dir([behavRoot filesep '*' expPrefix{2} '*']));
          
            [nRuns,~]=size(behavFileList);
%                 
%                 fn2=dir([behavRoot,subName,'_' expName, '*_run*_', taskStr, '_output.mat']);
%                 [nRuns2,~]=size(fn2);
%                 
%                 if nRuns1~=nRuns2
%                     error('input and output file structures do not match')
%                 end

            fprintf('Session %s: found %.f runs\n',sessfolders(ff).name,nRuns);

            for rr=1:nRuns

                thisfn = [behavRoot filesep behavFileList(rr).name];
%                     fprintf('loading %s...\n',thisfn);
                load(thisfn);
                    
                % for some runs - the script was calculating accuracy wrong
                % (multiplying by 10 instead of 100). Adjust these now
                if p.acc<10
                    p.acc = p.acc*10;
                end
                
                allaccs = [allaccs;p.acc];
                
                runlist = cat(1,runlist,[ss,ff,rr]);

            end


        end
%         end
    end       
end

for ss=1:nSubj

    subAccs(ss,1) = mean(allaccs(runlist(:,1)==ss));

end

fprintf('mean accuracy on contrast task: %.2f +/- %.2f\n',mean(subAccs(:,1)),std(subAccs(:,1)));
