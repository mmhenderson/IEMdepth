% Process behavior for scanner subjects on localizer task

% Localizer runs are spread across multiple sessions in some subjects, in
% other subjects they are all within one session. In this script we will
% go through all folders from all sessions, count how many .mat
% localizer files are found, and process them all. 
% Print out the final value of accuracy on this task.

% Note that for 2 runs of subject AP142, the script was exited prematurely
% so the final accuracy was not saved. The data from this run was processed
% as usual, but it does not contribute to the accuracy measured here.
%%
clear 
close all

root = '/usr/local/serenceslab/maggie/IEMdepth/';

subj = {'AI','AP','BB','BC','BD','BJ','BM','BN','BO'};
    
subsIgnore = {'AI131'};

nSubj = length(subj);

expName = 'IEMdepth';
expPrefix = {'Localizer'};

allaccs = [];
runlist = [];

subAccs = zeros(nSubj,1);

for ss=1:nSubj

    sessfolders = dir([root filesep subj{ss} '*']);
    
    for ff=1:length(sessfolders)
        if ~any(strcmp(subsIgnore, sessfolders(ff).name))            

                subName = sessfolders(ff).name;

                behavRoot=[root subName filesep subName '_Behav' filesep];

                %count the runs, then loop over them
                fn1=dir([behavRoot,subName,'_' expName, '*_LocalizerRun*.mat']);
                [nRuns1,~]=size(fn1);


                fprintf('Subj %s: found %.f runs for %s\n',subName,nRuns1,expPrefix{1});

                for rr=1:nRuns1

                    thisfn = [behavRoot fn1(rr).name];
%                     fprintf('loading %s...\n',thisfn);
                    load(thisfn);
                    
                    if strcmp(subName,'AI141') || strcmp(subName,'AP141') 
                         % for this subject, behavioral .mat files had a different
                     % format. Need to load a separate file.
                        thisfn2 = thisfn(1:end-4);
                        thisfn2 = [thisfn2 '_output.mat'];
                        load(thisfn2);
                        
                        allaccs = cat(1,allaccs,r.accuracy);
                        runlist = cat(1,runlist,[ss,ff,rr]);
                    else
                        if ~isfield(p,'accuracy')
                            fprintf('accuracy was not saved for %s\n',fn1(rr).name)
                        else
                            allaccs= cat(1,allaccs,p.accuracy);
                            runlist = cat(1,runlist,[ss,ff,rr]);
                        end
                    end

                    
                    

                end

              
%             end
        end
    end       
end

for ss=1:nSubj
%     for tt=1:length(expPrefix)
        subAccs(ss,1) = mean(allaccs(runlist(:,1)==ss));
%     end
end

fprintf('mean accuracy on localizer task: %.2f +/- %.2f\n',mean(subAccs(:,1)),std(subAccs(:,1)));
