% Process behavior for all scanner subjects .

% This script will load all behavioral files in .mat format, collate all
% performance metrics, print out the average accuracy on each task, and
% make plots comparing performance on the two tasks. 

% For the fixation task only, will also calculate performance as a function
% of the depth position of the stimulus. Make a plot of this quantity, and
% also run a repeated-measures ANOVA for the effect of depth position.

% Note that only the Fixation Task was used for further analysis, Stimulus
% Task is just shown here for the sake of completeness.
%% 
clear 
close all

root = '/usr/local/serenceslab/maggie/IEMdepth/';

subj = {'AI','AP','BB','BC','BD','BJ','BM','BN','BO'};

% list all the sessions and runs that we want to look at
sessionsAll={{'141','143'},...
            {'141','143','144'},...                
            {'141','142','143'},...
            {'141'},...
            {'141'},...
            {'141'},...
            {'141','142'},...
            {'141','142'},...
            {'141','142'}}; 

% list runs to process:
%this includes every run, except the AI141 runs (3,4) with the
%reversed grid scaling because the positions won't match with
%our model, and one run where AP141 coughed

goodRunsAll={{[1:2,5],[1:8]},...
            {[1:4,6:8],[1:6],[1:7]},...
            {[1:7],[1:6],[1:8]},...
            {[1:7]},...
            {[1:8]},...
            {[1:8]},...
            {[1:6],[1:6]},...
            {[1:7],[1:7]},...
            {[1:7],[1:9]}};
        
nSubj = length(subj);

expName = 'IEMdepth';
expPrefix = {'Stimulus','Fixation'};
nCond = length(expPrefix);

%% initialize some arrays

allaccs = [];
runlist = [];

allRT = [];
runlistRT = [];

dlims = [-1,5];
acclims = [0,100];

alld = [];

taskAccs = zeros(nSubj,2);


% this will store the number of times a subject is correct, incorrect, or
% fails to respond for fixation tasks, as a function of the depth position
% of the stimulus.
n_each = zeros(nSubj,6,3);

n_runs = zeros(nSubj,1);
n_sess = zeros(nSubj,1);


%% load data for each subject

for ss=1:nSubj
    subName=subj{ss};
    
    theseSessions=sessionsAll{ss};
    theseGoodRuns=goodRunsAll{ss};
   
    % loop over sessions
    for ff=1:length(theseSessions)
        
        subName = [subj{ss} theseSessions{ff}];
           
        folder = [root subName filesep subName '_Behav/'];
        if ~isdir(folder)
            error('folder for session %s not found',subName);
        end
        
        % loop over tasks
        for tt=1:length(expPrefix)
            
            taskStr = sprintf('task%01.f',tt);    

            fn1=dir([folder, subName '_' expName, '*_run*_', taskStr, '.mat']);
            [nRuns1,~]=size(fn1);

            if nRuns1==0
                error('cannot find behavioral files, check your path')
            end
            
            % loop over runs
            for rr=1:nRuns1

                if ~ismember(rr,theseGoodRuns{ff})
                    fprintf('Skipping run %d for %s (experimental error)\n',rr,subName);
                    continue
                end
                    
                thisfn = [folder filesep fn1(rr).name];

                load(thisfn);

                if strcmp(subName,'AI141') || strcmp(subName,'AP141')
                    % for this subject, behavioral .mat files had a different
                    % format. Need to load a separate file.
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


                %% get accuracy, RT, and d' for this run 
                if tt==1
                    reallabs = p.hasRotation;
                else
                    reallabs = p.hasFixChange;
                end

                predlabs_orig = predlabs;
                reallabs_orig = reallabs;
                
                % for dprime calculation - want to remove any no-response
                % trials.
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

                %% look at behavior on the fixation task, separated by the depth position of the stimulus.
                
                if tt==2
                    
                    unpos = unique(p.stimLocs(~isnan(p.hasFixChange),3));
                    nZPos = length(unpos);
                    nEach = p.nTargTrials/nZPos;
                    for zz =1:nZPos
                       
                        thispos = p.stimLocs(:,3)==unpos(zz) & ~isnan(p.hasFixChange);
                        
                        if sum(thispos)~=nEach
                            error('wrong number of trials at this z position')
                        end
                        
                        % 1 is correct, 2 is incorect, 3 is no response
                        n_each(ss,zz,1) = n_each(ss,zz,1) + sum(~isnan(predlabs_orig(thispos)) & predlabs_orig(thispos)==reallabs_orig(thispos));
                        n_each(ss,zz,2) = n_each(ss,zz,2) + sum(~isnan(predlabs_orig(thispos)) & predlabs_orig(thispos)~=reallabs_orig(thispos));
                        n_each(ss,zz,3) = n_each(ss,zz,3) + sum(isnan(predlabs_orig(thispos)));
                        
                    end
                    
                end
                clear p
                clear r
            end

              
        end
    end
end       

%% Print out the overall accuracy on each task
for ss=1:nSubj
    for tt=1:length(expPrefix)
        taskAccs(ss,tt) = mean(allaccs(runlist(:,1)==ss & runlist(:,2)==tt));
    end
end

fprintf('mean accuracy on target task: %.2f +/- %.2f\n',mean(taskAccs(:,1)),std(taskAccs(:,1)));

fprintf('mean accuracy on fixation task: %.2f +/- %.2f\n',mean(taskAccs(:,2)),std(taskAccs(:,2)));

%% plot subject d-prime

figure;hold all;

title(sprintf('D-prime'));

for su = 1:length(subj)
    
    group1 = alld(runlist(:,2)==1 & runlist(:,1)==su);
    group2 = alld(runlist(:,2)==2 & runlist(:,1)==su);

    plot([mean(group1),mean(group2)],'-o');
    set(gca, 'XLim',[.5,2.5],'XTick', 1:2, 'YLim',dlims,...
                    'XTickLabel', expPrefix,'XTickLabelRotation',90);

end

legend(subj,'Location','EastOutside');

%% plot subject accuracy

figure;hold all;

title(sprintf('Accuracy'));

for su = 1:length(subj)
    
    group1 = allaccs(runlist(:,2)==1 & runlist(:,1)==su);
    group2 = allaccs(runlist(:,2)==2 & runlist(:,1)==su);

    plot([mean(group1),mean(group2)],'-o');
    set(gca, 'XLim',[.5,2.5],'XTick', 1:2, 'YLim',acclims,...
                    'XTickLabel', expPrefix,'XTickLabelRotation',90);

end

legend(subj,'Location','EastOutside');

%% plot subject accuracy on fixation task only.

ax = figure;hold all;

cm = plasma(10);

% title(sprintf('Accuracy'));

vals = zeros(length(subj),1);

for su = 1:length(subj)
    
%     group1 = allaccs(runlist(:,2)==1 & runlist(:,1)==su);
    group2 = allaccs(runlist(:,2)==2 & runlist(:,1)==su)./100;

    vals(su) = mean(group2);
    
    scatter(1:2,[mean(group2), mean(group2)],[],cm(su,:),'filled');
    set(gca, 'XLim',[0.5,2.5],'XTick', 1, 'YLim',[0.5,1],...
                    'XTickLabel', expPrefix{2},'XTickLabelRotation',90);
%     alpha(0.25);

end

errorbar(1,mean(vals),std(vals)/sqrt(length(subj)),'Color','k','LineWidth',2);
ylabel('Proportion Accuracy')
legend(subj,'Location','EastOutside');

saveas(gcf,'BehavFixTaskAvg','epsc')

%% Fixation task only: plot the accuracy as a function of depth position.

ax = figure;hold all;

cm = plasma(10);

pct_each = n_each./repmat(sum(n_each,3),1,1,3);

for ss=1:nSubj
   plot(1:6, squeeze(pct_each(ss,:,1)),'Color',[cm(ss,:),0.25],'LineWidth',2);hold on
%    alpha(gca,0.25)
end

meanvals = mean(pct_each(:,:,1),1);
semvals = std(pct_each(:,:,1),[],1)./sqrt(nSubj-1);
errorbar(1:6, meanvals,semvals,'Color','k','LineWidth',2)
% title('Percent correct (+/- SEM)')
% legend({'Correct','Incorrect','No Response'},'Location','EastOutside')
set(gca, 'XTick',1:6, 'XTickLabels',unpos)
xlabel('Disparity (arcmin)')
ylabel('Proportion Accuracy')
xlim([0,7])
ylim([.50,1])

saveas(gcf,'BehavByDepthPos','epsc')

%% Fixation task only: RM anova with factor of depth position
ii=0;
data = squeeze(pct_each(:,:,1));

% Create a table storing the respones
varNames = {'Y1','Y2','Y3','Y4','Y5','Y6'};

t = array2table(data,'VariableNames',varNames);
% Create a table reflecting the within subject factors 'TestCond', 'Attention', and 'TMS' and their levels
factorNames = {'DepthPosition'};

strs = cell(6,1);
for ii=1:6
    strs{ii} = sprintf('Position %d', ii);
end

within = table(strs,'VariableNames',factorNames);


% fit the repeated measures model
rm = fitrm(t,'Y1-Y6~1','WithinDesign',within);


mauchly_tbl = mauchly(rm);


% run my repeated measures anova here
[ranovatbl] = ranova(rm, 'WithinModel','DepthPosition')
