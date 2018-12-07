function doPrt_IEMdepth_ZPred(Subjects)
% makes PRT with only valid trials, sorted by trial numbep.  A
% corresponding matlab struct will contain all relevant info about trial.

% starttime corresponds to beginning of delay period, endtime corresponds
% to end of delay period

% made for AI51, pilot of wm_delay_split.m on 11/25/12
% TCS, 11/25/12

Subjects={'BM141','BM142'};


root = '/usr/local/serenceslab/maggie/IEMdepth/';
expFileName = ['IEMdepth'];
expPrefix = 'Depth';  % NEW: takes form wmdelay%1.f%1.f - run, subrun

expSuffix = {'Stim','Fix'};

%nLoc = 8; % number of horizontal and vertical locations


%total_trials = 96;
%num_subruns = 3;

%numConds = total_trials/num_subruns;

%numConds = 96/4nLoc*nLoc; % 36 relevant trials
vOrMS = 0;
TR = 2.00;
%nTRs = 160;
targetDur = 0.5; % s, how long target is shown before stimulus chex
nDummyTRs = 0;

% parameters (for SDM)
%params.hshape  = {'twogamma'};
modelFIR = 0;
params.hpttp  = 5; %    HRF time to peak for positive response (5)
params.hnttp = 15; %    HRF5 time to peak for negative response (15)
params.hpnr = 6; %      HRF positive/negative ratio (6)
params.hons = 0; %      HRF onset (0)
params.hpdsp = 1;%     HRF positive response dispersion (1)
params.hndsp = 1;%     HRF negative response dispersion (1)
params.ndcreg = 10;%    deconvolution number of lags (regressors, 8)
params.nderiv = [];%    add derivatives to HRF (1xN double, [])
%params.nvol = numVols;  %      number of volumes (data points, 200)
params.ortho = 0;%     if given and true, orthogonalize derivatives
params.params = [];%    1xN struct for parametrical designs (empty)
%params.cond = 1:numConds;%    condition name or number
params.name = [];%    parameter name
%params.opts = 0;%    optional settings for parameter
%params.comp = '';%  compress parameter, 'log', 'sqrt'
params.norm = 0;%  z-normalization of parameter (otherwise mean removed)
params.ortho = 0;% normalize convolved regressor against design
%params.pval = 1;%    parameter values (must match number of onsets)
params.prtr = TR*1000;%      TR (in ms, 2000)
params.rcond = 0;%     Conditions to remove (rest, 1)
% params.regi      VxR regressors of interest or RTC (not orthogonalized)
% params.regni     VxR regressors of no interest or RTC (orthogonalized)
params.rnorm = 1;%     Normalization (default: normalize HRF plateau to 1)
params.tshift = 0;%    Temporally shift convoluted regressor (in ms, 0)
if modelFIR
    params.type = 'fir';%      HRF or FIR model ({'hrf'}, 'fir')
else
    params.type = 'hrf';%      HRF or FIR model ({'hrf'}, 'fir')
end

params.nvol = 150;%ceil(p.expDur/p.TR) - p.nTRsWait;

for s=1:length(Subjects)
    
    %% find the number of runs 
    sn=char(Subjects{s});
    
    fn=dir([root, sn, '/', sn, '_Behav', '/',sn,'_' expFileName, '_', '*task1.mat']);
    [runs,c]=size(fn);
    
    nRuns = runs;
 
    fprintf('Subj %s: %.f full runs found\n',Subjects{s},nRuns);
    
    %% loop over runs
    for rr=1:nRuns        
        
        %% loop over subruns
        for sr = 1:length(expSuffix)
            fn3 = sprintf( '%s%s/%s_Behav/%s_%s_sess01_run%02.f_task%i.mat' , root, sn, sn, sn,expFileName, rr, sr );
            load(fn3);disp(sprintf('loading %s...',fn3));
            
            %load the file "output" if it was one of our first sessions
            if strcmp(sn,'AP141') || strcmp(sn,'AI141')
                fn2 = sprintf( '%s%s/%s_Behav/%s_%s_sess%02.f_run%02.f_task%i_output.mat' , root, sn, sn, sn,expFileName, str2num(sn(end)), rr, sr );
                load(fn2);
                p.stimStart=r.stimStart;
                p.stimEnd=r.stimEnd;
                p.resp=r.resp;
                p.startExp=r.startExp;
            end
            
            %% correct a few weird runs
            
             
            if strcmpi(Subjects{s},'AI141') & rr == 5 & sr == 1
                p.startExp = p.startExp-2;
                params.nvol=150;
            elseif strcmpi(Subjects{s},'BM141') && rr==2 && sr==1
                % task 1 run 2: cut off the first 2 volumes because of high
                % signal 
                TR=2.000;
                p.startExp=p.startExp+2*TR;
                params.nvol=148;
            else
                params.nvol=150;
            end
           
            nullTrials = boolean(isnan(p.stimStart));
                    
            %idx is the variable that gets used below to define events 
            if strcmpi(Subjects{s},'BM141') & rr==1 & sr==2
                %early trigger - remove the first run
                p.startExp = p.startExp+10;            
                early = p.stimStart-p.startExp<0;               
                idx = find(~nullTrials & ~early);
                
            elseif strcmpi(Subjects{s}, 'BM141') & rr==6 & sr==1
                %goggles cut out for a couple of trials
                cutout = ~isnan(p.hasRotation) & isnan(p.resp);
                idx = find(~nullTrials & ~cutout);                
            else              
                idx = find(~nullTrials);

            end
            
            %% define the start/end times
            
            stimStart = p.stimStart-p.startExp;
            stimEnd = p.stimEnd-p.startExp;
            
            stimStart=stimStart(idx);
            stimEnd=stimEnd(idx);
            
            %% define the preds
            
            stimLocsAll=round(p.stimLocs(idx,3),4);
         
            unlocs=unique(stimLocsAll);
            %store the X position labels of each predictor
            stim.posLabs=unlocs;
            
            %store the X position labels of each predictor
            stim.posLabs=unlocs;
            
            numConds = length(unlocs);
            params.cond = 1:numConds;

            labels = {};
                               
            %% loop over the preds and assign events to each one 
            for i = 1:numConds
                
                theseinds=stimLocsAll==unlocs(i);
                
                if vOrMS
%                     evt_times(i).start = round( (p.delayStart(theseinds) )./TR)+1;
%                     evt_times(i).end = round( (p.delayEnd(theseinds) )./TR)+1;
                else
                    evt_times(i).start = round( (stimStart(theseinds)) * 1000);
                    evt_times(i).end = round( (stimEnd(theseinds)) * 1000);
                end
                
                labels{i} = sprintf('X=%.02f',unlocs(i));
            end
            
            %% save prt, sdm, mat
            
            %make the output directory, if it doesn't exist            
            fn_root =  [expPrefix, expSuffix{sr} num2str(rr)];
            
            runDir=[root, sn, '/' fn_root];% expPrefix, num2str(r), num2str(sr)];
            if ~exist(runDir)
                mkdir(runDir);
            end
            
            prtName = '';
            prtExt='_ZPred';
            prtName = [runDir, '/', sn, '_',fn_root prtExt '.prt'];
            
            %write the prt hdr info
            col = hsv(numConds);
            writePRT(prtName, expFileName, numConds, vOrMS, labels, evt_times, col,params);
            
            % save out the mat file...
            matName = [runDir, '/', sn, '_', fn_root prtExt '.mat'];
            %clear p;
            %p = pToSave;
            save(matName,'stim');
            clear r p;
        end
        
    end
end
%% plot one of the SDMs, just to visualize the way we set it up


lastSDM=[runDir, '/', sn, '_',fn_root prtExt '.sdm'];
thissdm=BVQXfile(lastSDM);
figure;plot(thissdm.SDMMatrix);
% legend(thissdm.PredictorNames)
title(sprintf('time course for each of %d predictors',size(thissdm.PredictorNames,2)));


%-------------------------------------------------------------------------
%support functions

function writePRT(prtName, expName, numConds, vOrMS, CondName, evt_times, col,params)
% use BVQXfile to write out the file
prt = BVQXfile('new:prt');  % make a new prt object, then we'll fill it up
if vOrMS
    prt.ResolutionOfTime = 'Volumes';
else
    prt.ResolutionOfTime = 'msec';
end
prt.Experiment = expName;

prt.NrOfConditions = numConds;
for i=1:numConds
    prt.Cond(i).ConditionName{1} = CondName{i};
    prt.Cond(i).NrOfOnOffsets = numel(evt_times(i).start);
    prt.Cond(i).OnOffsets = [evt_times(i).start, evt_times(i).end];
    prt.Cond(i).Weights = 1;
    prt.Cond(i).Color = [round(col(i,1)*255),round(col(i,2)*255),round(col(i,3)*255)];
end
prt.SaveAs(prtName);
sdm = prt.CreateSDM(params);
sdmName = [prtName(1:end-3) 'sdm'];
sdm.SaveAs(sdmName);

sdm.ClearObject;
prt.ClearObject;
return