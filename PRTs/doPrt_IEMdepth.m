function doPrt_IEMdepth(Subjects)
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
% TR = 2.00;
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
% params.prtr = TR*1000;%      TR (in ms, 2000)
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

% params.nvol = 150;%ceil(p.expDur/p.TR) - p.nTRsWait;

for s=1:length(Subjects)
    disp('doPrt_IEMdepth.m');
    sn=char(Subjects{s});
    
    %searchStr = [root, sn, '/', sn, '_Behav', '/', expName, '_',sn, '*_blk*.mat']
%     fn=dir([root, sn, '/', sn, '_Behav', '/',sn,'_' expFileName, '_', '*task1_output.mat'])
    fn=dir([root, sn, '/', sn, '_Behav', '/',sn,'_' expFileName, '_', '*task1.mat']);
    [runs,c]=size(fn);
    
    nRuns = runs;
    
      
    
    fprintf('Subj %s: %.f full runs found\n',Subjects{s},nRuns);
    
    %rOrder=sortBehavFiles(fn, 'sess1_blkNum');
    
    
    
    
    for rr=1:nRuns
        
        
        
        % for each full run, first look at all subruns
        for sr = 1:length(expSuffix)
%             fn2 = sprintf( '%s%s/%s_Behav/%s_%s_sess%02.f_run%02.f_task%i_output.mat' , root, sn, sn, sn,expFileName, str2num(sn(end)), rr, sr );
%             fn3 = sprintf( '%s%s/%s_Behav/%s_%s_sess%02.f_run%02.f_task%i.mat' , root, sn, sn, sn,expFileName, str2num(sn(end)), rr, sr );
            fn3 = sprintf( '%s%s/%s_Behav/%s_%s_sess01_run%02.f_task%i.mat' , root, sn, sn, sn,expFileName, rr, sr );
            disp(sprintf('loading %s...',fn3));
           
%             load(fn2);
            load(fn3);
             fprintf('acc on %d task %d: %.2f\n',rr,sr,p.accuracy)
            allgrid = [];
            for ii = 1:length(p.grid)
                allgrid = [allgrid; p.grid{ii}.'];
            end
                
            
            
%            numVols = nTRs-nDummyTRs;
            nullTrials = isnan(p.stimStart);
            nTotalTrials = sum(~isnan(p.stimStart));
            
            
            if rr == 1 && sr == 1
            
                % save a AI51_wmdelay_params.mat file in /AI51
                disp('saving wmdelay params...');
                
                IEMdepth_params.rad = p.spheresize;
                IEMdepth_params.grid = allgrid;
                %IEMdepth_params.jitterRadDeg = p.jitterRadDeg;
                %wmDrop_params.usedScreenSizeDeg = p.usedScreenSizeDeg;
                
                fn = sprintf('%s%s/%s_IEMdepth_params.mat',root,sn,sn);
                save(fn,'IEMdepth_params');
                clear IEMdepth_params fn;
            end
            
            if strcmpi(Subjects{s},'AI141') && rr == 5 && sr == 1
                p.startExp = p.startExp-2;
                TR=2.000;
%             elseif strcmpi(Subjects{s},'BM141') && rr==2 && sr==1
%                 % task 1 run 2: cut off the first 2 volumes because of high
%                 % signal 
%                 TR=2.000;
%                 p.startExp=p.startExp+2*TR;
%                
            elseif strcmp(sn,'BJ141')
                TR=2.250;
                params.nvol=134;
                            
            else
                TR=2.000;
                params.nvol=150;
%                 params.prtr=TR*1000;
            end
                
            
            
            
            params.prtr=TR*1000;
            
            %fn2 = dir([root, sn, '/', sn, '_Behav', '/', expFileName, '_',sn,'*run0' num2str(r) '_subsamp*.mat'])
            
            %curRun=rOrder(r);
            %fprintf('Crunching Run:\t%d\n', r)
            d = [];
            %load([root, sn, '/', sn, '_Behav', '/', fn(r).name])
            
            
            numConds = sum(nullTrials~=1);
            params.cond = 1:numConds;
            
            idx = find(nullTrials~=1);
            
            
            %stim.correct = p.correct(p.wmResp(idx)-1);
            %stim.correct = p.correct(idx);
            stim.resp = p.resp(idx);
            %stim.correct(~isnan(stim.correct)) = ~stim.correct(~isnan(stim.correct));
            
            % THIS CORRECTS FROM SCR COORDS TO VIS FIELD COORDS
            % in dva, cartesian axes (quadrant 1 is +x,+y)
            % 0,0 is at origin
            
            
            
            stim.stimLocs = p.stimLocs(idx,:);
            %
            %         stim.xLocDeg    = (p.stimLocsX(idx) - p.nLoc/2 - .5 +p.staggered*.25) * p.radDeg; % in DVA
            %         stim.yLocDeg    = (p.stimLocsY(idx) - p.nLoc/2 - .5 +p.staggered*.25) * p.radDeg;
            %         % location indices ( 1,1 at upper-left, nLoc,nLoc at bottom-right)
            %         stim.xLocIdx    = p.stimLocsX(idx);
            %         stim.yLocIdx    = p.stimLocsY(idx);
            %         %stim.staggered  = p.staggered * ones(size(idx));
            %         stim.cond       = p.cond * ones(size(idx));
            %
            
         
            
            

            p.stimStart = p.stimStart-p.startExp;
            p.stimEnd = p.stimEnd-p.startExp;
            
            %pToSave.radDeg = p.radDeg; % radius around which points are drawn, w/ jitter
            %pToSave.usedScreenSizeDeg = p.usedScreenSizeDeg;
            %pToSave.jitterRadDeg = p.jitterRadDeg;
            %pToSave.ppd = p.ppd;
            
            %pToSave.conditions = p.conditions;
            
            %st = cell(10,1);
            %et = cell(10,1);
            st_stim = [];
            et_stim = [];
            labels = {};
            
            
            
            % do the assignment here, and not that i'm making the assignment of
            % the times to a structure - this will make it easier to write out
            % later, particularly if there an uneven number of events in a
            % condition
            for i = 1:length(idx)
                if vOrMS
                    evt_times(i).start = round( (p.delayStart(idx(i)) )./TR)+1;
                    evt_times(i).end = round( (p.delayEnd(idx(i)) )./TR)+1;
                else
                    evt_times(i).start = round( (p.stimStart(idx(i))) * 1000);
                    evt_times(i).end = round( (p.stimEnd(idx(i))) * 1000);
                end
                
                labels{i} = sprintf('Trial_%02.f',i);
            end
            
            % make the output directory, if it doesn't exist
            
            fn_root =  [expPrefix, expSuffix{sr} num2str(rr)];
            
            runDir=[root, sn, '/' fn_root];% expPrefix, num2str(r), num2str(sr)];
            if ~exist(runDir)
                mkdir(runDir);
            end
            
            prtName = '';
            prtName = [runDir, '/', sn, '_',fn_root '.prt'];
            
            %write the prt hdr info
            col = hsv(numConds);
            writePRT(prtName, expFileName, numConds, vOrMS, labels, evt_times, col,params);
            
            % save out the mat file...
            matName = [runDir, '/', sn, '_', fn_root '.mat'];
            %clear p;
            %p = pToSave;
            save(matName,'stim');
            clear r p;
        end
        
    end
end

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