function IEMdepth_getChannelResp_1D(subj,varargin)
% compute channel responses on each trial based on voxel activation
% patterns in each ROI

% 1D version of IEMdepth_getChannelResp

% implements IEMdepth_makeX_1D to get design matrix

% loads data for each subject from IEMdepth_trialData, saves to
% IEMdepth_chanResp

% 12/15/2015, adapted from TCS by MMH 


if nargin < 1    
    subj = {'AP','AI','BA'};    
    VOIs = {'V1','V2','V3','V3A','V4','IPS0','IPS1','IPS2','IPS3'};
elseif nargin == 1
    VOIs = {'V1','V2','V3','V3A','V4','IPS0','IPS1','IPS2','IPS3'};
end

root = '/usr/local/serenceslab/vy/IEMdepth/';
% root = 'Y:/Maggie/IEMdepth/';
%root='/usr/local/serenceslab/maggie/IEMdepth/';

%% Parameters for the basis functions

% x by z - number of basis functions, each spanning [-1 1] uniformly from
% "fixation" position/depth
grid_size = 6;

% ratio of basis function (filter) size to spacing - 1.1 is good for
% triangular grid, 1.25 for square (see Sprague & Serences, 2013, Supp
% material)
filt_scale = 1.25; % FWHM: spacing ratio need to multiply by 1/rad2fwhm(1) to convert to cos7 size constant

stimSize=0.5454;

res=251;

%whethere we want to use average or beta
usingA=1;
usingB=0;

%% Loop over subjects
for ss = 1:size(subj,2)  
    
    %load the trialData for this subject
    trialData_fn = sprintf('%sIEMdepth_trialData/%s_all.mat',root,subj{ss});
    fprintf('Loading trial data from: %s...\n',trialData_fn);
    load(trialData_fn);

    for vv = 1:length(VOIs) 
        
        %get the relevant info from trialData
        thisTrialData=trialDataAll{vv};
        stimLocsAll=thisTrialData.stimLocsAll;
        aAll=thisTrialData.aAll;
        bAll=thisTrialData.bAll;
        runsAll=thisTrialData.runsAll;
        condAll=thisTrialData.condAll;
        
        % get rid of the (now-irrelevant) 2nd dimension of stimLocs
        stimLocs = stimLocsAll(:,3);

        %% Get design matrix, normalize
        
        fprintf('generating design matrix for subject %s, area %s\n', subj{ss}, VOIs{vv});

        [trnX,basis_set,basis_grid]= IEMdepth_makeX_1D(stimLocs,stimSize,grid_size,filt_scale,res);

        % normalize design matrix to be max 1
        trnX = trnX/max(trnX(:));
        %% Get channel responses, loop over all runs to cross-validate 

        chanRespAll = nan(size(aAll,1),size(trnX,2));

        runs = unique(runsAll);

            % cross-validate across runs
            for rr=1:length(runs)
                
                % figure out which trials we're using for model estimation
                % (training) and to compute channel responses (testing)
                trnidx = runsAll~=rr;
                tstidx = ~trnidx;

                if usingA
                    r_trn = aAll(trnidx,:);
                    r_tst = aAll(tstidx,:);
                elseif usingB
                    r_trn = bAll(trnidx,:);
                    r_tst = bAll(tstidx,:);
                end
                
                % compute weight matrix (nChannels x nVoxels)
                % nTrialsTraining x nChannels \ nTrialsTraining x nVoxels
                % w is solution to trnX(trnidx,:)*w = r_trn
                
                w = trnX(trnidx,:)\r_trn;   


                %% use estimated weight matrix to compute channel responses using testing data

                fprintf('Run %d : computing channel weights and responses\n', rr);

                %inv(w*w)'*w is also nChannels x nVoxels
                %a_tst is nTrialsTest x nVoxels
                %get nChannels x nTrials
                chan_resp = inv(w*w')*w*r_tst';

                % replace the channel response mean w/ the mean across voxels
                % on that trial
                % with (:) these are averages over all the values in the matrices
                chan_resp = chan_resp - nanmean(chan_resp(:)) + nanmean(r_tst(:));

                % fill in the relevant rows of chan_resp_all w/ those computed
                % on this CV iteration
                chanRespAll(tstidx,:) = chan_resp.';

            end % end cross-validation loop
            
            fn2s = sprintf('%sIEMdepth_chanResp/%s_%s_channelResp1D_trnAllCondAllSess.mat',root,subj{ss},VOIs{vv});
            fprintf('saving to %s...\n\n',fn2s);
            save(fn2s,'chanRespAll','basis_set','grid_size','basis_grid','condAll','stimLocsAll','runsAll');

         clear chanRespAll thisTrialData stimLocsAll aAll bAll condAll runsAll basis_set basis_gridx basis_gridz 
    end % end VOI loop
end % end subj loop

end % end fcn


