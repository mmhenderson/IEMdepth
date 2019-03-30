function IEMdepth_getChannelResp1D(subj,useAll,absStat,sepconds,usestereox)
% Make 1D reconstructions
% Edited VAV 11/27/2018 - take correct trialData structure

%% set up the parameters

if nargin<5
    useAll=1;
    absStat=1;
    sepconds=1;
    usestereox=1;
    subj = {'AI','AP','BB','BC','BD','BJ','BM','BN','BO'};
end

if useAll
    sessStr='_allSess';
else
    sessStr='_allGoodRuns';
end

% This is the directory where the folders "IEMdepth_trialData" and
% "IEMdepth_chanResp" are located.
root='/usr/local/serenceslab/maggie/IEMdepth/';

locSignStrs={'posVoxOnly','allVoxAbs'};
locSignStr=locSignStrs{absStat+1};

stereostrs={'_screenx','_stereox'};
stereostr=stereostrs{usestereox+1};

sepcondstrs={'_trnAllCond','_trnSepCond'};
sepcondstr=sepcondstrs{sepconds+1};

%% Parameters for the basis functions

% ratio of basis function (filter) size to spacing - 1.1 is good for
% triangular grid, 1.25 for square (see Sprague & Serences, 2013, Supp
% material)
filt_scale = 1.25; % FWHM: spacing ratio need to multiply by 1/rad2fwhm(1) to convert to cos7 size constant
stimSize = 0.5454;
basisDim = linspace(-2,2,6);
basisFWHM = (basisDim(end)-basisDim(end-1)) * filt_scale;

% whether we want to use average of timepoints or beta values from GLM
usingA = 1;   % if this is set to 0, will use the trial betas

% number of total conditions
nconds = 2;
res = 101;

% numVoxUse=Inf;
% termind=2;      % this uses the top F scores from this term in the ANOVA (x = 1, z = 2, int = 3)

%% Loop over subjects
for ss = 1:size(subj,2)  
    
    %load the trialData for this subject
    trialData_fn = sprintf('%sIEMdepth_trialData/%s_allROIs_%s%s.mat',...
        root,subj{ss},locSignStr,sessStr);
    fprintf('Loading trial data from: %s...\n',trialData_fn);
    load(trialData_fn);
        
    for vv = 1:length(VOIs)
        
        %get the relevant info from trialData
        thisTrialData=trialDataAll(vv);
        stimLocsAll=thisTrialData.stimLocsAll;
        aAll=thisTrialData.aAll;
        bAll=thisTrialData.bAll;
        runsAll=thisTrialData.runsAll;
        condAll=thisTrialData.condAll;
        
        % get rid of the (now-irrelevant) 2nd dimension of stimLocs
        stimLocsAll(:,2)= [];

        if usestereox
            %adjust all the coordinates to their percieved x positon, train
            %based on adjusted positions
            % undo the depth transformation on the stim locations
            % first sort by z location
            [zloc1d,~,zi] = unique(stimLocsAll(:,2));
            % this is the scaling factor...is this saved somewhere? I just
            % pulled it from MMH script
            sclen = fliplr(linspace(1,1.33,length(zloc1d)));
            for ii = 1:length(zloc1d)
                stimLocsAll(zi==ii,1) = stimLocsAll(zi==ii,1) / sclen(ii);
            end
            
        end
        
        %round all the stimulus location values
        stimLocsAll=round(stimLocsAll,4);
        
        % do a 1-dimensional IEM for every dimension of interest!
        for recdim = 1:2
            stimLocs = stimLocsAll(:,recdim);

            %% Get design matrix, normalize

            fprintf('generating design matrix for subject %s, area %s\n',...
                subj{ss}, VOIs{vv});

            [trnX,basis_set,basis_grid] = IEMdepth_makeX1_1D(stimLocs,...
                stimSize,basisDim,basisFWHM,res);
            %% Get channel responses, loop over all runs to cross-validate 

            if sepconds
                %train separately within each condition
                chanRespAll = nan(nconds,size(aAll,1),size(trnX,2));
                runs = unique(runsAll);
    
                for cc=1:nconds
                    % cross-validate across runs
                    for rr=1:length(runs)

                        % figure out which trials we're using for model estimation
                        % (training) and to compute channel responses (testing)
                        trnidx = runsAll~=rr & condAll==cc;
                        tstidx = runsAll==rr & condAll==cc;

                        if usingA
                            r_trn = aAll(trnidx,:);
                            r_tst = aAll(tstidx,:);
                        else
                            r_trn = bAll(trnidx,:);
                            r_tst = bAll(tstidx,:);
                        end

                        % normalize design matrix to be max 1
                        tX = trnX(trnidx,:)/max(max(trnX(trnidx,:)));
                        % compute weight matrix (nChannels x nVoxels)
                        % nTrialsTraining x nChannels \ nTrialsTraining x nVoxels
                        % w is solution to trnX(trnidx,:)*w = r_trn
                        w = tX\r_trn;   

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
                        chanRespAll(cc,tstidx,:) = chan_resp.';

                    end % end cross-validation loop
                end  %end condition loop
            else
                %train both conditions together
                chanRespAll = nan(size(aAll,1),size(trnX,2));
                runs = unique(runsAll);
                % cross-validate across runs
                for rr=1:length(runs)

                    % figure out which trials we're using for model estimation
                    % (training) and to compute channel responses (testing)
                    trnidx = runsAll~=rr;
                    tstidx = runsAll==rr;

                    if usingA
                        r_trn = aAll(trnidx,:);
                        r_tst = aAll(tstidx,:);
                    else
                        r_trn = bAll(trnidx,:);
                        r_tst = bAll(tstidx,:);
                    end

                    % normalize design matrix to be max 1
                    tX = trnX(trnidx,:)/max(max(trnX(trnidx,:)));
                    % compute weight matrix (nChannels x nVoxels)
                    % nTrialsTraining x nChannels \ nTrialsTraining x nVoxels
                    % w is solution to trnX(trnidx,:)*w = r_trn
                    w = tX\r_trn;   

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
            end
            
            recon1D(recdim).chanRespAll = chanRespAll;
            recon1D(recdim).basis_set = basis_set;
            recon1D(recdim).basis_grid = basis_grid;
            recon1D(recdim).condAll = condAll;
            % note: this is the transformed (stereo) x coord, not absolute
            % x coord
            recon1D(recdim).stimLocsAll = stimLocsAll;
            recon1D(recdim).runsAll = runsAll;
            recon1D(recdim).res = res;
            
            clear basis_set basis_grid stimLocs tX
        end
        clear chanRespAll thisTrialData stimLocsAll aAll bAll condAll runsAll
        
        %name of the file indicates how model was trained
        fn2s = sprintf('%sIEMdepth_chanResp/%s_%s_%s_channelResp_1D%s%s%s.mat',...
                root,subj{ss},locSignStr,VOIs{vv},sepcondstr,sessStr,stereostr);

        fprintf('saving to %s...\n\n',fn2s);
        save(fn2s,'recon1D');

        
    end % end VOI loop
end % end subj loop