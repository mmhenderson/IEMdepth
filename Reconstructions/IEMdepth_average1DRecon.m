function IEMdepth_average1DRecon(subj,varargin)
% load the channel responses & averages reconstructions that are at the
% same locations
% VAV 2/23/2016

if nargin < 1    
    subj = {'AI','AP','BB','BC','BD','BJ','BM','BN','BO'};
    VOIs = {'V1','V2','V3','V4','V3A','V3B','IPS0','LO1','LO2'};
elseif nargin == 1
    VOIs = {'V1','V2','V3','V4','V3A','V3B','IPS0','LO1','LO2'};
end

%%
% This is the directory where the folders "IEMdepth_trialData" and
% "IEMdepth_chanResp" are located.
root = '/usr/local/serenceslab/maggie/IEMdepth/';

sessi = 1;
sepconds = 1;
absStat = 1;
usestereox = 1;

sessStrs = {'_allGoodRuns','_allSess'};
sessStr = sessStrs{sessi+1};

sepcondstrs = {'_trnAllCond','_trnSepCond','_trnFixTstStim'};
sepcondstr = sepcondstrs{sepconds+1};

locSignStrs = {'posVoxOnly','allVoxAbs'};
locSignStr = locSignStrs{absStat+1};

stereostrs = {'_screenx','_stereox'};
stereostr = stereostrs{usestereox+1};


%%
for ss = 1:length(subj)
    for vv = 1:length(VOIs)
        
        fprintf('NaNs for %s %s: ',subj{ss},VOIs{vv});
        
        % BM_allVoxAbs_allDors_channelResp_1D_best100000ZVox_trnFixTstStim_allSess_stereox.mat
        % load the 1D recon file
        fn = sprintf('%sIEMdepth_chanResp/%s_%s_%s_channelResp_1D%s%s%s.mat',...
            root,subj{ss},locSignStr,VOIs{vv},sepcondstr,sessStr,stereostr);
        load(fn);
        
        % now find the unique stimulus positions
        % first get rid of y position, if it exists
        if size(recon1D(1).stimLocsAll,2) == 3
            recon1D(1).stimLocsAll(:,2) = [];
        end
        % want to average mirrored stimulus setups, so flip half these
        % stimuli about the x axis
        
        % just copied this code from MMH IEMdepth_plotChannelResp
        % makes the x,z position grid
        gridLocs = round(gridzigz2(6,6,1)*1.5,4);
        % then find all the locations we're averaging to
        allx = unique(gridLocs(:,1));
        uniqueLocs{1} = allx(1:2:end);    % only take every other x
        uniqueLocs{2} = unique(gridLocs(:,2));
        
        for bdim = 1:2
            for cc = 1:2
%                 for cc = 1  % VAV 5/30/17 - hacky!! this is for trn fix tst stim
                % turn the chan responses into recons
                % first get rid of the NaNs -- these happen because
                % the test trials on each run are not the right
                % condition
                ci = recon1D(1).condAll == cc;
                reconAll = squeeze(recon1D(bdim).chanRespAll(cc,ci,:)) * recon1D(bdim).basis_set';
                % get stimulus positions
                slocs = round(recon1D(bdim).stimLocsAll(ci,:),4);
                
                % now get rid of the NaNs
                nrow = any(isnan(reconAll),2);
                fprintf('%i ', sum(nrow));
                reconAll(nrow,:) = [];
                slocs(nrow,:) = [];
                
                ar = nan(size(recon1D(bdim).basis_set,2),size(reconAll,2));
                
                if isempty(reconAll)
                    fprintf('%s %s dim %i cond %i has no recons!\n',...
                        subj{ss},VOIs{vv},bdim,cc);
                    continue;
                end
                
                % now average across chanResp for same locs
                if bdim == 1
                    for ii = 1:size(uniqueLocs{bdim},1)
                        thisi = ismember(slocs(:,1),uniqueLocs{bdim}(ii,:));
                        flipi = ismember(slocs(:,1),uniqueLocs{bdim}(ii,:).*[-1 1]);

                        fliprecs = flipdim(reconAll(flipi,:),2);
                        avgthis = cat(1,fliprecs,reconAll(thisi,:));

                        % need to nanmean because when averaging per condition
                        % we left a bunch of NaNs where the trial was not in
                        % that condition
                        ar(ii,:) = nanmean(avgthis);
                        clear thisi flipi fliprecs avgthis
                    end
                else
                    for ii = 1:size(uniqueLocs{bdim},1)
                        thisi = ismember(slocs(:,2),uniqueLocs{bdim}(ii,:));
                        ar(ii,:) = nanmean(reconAll(thisi,:));
                    end
                end
                avgRecons(:,:,cc) = ar;
            end % end condition loop
            
            if ~exist('avgRecons','var')
                continue;
            else
                avgRec_1D(bdim).avgRecons = avgRecons;
                avgRec_1D(bdim).stimLocs = uniqueLocs{bdim};
                avgRec_1D(bdim).res = recon1D.res;                
                avgRec_1D(bdim).basis_set = recon1D.basis_set;
                avgRec_1D(bdim).basis_grid = recon1D.basis_grid;
                clear avgRecons
            end
            
        end     % end 1D loop           

        fns = sprintf('%sIEMdepth_chanResp/%s_%s_avgRecons1D_%s%s%s%s.mat',...
            root,subj{ss},VOIs{vv},locSignStr,...
            sepcondstr,sessStr,stereostr);
        
        if ~exist('avgRec_1D','var')
            fprintf('No data. Not saving %s...\n', fns);
        	continue;
        else
            save(fns,'avgRec_1D');
            disp(['Saved ' subj{ss} ' ' VOIs{vv} ' to ' fns]);
        end
        
    end % end VOI loop
end % end subj loop