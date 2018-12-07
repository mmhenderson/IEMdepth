
function IEMdepth_concatMRs()
% TCS 11/4/2015
%
% concatenates LH and RH of MR files, saves them as Bilat files

%this version MMH 8/19/16 works with output of
%IEMdepth_compileData_absStat, which calls compVOI2 (newer version)
%% set up subs, VOIs

subj = {'AI141','AI143','AP141','AP143','AP144','BB141','BB142','BB143',...
        'BC141','BD141','BJ141','BM141','BM142','BN141','BN142','BO141','BO142'};

% subj={};

% VOIs={'V1','V2','V3','V4','V3A','V3B','IPS0','IPS1','IPS2','IPS3','IPS1-3','allDors'};
VOIs = {'LO1','LO2'};

vrange = 1:length(VOIs);
%%

prefix = {'DepthFix','DepthStim'};
% prefix={'DepthStim'};

locstr = 'DepthLocalizer';

mr.absStat=1; %%take the absolute value of the localizer voxels, include negative scores too

locSignStrs={'posVoxOnly','allVoxAbs'};
locSignStr=locSignStrs{mr.absStat+1};


mrstr = 'VT1_IntInf_Norm2_RmMo1_Vol_1_';

truncOutliers = 0;
normalize = 0;
vKeep = inf;
keep = inf;
rawData = 3;
thresh = '05';

%nTRs = 14;
whichTRs_trn = [3 4];

hemis = {'LH','RH'};

root = '/usr/local/serenceslab/maggie/IEMdepth/';
%root where the prt is
% prtroot='/usr/local/serenceslab/Chaipat/fMRI/IEMdepth/';
prtroot=root;

for ss = 1:size(subj,2)   % loop over subjects
    sn=char(subj{1,ss});
    
    %loop over VOIs
    vcnt = 1;
    
    for vv=vrange
        % build a struct array of mr's - one for each cond (for each VOI)
        mrfn = cell(2,1);
        for c = 1:length(prefix)   % TRN & TEST
            
            for hh = 1:2
                
%                 myfs = [root 'IEMdepth_mrStructs/' subj{ss} '_' locstr '_' mrstr prefix{c} '_' hemis{hh} '-'  VOIs{vv} '.mat'];
%                 myf = dir(myfs);
                
                fn = sprintf('%sIEMdepth_mrStructs/%s_%s_%s_%s%s_%s-%s.mat',root,subj{ss},locstr,locSignStr,mrstr,prefix{c},hemis{hh},VOIs{vv});
                load(fn);
                
%                 load([root 'IEMdepth_mrStructs/' myf(1).name]);
%                 clear myf myfs;
                
                fprintf('Processing data from %s, %s numVox: %d\n', char(mr.voiNames{1}),[prefix{c}],size(mr.voiData(1).mod,2))
                
                % select the data type
                if rawData==1
                    tmpdata{hh} = mr.voiData(1).rd;
                    
                elseif rawData==2
                    tmpdata{hh} = mr.voiData(1).nd;
                    
                elseif rawData == 3
                    tmpdata{hh} = mr.voiData(1).mod;
                    
                end
                if hh == 1
                    clear mr;
                end
            end
            % now for RH...
            
            mr.data = [tmpdata{1} tmpdata{2}];
            
            
            mr.truncOutliers = truncOutliers;
            mr.rawData = rawData;
            mr.normalize = normalize;
            mr.thresh = thresh;
            %mr.numDirs = numDirs;
            
            
            mr.root = prtroot;
            mr.locstr = locstr;
            mr.thisVOI = ['Bilat-' VOIs{vv}];
            mrs(c) = mr;
            mr.cond = prefix{c};
            
            % save bilateral MR
            mrfn{c} = sprintf('%sIEMdepth_mrStructs/%s_%s_%s_%s%s_Bilat-%s.mat',root,subj{ss},locstr,locSignStr,mrstr,mr.cond,VOIs{vv});
            fprintf('saving bilateral mr struct to %s\n\n',mrfn{c});
            save(mrfn{c},'mr');
            
            
            %% extract timecourses from testing data
            
            % also get betas
            fprintf('using GLM to get betas...\n');
            
            [b_resp] = IEMdepth_extractSignal_betas(mrs(c));
                        
            [a_resp, scans, stimLocs] = IEMdepth_extractSignal_avg(mrs(c),whichTRs_trn);
            
            [betasX,predlabsX,runlabsX] = IEMdepth_extractSignal_betas_XPred(mrs(c));
            
            [betasZ,predlabsZ,runlabsZ] = IEMdepth_extractSignal_betas_ZPred(mrs(c));
            
            fn2s = sprintf('%sIEMdepth_trialData/%s_%s_%s_%s.mat',root,subj{ss},locSignStr,VOIs{vv},mr.cond);
            
            fprintf('Saving extracted data to: %s...\n',fn2s);
            save(fn2s,'mrfn','locstr','a_resp','scans','stimLocs','b_resp','whichTRs_trn',...
                    'betasX','predlabsX','runlabsX',...
                    'betasZ','predlabsZ','runlabsZ');

        end % end cond loop
        
        
    end
end
return;