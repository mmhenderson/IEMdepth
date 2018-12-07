
%plot distribution of decoding performance versus disparity
%difference in each ROI, compare them to zero

% MMH 10/3/17

%% define subjects and flags for what to do

addpath(genpath('/usr/local/serenceslab/maggie/mFiles/Plotting tools/Violinplot-Matlab-master/'))
rmpath(genpath('/mnt/neurocube/local/Rosanne/GeneralUseScripts/'))

subj = {'AI','AP','BB','BC','BD','BJ','BO','BN','BM'};
VOIs = {'V1','V2','V3','V4','V3A','V3B','IPS0','IPS1','IPS2','IPS3','LO1','LO2'};
nSubj = length(subj);
% VOIs={'V1','V2','V3','V4','V3A','V3B','IPS0','IPS1','IPS2','IPS3','IPS1-3','allDors'};

% vuse = 1:nVOIs;

vuse = [1:7,11:12];
nVOIs=length(vuse);

plotSig=1;

% numVoxUse=150;
% voxelStr=sprintf('take%dZVox',numVoxUse);
voxelStr = 'allVox';
classStr = 'svmtrain';
kernelStr='linear';
predStr = 'predA';
subMeanStr = 'noSubMean';
% subMeanStr = 'subMean2';
% tirange=[2];

typestr='Z2_oneVsOne';
% titlestrs='Slope of dprime versus disp';

% ymaxes=[1,1];
% drange=[-.2,0.6];
% 
% chanceVals=[1/2];

dimStrs = {'X position','Z position'};
nDim= length(dimStrs);

condStrs = {'trainStim','trainFixat'};
conduse = 2;

nIter=1000;

close all

labs = {'Across Fixation', 'Within Each Side of Fixation'};

%% set up file info, other params

root='/usr/local/serenceslab/maggie/IEMdepth/';
folder='IEMdepth_classif';

%% load the slope data

fnout=sprintf('%s%s/allSubj_allROIs_%s_%s_slopeByDisp_acrossWithin.mat',...
                    root,folder,typestr,voxelStr);
         
load(fnout);

%% make plot

for xx=1:2

figure;hold all;

title([labs{xx} ': Slope'])

ylabel('Slope');

vp = violinplot(squeeze(allslopes(:,xx,:))', VOIs,'ShowData',false, 'ShowMean',false, 'ViolinColor',[0.8,0.8,0.8])

set(gca,'XTickLabel',VOIs(vuse))

line(get(gca, 'XLim'),[0,0]);
ylim([-0.03,0.03])

end
