function IEMdepth_plotSlopeDist_Z6()

%plot confidence intervals of the slopes/intercept of decoding performance versus disparity
%difference in each ROI, compare them to zero

% MMH 10/3/17

%% define subjects and flags for what to do


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

%% set up file info, other params

root='/usr/local/serenceslab/maggie/IEMdepth/';
folder='IEMdepth_classif';

saveFig = 1;
figFolder='IEMdepth_figs';
ext='epsc';

sigLevels=[0.05,0.01];

horspacer=0.147;
verspacerbig = 0.01;
verspacersmall = 0.001;
markersize = 3;

    % arrays to store acc and d' 
    accs_allsub=nan(nVOIs,nSubj,nDim);
    accsrand_allsub=nan(nVOIs,nSubj,nDim,nIter);
    d_allsub=nan(nVOIs,nSubj,nDim);
    drand_allsub=nan(nVOIs,nSubj,nDim,nIter);
     
    % p vals for significance of decoding in each condition alone
    pValsAcc_allsub=nan(nVOIs,nDim);
    pValsD_allsub = nan(nVOIs,nDim);


%% loop over subs
for ss=1:nSubj   
        
    fns=sprintf('%s%s/%s_allROIs_%s_%s_%s_%s_%s_%s_%s.mat',...
                    root,folder,subj{ss},typestr,condStrs{conduse},voxelStr,predStr,classStr,kernelStr,subMeanStr);

    load(fns);

    for vv=1:length(vuse)
        
        accs_allsub(vv,ss,2) = mean(classStruct(vuse(vv)).accReal,1);
        accsrand_allsub(vv,ss,2,:) = mean(classStruct(vuse(vv)).accRand,1);
        
        d_allsub(vv,ss,2) =  mean(classStruct(vuse(vv)).dReal,1);
        drand_allsub(vv,ss,2,:) = mean(classStruct(vuse(vv)).dRand,1);
        
    end

end

%% load the slope data

fnout=sprintf('%s%s/allSubj_allROIs_%s_%s_slopeByDisp.mat',...
                    root,folder,typestr,voxelStr);
         
load(fnout);
                
medslopes = prctile(allslopes,50,2);
lbslopes = prctile(allslopes,2.5,2);
ubslopes = prctile(allslopes,97.5,2);

medint = prctile(allint,50,2);
lbint = prctile(allint,2.5,2);
ubint = prctile(allint,97.5,2);

%% Compare each slope value to zero, and FDR correct

fdrThresh = zeros(2,2);
fdrLabs = {'slope: q=0.05','int: q=0.05';'slope: q=0.01','int q=0.01'};

pValsSlope = min([mean(allslopes<0,2),mean(allslopes>0,2)],[],2);     

isSigSlope = zeros(nVOIs,2);

for aa=1:length(sigLevels)
    alpha=sigLevels(aa);

    [p_fdr, p_masked] = fdr(pValsSlope, alpha);
    isSigSlope(:,aa)=p_masked; 

    fdrThresh(aa,1) = p_fdr;
end

%% Compare each int value to zero, and FDR correct

pValsInt = min([mean(allint<0,2),mean(allint>0,2)],[],2);     

isSigInt = zeros(nVOIs,2);

for aa=1:length(sigLevels)
    alpha=sigLevels(aa);

    [p_fdr, p_masked] = fdr(pValsInt, alpha);
    isSigInt(:,aa)=p_masked; 

    fdrThresh(aa,2) = p_fdr;
end

%% make plot

figure;hold all;


%% slope

    subplot(2,1,1);hold all;

    title('Slope')

    ylabel('Slope');

    barvals = medslopes;
    errorvals = [medslopes-lbslopes,ubslopes-medslopes];

    colormap('jet')
    %     bar(barvals,'EdgeColor');
    errorbar((1:length(barvals))',barvals,errorvals(:,1),errorvals(:,2),'Marker','none',...
                'LineStyle','none','LineWidth',1,'Color',[0,0,0]);

    xlim([0,nVOIs+1]);
    ylim([-0.005,0.015]);
    set(gca, 'XTick', 1:nVOIs, ...
                    'XTickLabel', VOIs(vuse),'XTickLabelRotation',90);

    currxlim = get(gca,'XLim');
    line(currxlim,[0,0],'Color','k');
    % ylim([-0.01,0.02])

    astLocsEach=nan(nVOIs,2);
    for vv=1:nVOIs;

        if isSigSlope(vv,1)
            astLocsEach(vv,1)=barvals(vv)+errorvals(vv,2)+verspacersmall;
        end
        if isSigSlope(vv,2)
            astLocsEach(vv,2) = astLocsEach(vv,1);                    
        end

    end


    plot((1:nVOIs),astLocsEach(:,1),'o','Color','k','MarkerSize',markersize)
    plot((1:nVOIs),astLocsEach(:,2),'o','Color','k','MarkerFaceColor','k','MarkerSize',markersize)

%% int

    subplot(2,1,2);hold all;

    title('Intercept')

    ylabel('Intercept');

    barvals = medint;
    errorvals = [medint-lbint,ubint-medint];

    colormap('jet')
    %     bar(barvals,'EdgeColor');
    errorbar((1:length(barvals))',barvals,errorvals(:,1),errorvals(:,2),'Marker','none',...
                'LineStyle','none','LineWidth',1,'Color',[0,0,0]);
    xlim([0,nVOIs+1]);
    ylim([-.35,.25])
    set(gca, 'XTick', 1:nVOIs, ...
                    'XTickLabel', VOIs(vuse),'XTickLabelRotation',90);

    currxlim = get(gca,'XLim');
    line(currxlim,[0,0],'Color','k');
    % ylim([-0.01,0.02])

    astLocsEach=nan(nVOIs,2);
    for vv=1:nVOIs;

        if isSigInt(vv,1)
            astLocsEach(vv,1)=barvals(vv)+errorvals(vv,2)+verspacerbig;
        end
        if isSigInt(vv,2)
            astLocsEach(vv,2) = astLocsEach(vv,1);                    
        end

    end

    plot((1:nVOIs),astLocsEach(:,1),'o','Color','k','MarkerSize',markersize)
    plot((1:nVOIs),astLocsEach(:,2),'o','Color','k','MarkerFaceColor','k','MarkerSize',markersize)


    if saveFig
        fnFig = [root figFolder filesep 'SlopeInt_dprimeVsDisp_' typestr '_' voxelStr];

        fprintf('saving figure to %s...\n',fnFig);
        saveas(gcf,fnFig,ext);
    end

    %% tilefigs

if length(get(0,'Children')) > 1
    tilefigs();
end




