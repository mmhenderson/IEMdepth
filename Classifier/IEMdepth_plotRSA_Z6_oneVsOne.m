function IEMdepth_plotRSA_Z6_oneVsOne()

% make table of classifier performance for 6Z one vs one

% MMH 9/14/17

%% define subjects and flags for what to do
close all

addpath(genpath('/usr/local/serenceslab/maggie/mFiles/Plotting tools/'))

subj = {'AI','AP','BB','BC','BD','BJ','BM','BN','BO'};
vuse=[1:7,11:12];
VOIs={'V1','V2','V3','V4','V3A','V3B','IPS0','IPS1','IPS2','IPS3','LO1','LO2'};
for ss=1:length(subj)
    subj{ss} = [subj{ss} '_allROIs'];
end
nSubj=length(subj);
nVOIs=length(vuse);

saveFig = 1;
figFolder='IEMdepth_figs';
ext='epsc';

typestr = 'Z2_oneVsOne';

%parameters for the classifier
kernelStr='linear';
% functStr='svmtrain (-t 0 -q)';
classStr='svmtrain';
% subMeanStr = 'subMean2';
subMeanStr = 'noSubMean';
usingA=1;
predStrs={'predB','predA'};
predStr=predStrs{usingA+1};
voxelStr = 'allVox';
% voxelStr = 'take150ZVox';
condStrs = {'trainStim','trainFixat'};
conduse = 2;

plotConf = 0;
plotD = 1;

sigLevels = [0.05,0.01];

%% set up file info, other params

root='/usr/local/serenceslab/maggie/IEMdepth/';

% nIter=1000;

folder='IEMdepth_classif';
 
groupStrs = {'38.6','25.6','11.1','-5.2','-23.6','-44.6'};
% zPosList = linspace(-9.8,9.8,6);

bigMat = zeros(nSubj,length(vuse),6,6);

confMat = zeros(nSubj,length(vuse),6,6);
nPoss = zeros(6,6);

zDispList = [38.592,25.612,11.124,-5.152,-23.567,-44.573];

posPairList = combnk(1:6,2);
dispPairList = reshape(zDispList(posPairList(:)),15,2);
dispDiffList = dispPairList(:,1)-dispPairList(:,2);
[distListSort,distOrder] = sort(dispDiffList,'ascend');

undist = unique(dispDiffList);
nDist = length(undist);
    
% accs_allsub = zeros(nVOIs,nSubj,nDist);
% accsrand_allsub = zeros(nVOIs,nSubj,nDist,nIter);

nIter = 1000;

dMat_allsub = zeros(nSubj,nVOIs,nDist);
dMatrand_allsub = zeros(nSubj,nVOIs,nDist,nIter);


%% loop over subs
for ss=1:nSubj   
    
    fns=sprintf('%s%s/%s_%s_%s_%s_%s_%s_%s_%s.mat',...
                    root,folder,subj{ss},typestr,condStrs{conduse},voxelStr,predStr,classStr,kernelStr,subMeanStr);

    load(fns);
     
    posPairList = classStruct(1).posPairList;
    
    for vv=1:length(vuse)
        
        
        for pp=1:length(posPairList)

            gg1 = posPairList(pp,1);
            gg2 = posPairList(pp,2);
            bigMat(ss,vv,gg1,gg2) = classStruct(vuse(vv)).dReal(pp);
            bigMat(ss,vv,gg2,gg1) = classStruct(vuse(vv)).dReal(pp);
            
            confMat(ss,vv,gg1,gg1) = confMat(ss,vv,gg1,gg1) + sum(classStruct(vuse(vv)).realLabs(pp,:)==gg1 & classStruct(vuse(vv)).predLabs(pp,:)==gg1);
            confMat(ss,vv,gg1,gg2) = confMat(ss,vv,gg1,gg2) + sum(classStruct(vuse(vv)).realLabs(pp,:)==gg1 & classStruct(vuse(vv)).predLabs(pp,:)==gg2);
            confMat(ss,vv,gg2,gg1) = confMat(ss,vv,gg2,gg1) + sum(classStruct(vuse(vv)).realLabs(pp,:)==gg2 & classStruct(vuse(vv)).predLabs(pp,:)==gg1);
            confMat(ss,vv,gg2,gg2) = confMat(ss,vv,gg2,gg2) + sum(classStruct(vuse(vv)).realLabs(pp,:)==gg2 & classStruct(vuse(vv)).predLabs(pp,:)==gg2);

            
            if ss==1 && vv==1
                nTot = size(classStruct(vuse(vv)).realLabs,2);
                nPoss(gg1,gg1) = nPoss(gg1,gg1) + nTot;
                nPoss(gg1,gg2) = nPoss(gg1,gg2) + nTot;
                nPoss(gg2,gg1) = nPoss(gg2,gg1) + nTot;
                nPoss(gg2,gg2) = nPoss(gg2,gg2) + nTot;
                
            end
            
                        
%             theseaccs = classStruct(vuse(vv)).accReal(distOrder(dd),:);
%             accs_allsub(vv,ss,dd) = mean(theseaccs);
%             
%             theseaccs = classStruct(vuse(vv)).accRand(distOrder(dd),:);
%             accsrand_allsub(vv,ss,dd,:) = mean(theseaccs,1);
            
            thesed = classStruct(vuse(vv)).dReal(pp,:);
            dMat_allsub(ss,vv,pp) = thesed;
            
            thesed = classStruct(vuse(vv)).dRand(pp,:);
            dMatrand_allsub(ss,vv,pp,:) = thesed;
            
        end
       
        
    end
   
end


%% do significance test

pValsD_allsub = zeros(nVOIs,nDist);

for vv=1:nVOIs
    for dd=1:nDist

        %% all trials
        realD = squeeze(dMat_allsub(:,vv,dd));
        nullD = squeeze(dMatrand_allsub(:,vv,dd,:));

        pValsD_allsub(vv,dd) = 2*min([mean(mean(realD)<mean(nullD)),mean(mean(realD)>mean(nullD))]);


    end
end

%% FDR correct

isSigD = zeros(size(pValsD_allsub,1),size(pValsD_allsub,2),length(sigLevels));

for aa=1:length(sigLevels)

    alpha=sigLevels(aa);

    [p_fdr, p_masked] = fdr( pValsD_allsub, alpha);
    isSigD(:,:,aa)=p_masked;
    
end

%% put the significance values into matrix format

sigMat = zeros(nVOIs, 6, 6);

for vv=1:nVOIs
    
    for pp=1:nDist

        gg1 = posPairList(pp,1);
        gg2 = posPairList(pp,2);

        if isSigD(vv,pp,1)
            sigMat(vv,gg1,gg2) = 1;
%             sigMat(vv,gg2,gg1) = 1;
        end
        if isSigD(vv,pp,2)
            sigMat(vv,gg1,gg2) = 2;
%             sigMat(vv,gg2,gg1) = 2;
        end
        
    end
end
%%

if plotD

ax = [];
nGroups = 6;

figure('units','normalized','outerposition',[0,0,1,1]);
hold all;

for vv=1:nVOIs
    
    thisMat = squeeze(mean(bigMat(:,vv,:,:),1));

    for ii=1:6
        for jj=ii:6
            thisMat(ii,jj) =nan;
        end       
    end
    
    subplot(3,3,vv);hold all
   
    sanePColor([thisMat]);
    axis equal;

    xlim([.5, size(thisMat,1)+.5])
    ylim([.5, size(thisMat,1)+.5])


    set(gca,'XTick', 1:nGroups, ...
                'XTickLabel', groupStrs,'XTickLabelRotation',90);
    set(gca,'YTick', 1:nGroups, ...
        'YTickLabel', groupStrs);
    xlabel('Disparity (arcmin)');
    ylabel('Disparity (arcmin)');
    
    set(gca,'YDir','reverse')
    
    %% add significance stars

    [x,y] = ind2sub([6,6], find(squeeze(sigMat(vv,:,:))==1));    
    plot(x, y, 'ko')
    
    [x,y] = ind2sub([6,6], find(squeeze(sigMat(vv,:,:))==2));   
    plot(x, y, 'k.', 'MarkerSize',10)
    
    %%
    title(sprintf('%s: Z discriminability (d'')',VOIs{vuse(vv)}));
    colormap(viridis)
    colorbar();
    
    ax = [ax, gca];
end

match_clim(ax);

% set(gcf, 'PaperPositionMode','auto')

fnFig = [root figFolder filesep 'RSA_6Z_oneVsOne'];
        
fprintf('saving figure to %s...\n',fnFig);
% saveas(gcf,fnFig,ext); 
print(fnFig,'-depsc','-r0')

end

%%

if plotConf
    
ax = [];
nGroups = 6;

for vv=1:nVOIs
    thisMat = squeeze(mean(confMat(:,vv,:,:),1));
    thisMat = thisMat./nPoss;
%             thisMat = thisMat./(numOccs./5.*repmat(sum(thisMat,2),1,size(thisMat,2)));
%             thisMat = thisMat./repmat(sum(thisMat,2),1,size(thisMat,2));
    figure;
    hold all;
    imagesc(thisMat);
    axis equal;

    xlim([.5, size(thisMat,1)+.5])
    ylim([.5, size(thisMat,1)+.5])


    set(gca,'XTick', 1:nGroups, ...
                'XTickLabel', groupStrs,'XTickLabelRotation',90);
    set(gca,'YTick', 1:nGroups, ...
        'YTickLabel', groupStrs);
    xlabel('Classified Disparity (arcmin)');
    ylabel('Actual Disparity (arcmin)');

    title(sprintf('%s: confusion probability',VOIs{vuse(vv)}));
    colorbar();
    
    ax = [ax, gca];
end

match_clim(ax);

end
%% tilefigs

if length(get(0,'Children')) > 1
    tilefigs();
end




