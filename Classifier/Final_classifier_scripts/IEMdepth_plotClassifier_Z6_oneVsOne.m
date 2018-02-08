function IEMdepth_plotClassifier_Z6_oneVsOne()

% plot the performance at Z classification, for a one-vs-one pairwise
% decoding scheme
% there are 6 total positions in Z, we did each of 15 possible two-way
% classifications, then sorted the results according to the magnitude of
% the disparity difference between the compared points. 

% using all voxels, not using IPS2-3

% MMH 9/14/17

%% define subjects and flags for what to do
% close all

subj = {'AI','AP','BB','BC','BD','BJ','BM','BN','BO'};
vuse=[1:8,11:12];
VOIs={'V1','V2','V3','V4','V3A','V3B','IPS0','IPS1','IPS2','IPS3','LO1','LO2'};

nSubj=length(subj);
nVOIs=length(vuse);

saveFig = 1;
figFolder='IEMdepth_figs';
ext='epsc';

tirange=[1];

% typestrs={};
% typestrs={'Six-Z','Six-X Stereo'};
titlestrs={'Confusion of 6-way classifier'};
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

close all


%% set up file info, other params

root='/usr/local/serenceslab/maggie/IEMdepth/';

nIter=1000;
plotSig=1;
chanceVal=1/2;
sigLevels = [0.05,0.01];
plotAcc=0;
plotD=1;
ymax=1.2;
drange = [-.5,1];

folder='IEMdepth_classif';

verspacerbig = 0.03;
markersize = 3;

zDispList = [-35.8,-22.7,-8,8.5,27.2,48.4];

posPairList = combnk(1:6,2);
dispPairList = reshape(zDispList(posPairList(:)),15,2);
dispDiffList = dispPairList(:,2)-dispPairList(:,1);
[distListSort,distOrder] = sort(dispDiffList,'ascend');

undist = unique(dispDiffList);
nDist = length(undist);
    
accs_allsub = zeros(nVOIs,nSubj,nDist);
accsrand_allsub = zeros(nVOIs,nSubj,nDist,nIter);

d_allsub = zeros(nVOIs,nSubj,nDist);
drand_allsub = zeros(nVOIs,nSubj,nDist,nIter);


%% loop over subs
for ss=1:nSubj   
    
    fns=sprintf('%s%s/%s_allROIs_%s_%s_%s_%s_%s_%s_%s.mat',...
                    root,folder,subj{ss},typestr,condStrs{conduse},voxelStr,predStr,classStr,kernelStr,subMeanStr);

    load(fns);
            
    for vv=1:nVOIs
        
        for dd=1:nDist
            
            theseaccs = classStruct(vuse(vv)).accReal(distOrder(dd),:);
            accs_allsub(vv,ss,dd) = mean(theseaccs);
            
            theseaccs = classStruct(vuse(vv)).accRand(distOrder(dd),:);
            accsrand_allsub(vv,ss,dd,:) = mean(theseaccs,1);
            
            thesed = classStruct(vuse(vv)).dReal(distOrder(dd),:);
            d_allsub(vv,ss,dd) = mean(thesed);
            
            thesed = classStruct(vuse(vv)).dRand(distOrder(dd),:);
            drand_allsub(vv,ss,dd,:) = mean(thesed,1);
        end
        
    end
   
end

%% compute p vals for each roi, each cond separately 
tic
if plotSig
for vv=1:nVOIs
    for dd=1:nDist

        %% all trials
        realAccs = squeeze(accs_allsub(vv,:,dd));
        nullAccs = squeeze(accsrand_allsub(vv,:,dd,:));
%         
%         if length(subj)>1
%             realTAcc = get_tscore_nans(realAccs,chanceVal);
%             nullTAcc = get_tscore_nans(nullAccs',chanceVal);
%         else
%             realTAcc = realAccs;
%             nullTAcc = nullAccs';
%         end
        
        realD = squeeze(d_allsub(vv,:,dd));
        nullD = squeeze(drand_allsub(vv,:,dd,:));
        
%         if length(subj)>1
%             realTD = get_tscore_nans(realD,0);
%             nullTD = get_tscore_nans(nullD',0);
%         else
%             realTD = realD;
%             nullTD = nullD';
%         end

%         if length(nullTAcc)~=nIter 
%             error('wrong number of values in null distrib')
%         end
%          
%         pValsAcc_allsub(vv,dd) = 2*min([mean(realTAcc<nullTAcc),mean(realTAcc>nullTAcc)]);
%         pValsD_allsub(vv,cc) = 2*min([mean(realTD<nullTD),mean(realTD>nullTD)]);
        
        pValsAcc_allsub(vv,dd) = 2*min([mean(mean(realAccs)<mean(nullAccs)),mean(mean(realAccs)>mean(nullAccs))]);
        pValsD_allsub(vv,dd) = 2*min([mean(mean(realD)<mean(nullD)),mean(mean(realD)>mean(nullD))]);


    end
end
toc

%% FDR correction
isSigAcc = zeros(size(pValsAcc_allsub,1),size(pValsAcc_allsub,2),length(sigLevels));
isSigD = zeros(size(pValsAcc_allsub,1),size(pValsAcc_allsub,2),length(sigLevels));

for aa=1:length(sigLevels)
%     for dd=1:nDist
        alpha=sigLevels(aa);
    % 
    %     [p_fdr, p_masked] = fdr(pValsAcc_allsub_condDiff, alpha);
    %     isSigAcc_condDiff(:,:,aa) = p_masked;

    %     [p_fdr, p_masked] = fdr(pValsD_allsub_condDiff, alpha);
    %     isSigD_condDiff(:,:,aa) = p_masked;

        [p_fdr, p_masked] = fdr( pValsAcc_allsub(:,:), alpha);
        isSigAcc(:,:,aa)=p_masked; 

        [p_fdr, p_masked] = fdr( pValsD_allsub, alpha);
        isSigD(:,:,aa)=p_masked; 
%     end
end

end

%% load slopes

fnout=sprintf('%s%s/allSubj_allROIs_%s_%s_slopeByDisp.mat',...
                    root,folder,typestr,voxelStr);

load(fnout)

medlines = zeros(nVOIs,2);
% lblines = zeros(nVOIs,2);
% ublines = zeros(nVOIs,2);

for vv=1:nVOIs
    
    diff=abs(allslopes(vv,:)-prctile(allslopes(vv,:),50));
    thisind = find(diff==min(diff),1);
    medlines(vv,:) = [allslopes(vv,thisind),allint(vv,thisind)];
    
%     diff=abs(allslopes(vv,:)-prctile(allslopes(vv,:),2.5));
%     thisind = find(diff==min(diff),1);
%     lblines(vv,:) = [allslopes(vv,thisind),allint(vv,thisind)];
%     
%     diff=abs(allslopes(vv,:)-prctile(allslopes(vv,:),97.5));
%     thisind = find(diff==min(diff),1);
%     ublines(vv,:) = [allslopes(vv,thisind),allint(vv,thisind)];
    
end

% medlines = [mean(allslopes,2),zeros(nVOIs,1)];


%% make plot for all subs
if plotD
    figure;hold all;

    ii=0;    
    for vv=1:nVOIs
        ii=ii+1;
        ax=subplot(5,2,ii);hold all;
        
        title(sprintf('%s',VOIs{vuse(vv)}))

        
        if nSubj>1
            barErrsAll=nanstd(squeeze(d_allsub(vv,:,:)),[],1)./sqrt(nSubj);
            barMeansAll=nanmean(squeeze(d_allsub(vv,:,:)),1);

        else
            
            barMeansAll=squeeze(d_allsub(vv,:,:));
            barErrsAll=nan(size(barMeansAll));

        end

%         colormap('jet')
%         bar(barMeansAll);
        h1=errorbar(ax,distListSort,barMeansAll(:),barErrsAll(:),'LineWidth',1,'Color',[0,0,0]);
        
        set(gca,'YLim', [drange], 'XTick', distListSort, ...
                        'XTickLabel', distListSort,'XTickLabelRotation',90);
        xlabel('Disparity difference (in arcmin)')
        ylabel('d-prime');
        
        currxlim = get(gca,'XLim');
        line(currxlim,[0,0],'Color','k');
        
        xdat = (currxlim(1):currxlim(2))';
        ydat = medlines(vv,2)+medlines(vv,1)*xdat;

        line(xdat,ydat);
                
        if plotSig
            astLocsEach=nan(nDist,2);
            for dd=1:nDist
    %             for cc=1:nDim
                    if isSigAcc(vv,dd,1)
                        if barMeansAll(dd)>0
                            if nSubj>1
                                astLocsEach(dd,1)=barMeansAll(dd)+barErrsAll(dd)+verspacerbig;
                            else
                                astLocsEach(dd,1)=barMeansAll(dd)+verspacerbig;
                            end
                        else
                            if nSubj>1
                                astLocsEach(dd,1)=barMeansAll(dd)-barErrsAll(dd)-verspacerbig;
                            else
                                astLocsEach(dd,1)=barMeansAll(dd)-verspacerbig;
                            end
                        end
                    end
                    if isSigAcc(vv,dd,2)
                        astLocsEach(dd,2) = astLocsEach(dd,1);                    
                    end
    %             end
            end


            plot(distListSort,astLocsEach(:,1),'o','Color','k','MarkerSize',markersize)
            plot(distListSort,astLocsEach(:,2),'o','Color','k','MarkerFaceColor','k','MarkerSize',markersize)

            uistack(h1,'top')
            
        end
    end
   
    if saveFig
        fnFig = [root figFolder filesep 'Dprime_2Z_oneVsOne_byDisp_wLine_' voxelStr];
        
        fprintf('saving figure to %s...\n',fnFig);
        saveas(gcf,fnFig,ext);
    end

    %% tilefigs

if length(get(0,'Children')) > 1
    tilefigs();
end


end

