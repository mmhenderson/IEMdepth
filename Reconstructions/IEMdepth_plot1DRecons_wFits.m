function IEMdepth_plot1DRecons_wFits(varargin)

if nargin < 1    
    subj = {'AI','AP','BB','BC','BD','BJ','BM','BN','BO'};
%     plotVOIs = {'V3A','V3B'};
    plotVOIs = {'V1','V2','V3','V4','V3A','V3B','IPS0','LO1','LO2'};
elseif nargin == 1
    subj = varargin{1};
    plotVOIs = {'V1','V2','V3','V4','V3A','V3B','IPS0','LO1','LO2'};
end

root = '/usr/local/serenceslab/maggie/IEMdepth/';
dimlabels = {'X','Z'};
% clabels = {'stim','fix'};
% cm = lines(6);
% cplotstyles = {'--',':'};

sessi = 2;
sepconds = 2;
absStat = 2;
usestereox = 2;

sessStrs = {'_allGoodRuns','_allSess'};
sessStr = sessStrs{sessi};
sepcondstrs = {'_trnAllCond','_trnSepCond','_trnFixTstStim'};
sepcondstr = sepcondstrs{sepconds};
locSignStrs = {'posVoxOnly','allVoxAbs'};
locSignStr = locSignStrs{absStat};
stereostrs = {'_screenx','_stereox'};
stereostr = stereostrs{usestereox};

plotCond = 2;

% colors from colorbrewer
col1 = [8,88,158; 43,140,190; 78,179,211; 123,204,196; 168,221,181; 204,235,197] ./ 255;
col2 = [110,1,107; 136,65,157; 140,107,177; 140,150,198; 158,188,218; 191,211,230] ./ 255;
allc = cat(3, col1, col2);

%%
for ss = 1:length(subj)
    
    fnf = sprintf('%sIEMdepth_reconFits_vy/%s_%s_fits_1D%s%s%s.mat',...
            root,subj{ss},locSignStr,sepcondstr,sessStr,stereostr);
    fprintf('loading %s...\n\n',fnf);
    load(fnf);
    
    voin = find(ismember(VOIs,plotVOIs));    
    
    for dim = 1:2
%     for dim = 1
        cm = allc(:,:,dim);
        figure;
        for vx = 1:length(plotVOIs)
        
            % load the 1D recon file
            if sepconds == 2
                fn = sprintf('%sIEMdepth_chanResp/%s_%s_avgRecons1D_%s%s%s%s.mat',...
                    root,subj{ss},VOIs{voin(vx)},locSignStr,...
                    sepcondstr,sessStr,stereostr);
            else
                fn = sprintf('%sIEMdepth_chanResp/%s_%s_avgRecons1D.mat',...
                    root,subj{ss},VOIs{voin(vx)});
            end
            load(fn);

            %%
            for i = 1:size(avgRec_1D(dim).avgRecons,1)
                subplot(2,ceil(length(voin)/2),vx); hold all;

                hl = line(repmat(avgRec_1D(dim).stimLocs(i),1,2),...
                    [-0.5 1.5],'Color',[cm(i,:)]);
                hl.LineWidth = 1;
                plot(avgRec_1D(dim).basis_grid,...
                    avgRec_1D(dim).avgRecons(i,:,plotCond),...
                    'Color',[cm(i,:)],'LineWidth',1);
                hold on;
                % fitdat params: ndim x ncond x nrec x npts
                plot(avgRec_1D(dim).basis_grid,...
                    squeeze(fitdat(voin(vx)).allFits(dim,plotCond,i,:)),...
                    'Color',cm(i,:),'LineStyle',':','LineWidth',2);
                ylim([-0.5 1.5]);
                xlim([-3 3]);
            end
            % plot label stuff
            ax = gca;
%             if vx == 1
%                 title(dimlabels{dim});
%             else
            if vx == length(VOIs)
                xlabel('Position (openGL coords)');
                legend([{'stimulus location','recon','fit'}]);
            elseif vx == length(VOIs)
                ax.YTickLabel = {};
            end
            if vx ~= length(VOIs)
%                 ylabel(VOIs{voin(vx)}); 
                ax.XTickLabel = {};
            elseif vx ~= length(VOIs)
                ax.XTickLabel = {};
                ax.YTickLabel = {};
            end
            
            title(VOIs{voin(vx)});
%             pp = pp + 1;
        end % end VOI loop
        suptitle(sprintf('%s %s 1D recons, averaged',subj{ss},dimlabels{dim}));
        prepFigForExport;
        
    end % end dim loop
    
end % end sub loop