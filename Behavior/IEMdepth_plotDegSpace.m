
%% make figures in degree space

load('IEMdepth_allGridConversions.mat')

figFolder='IEMdepth_figs';
ext='epsc';

% degRange=12.5;
degRange = 14;

figure;hold all;
% subplot(1,2,1)
scatter(stimLocsDeg_odd(:,1),stimLocsDeg_odd(:,2),'.','k');
minz = 1;
viscircles([stimLocsDeg_odd(minz,1),stimLocsDeg_odd(minz,2)],sphereRadDeg(1));
maxz = 36;
viscircles([stimLocsDeg_odd(maxz,1),stimLocsDeg_odd(maxz,2)],sphereRadDeg(6));
title('Grid 1')
xlim([-degRange,degRange]);
ylim([-degRange,degRange]);
xlabel(sprintf('Horizontal axis (%c)',char(176)));
ylabel(sprintf('Depth axis (%c)',char(176)));
axis square;

fnFig = [root figFolder filesep 'Grid1_deg'];
fprintf('saving figure to %s...\n',fnFig);
saveas(gcf,fnFig,ext);


figure;hold all;
% subplot(1,2,2);
scatter(stimLocsDeg_even(:,1),stimLocsDeg_even(:,2),'.','k');
minz = 1;
viscircles([stimLocsDeg_even(minz,1),stimLocsDeg_even(minz,2)],sphereRadDeg(1));
maxz = 36;
viscircles([stimLocsDeg_even(maxz,1),stimLocsDeg_even(maxz,2)],sphereRadDeg(6));
title('Grid 2');
xlim([-degRange,degRange]);
ylim([-degRange,degRange]);
xlabel(sprintf('Horizontal axis (%c)',char(176)));
ylabel(sprintf('Depth axis (%c)',char(176)));
axis square;

fnFig = [root figFolder filesep 'Grid2_deg'];
fprintf('saving figure to %s...\n',fnFig);
saveas(gcf,fnFig,ext);


