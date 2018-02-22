
%% Plot the degrees versus disparity relation (negative is far)

load('IEMdepth_allGridConversions.mat')

figFolder='IEMdepth_figs';
ext='epsc';

figure;hold all;
scatter(stimLocsDeg_odd(:,1),zLocsArcMin_odd,'.','k');
title('DVA vs Disparity - grid 1')
xlim([-screenHeightDeg/2,screenHeightDeg/2]);
ylim([-50,50]);
xlabel('Left-Right axis (degrees)');
ylabel('Front-back axis (arcmin)');


fnFig = [root figFolder filesep 'Grid1_DispVsDeg'];
fprintf('saving figure to %s...\n',fnFig);
saveas(gcf,fnFig,ext);

figure;hold all;
scatter(stimLocsDeg_even(:,1),zLocsArcMin_even,'.','k');
title('DVA vs Disparity - grid 2')
xlim([-screenHeightDeg/2,screenHeightDeg/2]);
ylim([-50,50]);
xlabel('Left-Right axis (degrees)');
ylabel('Front-back axis (arcmin)');

fnFig = [root figFolder filesep 'Grid2_DispVsDeg'];
fprintf('saving figure to %s...\n',fnFig);
saveas(gcf,fnFig,ext);