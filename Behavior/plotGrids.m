% Make figures for the opengl/disparity space

clear 
close all

load('IEMdepth_allGridConversions.mat')

% figFolder='IEMdepth_figs';
% ext='epsc';


%% plot diagram from above

ylimsOGL = [-2.5,11];
  
figure;hold all;

minz = find(stimLocsOGL_odd(:,2)==-1.5,1);
viscircles([stimLocsOGL_odd(minz,1),stimLocsOGL_odd(minz,2)],sphereRadsOGL(1));
maxz = find(stimLocsOGL_odd(:,2)==1.5,1);
viscircles([stimLocsOGL_odd(maxz,1),stimLocsOGL_odd(maxz,2)],sphereRadsOGL(6));

for ll=1:12;    
    pt1 = [0,10];
    pt2 = stimLocsOGL_odd(ll,:);
    line([pt1(1),pt2(1)],[pt1(2),pt2(2)],'LineStyle',':');
end

xpts = [-.4,.4];
ypts = [10,10];
scatter(xpts,ypts,'k');

pt1 = [.4,10];
pt2 = [0,0];
line([pt1(1),pt2(1)],[pt1(2),pt2(2)],'Color','k');

pt1 = [-.4,10];
pt2 = [0,0];
line([pt1(1),pt2(1)],[pt1(2),pt2(2)],'Color','k');


scatter(stimLocsOGL_odd(:,1),stimLocsOGL_odd(:,2),'.','k');
xpts = [0];
ypts = [10];
scatter(xpts,ypts,'.','k');
title('OpenGL space - grid 1')
% xlim([-screenHeightOGL(1)/2,screenHeightOGL(1)/2]);
ylim(ylimsOGL);
xlabel(sprintf('Horizontal axis (openGL units)'))
ylabel(sprintf('Depth axis (openGL units)'))
axis equal 
line(get(gca,'XLim'),[0,0],'LineStyle','--')
scatter(0,0,'x','k')


% fnFig = [root figFolder filesep 'Grid1_OGL_fulldiagram'];
% fprintf('saving figure to %s...\n',fnFig);
% saveas(gcf,fnFig,ext);

%% plot disparity versus x-location 


figure;hold all;
scatter(stimLocsDeg_odd(:,1),zLocsArcMin_odd,'.','k');
title('Disparity vs Deg - grid 1')
xlim([-screenHeightDeg/2,screenHeightDeg/2]);
ylim([-50,50]);
xlabel('Left-Right axis (degrees)');
ylabel('Front-back axis (arcmin)');
plot(0,0,'square','Color','k')

% fnFig = [root figFolder filesep 'Grid1_DispVsDeg'];
% fprintf('saving figure to %s...\n',fnFig);
% saveas(gcf,fnFig,ext);

figure;hold all;
scatter(stimLocsDeg_even(:,1),zLocsArcMin_even,'.','k');
title('Disparity vs Deg - grid 2')
xlim([-screenHeightDeg/2,screenHeightDeg/2]);
ylim([-50,50]);
xlabel('Left-Right axis (degrees)');
ylabel('Front-back axis (arcmin)');
plot(0,0,'square','Color','k')
% 
% fnFig = [root figFolder filesep 'Grid2_DispVsDeg'];
% fprintf('saving figure to %s...\n',fnFig);
% saveas(gcf,fnFig,ext);

%% plot Z-position versus X-location

OGLRange = 2;

figure;hold all;

scatter(stimLocsDeg_odd(:,1),stimLocsOGL_odd(:,2),'.','k');

title('Z-Position vs Deg - Grid 1')
xlim([-screenHeightDeg/2,screenHeightDeg/2]);
ylim([-OGLRange,OGLRange]);
xlabel('Horizontal axis (Degrees)');
ylabel('Depth axis (OpenGL)');
axis square;
plot(0,0,'square','Color','k')
% 
% fnFig = [root figFolder filesep 'Grid1_OGLZ_DegX'];
% fprintf('saving figure to %s...\n',fnFig);
% saveas(gcf,fnFig,ext);


figure;hold all;

scatter(stimLocsDeg_even(:,1),stimLocsOGL_even(:,2),'.','k');

title('Z-Position vs Deg - Grid 2')
xlim([-screenHeightDeg/2,screenHeightDeg/2]);
ylim([-OGLRange,OGLRange]);
xlabel('Horizontal axis (Degrees)');
ylabel('Depth axis (OpenGL)');
axis square;
plot(0,0,'square','Color','k')

% fnFig = [root figFolder filesep 'Grid2_OGLZ_DegX'];
% fprintf('saving figure to %s...\n',fnFig);
% saveas(gcf,fnFig,ext);
