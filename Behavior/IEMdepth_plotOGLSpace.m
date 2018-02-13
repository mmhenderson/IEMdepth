
%% make figure for the opengl space


load('IEMdepth_allGridConversions.mat')

figFolder='IEMdepth_figs';
ext='epsc';

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


fnFig = [root figFolder filesep 'Grid1_OGL_fulldiagram'];
fprintf('saving figure to %s...\n',fnFig);
saveas(gcf,fnFig,ext);


%%

OGLRange =3;


figure;hold all;
% subplot(1,2,1)
scatter(stimLocsOGL_odd(:,1),stimLocsOGL_odd(:,2),'.','k');
minz = 1;
viscircles([stimLocsOGL_odd(minz,1),stimLocsOGL_odd(minz,2)],sphereRadsOGL(1));
maxz = 36;
viscircles([stimLocsOGL_odd(maxz,1),stimLocsOGL_odd(maxz,2)],sphereRadsOGL(6));
title('Grid 1')
xlim([-OGLRange,OGLRange]);
ylim([-OGLRange,OGLRange]);
xlabel('Horizontal axis (OpenGL)');
ylabel('Depth axis (OpenGL)');
axis square;

fnFig = [root figFolder filesep 'Grid1_OGL'];
fprintf('saving figure to %s...\n',fnFig);
saveas(gcf,fnFig,ext);


figure;hold all;
% subplot(1,2,2);
scatter(stimLocsOGL_even(:,1),stimLocsOGL_even(:,2),'.','k');
minz = 1;
viscircles([stimLocsOGL_even(minz,1),stimLocsOGL_even(minz,2)],sphereRadsOGL(1));
maxz = 36;
viscircles([stimLocsOGL_even(maxz,1),stimLocsOGL_even(maxz,2)],sphereRadsOGL(6));
title('Grid 2');
xlim([-OGLRange,OGLRange]);
ylim([-OGLRange,OGLRange]);
xlabel('Horizontal axis (OpenGL)');
ylabel('Depth axis (OpenGL)');
axis square;

fnFig = [root figFolder filesep 'Grid2_OGL'];
fprintf('saving figure to %s...\n',fnFig);
saveas(gcf,fnFig,ext);