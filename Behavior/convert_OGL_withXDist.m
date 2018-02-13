%% convert the units of openGL space to degrees visual angle and arcmin disparity
% all drawing in the psychtoolbox script was done in openGL units -
% so all the values saved by the script (stimLocs, sphereSize), are in
% these units. Convert them back to units that are interpretable for
% making figure axes. 

clear

close all

root='/usr/local/serenceslab/maggie/IEMdepth/';
saveFig = 1;
figFolder='IEMdepth_figs';
ext='epsc';

makeDegFigs=1;
makeOGLFigs=1;

%% load some behavioral files to get the exact params
%load an even and an odd run
fnbehav = '/mnt/neurocube/local/serenceslab/maggie/IEMdepth/BB143/BB143_Behav/BB143_IEMdepth_sess01_run01_task1.mat';
load(fnbehav)
pOdd=p;

fnbehav = '/mnt/neurocube/local/serenceslab/maggie/IEMdepth/BB143/BB143_Behav/BB143_IEMdepth_sess01_run02_task1.mat';
load(fnbehav)
pEven=p;

% load a localizer run
fnloc = '/mnt/neurocube/local/serenceslab/maggie/IEMdepth/BJ141/BJ141_Behav/BJ141_IEMdepth_sess01_LocalizerRun02.mat';
load(fnloc)
pLoc = p;

locLimsFrontOGL = [ -pLoc.quadLims{3}(1,2),pLoc.quadLims{3}(1,2);
                    pLoc.quadLims{3}(2,1),pLoc.quadLims{3}(2,2)];
locLimsBackOGL = [ -pLoc.quadLims{2}(1,2),pLoc.quadLims{2}(1,2);
                    pLoc.quadLims{2}(2,1),pLoc.quadLims{2}(2,2)];
                  
%%

screenSizePix = pOdd.sRect(3:4);
pixHalfX = screenSizePix(1)/2;
pixHalfY = screenSizePix(2)/2;

% this is half of the FOV, which is 25 total
% NOTE that the conversions here from openGL to degrees
% do NOT actually depend on the FOV in
% degrees. Degrees visual angle of each openGL position is only dependent
% on the location of the camera and the point in openGL space. 
% The FOV does determine the total range of position that are plot-able in 
% openGL space (i.e how much of the screen we sampled)
degHalfY = 25/2;
% the field of the view is set vertically (y) in degrees, at 25 degrees to
% each side (50 degrees total)
pixPerDeg = pixHalfY/degHalfY;
% approximate the FOV in x:
degHalfX = degHalfY*screenSizePix(1)/screenSizePix(2);

% screen size and pixel conversion will both change as a function of where
% we are in depth
oglPerPix = zeros(6,1);
screenSizeOGL = zeros(6,2);

% sort these by ascending X, and then ascending Z
stimLocsOGL_odd = pOdd.stimLocs((pOdd.stimLocs(:,1)~=0),[1,3]);
[~,sortorder] = sort(stimLocsOGL_odd(:,1),'ascend');
stimLocsOGL_odd = stimLocsOGL_odd(sortorder,:);
[~,sortorder] = sort(stimLocsOGL_odd(:,2),'ascend');
stimLocsOGL_odd = stimLocsOGL_odd(sortorder,:);

stimLocsOGL_even = pEven.stimLocs((pEven.stimLocs(:,1)~=0),[1,3]);
[~,sortorder] = sort(stimLocsOGL_even(:,1),'ascend');
stimLocsOGL_even = stimLocsOGL_even(sortorder,:);
[~,sortorder] = sort(stimLocsOGL_even(:,2),'ascend');
stimLocsOGL_even = stimLocsOGL_even(sortorder,:);

% z locations in openGL space, listed back to front
zLocsOGL = unique(stimLocsOGL_odd(:,2));

% distance from camera to center (fixation plane)
zDistFix = pOdd.zloccamera;
% this is the distance from the camera (viewer) to the point, going back to front  
zDistOGL = zDistFix - unique(stimLocsOGL_odd(:,2));

% distance from each eye to the center- hard coded into the script
eyeToCenterOGL = 0.4;

% calculate the angle of vergence at the fixation plane - using distance to
% screen an interocular distance
fixPlaneDeg = 2*rad2deg(atan(eyeToCenterOGL/zDistFix));
fixPlaneArcMin = fixPlaneDeg*60;

% we will calculate the disparity angle for each stimulus pos
zLocsArcMin_odd = zeros(size(stimLocsOGL_odd));
zLocsArcMin_even = zeros(size(stimLocsOGL_even));

% sphere sizes, which scale with distance
% listed back to front
sphereRadsOGL = unique(pOdd.stimScales(~isnan(pOdd.stimScales)));
sphereRadsOGL = sphereRadsOGL(end:-1:1);

% this will store the [x,z] position of each stimulus
% column 1 is X in degrees visual angle
% column 2 is "degrees in z", this is not really a measurable thing, it is
% just the same as how we convert OGL to degrees in Z (rough estimate of distance). 
stimLocsDeg_odd = zeros(size(stimLocsOGL_odd,1),2);
stimLocsDeg_even = zeros(size(stimLocsOGL_even,1),2);

% this is a single vector describing the Z disparity of each of the 36
% positions.
zLocsArcMin_odd = zeros(size(stimLocsOGL_even,1),1);
zLocsArcMin_even = zeros(size(stimLocsOGL_odd,1),1);

% store the screen height and conversion for each Z plane (assuming X=0)
screenHeightOGL = zeros(size(zDistOGL));
oglPerDeg = zeros(size(zDistOGL));

% loop over z positions, back to front
for dd=1:length(zDistOGL)

    % convert to opengl units (based on the z distance to the plane, and
    % the vertical angle in degrees)
    oglHalfY = (zDistOGL(dd))*tan(deg2rad(degHalfY));
    screenHeightOGL(dd) = oglHalfY*2;

     % calculate the conversion units
    oglPerDeg(dd) = oglHalfY/pixHalfY*pixPerDeg;

    % which indices have these Z positions? same in the even and odd grid
    % matrices
    theseinds = find(stimLocsOGL_odd(:,2)==zLocsOGL(dd));
    unX_odd = stimLocsOGL_odd(theseinds,1);
    unX_even = stimLocsOGL_even(theseinds,1);
     
    % get the adjusted x positions for these locations    
    stimLocsDeg_odd(theseinds,1) = unX_odd./oglPerDeg(dd);
    stimLocsDeg_even(theseinds,1) = unX_even./oglPerDeg(dd);

    % loop over the X positions, get actual disparity as function of both X
    % and Z   
    for xx=1:numel(unX_odd)
    
        % using the actual x,z coordinates - calculate the distance from
        % observer to the point (draw a horopter)
        actualZDist = sqrt(unX_odd(xx)^2 + zDistOGL(dd)^2);
        % get the disparity of this (x,z) position
        zDistDeg = 2*rad2deg(atan(eyeToCenterOGL/actualZDist));
        zLocsArcMin_odd(theseinds(xx)) = zDistDeg*60 - fixPlaneArcMin;

        % repeat for the even grid
        actualZDist = sqrt(unX_even(xx)^2 + zDistOGL(dd)^2);
        % get the disparity of this (x,z) position
        zDistDeg = 2*rad2deg(atan(eyeToCenterOGL/actualZDist));
        zLocsArcMin_even(theseinds(xx)) = zDistDeg*60 - fixPlaneArcMin;

    end
end

locLimsFrontDeg = locLimsFrontOGL/oglPerDeg(6);
locLimsBackDeg = locLimsBackOGL/oglPerDeg(1);

% a dumb way of getting a "Z distance" in degree units
stimLocsDeg_odd(:,2) = stimLocsOGL_odd(:,2)./mean(oglPerDeg);
stimLocsDeg_even(:,2) = stimLocsOGL_even(:,2)./mean(oglPerDeg);

% convert sphere sizes 
sphereRadDeg = sphereRadsOGL./oglPerDeg;

% sanity check - the height is the same number of degrees everywhere
screenHeightDeg = screenHeightOGL./oglPerDeg;
screenHeightDeg = screenHeightDeg(1);

degRange=12.5;
ylimsOGL = [-2.5,11];

%% Plot the degrees versus disparity relation (negative is far)

figure;hold all;
scatter(stimLocsDeg_odd(:,1),zLocsArcMin_odd);
title('DVA vs Disparity - grid 1')
xlim([-screenHeightDeg/2,screenHeightDeg/2]);
ylim([-50,50]);
xlabel('Left-Right axis (degrees)');
ylabel('Back-Front axis (arcmin)');

figure;hold all;
scatter(stimLocsDeg_even(:,1),zLocsArcMin_even);
title('DVA vs Disparity - grid 2')
xlim([-screenHeightDeg/2,screenHeightDeg/2]);
ylim([-50,50]);
xlabel('Left-Right axis (degrees)');
ylabel('Back-Front axis (arcmin)');


%% make figure for the opengl space
if makeOGLFigs
    
figure;hold all;

minz = find(stimLocsOGL_odd(:,2)==-1.5,1);
viscircles([stimLocsOGL_odd(minz,1),stimLocsOGL_odd(minz,2)],sphereRadsOGL(1))
maxz = find(stimLocsOGL_odd(:,2)==1.5,1);
viscircles([stimLocsOGL_odd(maxz,1),stimLocsOGL_odd(maxz,2)],sphereRadsOGL(6))

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

if saveFig
    fnFig = [root figFolder filesep 'Grid1_OGL'];
    fprintf('saving figure to %s...\n',fnFig);
    saveas(gcf,fnFig,ext);
end

end
%%

if makeDegFigs

figure;hold all;
% subplot(1,2,1)
scatter(stimLocsDeg_odd(:,1),stimLocsDeg_odd(:,2),'.','k');
minz = find(stimLocsDeg_odd(:,2)==min(stimLocsDeg_odd(:,2)),1);
viscircles([stimLocsDeg_odd(minz,1),stimLocsDeg_odd(minz,2)],sphereRadDeg(1))
maxz = find(stimLocsDeg_odd(:,2)==max(stimLocsDeg_odd(:,2)),1);
viscircles([stimLocsDeg_odd(maxz,1),stimLocsDeg_odd(maxz,2)],sphereRadDeg(6))
title('Grid 1')
xlim([-degRange,degRange]);
ylim([-degRange,degRange]);
xlabel(sprintf('Horizontal axis (%c)',char(176)))
ylabel(sprintf('Depth axis (%c)',char(176)))
axis square;
if saveFig
    fnFig = [root figFolder filesep 'Grid1_deg'];
    fprintf('saving figure to %s...\n',fnFig);
    saveas(gcf,fnFig,ext);
end
% markersize=1;
figure;hold all;
% subplot(1,2,2);
scatter(stimLocsDeg_even(:,1),stimLocsDeg_even(:,2),'.','k');
minz = find(stimLocsDeg_even(:,2)==min(stimLocsDeg_even(:,2)),1);
viscircles([stimLocsDeg_even(minz,1),stimLocsDeg_even(minz,2)],sphereRadDeg(1))
maxz = find(stimLocsDeg_even(:,2)==max(stimLocsDeg_even(:,2)),1);
viscircles([stimLocsDeg_even(maxz,1),stimLocsDeg_even(maxz,2)],sphereRadDeg(6))
title('Grid 2')
xlim([-degRange,degRange]);
ylim([-degRange,degRange]);
xlabel(sprintf('Horizontal axis (%c)',char(176)))
ylabel(sprintf('Depth axis (%c)',char(176)))
axis square;
if saveFig
    fnFig = [root figFolder filesep 'Grid2_deg'];
    fprintf('saving figure to %s...\n',fnFig);
    saveas(gcf,fnFig,ext);
end

end



% tilefigs()