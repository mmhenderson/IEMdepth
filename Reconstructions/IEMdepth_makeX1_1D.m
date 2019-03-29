function [X,basis_set,gridPts] = IEMdepth_makeX1_1D(stimLocs,stimSize,basisDim,fwhm,res)
% VAV 2/24/2016 adapted from IEMdepth_makeX_1D, but does not include basis
% set definition (which should be defined in the script calling this one &
% input into the parameter "basisDim")


% FOV will be defined by the basis set - for computing encoding model, at
% least
%
%
% note: does not!!! normalize X, need to do that wherever this is called
% (channelRespAmp)

% go from FWHM to rad
filt_size =  (1/rad2fwhm(1)) * fwhm;

% generate x, y grid for make_stim_mask and build_basis
% only need to go to the edge of the most distant stimulus x/y (as all
% other stim_mask values will be zero)
fov = 2*max(stimLocs(:)) + 4*stimSize;

% maybe add checks here to make sure spacing == fwhm

gridPts = linspace(-fov/2,fov/2,res);

% build our basis set (nPixels x nChannels)
basis_set = cos_basis_1D(basisDim,gridPts,filt_size); %pixels 

% get stimulus mask for each actual position (nPixels x nTrials) 
stim_mask = make_stim_mask_1D(stimLocs,stimSize,gridPts);

%get design matrix - each element gives weight of channel (column)
%to the actual stimulus location on each particular trial (row)
%nTrials x nChannels (length(stimLocs) x numel(rfPtsX))
X = stim_mask' * basis_set;

% intentionally not normalizing X here, it should be normalized across all
% runs, and this will work on a single run...

return