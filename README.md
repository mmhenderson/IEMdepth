# IEMdepth
We presented subjects with stereoscopic sphere stimuli, composed of multicolored flickering dots. Stimuli were presented at various positions evenly spaced along the horizontal (left-right) spatial axis, and the depth (near-far) spatial axis. We used an inverted encoding model to reconstruct the locations of stimuli along each of these spatial axes.

#### Contents:

## Reconstructions 
- Contains all code required to generate channel response functions (IEMdepth_getChannelResp1D.m), average these reconstructions over shared locations (IEMdepth_average1DRecon.m), and perform curve fitting to the reconstructions (IEMdepth_fitRecons1D.m). 
- Also contains code to plot reconstructions and analyze their accuracy over positions and ROIs.
## Classifier
- Contains all code needed to run a linear decoding analysis using libSVM 3.1 (libSVM is not included here but can be downloaded from https://www.csie.ntu.edu.tw/~cjlin/libsvm/)
## Behavior
- Contains all code needed to analyze subject behavior on our tasks, as well as visualize the grid used to present stimuli in the experiment. 