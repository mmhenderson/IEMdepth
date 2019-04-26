#setwd("D:R_scripts")
library(lme4)
library(tidyverse)
require(sjPlot)
library(lsmeans)

# read in my data
fn2 = "/mnt/neurocube/local/serenceslab/maggie/mFiles/IEMdepth/MixedModels/recon_err_tbl.txt"
fdat = read_csv(fn2)

# label the ROIs
fdat$ROI[fdat$ROI == 1] <- 'V1'
fdat$ROI[fdat$ROI == 2] <- 'V2'
fdat$ROI[fdat$ROI == 3] <- 'V3'
fdat$ROI[fdat$ROI == 4] <- 'V4'
fdat$ROI[fdat$ROI == 5] <- 'V3A'
fdat$ROI[fdat$ROI == 6] <- 'V3B'
fdat$ROI[fdat$ROI == 7] <- 'IPS0'
fdat$ROI[fdat$ROI == 8] <- 'LO1'
fdat$ROI[fdat$ROI == 9] <- 'LO2'

# make the subject variable categorical
fdat$ROI <- factor(fdat$ROI)
fdat$subject <- factor(fdat$subject)
fdat$position = factor(fdat$position)

fm0 = lmer(err~1+ (1|subject), data=fdat, REML=FALSE)
with(fm5@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)
fm1 = lmer(err~position + (1|subject), data=fdat, REML=FALSE)
fm2 = lmer(err~position + (1|subject) + (1|position:subject),  data=fdat, REML=FALSE)
fm3 = lmer(err~position + ROI + (1|subject) + (1|position:subject),  data=fdat, REML=FALSE)
fm4 = lmer(err~position + ROI + (1|subject) + (1|position:subject) + (1|ROI:subject), data=fdat, REML=FALSE)
fm5 = lmer(err~position*ROI + (1|subject)+ (1|position:subject) + (1|ROI:subject), data=fdat, REML=FALSE)
with(fm1@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)
with(fm2@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)
with(fm3@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)
with(fm4@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)
with(fm5@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)

anova(fm0,fm1,fm2,fm3,fm4,fm5)

# pairwise comparisons, correcting with the Tukey method.

lsmeans(fm3, pairwise~ROI, adjust='tukey')
lsmeans(fm1, pairwise~position, adjust='tukey')
