#setwd("D:R_scripts")
library(lme4)
library(tidyverse)
require(sjPlot)
library(lsmeans)
# Read in the classifier accuracies.
fn = '/mnt/neurocube/local/serenceslab/maggie/mFiles/IEMdepth/MixedModels/dprime_classifier_tbl.txt'

dat <- read.csv(fn)

# label the ROIs
dat$ROI[dat$ROI == 1] <- 'V1'
dat$ROI[dat$ROI == 2] <- 'V2'
dat$ROI[dat$ROI == 3] <- 'V3'
dat$ROI[dat$ROI == 4] <- 'V4'
dat$ROI[dat$ROI == 5] <- 'V3A'
dat$ROI[dat$ROI == 6] <- 'V3B'
dat$ROI[dat$ROI == 7] <- 'IPS0'
dat$ROI[dat$ROI == 8] <- 'LO1'
dat$ROI[dat$ROI == 9] <- 'LO2'

# make the subject variable categorical
dat$ROI <- factor(dat$ROI)
dat$subject <- factor(dat$subject)

# Each column in dat is: 'dprime','ROI','disparity','subject'

# run the mixed effects models
lm0 = lmer(dprime~1 + (1|subject), data=dat, REML=FALSE)
lm1 = lmer(dprime~disparity + (1|subject), data=dat, REML=FALSE)
lm2 = lmer(dprime~disparity + ROI + (1|subject), data=dat, REML=FALSE)
lm3 = lmer(dprime~disparity + ROI + (1|subject) + (1|ROI:subject),  data=dat, REML=FALSE)
lm4 = lmer(dprime~disparity + ROI + (1|subject) + (1|ROI:subject) + (disparity|subject), data=dat, REML=FALSE)
lm5 = lmer(dprime~disparity*ROI + (1|subject) + (1|ROI:subject) + (disparity|subject), data=dat, REML=FALSE)
anova(lm0,lm1,lm2,lm3,lm4,lm5)

# pairwise comparisons bw all ROIs. this is doing a bunch of paired t-tests 
#(after averaging across position), and then correcting with the Tukey method.
lsmeans(lm5, pairwise~ROI, adjust='tukey')

##################################################################
# Now the reconstruction error

fn2 = "/mnt/neurocube/local/serenceslab/maggie/mFiles/IEMdepth/MixedModels/recon_err_tbl_new.txt"
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
# fdat$err = abs(fdat$err)
fdat$position = factor(fdat$position)

fm_randsubs = lmer(err~(1|subject), data=fdat, REML=FALSE)
fm_randsubpos = lmer(err~(1|subject) + (1|position:subject), data = fdat, REML=FALSE)
fm_randsubroi = lmer(err~(1|subject) + (1|ROI:subject), data=fdat, REML=FALSE)
fm_allrand = lmer(err~(1|subject) + (1|position:subject) + (1|ROI:subject),  data=fdat, REML=FALSE)
# this last model gives an error - this is "spurious", can validate it with this line. It returns TRUE meaning 
# convergence was actually a success.
with(fm_allrand@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)

fm_pos = lmer(err~position + (1|subject) + (1|position:subject) + (1|ROI:subject),  data=fdat, REML=FALSE)
fm_roi = lmer(err~ROI + (1|subject) + (1|position:subject) + (1|ROI:subject),  data=fdat, REML=FALSE)
fm_add = lmer(err~ROI+position + (1|subject) + (1|position:subject) + (1|ROI:subject),  data=fdat, REML=FALSE)
fm_int = lmer(err~ROI*position + (1|subject) + (1|position:subject) + (1|ROI:subject),  data=fdat, REML=FALSE)
# this last model gives an error - this is "spurious", can validate it with this line. It returns TRUE meaning 
# convergence was actually a success.
with(fm_int@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)
# 
# # is there a random position/subject interaction?
# anova(fm_randsubpos, fm_randsubs)
# anova(fm_allrand, fm_randsubroi)
# 
# # is there a random ROI/subject interaction?
# anova(fm_randsubroi, fm_randsubs)
# anova(fm_allrand, fm_randsubpos)

# is position significant? (these give same answer)
anova(fm_pos, fm_allrand) 
anova(fm_add, fm_roi)

# is roi significant? (these give same answer)
anova(fm_roi, fm_allrand)
anova(fm_add, fm_pos) 

# is interaction significant?
anova(fm_int, fm_add)

#  pos/subject only
fm1 = lmer(err~(1|subject) + (1|position:subject),  data=fdat, REML=FALSE)
fm2 = lmer(err~position + (1|subject) + (1|position:subject), data=fdat, REML=FALSE)
fm3 = lmer(err~position + ROI + (1|subject) + (1|position:subject), data=fdat, REML=FALSE)
fm4 = lmer(err~position*ROI + (1|subject)+ (1|position:subject), data=fdat, REML=FALSE)
anova(fm1,fm2,fm3,fm4)

# roi/subject only
fm1 = lmer(err~(1|subject) + (1|ROI:subject),  data=fdat, REML=FALSE)
fm2 = lmer(err~position + (1|subject) + (1|ROI:subject), data=fdat, REML=FALSE)
fm3 = lmer(err~position + ROI + (1|subject) + (1|ROI:subject), data=fdat, REML=FALSE)
fm4 = lmer(err~position*ROI + (1|subject)+ (1|ROI:subject), data=fdat, REML=FALSE)
anova(fm1,fm2,fm3,fm4)

# all rand effects
fm1 = lmer(err~(1|subject) + (1|ROI:subject) + (1|position:subject),  data=fdat, REML=FALSE)
fm2 = lmer(err~position + (1|subject) + (1|ROI:subject) + (1|position:subject), data=fdat, REML=FALSE)
fm3 = lmer(err~position + ROI + (1|subject) + (1|ROI:subject) + (1|position:subject), data=fdat, REML=FALSE)
fm4 = lmer(err~position*ROI + (1|subject)+ (1|ROI:subject) + (1|position:subject), data=fdat, REML=FALSE)
anova(fm1,fm2,fm3,fm4)
lsmeans(fm4, pairwise~position, adjust="Tukey")

# just subject
fm1 = lmer(err~(1|subject),  data=fdat, REML=FALSE)
fm2 = lmer(err~ROI + (1|subject), data=fdat, REML=FALSE)
fm3 = lmer(err~position + ROI + (1|subject), data=fdat, REML=FALSE)
fm4 = lmer(err~position*ROI + (1|subject), data=fdat, REML=FALSE)
anova(fm1,fm2,fm3,fm4)

# the sequential model that we are using here
fm0 = lmer(err~1+ (1|subject), data=fdat, REML=FALSE)
fm1 = lmer(err~position + (1|subject), data=fdat, REML=FALSE)
fm2 = lmer(err~position + (1|subject) + (1|position:subject),  data=fdat, REML=FALSE)
fm3 = lmer(err~position + ROI + (1|subject) + (1|position:subject),  data=fdat, REML=FALSE)
fm4 = lmer(err~position + ROI + (1|subject) + (1|position:subject) + (1|ROI:subject), data=fdat, REML=FALSE)
fm5 = lmer(err~position*ROI + (1|subject)+ (1|position:subject) + (1|ROI:subject), data=fdat, REML=FALSE)
anova(fm0,fm1,fm2,fm3,fm4,fm5)

# fm0 = lmer(err~1 + (1|subject), data=fdat, REML=FALSE)
# fm1 = lmer(err~position + (1|subject), data=fdat, REML=FALSE)
# fm2 = lmer(err~position + ROI + (1|subject) ,  data=fdat, REML=FALSE)
# fm3 = lmer(err~position + ROI + (1|subject) + (1|ROI:subject),  data=fdat, REML=FALSE)
# fm4 = lmer(err~position + ROI + (1|subject) + (1|ROI:subject) + (1|position:subject), data=fdat, REML=FALSE)
# fm5 = lmer(err~position*ROI + (1|subject)+ (1|position:subject) + (1|ROI:subject), data=fdat, REML=FALSE)
# anova(fm0,fm1,fm2,fm3,fm4,fm5)

# pairwise comparisons bw all ROIs. this is doing a bunch of paired t-tests 
#(after averaging across position), and then correcting with the Tukey method.
lsmeans(fm5, pairwise~ROI, adjust='tukey')

# pairwise comparisons bw all ROIs. this is doing a bunch of paired t-tests 
#(after averaging across position), and then correcting with the Tukey method.
# lsmeans(fm_int, pairwise~ROI, adjust='tukey')
lsmeans(fm_int, pairwise~position, adjust='tukey')

summary(aov(err~position*ROI + Error(1/(subject*position*ROI)), data=fdat))
