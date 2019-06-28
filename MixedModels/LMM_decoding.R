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
lsmeans(lm2, pairwise~ROI, adjust='tukey')

