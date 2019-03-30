setwd("D:R_scripts")
library(lme4)
library(tidyverse)
require(sjPlot)

# Read in the classifier accuracies.
fn = '/mnt/neurocube/local/serenceslab/maggie/mFiles/IEMdepth/Mixed models/dprime_classifier_tbl.txt'

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

##################################################################
# Now the reconstruction error

fn2 = "/mnt/neurocube/local/serenceslab/maggie/mFiles/IEMdepth/Mixed models/fitErr2way.csv"
fdat = read_csv(fn2)
fdat = rename(fdat, V1_0=V1, V2_0=V2, V3_0=V3, V4_0=V4, V3A_0=V3A, V3B_0=V3B,
              IPS0_0=IPS0, LO1_0=LO1, LO2_0=LO2)

fitdat = gather(fdat, Cond, fitErr, V1_0:LO2_5, factor_key=TRUE)
fitdat2 = separate(fitdat, Cond, c("ROI","Position"), sep="_")
fitdat2$fitErr = abs(fitdat2$fitErr)
fitdat2$Position = factor(fitdat2$Position)

fm0 = lmer(fitErr~1 + (1|Subj), data=fitdat2, REML=FALSE)
fm1 = lmer(fitErr~Position + (1|Subj), data=fitdat2, REML=FALSE)
fm2 = lmer(fitErr~Position + ROI + (1|Subj) ,  data=fitdat2, REML=FALSE)
fm3 = lmer(fitErr~Position + ROI + (1|Subj) + (1|ROI:Subj),  data=fitdat2, REML=FALSE)
fm4 = lmer(fitErr~Position + ROI + (1|Subj) + (1|ROI:Subj) + (1|Position:Subj), data=fitdat2, REML=FALSE)
fm5 = lmer(fitErr~Position*ROI + (1|Subj)+ (1|Position:Subj) + (1|ROI:Subj), data=fitdat2, REML=FALSE)
anova(fm0,fm1,fm2,fm3,fm4,fm5)
