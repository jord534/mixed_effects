#setwd("/Users/jordanbonil/Documents/Mixed\ Effects\ Models/projet/")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# this is a test
#### PART 1 ###

# 1.1 #
#importing data :
data = read.table("arv_hiv.txt", header = TRUE)

# removing NA variables :
data = na.omit(data)
# normalising data :

par(mfrow=c(1,2))
qqnorm(data$CD4, main = "CD4 before normalisation")
data$CD4 = data$CD4^0.25
qqnorm(data$CD4, main = "CD4 after ^0.25")
# the normalisation for VL is unchanged

par(mfrow=c(1,2))
qqnorm(data$VL, main = "VL before normalisation")
data$VL = log(data$VL)/log(10)
qqnorm(data$VL, main = "qqplot for log10(VL)")
#the normalisation for VL is greatly improved
# but wee need to remove all points that are less than the 3rd quantile
VL_quantiles = qqnorm(data$VL, plot.it = FALSE)$x
mask = VL_quantiles > -1 # data points to keep
data = data[mask,]
par(mfrow=c(1,1))
qqnorm(data$VL, main = "qqplot for log10(VL) w/o non normalised pt")
qqline(data$VL, datax = FALSE, distribution = qnorm,
       probs = c(0.25, 0.75), qtype = 7)
data$dd4Tddl = as.numeric(data$RAN_GRP != 1)


# 1.2 #
summary(data)
# We obtain VL and CD4 values on 7 days for each patients
# and also the time difference between values


# in order to plot the average trajectory,
# on times where there is only one val, we print this val
# on times where there are multiple vals, we average the vals
par(mfrow=c(1,1))
# using aggregate(df$col_to_aggregate, list(df$col_to_group_by), FUN=mean)  ?
VL_avg = aggregate(data$VL, list(data$TIME), FUN=mean)
plot(VL_avg, type = "b", main = "Average trajectory for Viral Load wrt to Time", xlab = "days", ylab = "VL")

CD4_avg = aggregate(data$CD4, list(data$TIME), FUN=mean)
plot(CD4_avg, type = "b", main = "Average trajectory for CD4 concentration wrt to Time", xlab = "days", ylab = "cc/m3")

library(zoo)
plot(rollmean(CD4_avg$x, k=73), main = "Rolling Average of mean CD4 trajectory", xlab = "days")

plot(rollmean(VL_avg$x, k=73), main = "Rolling Average of mean VL trajectory", xlab = "days")

# 1.3 #

azt3tc_grp = subset(data, data$RAN_GRP != 2)
azt3tc_grp$TD2 = azt3tc_grp$TD^2
library(lme4)
cd4_mixed_RAN_GRP = lmer( azt3tc_grp$CD4 ~ azt3tc_grp$TD +  (1 | azt3tc_grp$RAN_GRP), data = azt3tc_grp)
summary(cd4_mixed_RAN_GRP)
confint(cd4_mixed_RAN_GRP)

cd4_mixed_dd4Tddl = lmer( azt3tc_grp$CD4 ~ azt3tc_grp$TD +  (1 | azt3tc_grp$dd4Tddl), data = azt3tc_grp)
summary(cd4_mixed_dd4Tddl)
