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
qqnorm(data$CD4, main = "CD4 after *0.25")
# the normalisation for VL is unchanged

par(mfrow=c(1,2))
qqnorm(data$VL, main = "VL before normalisation")
data$VL = log(data$VL)/log(10)
qqnorm(data$VL, main = "qqplot for log10(VL)")
#the normalisation for VL is greatly improved
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

