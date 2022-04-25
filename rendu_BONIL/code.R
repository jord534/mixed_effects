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
VL_avg = aggregate(data$VL, list(data$TD), FUN=mean)
plot(VL_avg, type = "b", main = "Average trajectory for Viral Load wrt to Time", xlab = "weeks", ylab = "VL")

CD4_avg = aggregate(data$CD4, list(data$TD), FUN=mean)
plot(CD4_avg, type = "b", main = "Average trajectory for CD4 concentration wrt to Time", xlab = "weeks", ylab = "cc/m3")


# 1.3 #

sub_grp = subset(data, data$RAN_GRP != 2)
sub_grp$TD2 = sub_grp$TD^2
library(lme4)

# trying to find the best model :
m1  <- lm(sub_grp$CD4 ~ sub_grp$TD2 , data=sub_grp) 
m2  <- lm(sub_grp$CD4 ~ sub_grp$TD + sub_grp$TD2 , data=sub_grp) 
m3  <- lm(sub_grp$CD4 ~ 1 + sub_grp$TD2:sub_grp$TD , data=sub_grp) 
m4  <- lm(sub_grp$CD4 ~ sub_grp$TD + sub_grp$TD2:sub_grp$TD , data=sub_grp) 
m5  <- lmer(sub_grp$CD4 ~ sub_grp$TD2 + (1|sub_grp$NUM_PAT) , data=sub_grp) 
m6  <- lmer(sub_grp$CD4 ~ sub_grp$TD + sub_grp$TD2 + (1|sub_grp$NUM_PAT) , data=sub_grp) 
m7  <- lmer(sub_grp$CD4 ~ 1 + sub_grp$TD2:sub_grp$TD + (1|sub_grp$NUM_PAT) , data=sub_grp) 
m8  <- lmer(sub_grp$CD4 ~ sub_grp$TD + sub_grp$TD2:sub_grp$TD + (1|sub_grp$NUM_PAT) , data=sub_grp) 
m9  <- lmer(sub_grp$CD4 ~ sub_grp$TD2 + (-1+sub_grp$TD2|sub_grp$NUM_PAT) , data=sub_grp) 
m10 <- lmer(sub_grp$CD4 ~ sub_grp$TD + sub_grp$TD2 + (-1+sub_grp$TD2|sub_grp$NUM_PAT) , data=sub_grp) 
m11 <- lmer(sub_grp$CD4 ~ 1 + sub_grp$TD2:sub_grp$TD + (-1+sub_grp$TD2|sub_grp$NUM_PAT) , data=sub_grp) 
m12 <- lmer(sub_grp$CD4 ~ sub_grp$TD + sub_grp$TD2:sub_grp$TD + (-1+sub_grp$TD2|sub_grp$NUM_PAT) , data=sub_grp) 
m13 <- lmer(sub_grp$CD4 ~ sub_grp$TD2 + (sub_grp$TD2|sub_grp$NUM_PAT) , data=sub_grp) 
m14 <- lmer(sub_grp$CD4 ~ sub_grp$TD + sub_grp$TD2 + (sub_grp$TD2|sub_grp$NUM_PAT) , data=sub_grp) 
m15 <- lmer(sub_grp$CD4 ~ 1 + sub_grp$TD2:sub_grp$TD + (sub_grp$TD2|sub_grp$NUM_PAT) , data=sub_grp) 
m16 <- lmer(sub_grp$CD4 ~ sub_grp$TD + sub_grp$TD2:sub_grp$TD + (sub_grp$TD2|sub_grp$NUM_PAT) , data=sub_grp) 
# Comparing BIC and AIC
AICtable = AIC(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16)
AICtable[which(AICtable$AIC == min(AICtable$AIC)),]
BICtable<-BIC(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16)
BICtable[which(BICtable$BIC==min(BICtable$BIC)),]

# m6 is designated by BIC as the best model
summary(m6)
sub_grp$pred.final <- fitted(m6)
library(ggplot2)
sub_grp_1 = subset(sub_grp, sub_grp$RAN_GRP == 1)
sub_grp_3 = subset(sub_grp, sub_grp$RAN_GRP == 3)

# for the AZT+3TC sub_grp_1
ggplot(data=sub_grp_1) + geom_point(aes(x=sub_grp_1$TD2,y=sub_grp_1$CD4), color="red", size=3) + 
  geom_line(aes(x=sub_grp_1$TD2,y=pred.final)) + facet_wrap(~sub_grp_1$NUM_PAT, ncol=5) 

# for the D4T+DDL sub_grp_3
ggplot(data=sub_grp_3) + geom_point(aes(x=sub_grp_3$TD2,y=sub_grp_3$CD4), color="red", size=3) + 
  geom_line(aes(x=sub_grp_3$TD2,y=pred.final)) + facet_wrap(~sub_grp_3$NUM_PAT, ncol=5) 

confint(m6)
confint(m6, method="boot")
library(lattice)
d = dotplot(ranef(m6, condVar=TRUE), strip=FALSE)
print(d[[1]])

# let's see which of the two groups has the better CD4 count reaction
# to treatment

par(mfrow=c(1,1))
# using aggregate(df$col_to_aggregate, list(df$col_to_group_by), FUN=mean)  ?
grp1_avg = aggregate(sub_grp_1$CD4, list(sub_grp_1$TD), FUN=mean)
plot(grp1_avg, type = "b", main = " CD4 count trajectories for grp 1 AZT+3TC", xlab = "weeks", ylab = "CD4 count")
grp1_avg_pred = aggregate(sub_grp_1$pred.final, list(sub_grp_1$TD), FUN=mean)
lines(grp1_avg_pred, type = "b", col = "red")

grp_3_avg = aggregate(sub_grp_3$CD4, list(sub_grp_3$TD), FUN = mean)
plot(grp_3_avg, type = "b", main = "CD4 count trajctories for grp 3 D4T+DDL", xlab = "weeks", ylab = "CD4 count")
grp_3_avg_pred = aggregate(sub_grp_3$pred.final, list(sub_grp_3$TD), FUN=mean)
lines(grp_3_avg_pred, type = "b", col = "blue")

# interpretation : CD4 count reaches a greater max value for the grp 1

#### PART2 ####

data = read.table("final_PK_data_exam2022.txt" , header = TRUE)
summary(data)
length(unique(data$ID)) # there are 20 unique patients
# we get multiple concentration levels at various time intervals
par(mfrow=c(1,1))
plot(data$TIME, data$CONC_ngmL, main = "Concentration wrt to time", xlab = "hours", ylab = "PK")

# 1.4 #

