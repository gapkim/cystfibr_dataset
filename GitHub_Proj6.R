##############################
### Initial Setup

knitr::opts_chunk$set(comment=NA, echo=FALSE, warning=FALSE, message=FALSE,
                      fig.align="center")
options(digits=4)

rm(list=ls())

##############################
### Get Data

library(ISwR)
data("cystfibr")

cat("\nNumber of NA's: ", sum(is.na(cystfibr)),"\n\n")

mydata = cystfibr
mydata$sex = factor(mydata$sex, levels=c(0,1), labels=c("male","female"))
str(mydata)

###################################
### EXPLORATORY DATA ANALYSIS
###################################

library(ggplot2)

# Figure 1: Relationship among age, height and weight by sex

ggplot(data=mydata) +
   geom_point(aes(x=age, y=height, size=weight, color=sex)) +
   geom_smooth(aes(x=age, y=height, color=sex), span=10, se=FALSE)

###################################
### ANALYSIS
###################################

###################################
### Principal Component Regression
###################################

# Perform Monte Carlo simulation to obtain distribution of ncomp that minimize MSEP

library(pls)
x = model.matrix(pemax ~ ., data=cystfibr)[,-1]
y = cystfibr$pemax

n_Monte = 100
ncomponent = rep(0, n_Monte)
pred.MSE_pcr = rep(NA, n_Monte)
MSEP.trainCV = matrix(NA, n_Monte, 9)

set.seed(3)
for(i in c(1:n_Monte)) {
   
   #Sampling for Monte Carlo simulation
   train = sample(c(1:nrow(cystfibr)), round(3*sqrt(nrow(cystfibr))-1))
   test = -train
   y.test = y[test]
   
   # PCR 5 fold CV on training dataset
   pcr.fit.train = pcr(pemax ~ ., data=cystfibr, subset=train,
                       scale=TRUE, validation="CV", segments=5)
   
   MSEP.train = MSEP(pcr.fit.train)
   ncomp.min = which.min(MSEP.train$val[1,,-1])
   
   MSEP.trainCV[i,] = MSEP.train$val[1,,-1] # MSE of CV error for each loop i
   ncomponent[i] = ncomp.min  # Optimal (lowest) number of components per loop i
   
   pcr.pred = predict(pcr.fit.train, x[test,], ncomp=ncomp.min)
   pred.MSE_pcr[i] = mean((pcr.pred - y.test)^2)  # Predicted MSE using test dataset per i
}

########################################################
### Get optimal number of components from training set

library(ggplot2)

freq.ncomp = table(ncomponent)

# Figure 2: Frequency of number of component selected as optimal from training dataset (100 Monte Carlo simulations)

ggplot() +
   geom_bar(aes(x=ncomponent)) +
   scale_x_continuous(labels=c(1:9), breaks=c(1:9)) +
   labs(x="Number of components in PCR", y="Frequency selected as optimal")

##############################
### Get MSEP of train dataset

library(reshape2)

# Mean of MSE CV in the training dataset
mean.MSEP.trainCV = apply(MSEP.trainCV, 2, mean)

# CV MSE of the training dataset
MSEP.trainCV_df = data.frame(MSEP.trainCV)
colnames(MSEP.trainCV_df) = c(1:9)

melt_MSEP.trainCV_df = melt(MSEP.trainCV_df, variable.name="ncomp", value.name="MSEP")

# Checkpoint (hidden) -------
freq.ncomp
mean.MSEP.trainCV
MSEP.trainCV_df
melt_MSEP.trainCV_df
#----------------------------

# Figure 3: PCR MSE of the training set from Monte Carlo CV simulation (sampling size=100)

ggplot(data=melt_MSEP.trainCV_df) +
#   geom_boxplot(aes(x=ncomp, y=MSEP)) +
   geom_jitter(aes(x=ncomp, y=MSEP, color=ncomp), width=0.25) +
   scale_y_continuous(limits=c(NA,5000)) +
   stat_summary(aes(x=ncomp, y=MSEP, group=1), 
                fun.y="mean", color="red", size=1, geom="line") +
   stat_summary(aes(x=ncomp, y=MSEP, group=1),
                fun.y="mean", color="red", size=3, geom="point") +
   scale_color_discrete(name="", breaks="") +
   labs(x="Number of components in PCR", y="Cross validated MSE of the training set")

###############################################
### Calculate prediction MSE using test dataset

index.max = which.max(freq.ncomp)  # Get the index of most frequent ncomp
id.pcr = as.integer(names(freq.ncomp)[index.max])  # Get optimal number of component

pred.pcr = pred.MSE_pcr[which(ncomponent == id.pcr)]
avg.pred.pcr = mean(pred.MSE_pcr[which(ncomponent == id.pcr)])

cat("Prediction MSE =", avg.pred.pcr, 
    "(where optimal number of PCR component is", id.pcr, ")")

index.max
names(freq.ncomp)[index.max]
length(pred.MSE_pcr)
length(pred.MSE_pcr[which(ncomponent == id.pcr)])
length(pred.MSE_pcr[which(ncomponent == names(freq.ncomp)[index.max])])
freq.ncomp

pcr.fit.full = pcr(pemax ~ ., data=cystfibr, scale=TRUE, ncomp=id.pcr)
summary(pcr.fit.full)

###########################################
### Partial Linear Squares Regression
###########################################

# Perform Monte Carlo simulation to obtain distribution of ncomp that minimize MSEP

ncomponent = rep(0, n_Monte)
pred.MSE_plsr = rep(NA, n_Monte)
MSEP.trainCV = matrix(NA, n_Monte, 9)

set.seed(3)
for(i in c(1:n_Monte)) {
   
   #Sampling for Monte Carlo simulation
   train = sample(c(1:nrow(cystfibr)), round(3*sqrt(nrow(cystfibr))-1))
   test = -train
   y.test = y[test]
   
   # plsr 5 fold CV on training dataset
   plsr.fit.train = plsr(pemax ~ ., data=cystfibr, subset=train,
                       scale=TRUE, validation="CV", segments=5)
   
   MSEP.train = MSEP(plsr.fit.train)
   ncomp.min = which.min(MSEP.train$val[1,,-1])
   
   MSEP.trainCV[i,] = MSEP.train$val[1,,-1] # MSE of CV error for each loop i
   ncomponent[i] = ncomp.min  # Optimal (lowest) number of components per loop i
   
   plsr.pred = predict(plsr.fit.train, x[test,], ncomp=ncomp.min)
   pred.MSE_plsr[i] = mean((plsr.pred - y.test)^2)  # Predicted MSE
}

########################################################
### Get optimal number of components from training set

freq.ncomp = table(ncomponent)

# Figure 4: Frequency of number of component selected as optimal from training dataset using PLSR (100 Monte Carlo simulations)

ggplot() +
   geom_bar(aes(x=ncomponent)) +
   scale_x_continuous(labels=c(1:9), breaks=c(1:9)) +
   labs(x="Number of components in PLSR", y="Frequency selected as optimal")

##############################
### Get MSEP of train dataset

# Mean of MSE CV in the training dataset
mean.MSEP.trainCV = apply(MSEP.trainCV, 2, mean)

# CV MSE of the training dataset
MSEP.trainCV_df = data.frame(MSEP.trainCV)
colnames(MSEP.trainCV_df) = c(1:9)

melt_MSEP.trainCV_df = melt(MSEP.trainCV_df, variable.name="ncomp", value.name="MSEP")

# Checkpoint (hidden) -------
freq.ncomp
mean.MSEP.trainCV
MSEP.trainCV_df
melt_MSEP.trainCV_df
#----------------------------

# Figure 5: PLSR MSE of the training set from Monte Carlo CV simulation (sampling size=100)

ggplot(data=melt_MSEP.trainCV_df) +
#   geom_boxplot(aes(x=ncomp, y=MSEP)) +
   geom_jitter(aes(x=ncomp, y=MSEP, color=ncomp), width=0.25) +
   scale_y_continuous(limits=c(NA,6000)) +
   stat_summary(aes(x=ncomp, y=MSEP, group=1), 
                fun.y="mean", color="red", size=1, geom="line") +
   stat_summary(aes(x=ncomp, y=MSEP, group=1),
                fun.y="mean", color="red", size=3, geom="point") +
   scale_color_discrete(name="", breaks="") +
   labs(x="Number of components in PLSR", y="Cross validated MSE of the training set")

###############################################
### Calculate prediction MSE using test dataset

index.max = which.max(freq.ncomp)  # Get the index of most frequent ncomp
id.plsr = as.integer(names(freq.ncomp)[index.max])  # Get optimal number of component

pred.plsr = pred.MSE_plsr[which(ncomponent == id.plsr)]
avg.pred.plsr = mean(pred.MSE_plsr[which(ncomponent == id.plsr)])

cat("Prediction MSE =", avg.pred.plsr, 
    "(where optimal number of PLSR component is", id.plsr, ")")

plsr.fit.full = plsr(pemax ~ ., data=cystfibr, scale=TRUE, ncomp=id.plsr)
summary(plsr.fit.full)

##################################################
### Best Subsets Regression (with Monte Carlo CV)
##################################################

library(leaps)

npred = ncol(cystfibr)-1  # Number of predictors

cv.errors = matrix(NA, n_Monte, npred)
train.errors = matrix(NA, n_Monte, npred)

set.seed(3)
for(j in c(1:n_Monte)) {

   #Sampling for Monte Carlo simulation
   train = sample(c(1:nrow(cystfibr)), round(3*sqrt(nrow(cystfibr))-1))
   test = -train

   best.subset = regsubsets(pemax ~ ., data=cystfibr[train,], nvmax=npred)
   train.errors[j,] = summary(best.subset)$rss / best.subset$nn
   
   x.test = model.matrix(pemax ~ ., data=cystfibr[test,])
   y.test = y[test]
   
   for(i in 1:npred) {
      coefi = coef(best.subset, id=i)
      pred = x.test[, names(coefi)] %*% coefi
      cv.errors[j,i] = mean((y.test - pred)^2)
   }
}


mean.cv.errors = apply(cv.errors, 2, mean)
mean.train.errors = apply(train.errors, 2, mean)

# Best subset optimal number of predictors
index.best = which.min(mean.cv.errors)  

# Test MSE of Best Subsets
pred.BS = cv.errors[,index.best]
test.MSE_bestsub = mean.cv.errors[index.best]

# Figure 6: Test and training MSE using Best subsets regression with Monte Carlo simulation (100 samplings)

ggplot() +
   geom_point(aes(x=c(1:npred),y=mean.cv.errors, color="Test MSE"), size=3) +
   geom_line(aes(x=c(1:npred),y=mean.cv.errors, color="Test MSE"), linetype="solid") +
   geom_point(aes(x=c(1:npred),y=mean.train.errors, color="Training MSE"), size=3) +
   geom_line(aes(x=c(1:npred),y=mean.train.errors, color="Training MSE"), linetype="solid") +
   scale_x_continuous(breaks=c(1:npred)) +
   labs(x="Number of Predictors", y="Mean Squared Error", color=c("Test/Training"))

#########################################################################
### Best subsets fit to full dataset with optimal number of predictors

reg.best = regsubsets(pemax ~ ., data=cystfibr, nvmax=npred)
coef(reg.best, index.best)

reg.best.lm = lm(pemax ~ weight, data=cystfibr)
summary_reg.best.lm = summary.lm(reg.best.lm)
summary_reg.best.lm$coef

##################################################
### Ridge Regression (with Monte Carlo CV)
##################################################

library(glmnet)

grid = exp(seq(8, -6, length=100))
n_Monte = 100
cv.errors = rep(NA, n_Monte)
bestlam_ridge = rep(NA, n_Monte)

set.seed(3)
for(j in c(1:n_Monte)) {

   #Sampling for Monte Carlo simulation
   train = sample(c(1:nrow(cystfibr)), round(3*sqrt(nrow(cystfibr))-1))
   test = -train
   y.test = y[test]
   
   ridge.cv.out = cv.glmnet(x[train,], y[train], alpha=0, nfolds=5, lambda=grid)
   bestlam_min_ridge = ridge.cv.out$lambda.min

   bestlam_ridge[j] = bestlam_min_ridge
   
   ridge.cv.pred = predict(ridge.cv.out, s=bestlam_ridge[j], newx=x[test,],
                           exact=T, x=x[train,], y=y[train])
   cv.errors[j] = mean((ridge.cv.pred - y.test)^2)
}

cv.output_ridge = data.frame("folds"=c(1:n_Monte), "lambda"=bestlam_ridge, "Test_MSE"=cv.errors)

avg.bestlam_ridge = mean(bestlam_ridge)

# Test MSE from Ridge Regression with lambda.min
pred.ridge = cv.errors
test.MSE_ridge.min = mean(cv.errors)

cat("Optimal lambda min =", avg.bestlam_ridge)

###########################################
### Fit Ridge Regression to full dataset

out_ridge = glmnet(x, y, alpha=0, lambda=grid)
coef.bestlam_ridge = coef(out_ridge, s=avg.bestlam_ridge, exact=T, x=x, y=y)[,1]

#Figure.7, fig.ap="Figure 15: Coefficients vs. log(lambda) for Ridge Regression with full dataset
plot(out_ridge, xvar="lambda", label=T,
     sub=paste("(Vertical line at lambda_min =", round(avg.bestlam_ridge, 1), ")"))
abline(v=log(avg.bestlam_ridge ))

###################################
### Lasso Rgression
###################################

grid = exp(seq(4, -5, length=100))

cv.errors = rep(NA, n_Monte)
bestlam_lasso = rep(NA, n_Monte)

set.seed(3)
for(j in c(1:n_Monte)) {

   #Sampling for Monte Carlo simulation
   train = sample(c(1:nrow(x)), round(3*sqrt(nrow(x))-1))
   test = -(train)
   y.test = y[test]
   
   lasso.cv.out = cv.glmnet(x[train,], y[train], alpha=1, nfolds=5, lambda=grid)
   bestlam_min_lasso = lasso.cv.out$lambda.min
   
   bestlam_lasso[j] = bestlam_min_lasso
   lasso.cv.pred = predict(lasso.cv.out, s=bestlam_lasso[j], newx=x[test,],
                           exact=T, x=x[train,], y=y[train])
   cv.errors[j] = mean((lasso.cv.pred - y.test)^2)
}

cv.output_lasso = data.frame("folds"=c(1:n_Monte), "lambda"=bestlam_lasso, "Test_MSE"=cv.errors)

avg.bestlam_lasso = mean(bestlam_lasso)
avg.cv.errors = mean(cv.errors)

# Test MSE from Lasso Regression with lambda.min
pred.lasso = cv.errors
test.MSE_lasso.min = avg.cv.errors

cat("Optimal lambda min =", avg.bestlam_lasso)

###########################################
### Fit Lasso Regression to full dataset

out_lasso = glmnet(x, y, alpha=1, lambda=grid)
coef.bestlam_lasso = coef(out_lasso, s=avg.bestlam_lasso, exact=T, x=x, y=y)[,1]

# Figure 8: Coefficients vs. log(lambda) for Ridge Regression with full dataset
plot(out_lasso, xvar="lambda", label=T,
     sub=paste("(Vertical line at lambda_min =", round(avg.bestlam_lasso, 1), ")"))
abline(v=log(avg.bestlam_lasso))

############################################
### DISCUSSIONS
############################################

MSE.all = list(pred.pcr, pred.plsr, pred.BS, pred.ridge, pred.lasso)
names(MSE.all) = c("PCR", "PLSR", "BS", "Ridge", "Lasso")

melt_test.MSE = melt(MSE.all)
colnames(melt_test.MSE) = c("Test.MSE", "Method")

melt_test.MSE = melt_test.MSE[c(2,1)]
melt_test.MSE$Method = factor(melt_test.MSE$Method, levels=c("PCR", "PLSR", "BS", "Ridge", "Lasso"))

# Figure 9: Comparison of Test MSE of different model selection

ggplot(data=melt_test.MSE, aes(x=Method, y=Test.MSE)) +
   coord_cartesian(ylim = c(400, 2000)) +
   stat_summary(fun.y=mean, fun.ymax=max, fun.ymin=min, geom="point", size=3) +
   stat_summary(fun.ymax=function(z) {quantile(z,0.75)}, 
                fun.ymin=function(z) {quantile(z,0.25)}, 
                geom="errorbar", width=0.2) +
   labs(x="Model selection method", y="Test MSE")

library(knitr)

# Table of model coefficients
coef.table = matrix(NA, 10, 3, dimnames=list(reg.best$xnames, c("BS", "Ridge", "Lasso")))

coef.table[names(coef(reg.best, index.best)), 1] = coef(reg.best, index.best)
coef.table[names(coef.bestlam_ridge), 2] = coef.bestlam_ridge
coef.table[names(coef.bestlam_lasso), 3] = coef.bestlam_lasso


coef.table[coef.table == 0] = NA
n.comp = apply(!is.na(coef.table), 2, sum)-1

coef.table = rbind(coef.table, n.comp)
rownames(coef.table)[11] = "Num.Comp"

# Table 1: Number of coefficients and values selected by the Best subsets, Ridge regression and Lasso regression

kable(coef.table, format="html", table.attr = "style='width:40%;'",
      caption="Table 1: Number of coefficients and values selected by the Best subsets, Ridge regression and Lasso regression")

coef.table1 = cbind(pcr.fit.full$loadings[, 1], plsr.fit.full$loadings[, 1])
colnames(coef.table1) = c("PCR-Comp1", "PLSR-Comp1")

# Table 2: Loadings of the first component for PCR and PLSR
kable(coef.table1, format="html", table.attr="style='width:40%;'",
      caption="Table 2: Loadings of the first component for PCR and PLSR")


temp.table = coef.table[-c(1, 11),]
temp.table[is.na(temp.table)] = 0

ss = sqrt(apply(temp.table^2, 2, sum)) # Magnitude of each method

temp.table[,1] = temp.table[,1] / ss[1]
temp.table[,2] = temp.table[,2] / ss[2]
temp.table[,3] = temp.table[,3] / ss[3]

# Table of spectral components
spectral.table = cbind(coef.table1, temp.table)

melt_spectral.table = melt(spectral.table)
colnames(melt_spectral.table) = c("Predictors", "Method", "Loadings")

# Figure 10: Spectral plot of coefficients (loadings) for different seletion methods
ggplot(data=melt_spectral.table, aes(group=Method, fill=Method)) +
   geom_col(aes(x=Predictors, y=Loadings), position="dodge") +
   labs(y="Loadings (normalized coefficients)")


fitted.table = cbind(pcr.fit.full$fitted.values[,,1],
      plsr.fit.full$fitted.values[,,1],
      reg.best.lm$fitted.values,
      predict(out_ridge, s=avg.bestlam_ridge, newx=x, exact=T, x=x, y=y),
      predict(out_lasso, s=avg.bestlam_lasso, newx=x, exact=T, x=x, y=y))
colnames(fitted.table) = c("PCR", "PLSR", "BS", "Ridge", "Lasso")

melt_fitted.table = melt(fitted.table)
colnames(melt_fitted.table) = c("Data", "Method", "Fitted.value")
melt_fitted.table = cbind(melt_fitted.table, cystfibr$pemax)
colnames(melt_fitted.table)[4] = "Obs.value"

# Figure 11: Comparison of fitted vs. observed values

ggplot(data=melt_fitted.table) +
   geom_point(aes(x=Obs.value, y=Fitted.value, color=Method), size=3) +
   geom_abline(slope=1, size=1, linetype=3) +
   annotate("text", x=142, y=155, label="Slope = 1") +
   geom_smooth(aes(x=Obs.value, y=Fitted.value, color=Method), method="lm", se=F) +
   labs(x="Observed values (pemax)", y="Fitted values (pemax)")

mean.pemax = mean(cystfibr$pemax)
TSS = sum((cystfibr$pemax - mean.pemax)^2)

RSS.pcr = sum((pcr.fit.full$fitted.values[,,1] - cystfibr$pemax)^2)
RSS.plsr = sum((plsr.fit.full$fitted.values[,,1] - cystfibr$pemax)^2)
RSS.BS = sum((reg.best.lm$fitted.values - cystfibr$pemax)^2)
RSS.ridge = sum((predict(out_ridge, s=avg.bestlam_ridge, newx=x, exact=T, x=x, y=y) - cystfibr$pemax)^2)
RSS.lasso = sum((predict(out_lasso, s=avg.bestlam_lasso, newx=x, exact=T, x=x, y=y) - cystfibr$pemax)^2)

R2 = cbind(1-RSS.pcr/TSS, 1-RSS.plsr/TSS, 1-RSS.BS/TSS, 1-RSS.ridge/TSS, 1-RSS.lasso/TSS)
colnames(R2) = c("PCR", "PLSR", "BS", "Ridge", "Lasso")

melt.R2 = melt(R2)[,-1]
colnames(melt.R2) = c("Method", "R2")

# Figure 12: R2 values of the model selection method with the training dataset

ggplot(data=melt.R2) +
   geom_col(aes(x=Method, y=R2*100)) +
   labs(x="Model selection method", y="R2 (%)")
