###############################
# AutoClaim dataset
###############################

#import package
library(MSTweedie)

#load data
data(AutoClaim)

#display head of dataset
head(AutoClaim)

#classify the policies by REVOLKED and whether there was a claim or not
table(AutoClaim$REVOLKED, AutoClaim$CLM_AMT5 > 0)

##############################
# MSTweedie
##############################

#import package
library(MSTweedie)

#load data
data(AutoClaim)

#fit the MSTweedie model with L1/Linf regularization
# y=1 sets CLM_AMT5 as the response, source=4 sets REVOLKED as the source index
fit <- MSTweedie(x = AutoClaim, y=1, source=4, reg='Linf')

##############################
# coef.MSTweedie
##############################

#import package
library(MSTweedie)

#load data
data(AutoClaim)

# performs 10-folds CV with L1/Linf regularization
fit <- MSTweedie(x = AutoClaim, y=1, source=4, reg='Linf')

# extract coefficients at 34th lambda
coef.MSTweedie(fit, s=34:36)

##############################
# predict.MSTweedie
##############################

#import package
library(MSTweedie)

#load data
data(AutoClaim)

# performs 10-folds CV with L1/Linf regularization
fit <- MSTweedie(x = AutoClaim, y=1, source=4, reg='Linf')

# predict first source at 34th lambda
head(predict.MSTweedie(fit, s=34L)[[1]])

##############################
# print.MSTweedie
##############################

#import package
library(MSTweedie)

#load data
data(AutoClaim)

# performs 10-folds CV with L1/Linf regularization
fit <- MSTweedie(x = AutoClaim, y=1, source=4, reg='Linf')

# prints number of selected variables along solution path
print.MSTweedie(fit)


##############################
# plot.MSTweedie
##############################

#import package
library(MSTweedie)

#load data
data(AutoClaim)

# performs 10-folds CV with L1/Linf regularization
fit <- MSTweedie(x = AutoClaim, y=1, source=4, reg='Linf')

# plot solution path of the norm of the coefficients
plot.MSTweedie(fit, type.coef='norm')

# plot solution path of the the coefficients
par(mfrow=c(2,1))
plot.MSTweedie(fit, type.coef='coef')

##############################
# kkt.check
##############################

#import package
library(MSTweedie)

#load data
data(AutoClaim)

# performs 10-folds CV with L1/Linf regularization
fit <- MSTweedie(x = AutoClaim, y=1, source=4, reg='Linf')

# plot the two kkt conditions
par(mfrow=c(2,1))
kkt.check(fit)

##############################
# deviance.MSTweedie
##############################

#import package
library(MSTweedie)

#load data
data(AutoClaim)

# performs 10-folds CV with L1/Linf regularization
fit <- MSTweedie(x = AutoClaim, y=1, source=4, reg='Linf')

# copmute deviance of the model at 34th lambda
deviance.MSTweedie(y= fit$y,
                   mu = predict.MSTweedie(fit, s=34L))

##############################
# cv.MSTweedie
##############################

#import package
library(MSTweedie)

#load data
data(AutoClaim)

# performs 10-folds CV with L1/Linf regularization
cv<-cv.MSTweedie(x = AutoClaim, y=1, source=4, reg='Linf')


##############################
# coef.cv.MSTweedie
##############################

#import package
library(MSTweedie)

#load data
data(AutoClaim)

# performs 10-folds CV with L1/Linf regularization
cv <- cv.MSTweedie(x = AutoClaim, y=1, source=4, reg='Linf')

# extract coefficients at lambda.1se
coef.cv.MSTweedie(cv)

##############################
# predict.cv.MSTweedie
##############################

#import package
library(MSTweedie)

#load data
data(AutoClaim)

# performs 10-folds CV with L1/Linf regularization
cv <- cv.MSTweedie(x = AutoClaim, y=1, source=4, reg='Linf')

# extract coefficients at lambda.1se
head(predict.cv.MSTweedie(cv)[[1]])

##############################
# plot.cv.MSTweedie
##############################

#import package
library(MSTweedie)

#load data
data(AutoClaim)

# performs 10-folds CV with L1/Linf regularization
cv <- cv.MSTweedie(x = AutoClaim, y=1, source=4, reg='Linf')

# plot CV deviance mean and std. err., and lambda.min, lambda.1se
plot.cv.MSTweedie(cv)
