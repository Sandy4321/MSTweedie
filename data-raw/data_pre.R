################################
##import data and preprocess
################################
#import and drop unused variables
library(cplm)
library(rrcovNA)
names(AutoClaim)

data <- subset(AutoClaim, select = -c(1,2,3,5,16,19,20,26,27))
#imputation kNN
dat <- impSeq(data)
dat <- as.data.frame(dat)

#scale
n <- nrow(dat)
dat$CLM_AMT5 <- dat[,"CLM_AMT5"]/1000
dat$BLUEBOOK <- log(dat[,"BLUEBOOK"])

table(dat$tasks)
# create dummy categorical variables
attach(dat)
##CAR_TYPE
ng <- max(CAR_TYPE)
CAR_TYPE_dummy <- matrix(NA,n,ng-1)
name <- rep(NA,ng-1)
for (i in 2:ng) {
   CAR_TYPE_dummy[,i-1] <- as.numeric(CAR_TYPE==i)
   name[i-1] <- paste("CAR_TYPE",i,sep="_")
}
colnames(CAR_TYPE_dummy) <- name
head(CAR_TYPE_dummy)
##JOBCLASS
ng <- max(JOBCLASS)
name <- rep(NA,ng-1) # the category "Unknown" has no observations
JOBCLASS_dummy <- matrix(NA,n,ng-1)
for (i in 2:ng) {
   JOBCLASS_dummy[,i-1] <- as.numeric(JOBCLASS==i)
   name[i-1] <- paste("JOBCLASS",i,sep="_")
}
colnames(JOBCLASS_dummy) <- name
head(JOBCLASS_dummy)
##MAX_EDUC
ng <- max(MAX_EDUC)
name <- rep(NA,ng-1)
MAX_EDUC_dummy <- matrix(NA,n,ng-1)
for (i in 2:ng) {
   MAX_EDUC_dummy[,i-1] <- as.numeric(MAX_EDUC==i)
   name[i-1] <- paste("MAX_EDUC",i,sep="_")
}
colnames(MAX_EDUC_dummy) <- name
head(MAX_EDUC_dummy)
##AGE_CAT
#age as categories
AGE_CAT <- cut(AGE, breaks = c(0,30,40,50,60,100), labels = F)

ng <- max(AGE_CAT)
name <- rep(NA,ng-1)
AGE_CAT_dummy <- matrix(NA,n,ng-1)
for (i in 2:ng) {
   AGE_CAT_dummy[,i-1] <- as.numeric(AGE_CAT==i)
   name[i-1] <- paste("AGE_CAT",i,sep="_")
}
colnames(AGE_CAT_dummy) <- name
head(AGE_CAT_dummy)


##recreate data frame
dat <- cbind(dat,CAR_TYPE_dummy,JOBCLASS_dummy,MAX_EDUC_dummy,AGE_CAT_dummy)
AutoClaim <- subset(dat, select = -c(6,8,12,17,18,20))
