# Brandon Wong Thesis Objective 3
# Exposure of interest: temperature
# Outcome of interest: stroke-related ED visits
# Potential confounder: humidity

# 1 Data preparation

# Upload packages
library ("dlnm")
library ("mixmeta")
library ('splines')
library ('survival')
library ("boot")
library ("lme4")
library ("dplyr")
library ("readxl")
library ("timeDate")
library ("lubridate")

# Set working directory
setwd("/Users/Brandon W/Documents/Western/THESIS/Analysis/FINAL ANALYSIS")

# Load dataset
df <- read.csv("stroke_pm25_temp_imputed.csv", stringsAsFactors = FALSE)

df  <- df %>% rename(stroke_type = stoke_type)                 # fix typo
keep_vars <- c(
  "case","match","fis_yr","visit_d","fsa","age_yr","male","stroke_type",
  paste0("temp_avg_lag",0:7), paste0("temp_max_lag",0:7),
  paste0("humidity_lag",0:7)
)
thesis_data <- df[, keep_vars]

## Attach health-region codes
fsa_hr <- read_xlsx("Canada_FSA_to_HR18_Correspondence.xlsx") %>% 
  rename(fsa = FSA)
final  <- merge(thesis_data, fsa_hr[,c("fsa","HR_UID")], by = "fsa")

## Add visit date and holiday flag
final$visit_date <- dmy(final$visit_d)
holiday_dates    <- as.Date(holidayTSX(unique(year(final$visit_date))))
final$stat_hol   <- as.integer(final$visit_date %in% holiday_dates)

# Add categorical age variable
final$age_cat <- cut(final$age_yr,  breaks = c(0, 65, 75, Inf),
                     labels = c(1, 2, 3),
                     right = TRUE, # Includes the rightmost break value in the interval
                     include.lowest = TRUE) # Includes the lowest value in the first interval

# Remove rows where age_cat is NA
final <- final[!is.na(final$age_cat), ]
# Merge HR_UID with age_cat
final$hr_uid_age_cat <- paste(final$HR_UID,
                              final$age_cat,
                              sep = "_")

colnames(final)


# ARRANGE THE DATA AS A LIST OF DATASETS
regions <- as.character(unique(final$hr_uid_age_cat))
data1list <- lapply(regions, function(region) final[final$hr_uid_age_cat==region,])
names(data1list) <- regions
m <- length(data1list)
m
regions

# METADATA FOR LOCATIONS
cities <- data.frame(
  city = regions)

# CREATE  META-PREDICTORS
# coerce age_cat to numeric 1,2,3
agemean <- sapply(data1list, function(x) {
  mean(as.numeric(as.character(x$age_cat)), na.rm = TRUE)
})

# HR_UID is already numeric—or if not, coerce it too
hruidmean <- sapply(data1list, function(x) {
  mean(as.numeric(as.character(x$HR_UID)), na.rm = TRUE)
})


# METADATA FOR LOCATIONS 

# IDENTIFY THE MIN AND MAX (RANGE) FOR EACH  VARIABLE, PER region 
# WE USE Lag0 (INDEX DATE) FOR SIMPLICITY SINCE WE ASSUME IT WON'T CHANGE CONSIDERABLY COMPARED 
# TO THE OTHER Lags. So you will need to do so in your dataset. 

# B - TEMPERATURE & humidity
ranges_temp <- t(sapply(data1list,function(x) range(x$temp_avg_lag0 ,na.rm=T)))
ranges_hum <- t(sapply(data1list,function(x) range(x$humidity_lag0,na.rm=T)))


bound <- colMeans(ranges_temp)
exposure <- "Temperature"


# Define the centering point, type and degrees of freedom:

cen <- 0
type <- "ns"
df <- 3

# Temperature is non-linear
cen_temp <- 20
type_temp <- "ns"
df_temp <- 3


#View(final)

#Change lag 0 to 6
temp <- as.matrix(final[,rep(9:15, each=1)])   
#temp[1:15,1:7]
templist <- lapply(regions, function(region) temp[final$hr_uid_age_cat==region,])
names(templist) <- regions

humid <- as.matrix(final[,rep(25:31, each=1)]) 
#humid[1:15,1:7]
humidlist <- lapply(regions, function(region) humid[final$hr_uid_age_cat==region,])
names(humidlist) <- regions

#lag space specification#
lag6 <- c(0,6)
lagknots6 <- logknots(6,df=3,int=T)
arglag6 <- list(fun="ns",knots=lagknots6)

#unconstrained lag#
arglag <- list(fun="integer")


# #lag 5#
# lag5 <- c(0,5)
# lagknots5 <- logknots(5,df=3,int=T)
# arglag5 <- list(fun="ns",knots=lagknots5)

# lag 7
# lag7 <- c(0,7)
# lagknots7 <- logknots(7,df=3,int=T)
# arglag7 <- list(fun="ns",knots=lagknots7)

#lag 1# Changed to 1
#lag1 <- c(0,1)
#lagknots1 <- logknots(1,df=3,int=T)
#arglag1 <- list(fun= "strata")

# Specification of exposure function
varfun <- "ns"
vardegree <- NULL
varper <- c(50,90) 

argvar3 <- list(fun="ns",cen=0)

# 2 Empty objects to store coeff-vcov and other city-specific info
coefall <- matrix(NA,m,3,dimnames=list(regions,paste("spl",seq(3),sep=""))) #3 because of 2 knots and 1 intercept#

vcovall <- vector("list", m)
names(vcovall) <- regions

# PREDICTOR SPECIFIC SUMMARIES
coef_spec <- matrix(NA,m,3,dimnames=list(regions,paste("spl",seq(3),sep=""))) #3 df for the lag space#
vcov_spec <- vector("list",m) 
names(vcov_spec) <- regions 



#------------------ Total stroke ----------------------------#

# Warning for prediction beyond boundaries suppressed
options(warn=-1)

system.time(
  for(i in seq(m)) {
    
    # PRINT ITERATION
    cat(i,"")
    
    # LOAD
    data <- data1list[[i]]
    
    # DEFINE THE CROSSBASIS
    argvar1 <- list(
      fun = varfun,
      knots = quantile(templist[[i]], varper / 100, na.rm = T),
      Bound = range(templist[[i]], na.rm = T) )
    
    
    # CREATE THE CROSSBASIS (Adjusting factors)
    cb.temp <- crossbasis(templist[[i]], lag=lag6, argvar=argvar1, arglag=arglag6)
    cb.humid <- crossbasis(humidlist[[i]], lag=lag6, argvar=argvar3, arglag=arglag6)
    
    
    # FIT THE MODEL (conditional logistic regression model here)
    
    model <- clogit(data$case ~  cb.temp   + cb.humid + strata(data$match))
    
    
    # REDUCTION TO OVERALL CUMULATIVE
    redall <- crossreduce(cb.temp,model,cen =mean(templist[[i]],na.rm = T))
    coefall[i,] <- coef(redall)
    vcovall[[i]] <- vcov(redall)
    
  })

# RESET WARNING
options(warn=0)

model


method <- "reml"

##############################################################################
####### - 2nd STAGE ANALYSIS
##############################################################################


###### MULTIVARIATE META-REGRESSION with BLUP #########

mvmeta <- mixmeta(coefall~ 1 + as.factor(agemean), random= ~ 1|hruidmean/agemean  ,vcovall,method="reml")
summary(mvmeta)


summary(mvmeta)

# FUNCTION FOR COMPUTING THE P-VALUE OF A WALD TEST 
fwald <- function(model,var) { 
  ind <- grep(var,names(coef(model))) 
  coef <- coef(model)[ind] 
  vcov <- vcov(model)[ind,ind] 
  waldstat <- coef%*%solve(vcov)%*%coef 
  df <- length(coef) 
  return(1-pchisq(waldstat,df)) 
}

# TEST THE EFFECTS
fwald(mvmeta,"agemean")
#fwald(mvmeta,"hruidmean")
#META PREDICTORS ARE NOT SIGNIFICANT


mvmeta <- mixmeta(coefall~ 1 + agemean, random= ~ 1|hruidmean/agemean  ,vcovall,method="reml")
summary(mvmeta)

### META-ANALYSIS METHOD (OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL) ###

mvall <- mixmeta(coefall~1,vcovall,method="reml") #cp, coef et vcov n'existent pas, overall cumulative are named coefall et vcovall in the loop

summary(mvall)

# OBTAIN BLUPS from meta regression including average mean and range of mean temperature#
blup <- blup(mvall, vcov = T)
#blupall <- blup(mvall, vcovall = T) #cp

#write.table(blupall,file="blupall.csv")
#data_temp <-  read.table("blupall.csv",header=TRUE,row.names=1)


#(mvmeta1 <- mixmeta(blup~1  ,vcov,method=c("reml"), data=data_temp))
#summary(mvmeta1)


### Re-centering
# GENERATE THE MATRIX FOR STORING THE RESULTS
minperccity <- mintempcity <- rep(NA,length(data1list))
names(mintempcity) <- names(minperccity) <- cities$city

# DEFINE MINIMUM RISK VALUES (MRT): EXCLUDE LOW AND VERY HOT TEMPERATURE
for(i in seq(length(data1list))) {
  data <- data1list[[i]]
  predvar <- quantile(templist[[i]],50:75/100,na.rm=T)#
  #predvar <- quantile(templist[[i]],1:99/100,na.rm=T)#
  #predvar <- quantile(templist[[i]],5:95/100,na.rm=T)# 
  # REDEFINE THE FUNCTION USING ALL THE ARGUMENTS (BOUNDARY KNOTS INCLUDED)
  argvar <- list(x=predvar,fun=varfun,
                 knots=quantile(templist[[i]],varper/100,na.rm=T),
                 Bound=range(templist[[i]],na.rm=T))
  bvar <- do.call(onebasis,argvar)
  minperccity[i] <- (50:75)[which.min((bvar%*%blup[[i]]$blup))]
  #minperccity[i] <- (1:99)[which.min((bvar%*%blup[[i]]$blup))]
  #minperccity[i] <- (90)[which.min((bvar%*%blup[[i]]$blup))]
  mintempcity[i] <- quantile(templist[[i]],minperccity[i]/100,na.rm=T)
}

#just for information, these are the minimum risk temperature percentiles and minimum risk absolute temperature
minperccity
mintempcity


# CITY-SPECIFIC FIRST-STAGE SUMMARIES

regall <- lapply(seq(nrow(coefall)),function(i) crosspred(bvar,coef=coefall[i,],
                                                          vcov=vcovall[[i]],model.link="log",cen = mintempcity[i] ))

# OVERALL CUMULATIVE SUMMARY ASSOCIATION FOR MAIN MODEL

cpall <- crosspred(bvar,coef=coef(mvall),vcov=vcov(mvall),model.link="log",
                   by=0.1,from=bound[1],to=bound[2], cen = mintempcity[i] )



#Exposure-response curve#
cpall <- crosspred(bvar,coef=coef(mvall),vcov=vcov(mvall),model.link="log",
                   by=0.1,from=bound[1],to=bound[2], cen = mintempcity[i])


plot(cpall,ci="lines",ylab="Odds Ratio",xlab="Temperature (\u00B0C)",col=2, lwd=2)
mtext(text=paste("Total stroke"),cex=1)



########### Run models again with MRT and obtain estimates at 99th percentile of temperature per city and obtain pooled effect ############

options(warn=-1)

system.time(
  for(i in seq(m)) {
    
    # PRINT ITERATION
    cat(i,"")
    
    # LOAD
    data <- data1list[[i]]
    
    # DEFINE THE CROSSBASIS
    argvar <- list(
      fun = varfun,
      knots = quantile(templist[[i]], varper / 100, na.rm = T),
      Bound = range(templist[[i]], na.rm = T) )
    
    # CREATE THE CROSSBASIS (Adjusting factors)
    cb.temp <- crossbasis(templist[[i]], lag=lag6, argvar=argvar, arglag=arglag6)
    cb.humid <- crossbasis(humidlist[[i]], lag=lag6, argvar=argvar3, arglag=arglag6)
    
    
    
    # FIT THE MODEL (conditional logistic regression model here)
    
    model <- clogit(data$case ~  cb.temp    + cb.humid + strata(data$match))
    # model <- clogit(data$Case ~  cb.temp  +  strata(data$match))
    
    # REDUCTION TO OVERALL CUMULATIVE
    redall <- crossreduce(cb.temp,model, cen = mintempcity[i])
    coefall[i,] <- coef(redall)
    vcovall[[i]] <- vcov(redall)
    
    
    # REDUCTION TO SPECIFIC PREDICTORS
    
    red_spec <- crossreduce(cb.temp,model,type="var",value= quantile(templist[[i]], 99/100,na.rm=T), cen = mintempcity[i])
    coef_spec[i,] <- coef(red_spec)
    vcov_spec[[i]] <- vcov(red_spec)
    
  })

# RESET WARNING
options(warn=0)


#mvall <- mixmeta(coef~1,vcov,method=method)
mvall <- mixmeta(coefall~1,vcovall,method="reml") #cp


# CREATE BASES FOR PREDICTION

xvar <- seq(bound[1],bound[2],by=0.1)
bvar <- do.call("onebasis",c(list(x=xvar),attr(cb.temp,"argvar")))
xlag <- 0:60/10 #Changed to lag 6
blag <- do.call("onebasis",c(list(x=xlag),attr(cb.temp,"arglag")))


###### MULTIVARIATE META-REGRESSION with BLUP #########

mvmeta <- mixmeta(coefall~ 1 + as.factor(agemean), random= ~ 1|hruidmean/agemean  ,vcovall,method="reml")
summary(mvmeta)


# FUNCTION FOR COMPUTING THE P-VALUE OF A WALD TEST 
fwald <- function(model,var) { 
  ind <- grep(var,names(coef(model))) 
  coef <- coef(model)[ind] 
  vcov <- vcov(model)[ind,ind] 
  waldstat <- coef%*%solve(vcov)%*%coef 
  df <- length(coef) 
  return(1-pchisq(waldstat,df)) 
}

# TEST THE EFFECTS
fwald(mvmeta,"agemean")


###### MULTIVARIATE META-REGRESSION with BLUP with agemean as an effect modifier #########

#age cat 1#
mvmeta <- mixmeta(coefall~ 1  ,vcovall,method="reml", subset = agemean==1)
summary(mvmeta)

# CITY-SPECIFIC FIRST-STAGE SUMMARIES

regall <- lapply(seq(nrow(coefall)),function(i) crosspred(bvar,coef=coefall[i,],
                                                          vcov=vcovall[[i]],model.link="log",cen = mintempcity[i] ))

# OVERALL CUMULATIVE SUMMARY ASSOCIATION FOR MAIN MODEL


cpall <- crosspred(bvar,coef=coef(mvmeta),vcov=vcov(mvmeta),model.link="log",
                   by=0.1,from=bound[1],to=bound[2], cen = mintempcity[i], at=quantile(templist[[i]], 99/100,na.rm=T))

# OVERALL EFFECTS OF TEMPERATURE over 0 to 6 days#
round(with(cpall,cbind(allRRfit,allRRlow,allRRhigh)[,]),3)

# #Exposure-response curve#
# cpall <- crosspred(bvar,coef=coef(mvall),vcov=vcov(mvall),model.link="log",
#                    by=0.1,from=bound[1],to=bound[2], cen = mintempcity[i])
# 
# 
# plot(cpall,ci="lines",ylab="Odds Ratio",xlab="Temperature (\u00B0C)",col=2, lwd=2)
# mtext(text=paste("<65"),cex=1)

#age cat 2#

mvmeta <- mixmeta(coefall~ 1  ,vcovall,method="reml", subset = agemean==2)
summary(mvmeta)

# CITY-SPECIFIC FIRST-STAGE SUMMARIES

regall <- lapply(seq(nrow(coefall)),function(i) crosspred(bvar,coef=coefall[i,],
                                                          vcov=vcovall[[i]],model.link="log",cen = mintempcity[i] ))

# OVERALL CUMULATIVE SUMMARY ASSOCIATION FOR MAIN MODEL


cpall <- crosspred(bvar,coef=coef(mvmeta),vcov=vcov(mvmeta),model.link="log",
                   by=0.1,from=bound[1],to=bound[2], cen = mintempcity[i], at=quantile(templist[[i]], 99/100,na.rm=T))

# OVERALL EFFECTS OF TEMPERATURE over 0 to 6 days#
round(with(cpall,cbind(allRRfit,allRRlow,allRRhigh)[,]),3)

# #Exposure-response curve#
# cpall <- crosspred(bvar,coef=coef(mvall),vcov=vcov(mvall),model.link="log",
#                    by=0.1,from=bound[1],to=bound[2], cen = mintempcity[i])
# 
# 
# plot(cpall,ci="lines",ylab="Odds Ratio",xlab="Temperature (\u00B0C)",col=2, lwd=2)
# mtext(text=paste("65-74"),cex=1)


#age cat 3#

mvmeta <- mixmeta(coefall~ 1  ,vcovall,method="reml", subset = agemean==3)
summary(mvmeta)

# CITY-SPECIFIC FIRST-STAGE SUMMARIES

regall <- lapply(seq(nrow(coefall)),function(i) crosspred(bvar,coef=coefall[i,],
                                                          vcov=vcovall[[i]],model.link="log",cen = mintempcity[i] ))

# OVERALL CUMULATIVE SUMMARY ASSOCIATION FOR MAIN MODEL


cpall <- crosspred(bvar,coef=coef(mvmeta),vcov=vcov(mvmeta),model.link="log",
                   by=0.1,from=bound[1],to=bound[2], cen = mintempcity[i], at=quantile(templist[[i]], 99/100,na.rm=T))


# OVERALL EFFECTS OF TEMPERATURE over 0 to 6 days#
round(with(cpall,cbind(allRRfit,allRRlow,allRRhigh)[,]),3)

# #Exposure-response curve#
# cpall <- crosspred(bvar,coef=coef(mvall),vcov=vcov(mvall),model.link="log",
#                    by=0.1,from=bound[1],to=bound[2], cen = mintempcity[i])
# 
# 
# plot(cpall,ci="lines",ylab="Odds Ratio",xlab="Temperature (\u00B0C)",col=2, lwd=2)
# mtext(text=paste(">75"),cex=1)





###################### sex as effect modifier ###############################


# Set working directory
setwd("/Users/Brandon W/Documents/Western/THESIS/Analysis/FINAL ANALYSIS")

# Load dataset
df <- read.csv("stroke_pm25_temp_imputed.csv", stringsAsFactors = FALSE)

df  <- df %>% rename(stroke_type = stoke_type)                 # fix typo
keep_vars <- c(
  "case","match","fis_yr","visit_d","fsa","age_yr","male","stroke_type",
  paste0("temp_avg_lag",0:7), paste0("temp_max_lag",0:7),
  paste0("humidity_lag",0:7)
)
thesis_data <- df[, keep_vars]

## Attach health-region codes
fsa_hr <- read_xlsx("Canada_FSA_to_HR18_Correspondence.xlsx") %>% 
  rename(fsa = FSA)
final  <- merge(thesis_data, fsa_hr[,c("fsa","HR_UID")], by = "fsa")

## Add visit date and holiday flag
final$visit_date <- dmy(final$visit_d)
holiday_dates    <- as.Date(holidayTSX(unique(year(final$visit_date))))
final$stat_hol   <- as.integer(final$visit_date %in% holiday_dates)

# Add categorical age variable
final$age_cat <- cut(final$age_yr,  breaks = c(0, 65, 75, Inf),
                     labels = c(1, 2, 3),
                     right = TRUE, # Includes the rightmost break value in the interval
                     include.lowest = TRUE) # Includes the lowest value in the first interval

# Remove rows where male is NA
final <- final[!is.na(final$male), ]
# Merge HR_UID with age_cat
final$hr_uid_male <- paste(final$HR_UID,
                              final$male,
                              sep = "_")

table(final$hr_uid_male)

class_data <- data.frame(student=c('3526_0', '3526_1', '3527_0', '3527_1', '3530_0', '3530_1', '3533_0', '3533_1', '3534_0', '3534_1', '3535_0', '3535_1', '3536_0', '3536_1', '3537_0',
                                   '3537_1', '3538_0', '3538_1', '3539_0', '3539_1', '3540_0', '3540_1', '3541_0', '3541_1', '3542_0', '3542_1', '3543_0', '3543_1', '3544_0', '3544_1',
                                   '3546_0', '3546_1', '3547_0', '3547_1', '3549_0', '3549_1', '3551_0', '3551_1', '3553_0', '3553_1', '3554_0', '3554_1', '3555_0', '3555_1', '3556_0',
                                   '3556_1', '3557_0', '3557_1', '3558_0', '3558_1', '3560_0', '3560_1', '3561_0', '3561_1', '3562_0', '3562_1', '3563_0', '3563_1', '3565_0', '3565_1',
                                   '3566_0', '3566_1', '3568_0', '3568_1', '3570_0', '3570_1', '3575_0', '3575_1', '3595_0', '3595_1'))

# ARRANGE THE DATA AS A LIST OF DATASETS
regions <- as.character(unique(class_data$student))
regions


data1list <- lapply(regions, function(region) final[final$hr_uid_male==region,])
names(data1list) <- regions
m <- length(data1list)
m

table(regions)
regions
summary(regions)


# METADATA FOR LOCATIONS
cities <- data.frame(
  city = regions)


# CREATE  META-PREDICTORS
gendermean <- sapply(data1list,function(x) mean(x$male,na.rm=T))
hruidmean <- sapply(data1list,function(x) mean(x$HR_UID,na.rm=T))


# METADATA FOR LOCATIONS 

# IDENTIFY THE MIN AND MAX (RANGE) FOR EACH  VARIABLE, PER region 
# WE USE Lag0 (INDEX DATE) FOR SIMPLICITY SINCE WE ASSUME IT WON'T CHANGE CONSIDERABLY COMPARED 
# TO THE OTHER Lags. So you will need to do so in your dataset. 

# B - TEMPERATURE & humidity
ranges_temp <- t(sapply(data1list,function(x) range(x$temp_avg_lag0 ,na.rm=T)))
ranges_hum <- t(sapply(data1list,function(x) range(x$humidity_lag0,na.rm=T)))

bound <- colMeans(ranges_temp)
exposure <- "Temperature"

# DEFINE THE CENTERING POINT, TYPE AND DEGREES OF FREEDOM: 
#  POLLUANTS ARE CONSIDERED LINEAR

cen <- 0
type <- "ns"
df <- 3

# TEMPERATURE IS CONSIDERED NON-LINEAR AND WILL USE NS (NATURAL SPLINES) (you can use this specification if we want to test non linear associations of air pollution)
cen_temp <- 20
type_temp <- "ns"
df_temp <- 3



#Change lag 0 to 6
temp <- as.matrix(final[,rep(9:15, each=1)])   
#temp[1:15,1:7]
templist <- lapply(regions, function(region) temp[final$hr_uid_male==region,])
names(templist) <- regions

humid <- as.matrix(final[,rep(25:31, each=1)]) 
#humid[1:15,1:7]
humidlist <- lapply(regions, function(region) humid[final$hr_uid_male==region,])
names(humidlist) <- regions

#lag space specification#
lag6 <- c(0,6)
lagknots6 <- logknots(6,df=3,int=T)
arglag6 <- list(fun="ns",knots=lagknots6)

#unconstrained lag#
arglag <- list(fun="integer")


# #lag 5#
# lag5 <- c(0,5)
# lagknots5 <- logknots(5,df=3,int=T)
# arglag5 <- list(fun="ns",knots=lagknots5)

# lag 7
# lag7 <- c(0,7)
# lagknots7 <- logknots(7,df=3,int=T)
# arglag7 <- list(fun="ns",knots=lagknots7)

#lag 1# Changed to 1
#lag1 <- c(0,1)
#lagknots1 <- logknots(1,df=3,int=T)
#arglag1 <- list(fun= "strata")

# Specification of exposure function
varfun <- "ns"
vardegree <- NULL
varper <- c(50,90) 

argvar3 <- list(fun="ns",cen=0)


# 2 Empty objects to store coeff-vcov and other city-specific info
coefall <- matrix(NA,m,3,dimnames=list(regions,paste("spl",seq(3),sep=""))) #3 because of 2 knots and 1 intercept#

vcovall <- vector("list", m)
names(vcovall) <- regions

# PREDICTOR SPECIFIC SUMMARIES
coef_spec <- matrix(NA,m,3,dimnames=list(regions,paste("spl",seq(3),sep=""))) #3 df for the lag space#
vcov_spec <- vector("list",m) 
names(vcov_spec) <- regions 


#------------------ Total stroke ----------------------------#

# Warning for prediction beyond boundaries suppressed
options(warn=-1)

system.time(
  for(i in seq(m)) {
    
    # PRINT ITERATION
    cat(i,"")
    
    # LOAD
    data <- data1list[[i]]
    
    # DEFINE THE CROSSBASIS
    argvar1 <- list(
      fun = varfun,
      knots = quantile(templist[[i]], varper / 100, na.rm = T),
      Bound = range(templist[[i]], na.rm = T) )
    
    
    # CREATE THE CROSSBASIS (Adjusting factors)
    cb.temp <- crossbasis(templist[[i]], lag=lag6, argvar=argvar1, arglag=arglag6)
    cb.humid <- crossbasis(humidlist[[i]], lag=lag6, argvar=argvar3, arglag=arglag6)
    
    
    # FIT THE MODEL (conditional logistic regression model here)
    
    model <- clogit(data$case ~  cb.temp   + cb.humid + strata(data$match))
    
    
    # REDUCTION TO OVERALL CUMULATIVE
    redall <- crossreduce(cb.temp,model,cen =mean(templist[[i]],na.rm = T))
    coefall[i,] <- coef(redall)
    vcovall[[i]] <- vcov(redall)
    
  })

# RESET WARNING
options(warn=0)

model


method <- "reml"


###### MULTIVARIATE META-REGRESSION with BLUP #########

mvmeta <- mixmeta(coefall~ 1 + as.factor(gendermean), random= ~ 1|hruidmean/gendermean  ,vcovall,method="reml")
summary(mvmeta)


# FUNCTION FOR COMPUTING THE P-VALUE OF A WALD TEST 
fwald <- function(model,var) { 
  ind <- grep(var,names(coef(model))) 
  coef <- coef(model)[ind] 
  vcov <- vcov(model)[ind,ind] 
  waldstat <- coef%*%solve(vcov)%*%coef 
  df <- length(coef) 
  return(1-pchisq(waldstat,df)) 
}

# TEST THE EFFECTS
fwald(mvmeta,"gendermean")
#fwald(mvmeta,"hruidmean")

mvmeta <- mixmeta(coefall~ 1 + gendermean, random= ~ 1|hruidmean/gendermean  ,vcovall,method="reml")
summary(mvmeta)

#META PREDICTORS ARE NOT SIGNIFICANT

### META-ANALYSIS METHOD (OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL) ###

mvall <- mixmeta(coefall~1,vcovall,method="reml") #cp, coef et vcov n'existent pas, overall cumulative are named coefall et vcovall in the loop

summary(mvall)


# OBTAIN BLUPS from meta regression including average mean and range of mean temperature#
blup <- blup(mvall, vcov = T)
#blupall <- blup(mvall, vcovall = T) #cp

#write.table(blupall,file="blupall.csv")
#data_temp <-  read.table("blupall.csv",header=TRUE,row.names=1)


#(mvmeta1 <- mixmeta(blup~1  ,vcov,method=c("reml"), data=data_temp))
#summary(mvmeta1)


### Re-centering
# GENERATE THE MATRIX FOR STORING THE RESULTS
minperccity <- mintempcity <- rep(NA,length(data1list))
names(mintempcity) <- names(minperccity) <- cities$city

# DEFINE MINIMUM RISK VALUES (MRT): EXCLUDE LOW AND VERY HOT TEMPERATURE
for(i in seq(length(data1list))) {
  data <- data1list[[i]]
  predvar <- quantile(templist[[i]],50:75/100,na.rm=T)#
  #predvar <- quantile(templist[[i]],1:99/100,na.rm=T)#
  #predvar <- quantile(templist[[i]],5:95/100,na.rm=T)# 
  # REDEFINE THE FUNCTION USING ALL THE ARGUMENTS (BOUNDARY KNOTS INCLUDED)
  argvar <- list(x=predvar,fun=varfun,
                 knots=quantile(templist[[i]],varper/100,na.rm=T),
                 Bound=range(templist[[i]],na.rm=T))
  bvar <- do.call(onebasis,argvar)
  minperccity[i] <- (50:75)[which.min((bvar%*%blup[[i]]$blup))]
  #minperccity[i] <- (1:99)[which.min((bvar%*%blup[[i]]$blup))]
  #minperccity[i] <- (90)[which.min((bvar%*%blup[[i]]$blup))]
  mintempcity[i] <- quantile(templist[[i]],minperccity[i]/100,na.rm=T)
}

#just for information, these are the minimum risk temperature percentiles and minimum risk absolute temperature
minperccity
mintempcity


# CITY-SPECIFIC FIRST-STAGE SUMMARIES

regall <- lapply(seq(nrow(coefall)),function(i) crosspred(bvar,coef=coefall[i,],
                                                          vcov=vcovall[[i]],model.link="log",cen = mintempcity[i] ))

# OVERALL CUMULATIVE SUMMARY ASSOCIATION FOR MAIN MODEL

cpall <- crosspred(bvar,coef=coef(mvall),vcov=vcov(mvall),model.link="log",
                   by=0.1,from=bound[1],to=bound[2], cen = mintempcity[i] )



#Exposure-response curve#
cpall <- crosspred(bvar,coef=coef(mvall),vcov=vcov(mvall),model.link="log",
                   by=0.1,from=bound[1],to=bound[2], cen = mintempcity[i])


plot(cpall,ci="lines",ylab="Odds Ratio",xlab="Temperature (\u00B0C)",col=2, lwd=2)
mtext(text=paste("Total stroke"),cex=1)



########### Run models again with MRT and obtain estimates at 99th percentile of temperature per city and obtain pooled effect ############

options(warn=-1)

system.time(
  for(i in seq(m)) {
    
    # PRINT ITERATION
    cat(i,"")
    
    # LOAD
    data <- data1list[[i]]
    
    # DEFINE THE CROSSBASIS
    argvar <- list(
      fun = varfun,
      knots = quantile(templist[[i]], varper / 100, na.rm = T),
      Bound = range(templist[[i]], na.rm = T) )
    
    # CREATE THE CROSSBASIS (Adjusting factors)
    cb.temp <- crossbasis(templist[[i]], lag=lag6, argvar=argvar, arglag=arglag6)
    cb.humid <- crossbasis(humidlist[[i]], lag=lag6, argvar=argvar3, arglag=arglag6)
    
    
    
    # FIT THE MODEL (conditional logistic regression model here)
    
    model <- clogit(data$case ~  cb.temp    + cb.humid + strata(data$match))
    # model <- clogit(data$Case ~  cb.temp  +  strata(data$match))
    
    # REDUCTION TO OVERALL CUMULATIVE
    redall <- crossreduce(cb.temp,model, cen = mintempcity[i])
    coefall[i,] <- coef(redall)
    vcovall[[i]] <- vcov(redall)
    
    
    # REDUCTION TO SPECIFIC PREDICTORS
    
    red_spec <- crossreduce(cb.temp,model,type="var",value= quantile(templist[[i]], 99/100,na.rm=T), cen = mintempcity[i])
    coef_spec[i,] <- coef(red_spec)
    vcov_spec[[i]] <- vcov(red_spec)
    
  })

# RESET WARNING
options(warn=0)

#mvall <- mixmeta(coef~1,vcov,method=method)
mvall <- mixmeta(coefall~1,vcovall,method="reml") #cp
summary(mvall)

# CREATE BASES FOR PREDICTION

xvar <- seq(bound[1],bound[2],by=0.1)
bvar <- do.call("onebasis",c(list(x=xvar),attr(cb.temp,"argvar")))
xlag <- 0:60/10 #Changed to lag 6
blag <- do.call("onebasis",c(list(x=xlag),attr(cb.temp,"arglag")))

###### MULTIVARIATE META-REGRESSION with BLUP #########

mvmeta <- mixmeta(coefall~ 1 + as.factor(gendermean), random= ~ 1|hruidmean/gendermean  ,vcovall,method="reml")
summary(mvmeta)


# FUNCTION FOR COMPUTING THE P-VALUE OF A WALD TEST 
fwald <- function(model,var) { 
  ind <- grep(var,names(coef(model))) 
  coef <- coef(model)[ind] 
  vcov <- vcov(model)[ind,ind] 
  waldstat <- coef%*%solve(vcov)%*%coef 
  df <- length(coef) 
  return(1-pchisq(waldstat,df)) 
}

# TEST THE EFFECTS
fwald(mvmeta,"gendermean")



###### MULTIVARIATE META-REGRESSION with BLUP, with gendermean as an effect modifier #########

#gender 1 men#
mvmeta <- mixmeta(coefall~ 1  ,vcovall,method="reml", subset = gendermean==1)
summary(mvmeta)

# CITY-SPECIFIC FIRST-STAGE SUMMARIES

regall <- lapply(seq(nrow(coefall)),function(i) crosspred(bvar,coef=coefall[i,],
                                                          vcov=vcovall[[i]],model.link="log",cen = mintempcity[i] ))

# OVERALL CUMULATIVE SUMMARY ASSOCIATION FOR MAIN MODEL


cpall <- crosspred(bvar,coef=coef(mvmeta),vcov=vcov(mvmeta),model.link="log",
                   by=0.1,from=bound[1],to=bound[2], cen = mintempcity[i], at=quantile(templist[[i]], 99/100,na.rm=T))

# OVERALL EFFECTS OF TEMPERATURE over 0 to 6 days#
round(with(cpall,cbind(allRRfit,allRRlow,allRRhigh)[,]),3)

#gender 0 females#

mvmeta <- mixmeta(coefall~ 1  ,vcovall,method="reml", subset = gendermean==0)
summary(mvmeta)

# CITY-SPECIFIC FIRST-STAGE SUMMARIES

regall <- lapply(seq(nrow(coefall)),function(i) crosspred(bvar,coef=coefall[i,],
                                                          vcov=vcovall[[i]],model.link="log",cen = mintempcity[i] ))

# OVERALL CUMULATIVE SUMMARY ASSOCIATION FOR MAIN MODEL


cpall <- crosspred(bvar,coef=coef(mvmeta),vcov=vcov(mvmeta),model.link="log",
                   by=0.1,from=bound[1],to=bound[2], cen = mintempcity[i], at=quantile(templist[[i]], 99/100,na.rm=T))

# OVERALL EFFECTS OF TEMPERATURE over 0 to 6 days#
round(with(cpall,cbind(allRRfit,allRRlow,allRRhigh)[,]),3)
