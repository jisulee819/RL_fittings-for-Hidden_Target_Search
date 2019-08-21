library(MASS)
library(readr)
library(mvtnorm)
library(corpcor)
library(Hmisc)
library(reshape2)
library(ggplot2)
library(writexl)
library(readxl)
#####################################
###### Reward Rate per Block ######
####################################

subjectpool <- c(1:7, 9:17,20)
dat <- NA
data <- NA

for (subjectID in subjectpool){
  
  if (subjectID <10){ 
    filename <- paste0("HT0",subjectID)
  } else {
    filename <- paste0("HT",subjectID)
  }
  dat <- read.csv(paste0("/Volumes/clmnlab/HT/behavior_data/raw_data/",filename,"/",filename,".csv"), head=FALSE)
  dat$block <- rep(1:36, each=10)
  for (i in 1:36){
      dat$reward.rate[(10*(i-1)+1):(10*(i))] <- sum(dat$V3[dat$block == i])/10
  }
  
  
  if (subjectID == 1){
    data <- dat
  } else{
    data <- rbind(data, dat)  
  }
}

head(data)
mean(data$reward.rate)

data <- data.frame(data)
rr <- data$reward.rate





ID <- rep(subjectpool, each=360)
data$V1 <- ID


write_xlsx(data, "/Volumes/clmnlab/HT/regressors/MAP_regressor/HT_reward_rate.xlsx")

