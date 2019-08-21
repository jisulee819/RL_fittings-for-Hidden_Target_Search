library(writexl)
library(readxl)

#####################################
###### RP and RPE Correlation ######
####################################

for (subjectID in c(1:7,9)){
  filename <- paste0("HT0",subjectID)
  regressor <- read_excel(paste0("/Users/jisu/Documents/Hidden Target/Behavioral Data/regressor/",filename,"_regressor.xlsx"))
  regressor <- cbind(regressor, subjectID)
    if (subjectID == 1){
      regressors <- regressor
    } 
    else {
     regressors <- rbind(regressors, regressor)  
    }
}
regressors <- regressors[,-1]



for (subjectID in c(1:7,9)){
  filename <- paste0("HT0",subjectID)
  raw <- read.csv(paste0("/Users/jisu/Documents/Hidden Target/Behavioral Data/raw_data/",filename,".csv"), header = FALSE)
  if (subjectID == 1){
    raws <- raw
  } 
  else {
    raws <- rbind(raws, raw)  
  }
}

dat <- cbind(regressors$data_opt_h.results.prob_reward, raws$V3)
colnames(dat) <- c("prob","reward")
pred.err <- dat[,2] - dat[,1]
dat <- cbind(dat, pred.err)
dat <- data.frame(dat)
cor(dat)

plot(x=dat$prob, y=dat$pred.err)
plot(x=dat$reward, y=dat$pred.err)



#####################################
###### New Error Ratio (Squared) ######
####################################

for (subjectID in c(1:7,9)){
  filename <- paste0("HT0",subjectID)
  fit <- read.csv(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename,"/",filename,"_fitted.csv"), header = TRUE)
  if (subjectID == 1){
    fits <- fit
  } 
  else if (subjectID > 4){
    fit <- fit[,1:12]
    fits <- rbind(fits, fit)
  }
  else {
    fits <- rbind(fits, fit)  
  }
}
fits <- fits[,-1]
colnames(fits) <- c("inferX", "inferY","migX","migY","mapX","mapY","prob","rew","type","err.ratio","fitting.err")
fits <- data.frame(fits)
fits$subjID <- rep(c(1:7,9), each=360)
fits$map.err <- ifelse(fits$type=="MAP", fits$fitting.err, fits$err.ratio * fits$fitting.err)
fits$mig.err <- ifelse(fits$type=="MIG", fits$fitting.err, fits$fitting.err/fits$err.ratio)
fits$squared.err <- (fits$err.ratio)^2
write_xlsx(fits, "/Users/jisu/Documents/Hidden Target/Behavioral Data/regressor/error.xlsx")













