library(writexl)

#####################################
###### New Error Ratio (Squared) ######
####################################
nsubject <- 10
subjectpool <- c(1:7, 9:nsubject)
behavdata <- NA

for (subjectID in subjectpool){

  if (subjectID <10){ 
    filename <- paste0("HT0",subjectID)
  }
  else {
    filename <- paste0("HT",subjectID)
  }
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




