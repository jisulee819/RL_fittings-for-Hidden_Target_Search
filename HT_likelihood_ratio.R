
#####################################
###### Regression with Error Ratio ######
####################################
library(writexl)
subjectpool <- c(1:22,25)
#subjectpool <- 1
#data <- matrix(NA, nrow=360, ncol=max(subjectpool))
data <- NA
for (subjectID in subjectpool){

  if (subjectID < 10){
    filename <- paste0("HT0",subjectID)
  }    else {filename <- paste0("HT",subjectID)
  }
  fit <- read.csv(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename,"/",filename,"_fitted.csv"), header = TRUE)
  
  if (subjectID == 1){
    data <- fit[,11]
  } else {
    data <- c(data, fit[,11])
  }
  
  ########## why NA ##########
  ########## why NA ##########
  crit1 <- quantile(data, c(.50), na.rm=TRUE) 
  crit2 <- quantile(data, c(1/3, 2/3), na.rm=TRUE) 
  
  label1 <- ifelse(data < crit1, 1, 2)
  label2 <- ifelse(data < crit2[1], 1, ifelse(data > crit2[2], 2, 0))
  dataset <- data.frame(data, label1, label2)
  ########## why NA ##########
  ########## why NA ##########
  
  #data[,subjectID] <- fit[,11]

}

ID <- rep(subjectpool, each=360)
dataset <- data.frame(cbind(ID, dataset))
dataset$trial <- rep(c(1:10))

#a <- NA
#for (i in 1:10){
#a[i] <- mean(data$data[data$trial == i], na.rm=TRUE)
#}
#a  #since log(probability) < 0, this index indicates MAP as it decreases

write_xlsx(dataset, "/Volumes/clmnlab/HT/regressors/MAP_regressor/HT_likelihood_ratio.xlsx")
write_xlsx(dataset, "/Users/jisu/desktop/HT_likelihood_ratio.xlsx")

