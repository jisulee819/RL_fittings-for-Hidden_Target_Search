library(writexl)

subjectpool <- c(1:22, 24, 25)
for (subjectID in subjectpool){
  
  if (subjectID < 10){
    filename <- paste0("HT0",subjectID)
  } else {filename <- paste0("HT",subjectID)
  }
  
  RP <-  readxl::read_excel(paste0("/Volumes/clmnlab/HT/regressors/MAP_regressor/",filename,"_regressor.xlsx"))
  names(RP) <- c("number", "RP", "label", "entropy")
  if (subjectID == 1){
    RPdata <- data.frame(RP)
  } else {
    RPdata <- rbind(RPdata, RP)
  }
}


RPdata$ID <- rep(subjectpool, each=360)
RPdata$ent.pr <- RPdata$RP*(log(1/RPdata$RP)) + (1-RPdata$RP)*log(1/(1-RPdata$RP))

#plot(x = RPdata$RP, y= RPdata$ent.pr)
RPdata$trial <- rep(c(1:10))
write_xlsx(RPdata, "/Volumes/clmnlab/HT/regressors/MAP_regressor/HT_entropy_rp.xlsx")
