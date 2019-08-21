likdata <-  readxl::read_excel(paste0("/Volumes/clmnlab/HT/regressors/MAP_regressor/HT_likelihood_ratio.xlsx"))

subjectpool <- c(1:18,20)
nsubject <- 19
for (subjectID in subjectpool){
  
  if (subjectID < 10){
    filename <- paste0("HT0",subjectID)
  } else {filename <- paste0("HT",subjectID)
  }
  
  dat <- read.csv(paste0("/Volumes/clmnlab/HT/behavior_data/raw_data/",filename,"/",filename,".csv"), header=FALSE)
  
  if (subjectID == 1){
    data <- data.frame(dat)
  } else {
    data <- rbind(data, dat)
  }
}

dim(data)
names(data) 
head(likdata)




tapply(data$V10, data$V3, mean)
data$ID <- rep(subjectpool, each=360)
data$tapplylabel <- (2*data$ID+data$V3-1)
y <- tapply(data$V10, data$tapplylabel, mean)

sub.mean <- matrix(NA, ncol=2, nrow=nsubject)
for (i in 1:nsubject){
  sub.mean[i,] <- y[(2*i-1):(2*i)]
}


data2 <- data.frame(cbind(sub.mean), "ID" = subjectpool)
names(data2) = c("m.Exploit", "m.Explore") 

t.test(data2$m.Exploit, data2$m.Explore, paired=TRUE)

#data:  data2$m.Exploit and data2$m.Explore
#t = 2.5914, df = 14, p-value = 0.02133
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  10.6245 112.6665
#sample estimates:
#  mean of the differences 
#61.64549 


# MAP likelihood / MIG likelihood



#####################################
########## Average in Subject ##########
#####################################


likdata$RT <- data$V10
likdata$tapplylabel <- (2*likdata$ID+likdata$label1-2)
y <- tapply(likdata$RT, likdata$tapplylabel, mean)

sub.mean <- matrix(NA, ncol=2, nrow=nsubject)
for (i in 1:nsubject){
  sub.mean[i,] <- y[(2*i-1):(2*i)]
}
  

data2 <- data.frame(cbind(sub.mean), "ID" = subjectpool)
names(data2) = c("m.Exploit", "m.Explore") 

t.test(data2$m.Exploit, data2$m.Explore, paired=TRUE)
#Paired t-test
#
#data:  data2$m.Exploit and data2$m.Explore
#t = 4.3673, df = 18, p-value = 0.0003714
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  64.8248 185.0108
#sample estimates:
#  mean of the differences 
#124.9178 

#########################################################################################################
data3 <- data.frame("Reaction Time" = c(data2$m.Exploit, data2$m.Explore))
data3$label <- rep(c("Exploit", "Explore"), each=nsubject)


ggplot(data3, aes(x=label, y=Reaction.Time)) +
  scale_y_continuous(name="Reaction Time") +#, limits=c(0, 1700), breaks = seq(0, 1700, 100)) +
  scale_x_discrete(name="Strategy") +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) +
  #geom_point(aes(x=label, y=mean, col="red")) +
  theme(legend.position = "none")



#####################################
########## Average in Block ##########
#####################################

likdata$block <- rep(1:36, each=10)
ttestvalues <- matrix(NA, ncol=2, nrow=nsubject+1)
for (i in subjectpool){
  subdata <- likdata[likdata$ID == i,]
  subdata$tapplylabel <- (2*subdata$block + subdata$label1 - 2)
  y <- tapply(subdata$RT, subdata$tapplylabel, mean)
  sub.mean <- matrix(NA, ncol=2, nrow=36)
    for (j in 1:36){
      sub.mean[j,] <- y[(2*j-1):(2*j)]
    }
  #sub.mean
  t <- t.test(sub.mean[,1], sub.mean[,2], paired=TRUE)
  #print(c(t$statistic, t$p.value))
  ttestvalues[i,] <- c(t$statistic, t$p.value)
}


write.csv(ttestvalues, "ttest.csv")
