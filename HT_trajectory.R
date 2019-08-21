
library(ggplot2)
#####################################################
###### Analysis of behavioral data: Trajectory ######
#####################################################


##################################
###### Trajectory Summation ######
##################################


data <- NA
nsubject <- 15

###################################
###### Trajectory Comparison ######
###################################
subjectpool <- c(5:18,20)
for (subjectID in subjectpool){
  
  if (subjectID < 10){
    filename <- paste0("HT0",subjectID)
  } else {filename <- paste0("HT",subjectID)
  }

  dat <- read.csv(paste0("/Volumes/clmnlab/HT/behavior_data/raw_data/",filename,"/",filename,"_trajdistsum.csv"))
                  
  if (subjectID == 5){
    data <- data.frame(dat)
  } else {
    data <- rbind(data, dat)
  }
}
dim(data)
data

##################################
###### Exploit / Explore ######
##################################



likdata <-  readxl::read_excel(paste0("/Volumes/clmnlab/HT/regressors/MAP_regressor/HT_likelihood_ratio.xlsx"))
likdata <- likdata[likdata$ID > 4,]
dim(likdata)

data <- cbind(likdata, data)
names(data) <- c("ID","like.ratio","label1","label2","trial","number","trajsum")
data$tapplylabel <- 2*data$ID + (data$label1-1)
data$label1 <- ifelse(data$label1==1, "Exploit", "Explore")

y <- tapply(data$trajsum, (data$tapplylabel), mean)
y2 <- tapply(data$trajsum, (data$tapplylabel), sd)

sub.mean <- matrix(NA, ncol=2, nrow=nsubject)
sub.sd <- matrix(NA, ncol=2, nrow=nsubject)

for (i in 1:nsubject){
  sub.mean[i,] <- y[(2*i-1):(2*i)]
  sub.sd[i,] <- y2[(2*i-1):(2*i)]
}

data2 <- data.frame(cbind(sub.mean, sub.sd), "ID" = subjectpool)
names(data2) = c("m.Exploit", "m.Explore", "sd.Exploit", "sd.Explore", "ID") 

t.test(data2$m.Exploit, data2$m.Explore)
#data:  data2$m.Exploit and data2$m.Explore
#t = 2.3365, df = 23.726, p-value = 0.02825
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  3.275133 53.133744
#sample estimates:
#  mean of x mean of y 
# 438.5432  410.3387 

#########################################################################################################
data3 <- data.frame("Trajectory" = c(data2$m.Exploit, data2$m.Explore))
data3$label <- rep(c("Exploit", "Explore"), each=nsubject)


ggplot(data3, aes(x=label, y=Trajectory)) +
  scale_y_continuous(name="Trajectory") +#, limits=c(0, 1700), breaks = seq(0, 1700, 100)) +
  scale_x_discrete(name="Strategy") +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) +
  #geom_point(aes(x=label, y=mean, col="red")) +
  theme(legend.position = "none")





##################################
###### RPE and Next trajectory ######
##################################

subjectpool <- c(5:18,20)
for (subjectID in subjectpool){
  
  if (subjectID < 10){
    filename <- paste0("HT0",subjectID)
  } else {filename <- paste0("HT",subjectID)
  }
  
  RP <-  readxl::read_excel(paste0("/Volumes/clmnlab/HT/regressors/MAP_regressor/",filename,"_regressor.xlsx"))
  dat <- read.csv(paste0("/Volumes/clmnlab/HT/behavior_data/raw_data/",filename,"/",filename,".csv"), header = FALSE)
  names(RP) <- c("number", "RP", "label", "entropy")
  if (subjectID == 5){
    RPdata <- data.frame(RP)
    data <- data.frame(dat)
  } else {
    RPdata <- rbind(RPdata, RP)
    data <- rbind(data, dat)
  }
}

head(data)
dim(data)
dim(RPdata)


RPEdata <- data.frame("Probability Reward" = RPdata$RP, "Reward" = data$V3)
RPEdata$RPE <- RPEdata$Reward - RPEdata$Probability.Reward
RPEdata$ID <- likdata$ID
RPEdata$trajsum <- trajdata$trajsum
RPEdata$trial <- rep(1:10)
RPEdata$label <- ifelse(RPEdata$Reward == 1, "smiley", "unsmiley")
RPEdata$label2 <- likdata$label1

cor.test(RPEdata$RPE, RPEdata$trajsum)
# -0.02135021, p-value = 0.1167

cor.test(RPEdata$Probability.Reward, RPEdata$trajsum)
# -0.07293565, p-value =  8.07e-08


ggplot(RPEdata, aes(x= RPE, y=trajsum, col=ID)) +
  geom_point()

ggplot(RPEdata, aes(x= Probability.Reward, y=trajsum, col=label)) +
  geom_point() +
  geom_line()

# RPE 1~9까지와 trajsum 2~10까지를 비교
RPE2 <- RPEdata$RPE[RPEdata$trial!=10]
label <- RPEdata$label[RPEdata$trial!=10]
label2 <- RPEdata$label2[RPEdata$trial!=1]
trajsum2 <- RPEdata$trajsum[RPEdata$trial!=1]
length(RPE2)
length(trajsum2)

  
data2 <- data.frame(cbind("prev RPE" = RPE2, "post Traj" = trajsum2, "label" = label, "strategy" = label2))

data2$prev.RPE <- as.numeric(paste(data2$prev.RPE))
data2$post.Traj <- as.numeric(paste(data2$post.Traj))
cor.test(data2$prev.RPE, data2$post.Traj)
#[1] -0.0404648, p-value = 0.004782


ggplot(data2, aes(x=prev.RPE, y=post.Traj, col=label)) +
  geom_point()


RPEdata2 <- data2


##################################
###### RPE and next strategy ######
##################################

dim(likdata)
dim(RPEdata)
RPE2 <- RPEdata$RPE[RPEdata$trial!=10]
label <- likdata$label1[likdata$trial!=1]
data2 <- data.frame(cbind("prev RPE" = as.numeric(paste(RPE2)), "post Strategy" = label))
data2$ID <- rep(subjectpool, each=324)
data2$tapplylabel <- 2*data2$ID + (data2$post.Strategy-1)
data2$post.Strategy <- ifelse(data2$post.Strategy==1, "Exploit", "Explore")
head(data2)


y <- tapply(data2$prev.RPE, (data2$tapplylabel), mean)
y2 <- tapply(data2$prev.RPE, (data2$tapplylabel), sd)

sub.mean <- matrix(NA, ncol=2, nrow=nsubject)
sub.sd <- matrix(NA, ncol=2, nrow=nsubject)

for (i in 1:nsubject){
  sub.mean[i,] <- y[(2*i-1):(2*i)]
  sub.sd[i,] <- y2[(2*i-1):(2*i)]
}

data2 <- data.frame(cbind(sub.mean, sub.sd), "ID" = subjectpool)
names(data2) = c("m.postExploit", "m.postExplore", "sd.postExploit", "sd.postExplore", "ID") 

t.test(data2$m.postExploit, data2$m.postExplore)  
#data:  data2$m.postExploit and data2$m.postExplore
#t = 15.694, df = 20.367, p-value = 7.633e-13
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.2482674 0.3242821
#sample estimates:
#  mean of x mean of y 
#0.3077937 0.0215189   
  

#########################################################################################################
data3 <- data.frame("Previous RPE" = c(data2$m.postExploit, data2$m.postExplore))
data3$label <- rep(c("Exploit", "Explore"), each=nsubject)


ggplot(data3, aes(x=label, y=Previous.RPE)) +
  scale_y_continuous(name="Previous RPE", limits=c(-1, 1), breaks = seq(-1, 1, 0.25)) +
  scale_x_discrete(name="Strategy") +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) +
  #geom_point(aes(x=label, y=mean, col="red")) +
  theme(legend.position = "none")

data <- data2

##################################
###### RP and present strategy ######
##################################

dim(likdata)
dim(RPEdata)
RP <- RPEdata$Probability.Reward
label <- likdata$label1
data2 <- data.frame(cbind("RP" = as.numeric(paste(RP)), "Strategy" = label))
data2$ID <- rep(subjectpool, each=360)
data2$tapplylabel <- 2*data2$ID + (data2$Strategy-1)

data2$Strategy <- ifelse(data2$Strategy==1, "Exploit", "Explore")


y <- tapply(data2$RP, (data2$tapplylabel), mean)
y2 <- tapply(data2$RP, (data2$tapplylabel), sd)

sub.mean <- matrix(NA, ncol=2, nrow=nsubject)
sub.sd <- matrix(NA, ncol=2, nrow=nsubject)

for (i in 1:nsubject){
  sub.mean[i,] <- y[(2*i-1):(2*i)]
  sub.sd[i,] <- y2[(2*i-1):(2*i)]
}

data2 <- data.frame(cbind(sub.mean, sub.sd), "ID" = subjectpool)
names(data2) = c("m.postExploit", "m.postExplore", "sd.postExploit", "sd.postExplore", "ID") 

t.test(data2$m.postExploit, data2$m.postExplore)  
#Welch Two Sample t-test
#
#data:  data2$m.postExploit and data2$m.postExplore
#t = 1.2285, df = 27.861, p-value = 0.2296
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.01077090  0.04302696
#sample estimates:
#  mean of x mean of y 
#0.462215  0.446087 

data3 <- data.frame("RP" = c(data2$m.postExploit, data2$m.postExplore))
data3$label <- rep(c("Exploit", "Explore"), each=nsubject)

ggplot(data3, aes(x=label, y=RP)) +
  scale_y_continuous(name="Probability of Reward", limits=c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_x_discrete(name="Strategy") +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) +
  #geom_point(aes(x=label, y=mean, col="red")) +
  theme(legend.position = "none")




##################################
###### Reward and next strategy ######
##################################


data2$prev.reward <- ifelse(data2$prev.RPE > 0, 1, 0)
data2
t.test(data2$prev.reward[data2$post.Strategy=="Exploit"], data2$prev.reward[data2$post.Strategy=="Explore"])
#t = 21.097, df = 4217.9, p-value < 2.2e-16
#0.7709414 0.4870490 

visu.data <- NA
visu.data <- c(mean(data2$prev.reward[data2$post.Strategy == "Exploit"]), mean(data2$prev.reward[data2$post.Strategy == "Explore"]))
visu.data2 <- c(sd(data2$prev.reward[data2$post.Strategy == "Exploit"]), sd(data2$prev.reward[data2$post.Strategy == "Explore"]))
visu.data <- data.frame("mean" = visu.data, "sd" = visu.data2, "label" = c("Exploit", "Explore"))

p <- ggplot(visu.data, aes(x=label, y=mean)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) 

p

#The number of Exploit/Explore is not equal since the first search is omitted

dim(RPEdata2)
RPEdata2$ID <- rep(subjectpool, each = 324)
RPEdata2$strategy <- ifelse(RPEdata2$strategy == 1, "Exploit", "Explore")

RPEdata2$tapplylabel1 <- ifelse(RPEdata2$label == "smiley" & RPEdata2$strategy == "Exploit", 1, 
                                ifelse(RPEdata2$label == "smiley" & RPEdata2$strategy == "Explore", 2, 
                                       ifelse(RPEdata2$label == "unsmiley" & RPEdata2$strategy == "Exploit", 3, 4))) 

RPEdata2$tapplylabel <- ifelse(RPEdata2$label == "smiley" & RPEdata2$strategy == "Exploit", "rewarded Exploit", 
                                ifelse(RPEdata2$label == "smiley" & RPEdata2$strategy == "Explore", "rewarded Explore",
                                       ifelse(RPEdata2$label == "unsmiley" & RPEdata2$strategy == "Exploit", "non-rewarded Exploit", "non-rewarded Explore"))) 


table(RPEdata2$tapplylabel1)

RPEdata2$tapplylabel2 <- 4*RPEdata2$ID + RPEdata2$tapplylabel1 - 1


ggplot(RPEdata2, aes(x=tapplylabel)) +
  scale_x_discrete(name = "Strategy") +
  geom_bar()




y <- tapply(RPEdata2$prev.RPE, (RPEdata2$tapplylabel2), mean)
#y2 <- tapply(RPEdata2$prev.RPE, (RPEdata2$tapplylabel2), sd)

sub.mean <- matrix(NA, ncol=4, nrow=nsubject)
#sub.sd <- matrix(NA, ncol=4, nrow=nsubject)

for (i in 1:nsubject){
  sub.mean[i,] <- y[(4*i-3):(4*i)]
#  sub.sd[i,] <- y2[(4*i-3):(4*i)]
}
RPEdata3 <- data.frame(sub.mean, c(subjectpool))
names(RPEdata3) = c("rewarded Exploit", "rewarded Explore", "non-rewarded Exploit", "non-rewarded Explore", "ID") 

RPEdata3 <- data.frame("mean" = c(RPEdata3[,1],RPEdata3[,2],RPEdata3[,3],RPEdata3[,4]), 
                       "ID" = rep(RPEdata3$ID), "group" = rep(names(RPEdata3)[1:4], each=nsubject))


#fit <- aov(mean ~ group, data=RPEdata3)

ggplot(RPEdata2, aes(x=group)) +
  geom_bar()


ggplot(RPEdata3, aes(x=group, y=mean)) +
  scale_y_continuous(name="Probability of Reward", limits=c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_x_discrete(name="Strategy") +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) +
  #geom_point(aes(x=label, y=mean, col="red")) +
  theme(legend.position = "none")


ggplot(data3, aes(x=label, y=RP)) +
  scale_y_continuous(name="Probability of Reward", limits=c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_x_discrete(name="Strategy") +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) +
  #geom_point(aes(x=label, y=mean, col="red")) +
  theme(legend.position = "none")








# Reinforcement -> RT
# duration modulator: trial start - click

