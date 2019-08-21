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
###### Trial - Normalized Entropy ######
####################################

nsubject <- 10
subjectpool <- c(1:7, 9:nsubject)
fit <- NA
fits <- NA

for (subjectID in subjectpool){
  
  if (subjectID <10){ 
    filename <- paste0("HT0",subjectID)
  } else {
    filename <- paste0("HT",subjectID)
  }
  fit <- read_xlsx(paste0("/Users/jisu/Documents/Hidden Target/Behavioral Data/regressor/",filename,"_regressor.xlsx"))
  fit <- fit[,-1]
  fit[,dim(fit)[2]+1] <- subjectID
  names(fit) <- c("prob_reward","infer_type","entropy","subjectID")
  
  if (subjectID == 1){
    fits <- fit
  } else{
    fits <- rbind(fits, fit)  
  }
}
fits$trial <- rep(1:10)


ent.mean <- NA
ent.sd <- NA
for (i in 1:10){
  ent.mean[i] <- mean(fits$entropy[fits$trial==i])
}

for (i in 1:10){
  ent.sd[i] <- sd(fits$entropy[fits$trial==i])
}
ent.sd
fits$ent.mean <- ent.mean
fits$ent.sd <- ent.sd
fits$ent.norm <- (fits$entropy - fits$ent.mean) / fits$ent.sd
fits$block <- rep(c(1:36), each=10)
fits <- data.frame(fits)
write_xlsx(fits, "/Users/jisu/Documents/Hidden Target/Behavioral Data/regressor/HT_entropy_normalization.xlsx")


ggplot(data = fits, aes(x = trial, y = entropy)) + 
  #scale_y_continuous(name="Reward Rate", limits=c(0, 25)) + 
  scale_x_continuous("Trial",limits=c(1,10),breaks=c(1:10)) +
  geom_point(aes(col=trial)) +
  geom_line(aes(x=trial, y=ent.mean)) +
  ggtitle("Entropy per Trial") + 
  theme(plot.title=element_text(face="bold", size=20, hjust=0.5))


ggplot(data = fits[fits$block == 1,], aes(x = trial, y = entropy, group=subjectID)) + 
  #scale_y_continuous(name="Reward Rate", limits=c(0, 25)) + 
  scale_x_continuous("Trial",limits=c(1,10),breaks=c(1:10)) +
  #geom_point(aes(col=trial)) +
  geom_line(data = fits[fits$block == 2,], aes(x = trial, y = entropy, group=subjectID, col=subjectID)) +
  geom_line(data = fits[fits$block == 3,], aes(x = trial, y = entropy, group=subjectID, col=subjectID)) +
  geom_line(data = fits[fits$block == 4,], aes(x = trial, y = entropy, group=subjectID, col=subjectID)) +
  geom_line(data = fits[fits$block == 5,], aes(x = trial, y = entropy, group=subjectID, col=subjectID)) +
  geom_line(data = fits[fits$block == 6,], aes(x = trial, y = entropy, group=subjectID, col=subjectID)) +
  geom_line(data = fits[fits$block == 7,], aes(x = trial, y = entropy, group=subjectID, col=subjectID)) +
  geom_line(data = fits[fits$block == 8,], aes(x = trial, y = entropy, group=subjectID, col=subjectID)) +
  geom_line(data = fits[fits$block == 9,], aes(x = trial, y = entropy, group=subjectID, col=subjectID)) +
  geom_line(data = fits[fits$block == 10,], aes(x = trial, y = entropy, group=subjectID, col=subjectID)) +
  ggtitle("Entropy per Trial(Block:1~10)") + 
  theme(plot.title=element_text(face="bold", size=20, hjust=0.5))


#####################################
###### Entropy Change Amount $ Ratio ######
####################################
fits2 <- data.frame(subjectID = rep(c(1:7, 9:nsubject), each = 9*36))
fits2$block <- rep(c(1:36), each=9, nsubject-1)
fits2$pretrial <- rep(c(1:9))

a = c(1:(360*(nsubject-1)))[-seq(10,(360*(nsubject-1)),10)]

b <- NA
c <- NA
for (i in a){
  b[i] <- fits$entropy[i+1] - fits$entropy[i]
  c[i] <- (fits$entropy[i+1] - fits$entropy[i])/(fits$entropy[i])
}

fits2$ent.diff <- na.omit(b)
fits2$ent.diff.ratio <- na.omit(c)

p <- ggplot(data = fits2) +
  scale_y_continuous("Entropy Difference",limits=c(-1,1),breaks=seq(-1,1,0.1)) 
for (i in 1:36){
  p <- p + geom_line(data = fits2[fits$block==i,], aes(x=pretrial, y=ent.diff, group=subjectID, col=subjectID))
}
p  

#p2 <- ggplot(data = fits2) +
#  scale_y_continuous("Entropy Difference",limits=c(-0.2,0.2),breaks=seq(-0.2,0.2,0.1)) 
#for (i in 1:36){
#  p2 <- p2 + geom_line(data = fits2[fits$block==i,], aes(x=pretrial, y=ent.diff.ratio, group=subjectID, col=subjectID))
#}
#p2  

# Since the correlation btw diff and diff ratio is about 95%, only use the diff data

write_xlsx(fits2, "/Users/jisu/Documents/Hidden Target/Behavioral Data/regressor/HT_entropy_diff.xlsx")
write_xlsx(fits2, "/Volumes/clmnlab/HT/regressors/MAP_regressor/HT_entropy_diff.xlsx")



