

library(ggplot2)
multiplot <- function(..., plotlist=NULL, file, cols=2, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

##############################
### Linear Trend per Trial ###
##############################


likdata <-  readxl::read_excel(paste0("/Volumes/clmnlab/HT/regressors/MAP_regressor/HT_likelihood_ratio.xlsx"))
likdata$label<- ifelse(likdata$label1 == 1, "Exploit", "Explore")
dim(likdata)



subjectpool <- c(1:20)
nsubject <- length(subjectpool)
for (subjectID in subjectpool){
  
  if (subjectID < 10){
    filename <- paste0("HT0",subjectID)
  } else {filename <- paste0("HT",subjectID)
  }
  
  dat <- read.csv(paste0("/Volumes/clmnlab/HT/behavior_data/raw_data/",filename,"/",filename,".csv"), header=FALSE)
  likdat <-  read.csv(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename,"/",filename,"_fitted.csv"))
  
  if (subjectID == 1){
    data <- data.frame(dat)
    likdata2 <- data.frame(likdat)
  } else {
    data <- rbind(data, dat)
    likdata2 <- rbind(likdata2, likdat)
  }
}

subjectpool <- c(5:20)
for (subjectID in subjectpool){
  
  if (subjectID < 10){
    filename <- paste0("HT0",subjectID)
  } else {filename <- paste0("HT",subjectID)
  }
  
  tjdat <- read.csv(paste0("/Volumes/clmnlab/HT/behavior_data/raw_data/",filename,"/",filename,"_trajdistsum.csv"))
  
  if (subjectID == 5){
    tjdata <- data.frame(tjdat)
  } else {
    tjdata <- rbind(tjdata, tjdat)
  }
}
falsetj <- cbind("X" = rep(1:360, 4), "x" = 0)
tjdata <- rbind(falsetj, tjdata)

subjectpool <- c(1:20)
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

RPdata$entropy2 <-  RPdata$RP*(log(1/RPdata$RP)) + (1-RPdata$RP)*log(1/(1-RPdata$RP))

visu.data <- data.frame("ID" = likdata$ID, "trial" = rep(c(1:10), 36*length(subjectpool)), "search" = rep(c(1:360), length(subjectpool)), "traj" = tjdata$x, "LR" = (likdata2$V12 - likdata2$V11), "RT" = data$V10, "label" = likdata$label, "rew" = data$V3, "ent" = RPdata$entropy2, "RP" = RPdata$RP, "RPE" = (data$V3 - RPdata$RP))



for (j in c(1:3, 5:8)){
  variablenames = c("subdata$traj", "subdata$LR", "subdata$RT", "subdata$label", "subdata$rew", "subdata$ent", "subdata$RP", "subdata$RPE")
  visualnames = c("Trajectory", "Likelihood Ratio", "Reaction Time", "Exploit / Explore", "Reward Rate", "Entropy", "Reward Probability", "Reward Prediction Error")
  
  y <- matrix(NA, ncol=2, nrow=(nsubject+1)*10)
  for (i in c(1:20)){
      subdata <- visu.data[visu.data$ID == i,]
      y[(10*i-9):(10*i),2] <- tapply(subdata[,j+3], subdata$trial, mean)
  }
  
  x <- tapply(visu.data[,j+3], visu.data$trial, mean)
  y[((nsubject+1)*10-9):((nsubject+1)*10), 2] <- x
  y[,1] <- rep(1:10)
  ID <- rep(c(subjectpool,0), each=10)
  y <- cbind(y, ID)
  
  ggdata <- data.frame(y)
  names(ggdata) <- c("Trial", "Value", "ID")
  
  
  for (id in c(1:20,0)){
    if (id != 0){
     p <- ggplot(data = ggdata[ggdata$ID == id,]) +
           scale_x_continuous(paste0("Subject",id), breaks = c(1:10)) +
           scale_y_continuous(visualnames[j])+#, breaks=seq(-1,1,0.25)) +
           #coord_cartesian(ylim = c(min(ggdata$Value),max(ggdata$Value))) +
           geom_point(aes(x=Trial, y=Value)) +
           geom_line(aes(x=Trial, y=Value))
    } else {
      p <- ggplot(data = ggdata[ggdata$ID == id,]) +
            scale_x_continuous("Average", breaks = c(1:10)) +
            scale_y_continuous(visualnames[j])+ #, breaks=seq(0,600,50)) +
            #coord_cartesian(ylim = c(min(ggdata$Value),max(ggdata$Value))) +
            geom_point(aes(x=Trial, y=Value), size = 3, col="#FF8888") +
            geom_line(aes(x=Trial, y=Value), size = 1, col="#FF8888")
    }
     
     assign(paste0("p",id), p)
  }
  

 #pdf(paste0("/Users/jisu/Documents/Hidden Target/Visualizations/",visualnames[j],"_indiv_crit.pdf"), width = 15, height = 10)
 #multiplot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p20, p0, cols=4)
 #dev.off()
 
  

 png(paste0("/Users/jisu/Documents/Hidden Target/Visualizations/",visualnames[j],"_indiv_crit.png"), width = 4000, height = 3200, res=220)
 multiplot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p0, cols=3) 
 dev.off()
 
  }




visu.data <- visu.data[,-7]


for (i in c(1:20,0)){
  if (i != 0){
    subdata <- visu.data[visu.data$ID == i,]
    subdata2 <- data.frame(subdata[,c(4:10)])
    subcor <- cor(subdata2)
    subcor[is.na(subcor)] <- 0
    print(subcor)
    subcormelt <- reshape2::melt(subcor)
    p <- ggplot(data = subcormelt, aes(x=Var1, y=Var2, fill=value)) + 
      geom_tile() +
      scale_x_discrete(paste0("Subject",i)) +
      scale_y_discrete("") + 
      scale_fill_gradient2(low = "#3B9AB2",
                           mid = "#EEEEEE", 
                           high = "#f55e4c",
                           limit = c(-0.5,1))
    
    assign(paste0("p",i), p)
  } else {
    subdata <- visu.data
    subdata2 <- data.frame(subdata[,c(4:10)])
    subcor <- cor(subdata2)
    subcormelt <- reshape2::melt(subcor)
    p <- ggplot(data = subcormelt, aes(x=Var1, y=Var2, fill=value)) + 
      geom_tile() +
      scale_x_discrete("Average") +
      scale_y_discrete("") + 
      scale_fill_gradient2(low = "#3B9AB2",
                           mid = "#EEEEEE", 
                           high = "#f55e4c",
                           limit = c(-0.5,1))
    
    assign(paste0("p",i), p)
  }
}

png(paste0("/Users/jisu/Documents/Hidden Target/Visualizations/Correlation matrix.png"), width = 2000, height = 1600, res=220)
multiplot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p0, cols=5)
dev.off()



for (i in c(1:20,0)){

subdata <- visu.data[visu.data$ID == i,]
subdata2 <- data.frame(subdata[,c(4:10)])
subdata2$search <- rep(1:360)
#subcor <- cor(subdata2)
#subcor[is.na(subcor)] <- 0
#print(subcor)

max.ratio <- max(subdata2$LR) / max(subdata2$ent)


p <- ggplot(subdata2) +
  coord_cartesian(ylim = c(-2, 2)) +
  #geom_point(aes(x=search, y=LR , col = "#3B9AB2")) +
  geom_line(aes(x=search,  y=LR , col = "#3B9AB2"))  +
  #geom_point(aes(x=search, y=ent * max.ratio, col = "#f55e4c")) +
  geom_line(aes(x=search,  y=ent * max.ratio, col = "#f55e4c")) +
  coord_cartesian(ylim = c(-4,4)) +
  scale_y_continuous("LR(Red)", 
                     sec.axis = (sec_axis(~ ./max.ratio, name = "Ent(Green)"))) +
  scale_x_continuous(paste0("Subject",i), breaks = seq(1,360,10)) +
  theme(legend.position="none")
assign(paste0("p",i), p)
}




png(paste0("/Users/jisu/Documents/Hidden Target/Visualizations/LRandENT.png"), width = 4000, height = 3000, res=220)
multiplot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, cols=2)
dev.off()






##############################
### Linear Trend per Search ###
##############################

#
#for (j in c(1)){
#  variablenames = c("subdata$traj", "subdata$LR", "subdata$RT", "subdata$label", "subdata$rew", "subdata$ent", "subdata$RP", "subdata$RPE")
#  visualnames = c("Trajectory", "Likelihood Ratio", "Reaction Time", "Exploit / Explore", "Reward Rate", "Entropy", "Reward Probability", "Reward Prediction Error")
#  y <- matrix(NA, ncol=2, nrow=200)
#  
#  for (i in subjectpool){
#    subdata <- visu.data[visu.data$ID == i,]
#  }
#  
#  for (id in c(1:18,20,0)){
#    if (id != 0){
#      p <- ggplot(data = subdata[subdata$ID == id,]) +
#        scale_x_continuous(paste0("Subject",id), breaks = seq(1, 360, 10)) +
#        scale_y_continuous(visualnames[j])+#, breaks=seq(-1,1,0.25)) +
#        #coord_cartesian(ylim = c(min(ggdata$Value),max(ggdata$Value))) +
#        geom_line(aes(x=search,  y=LR))
#    } else {
#      p <- ggplot(data = subdata[subdata$ID == id,]) +
#        scale_x_continuous("Average", breaks = seq(1, 360, 10)) +
#        scale_y_continuous(visualnames[j])+ #, breaks=seq(0,600,50)) +
#        #coord_cartesian(ylim = c(min(ggdata$Value),max(ggdata$Value))) +
#        #geom_point(aes(x=search, y=traj), size = 3, col="#FF8888") +
#        geom_line(aes(x=search,  y=LR), size = 1, col="#FF8888")
#    }
#
#    assign(paste0("p",id), p)
#  }
#  
#  
#  pdf(paste0("/Users/jisu/Documents/Hidden Target/Visualizations/",visualnames[j],"_indiv_crit.pdf"), width = 20, height = 10)
#  multiplot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p20, p0, cols=4)
#  dev.off()
#  
#  #png(paste0("/Users/jisu/Documents/Hidden Target/Visualizations/",visualnames[j],"_indiv_crit.png"), width = 4000, height = 3200, res=220)
#  #multiplot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p20, p0, cols=4)
#  #dev.off()
#  
#}
#
#data$ID <- rep(subjectpool, each = 360)
#subdata <- data[data$ID ==18, ]
#subdata <- subdata[subdata$V2 == 1, ]
#subdata$V6 <- subdata$V6/16
#subdata$V7 <- subdata$V7/16
#
#RPdata$ID <- rep(subjectpool,each =360)
#RPdata$trial <- rep(c(1:10))
#RPsub <- RPdata[(RPdata$ID == 18) & (RPdata$trial == 1),]
#
#######################
#### Average Trend ###
#######################
#
#
#subjectpool <- c(1:18,20)
#nsubject <- 19
#for (subjectID in subjectpool){
#  
#  if (subjectID < 10){
#    filename <- paste0("HT0",subjectID)
#  } else {filename <- paste0("HT",subjectID)
#  }
#  
#  dat <- read.csv(paste0("/Volumes/clmnlab/HT/behavior_data/raw_data/",filename,"/",filename,".csv"), header=FALSE)
#  
#  if (subjectID == 1){
#    data <- data.frame(dat)
#  } else {
#    data <- rbind(data, dat)
#  }
#}
#
#dim(data)
#names(data) 
#head(likdata)
#likdata$label<- ifelse(likdata$label1 == 1, "Exploit", "Explore")







subjectpool <- c(5:18,20)
for (subjectID in subjectpool){
  
  if (subjectID < 10){
    filename <- paste0("HT0",subjectID)
  } else {filename <- paste0("HT",subjectID)
  }
  
  tjdat <- read.csv(paste0("/Volumes/clmnlab/HT/behavior_data/raw_data/",filename,"/",filename,"_trajdistsum.csv"))
  
  if (subjectID == 5){
    tjdata <- data.frame(tjdat)
  } else {
    tjdata <- rbind(tjdata, tjdat)
  }
}


visu.data <- data.frame("ID" = likdata$ID, "label" = likdata$label, "RT" = data$V10, "trial" = rep(1:10))
visu.data$ones <- 1
RT.mean <- tapply(visu.data$RT, visu.data$trial, mean)

label.sum <- matrix(NA, ncol=2, nrow=10)
for (i in 1:10){
  label.sum[i,1] <- sum((visu.data$label == "Exploit") & (visu.data$trial == i))
  label.sum[i,2] <- sum((visu.data$label == "Explore") & (visu.data$trial == i))
}

data2 <- data.frame("RT" = RT.mean/1000, "Exploit" = label.sum[,1], "Explore" = label.sum[,2], "ratio" = (label.sum[,1] /684), "trial" = c(1:10))



ggplot(data2) +
  geom_bar(aes(x=trial, y=ratio), stat = "Identity", position=position_dodge()) +
  geom_point(aes(x=trial, y=RT -1.25), col="#FF3333", size=3) +#, stat = "Identity", position=position_dodge(), col="red") +
  geom_line(aes(x=trial, y=RT -1.25), col="#FF3333", size=0.5) +
  scale_x_continuous("Trial", breaks = c(1:10)) +
  scale_y_continuous("Ratio of Exploitation (gray)", limits = c(0,1), 
                     sec.axis = (sec_axis(~ . + 1.25, name = "Reaction Time (red)")))



visu.data2 <- visu.data[visu.data$ID > 4, ]
visu.data2 <- cbind(visu.data2, "traj" = tjdata$x)
RT.mean2 <- tapply(visu.data2$RT, visu.data2$trial, mean)
Traj.mean2 <- tapply(visu.data2$traj, visu.data2$trial, mean)

data3 <- data.frame("RT" = RT.mean, "trajectory" = Traj.mean2, "trial" = c(1:10))
max.ratio <-  max(data3$trajectory) / max(data3$RT)

ggplot(data3) +
  geom_point(aes(x=trial, y=RT * max.ratio), col="#FF8888", size=3) +#, stat = "Identity", position=position_dodge(), col="red") +
  geom_line(aes(x=trial,  y=RT * max.ratio), col="#FF8888", size=0.5) +
  geom_point(aes(x=trial, y=trajectory), col="#55CC99", size=3) +#, stat = "Identity", position=position_dodge(), col="red") +
  geom_line(aes(x=trial,  y=trajectory), col="#55CC99", size=0.5) +
  scale_x_continuous("Trial", breaks = c(1:10)) +
  scale_y_continuous("Length of Trajectory (green)", limits = c(350,450), 
                     sec.axis = (sec_axis(~ ./max.ratio, name = "Recation Time (red)"))) 
  
  

