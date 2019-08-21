library(matlib)
library(pracma)



subjectpool <- c(22, 23, 24, 25, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)#, 22, 23, 24, 25)
nsubject <- length(subjectpool)
for (subjectID in subjectpool){
  

    
    if (subjectID < 10){
      filename <- paste0("HT0",subjectID)
    } else {filename <- paste0("HT",subjectID)
    }
    dat <- read.csv(paste0("/Volumes/clmnlab/HT/behavior_data/raw_data/",filename,"/",filename,"_traj.csv"), header = FALSE)
    dim(dat)
    colnames(dat) <- rep(1:360, each=2)
    distdata <- matrix(NA,nrow=200, ncol=360)
    
    for (j in 1:360){
      for (i in 1:199){ 
        #print(c(i,j))
        tryCatch(
          {
            datset <- rbind(dat[i, (2*j-1):(2*j)], dat[(i+1),(2*j-1):(2*j)])
            distdata[i,j] <- dist(datset, method = "euclidean")
            #print(datset
              },
            error = function(e){
            #distdata[i,j] <- NA
            #print(c(i,j))
        }
        )
        sum.dist <- apply(distdata, 2, sum, na.rm=TRUE)
        write.csv(distdata, paste0("/Volumes/clmnlab/HT/behavior_data/raw_data/",filename,"/",filename,"_trajdist.csv"))
        write.csv(sum.dist, paste0("/Volumes/clmnlab/HT/behavior_data/raw_data/",filename,"/",filename,"_trajdistsum.csv"))
      }
    }
  
  dat <- read.csv(paste0("/Volumes/clmnlab/HT/behavior_data/raw_data/",filename,"/",filename,"_trajdistsum.csv"))
  traj <- dat[,2]
  
  dat <- read.csv(paste0("/Volumes/clmnlab/HT/behavior_data/raw_data/",filename,"/",filename,".csv"), header=FALSE)
  RT <- dat$V10

  rtraj <- (cbind(RT, traj))

  orth2 <- gramSchmidt(rtraj)
  orth2$Q[,1]%*%orth2$Q[,2]
  rtraj[,1] <- orth2$Q[,1]
  rtraj[,2] <- orth2$Q[,2]
  rtraj <- data.frame(rtraj)
  rtraj$ID <- rep(subjectID)
  

  writexl::write_xlsx(subdata, paste0("/Volumes/clmnlab/HT/regressors/RTandTRAJ/",filename,"_rtraj.xlsx"))

}


