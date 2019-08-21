
subjectID <- 1
filename <- paste0("HT0",subjectID)
dir.create(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename))
dir.create(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename,"/posterior"))
dir.create(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename,"/entropy difference"))
datasetname <- paste0(filename,"_mixed.fitted")
#assign(datasetname, matrix(NA,nrow=1,ncol=8))
pilot01_mixed.fitted <- matrix(NA,nrow=1,ncol=11)

for (r in 1:36){
  
  data <- modelfitting_MIX(12, 0.78, r, 1)
  fitted.data <- cbind(data$inference[,1],data$inference[,2],data$`mig history`[,1],data$`mig history`[,2],data$`map history`[,1],data$`map history`[,2],data$`prob reward history`,data$`reward history`,data$`inference type`, data$`error ratio`, data$`fitting error`)
  pilot01_mixed.fitted <- rbind(pilot01_mixed.fitted, fitted.data)
  #assign(datasetname, rbind(datasetname, fitted.data))
  
  pdf(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/HT01/posterior/Round",r,".pdf"), width = 13, height = 20)
  gg_multiheatmap_10_mixmodel(data$`prior history`)
  dev.off()
  pdf(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/HT01/entropy difference/Round",r,".pdf"), width = 13, height = 20)
  gg_multiheatmap_10(data$`entropy diff`,"blue")
  dev.off()
}
#assign(datasetname, datasetname[-1,])
pilot01_mixed.fitted <- pilot01_mixed.fitted[-1,]
write.csv(pilot01_mixed.fitted,"/Users/jisu/Documents/Hidden Target/Model Fit/HT01/HT01_fitted.csv")


subjectID <- 2
filename <- paste0("HT0",subjectID)
dir.create(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename))
dir.create(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename,"/posterior"))
dir.create(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename,"/entropy difference"))
datasetname <- paste0(filename,"_mixed.fitted")
#assign(datasetname, matrix(NA,nrow=1,ncol=8))
pilot02_mixed.fitted <- matrix(NA,nrow=1,ncol=11)

for (r in 1:36){
  
  data <- modelfitting_MIX(12, 0.78, r, 2)
  fitted.data <- cbind(data$inference[,1],data$inference[,2],data$`mig history`[,1],data$`mig history`[,2],data$`map history`[,1],data$`map history`[,2],data$`prob reward history`,data$`reward history`,data$`inference type`, data$`error ratio`, data$`fitting error`)
  pilot02_mixed.fitted <- rbind(pilot02_mixed.fitted, fitted.data)
  #assign(datasetname, rbind(datasetname, fitted.data))
  
  pdf(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/HT02/posterior/Round",r,".pdf"), width = 13, height = 20)
  gg_multiheatmap_10_mixmodel(data$`prior history`)
  dev.off()
  pdf(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/HT02/entropy difference/Round",r,".pdf"), width = 13, height = 20)
  gg_multiheatmap_10(data$`entropy diff`,"blue")
  dev.off()
}
#assign(datasetname, datasetname[-1,])
pilot02_mixed.fitted <- pilot02_mixed.fitted[-1,]
write.csv(pilot02_mixed.fitted,"/Users/jisu/Documents/Hidden Target/Model Fit/HT02/HT02_fitted.csv")


subjectID <- 3
filename <- paste0("HT0",subjectID)
dir.create(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename))
dir.create(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename,"/posterior"))
dir.create(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename,"/entropy difference"))
datasetname <- paste0(filename,"_mixed.fitted")
#assign(datasetname, matrix(NA,nrow=1,ncol=8))
pilot03_mixed.fitted <- matrix(NA,nrow=1,ncol=11)

for (r in 1:36){
  
  data <- modelfitting_MIX(12, 0.78, r, 3)
  fitted.data <- cbind(data$inference[,1],data$inference[,2],data$`mig history`[,1],data$`mig history`[,2],data$`map history`[,1],data$`map history`[,2],data$`prob reward history`,data$`reward history`,data$`inference type`, data$`error ratio`, data$`fitting error`)
  pilot03_mixed.fitted <- rbind(pilot03_mixed.fitted, fitted.data)
  #assign(datasetname, rbind(datasetname, fitted.data))
  
  pdf(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/HT03/posterior/Round",r,".pdf"), width = 13, height = 20)
  gg_multiheatmap_10_mixmodel(data$`prior history`)
  dev.off()
  pdf(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/HT03/entropy difference/Round",r,".pdf"), width = 13, height = 20)
  gg_multiheatmap_10(data$`entropy diff`,"blue")
  dev.off()
}
#assign(datasetname, datasetname[-1,])
pilot03_mixed.fitted <- pilot03_mixed.fitted[-1,]
write.csv(pilot03_mixed.fitted,"/Users/jisu/Documents/Hidden Target/Model Fit/HT03/HT03_fitted.csv")


subjectID <- 4
filename <- paste0("HT0",subjectID)
dir.create(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename))
dir.create(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename,"/posterior"))
dir.create(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename,"/entropy difference"))
datasetname <- paste0(filename,"_mixed.fitted")
#assign(datasetname, matrix(NA,nrow=1,ncol=8))
pilot04_mixed.fitted <- matrix(NA,nrow=1,ncol=11)

for (r in 1:36){
  
  data <- modelfitting_MIX(12, 0.78, r, 4)
  fitted.data <- cbind(data$inference[,1],data$inference[,2],data$`mig history`[,1],data$`mig history`[,2],data$`map history`[,1],data$`map history`[,2],data$`prob reward history`,data$`reward history`,data$`inference type`, data$`error ratio`, data$`fitting error`)
  pilot04_mixed.fitted <- rbind(pilot04_mixed.fitted, fitted.data)
  #assign(datasetname, rbind(datasetname, fitted.data))
  
  pdf(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/HT04/posterior/Round",r,".pdf"), width = 13, height = 20)
  gg_multiheatmap_10_mixmodel(data$`prior history`)
  dev.off()
  pdf(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/HT04/entropy difference/Round",r,".pdf"), width = 13, height = 20)
  gg_multiheatmap_10(data$`entropy diff`,"blue")
  dev.off()
}
#assign(datasetname, datasetname[-1,])
pilot04_mixed.fitted <- pilot04_mixed.fitted[-1,]
write.csv(pilot04_mixed.fitted,"/Users/jisu/Documents/Hidden Target/Model Fit/HT04/HT04_fitted.csv")



subjectID <- 5
filename <- paste0("HT0",subjectID)
dir.create(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename))
dir.create(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename,"/posterior"))
dir.create(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename,"/entropy difference"))
datasetname <- paste0(filename,"_mixed.fitted")
#assign(datasetname, matrix(NA,nrow=1,ncol=8))
pilot04_mixed.fitted <- matrix(NA,nrow=1,ncol=11)

for (r in 1:36){
  
  data <- modelfitting_MIX(12, 0.78, r, 5)
  fitted.data <- cbind(data$inference[,1],data$inference[,2],data$`mig history`[,1],data$`mig history`[,2],data$`map history`[,1],data$`map history`[,2],data$`prob reward history`,data$`reward history`,data$`inference type`, data$`error ratio`, data$`fitting error`)
  pilot04_mixed.fitted <- rbind(pilot04_mixed.fitted, fitted.data)
  #assign(datasetname, rbind(datasetname, fitted.data))
  
  pdf(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/HT04/posterior/Round",r,".pdf"), width = 13, height = 20)
  gg_multiheatmap_10_mixmodel(data$`prior history`)
  dev.off()
  pdf(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/HT04/entropy difference/Round",r,".pdf"), width = 13, height = 20)
  gg_multiheatmap_10(data$`entropy diff`,"blue")
  dev.off()
}
#assign(datasetname, datasetname[-1,])
pilot04_mixed.fitted <- pilot04_mixed.fitted[-1,]
write.csv(pilot04_mixed.fitted,"/Users/jisu/Documents/Hidden Target/Model Fit/HT04/HT04_fitted.csv")