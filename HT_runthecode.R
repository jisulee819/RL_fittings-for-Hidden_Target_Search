subjectID <- 2
filename <- paste0("pilot0",subjectID)
dir.create(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename))
dir.create(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename,"/posterior"))
dir.create(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename,"/entropy difference"))
datasetname <- paste0(filename,"_mixed.fitted")
#assign(datasetname, matrix(NA,nrow=1,ncol=8))
pilot02_mixed.fitted <- matrix(NA,nrow=1,ncol=11)

for (r in 1:36){
  
  data <- modelfitting_MIX(25/2,  5/6, 40/2, r, 4)
  fitted.data <- cbind(data$inference[,1],data$inference[,2],data$`mig history`[,1],data$`mig history`[,2],data$`map history`[,1],data$`map history`[,2],data$`prob reward history`,data$`reward history`,data$`inference type`, data$`error ratio`, data$`fitting error`)
  pilot02_mixed.fitted <- rbind(pilot02_mixed.fitted, fitted.data)
  if (r == 1){
    pilot02_mixed.fitted
  }
  #assign(datasetname, rbind(datasetname, fitted.data))
  
  pdf(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/pilot02/posterior/Round",r,".pdf"), width = 15, height = 10)
  gg_multiheatmap_12_mixmodel(data$`prior history`)
  dev.off()
  pdf(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/pilot02/entropy difference/Round",r,".pdf"), width = 15, height = 10)
  gg_multiheatmap_12(data$`entropy diff`,"blue")
  dev.off()
}
#assign(datasetname, datasetname[-1,])
pilot02_mixed.fitted <- pilot02_mixed.fitted[-1,]
write.csv(pilot02_mixed.fitted,"/Users/jisu/Documents/Hidden Target/Model Fit/pilot02/mixed_fitted_data.csv")


subjectID <- 3
filename <- paste0("pilot0",subjectID)
dir.create(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename))
dir.create(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename,"/posterior"))
dir.create(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename,"/entropy difference"))
datasetname <- paste0(filename,"_mixed.fitted")
#assign(datasetname, matrix(NA,nrow=1,ncol=8))
pilot03_mixed.fitted <- matrix(NA,nrow=1,ncol=11)

for (r in 1:36){
  
  data <- modelfitting_MIX(25/2,  5/6, 40/2, r, 4)
  fitted.data <- cbind(data$inference[,1],data$inference[,2],data$`mig history`[,1],data$`mig history`[,2],data$`map history`[,1],data$`map history`[,2],data$`prob reward history`,data$`reward history`,data$`inference type`, data$`error ratio`, data$`fitting error`)
  pilot03_mixed.fitted <- rbind(pilot03_mixed.fitted, fitted.data)
  #assign(datasetname, rbind(datasetname, fitted.data))
  
  pdf(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/pilot03/posterior/Round",r,".pdf"), width = 15, height = 10)
  gg_multiheatmap_12_mixmodel(data$`prior history`)
  dev.off()
  pdf(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/pilot03/entropy difference/Round",r,".pdf"), width = 15, height = 10)
  gg_multiheatmap_12(data$`entropy diff`,"blue")
  dev.off()
}
#assign(datasetname, datasetname[-1,])
pilot03_mixed.fitted <- pilot03_mixed.fitted[-1,]
write.csv(pilot03_mixed.fitted,"/Users/jisu/Documents/Hidden Target/Model Fit/pilot03/mixed_fitted_data.csv")


subjectID <- 4
filename <- paste0("pilot0",subjectID)
dir.create(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename))
dir.create(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename,"/posterior"))
dir.create(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename,"/entropy difference"))
datasetname <- paste0(filename,"_mixed.fitted")
#assign(datasetname, matrix(NA,nrow=1,ncol=8))
pilot04_mixed.fitted <- matrix(NA,nrow=1,ncol=11)

for (r in 1:36){
  
  data <- modelfitting_MIX(25/2,  5/6, 40/2, r, 4)
  fitted.data <- cbind(data$inference[,1],data$inference[,2],data$`mig history`[,1],data$`mig history`[,2],data$`map history`[,1],data$`map history`[,2],data$`prob reward history`,data$`reward history`,data$`inference type`, data$`error ratio`, data$`fitting error`)
  pilot04_mixed.fitted <- rbind(pilot04_mixed.fitted, fitted.data)
  #assign(datasetname, rbind(datasetname, fitted.data))
  
  pdf(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/pilot04/posterior/Round",r,".pdf"), width = 15, height = 10)
  gg_multiheatmap_12_mixmodel(data$`prior history`)
  dev.off()
  pdf(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/pilot04/entropy difference/Round",r,".pdf"), width = 15, height = 10)
  gg_multiheatmap_12(data$`entropy diff`,"blue")
  dev.off()
}
#assign(datasetname, datasetname[-1,])
pilot04_mixed.fitted <- pilot04_mixed.fitted[-1,]
write.csv(pilot04_mixed.fitted,"/Users/jisu/Documents/Hidden Target/Model Fit/pilot04/mixed_fitted_data.csv")


subjectID <- 5
filename <- paste0("pilot0",subjectID)
dir.create(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename))
dir.create(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename,"/posterior"))
dir.create(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename,"/entropy difference"))
datasetname <- paste0(filename,"_mixed.fitted")
#assign(datasetname, matrix(NA,nrow=1,ncol=8))
pilot05_mixed.fitted <- matrix(NA,nrow=1,ncol=11)

for (r in 1:36){
  
  data <- modelfitting_MIX(25/2,  5/6, 40/2, r, 4)
  fitted.data <- cbind(data$inference[,1],data$inference[,2],data$`mig history`[,1],data$`mig history`[,2],data$`map history`[,1],data$`map history`[,2],data$`prob reward history`,data$`reward history`,data$`inference type`, data$`error ratio`, data$`fitting error`)
  pilot05_mixed.fitted <- rbind(pilot05_mixed.fitted, fitted.data)
  #assign(datasetname, rbind(datasetname, fitted.data))
  
  pdf(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/pilot05/posterior/Round",r,".pdf"), width = 15, height = 10)
  gg_multiheatmap_12_mixmodel(data$`prior history`)
  dev.off()
  pdf(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/pilot05/entropy difference/Round",r,".pdf"), width = 15, height = 10)
  gg_multiheatmap_12(data$`entropy diff`,"blue")
  dev.off()
}
#assign(datasetname, datasetname[-1,])
pilot05_mixed.fitted <- pilot05_mixed.fitted[-1,]
write.csv(pilot05_mixed.fitted,"/Users/jisu/Documents/Hidden Target/Model Fit/pilot05/mixed_fitted_data.csv")

