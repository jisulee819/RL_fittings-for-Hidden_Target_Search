
# For fitting MAP/MIG and saving posterior heatmap
# Modified 190507 by JISU
# Modified 190823 by JISU (ver.13)


library(MASS)
library(readr)
library(mvtnorm)
library(corpcor)
library(Hmisc)
library(reshape2)
library(ggplot2)


#### Model Fitting : saving 1.inferences from the models 2.likelihood ####

nx=32
ny=24
n=10

modelfitting_MIX <- function(s,lambda, roundn, subjID){
  
  #### Prior Distribution ####
  # prior of the reward probability is analytically calculated in advance 
  pr <- read.csv("/Volumes/clmnlab/HT/behavior_data/raw_data/priorreduced2.csv")
  prior <<- data.matrix(pr[,-1])
  
  #### Settings ####
  sigma <- matrix(c(s^2, 0, 0, s^2),nrow=2)
  rho <- 0
  ratio <- 1024/nx
  
  subdata <- behavdata[((roundn-1)*10+1):(roundn*10),]
  treasure_location <- c((subdata$treasureX[10]/ratio)+nx/2,(subdata$treasureY[10]/ratio)+ny/2)
  
  
  prior_history <- prior
  entropy_diff_history <- rep(NA, ny)
  infer_history <- c(NA, NA)
  rew_history <- 0
  prob_reward_history <- NA
  likelihood <- NA
  lkh.ratio <- NA
  mig.history <- matrix(NA,nrow=10,ncol=2)
  migs.history <- matrix(0, nrow = 10*10, ncol=2)
  map.history <- matrix(NA,nrow=10,ncol=2)
  map.likelihood <- NA
  mig.likelihood <- NA
  z.map <- NA
  z.mig <- NA
  l.map <- NA
  l.mig <- NA
  f.map <- NA
  f.mig <- NA
  
  n = 10
  
  #### Iteration ####
  for (i in 1:10){
    print(paste0("Trial ",i))
    
    #### Radius and Reward ####
    if (i == 1){rad <- subdata$Radius[1]/16
    } else {rad <- subdata$Radius[i-1]/16}
    
    rew <- subdata$Smiley[i]
    rew_history[i] <- rew
    
    #### Truncate Inferences ####
    infer <- c((subdata$inferX[i]/ratio)+nx/2, (subdata$inferY[i]/ratio)+ny/2)
    if (infer[1] < 1){
      infer[1] = 1
    } else if (infer[1] > 32){
      infer[1] = 32
    }
    
    if (infer[2] < 1){
      infer[2] = 1
    } else if (infer[2] > 24) {
      infer[2] = 24
    }
    
    infer[1] <- round(infer[1])
    infer[2] <- round(infer[2])
    
    
    #### MIG fitting and saving ####
    mig_dataset <- inference_as(prior, sigma, rad)
    infer_model.migs <- mig_dataset$`inference mig`
    infer_model.mig <- infer_model.migs[1,]
    mig.history[i,] <- infer_model.mig
    migs.history[((i-1)*10+1):(i*10),] <- infer_model.migs
    
    entropy_diff_normalized <- mig_dataset$`entropy_difference`/sum(mig_dataset$`entropy_difference`)
    entropy_diff_history <- rbind(entropy_diff_history, entropy_diff_normalized)
    
    
    #### Probability of Reward ####
    
    checker <- matrix(NA, nrow=nx, ncol=ny)
    for(x in 1:nx){for(y in 1:ny){
        if (sqrt((x-infer[1])^2 + (y-infer[2])^2) < sqrt(sigma[1,1])){ 
          checker[x,y] <- TRUE
        } else {checker[x,y] <- FALSE}}}
      
    prob_reward <- sum(checker * prior)
    prob_reward_history[i] <- prob_reward  
    
    
    #### MAP fitting and saving ####
    infer_model.map <- inference_bs(prior) 
    map.history[i,] <- infer_model.map
    
    #migdist <- dist(rbind(infer, infer_model.mig))
    #mapdist <- dist(rbind(infer, infer_model.map))

    
    #### Bayesian Update of Prior ####
    posterior <- bayesian_update(prior, infer, sigma, rew)
    posterior <- posterior / sum(posterior)
    
    mat.prior <- melt(prior)
    mat.ent <- melt(entropy_diff_normalized)
    
    
    #### Calculate the Index of Model fitting ####
    # z.index: prior & entropy difference in the normalized space
    # l.index: prior & entropy difference linearly scaled btw [-1,1]
    # f.index: prior & entropy difference Fisher transformed
    
    z.prior <- (prior - mean(mat.prior[,3])) / sd(mat.prior[,3])
    z.ent <- (entropy_diff_normalized - mean(mat.ent[,3])) / sd(mat.ent[,3])
    
    z.map[i] <- z.prior[infer[1], infer[2]]
    z.mig[i] <- z.ent[infer[1], infer[2]]

    mat.z.prior <- melt(z.prior)
    mat.z.ent <- melt(z.ent)
    
    l.map[i] <- (2 / (max(mat.z.prior[,3]) - min(mat.z.prior[,3])))  * (z.map[i] - min(mat.z.prior[,3])) - 1
    l.mig[i] <- (2 / (max(mat.z.ent[,3]) - min(mat.z.ent[,3])))  * (z.mig[i] - min(mat.z.ent[,3])) - 1
    
    f.map[i] <- (1/2) * log((1+l.map[i])/(1-l.map[i])) 
    f.mig[i] <- (1/2) * log((1+l.mig[i])/(1-l.mig[i])) 
   
    map.lkh <- log(prior[infer[1],infer[2]]+0.0000001)
    mig.lkh <- log(entropy_diff_normalized[infer[1],infer[2]]+0.0000001)
    isEmpty <- function(x) {
    return(length(x)==0)
    }
    
    # index.lkh : prior & entropy difference in the original space
    
    map.likelihood[i] <- map.lkh
    mig.likelihood[i] <- mig.lkh
    
    if (!isEmpty(map.lkh - mig.lkh)){
    lkh.ratio[i] <- (map.lkh - mig.lkh) 
    } else {lkh.ratio[i] <- NA}
    
    prior <- posterior
    prior_history <- rbind(prior_history, data.matrix(prior))
    infer_history <- rbind(infer_history, infer)
    s <- s * lambda^(rew)
    sigma <- matrix(c(s^2, 0, 0, s^2),nrow=2) 
    
  }

  #### Save ####
  data = list('treasureX' = treasure_location[1],'treasureY' = treasure_location[2], "prior history" = prior_history, 'inference' = infer_history[-1,],
              'reward history' = rew_history, 'prob reward history' = prob_reward_history, "sigma" = s, "lambda" = lambda, 'entropy diff' = entropy_diff_history[-1,],
              'mig history' = mig.history, 'migs.history' = migs.history, 'map history' = map.history, 'lkh ratio' = lkh.ratio, 'map.likelihood' = map.likelihood, 'mig.likelihood' = mig.likelihood,
              'prior_history' = prior_history, 'entropy_diff_history' = entropy_diff_history, 'z.map' = z.map, 'z.mig' = z.mig, 'f.map' = f.map, 'f.mig' = f.mig)
  
  return(data)
}


#################################
#################################
#################################
######### MIG Inference #########
#################################
#################################
#################################


inference_as <- function(prior, sigma, radius){
  
  #### Setting ####
  entropy_prior <- NA
  entropy_post_w_rew <- matrix(NA, nrow=nx, ncol=ny)
  entropy_post_wo_rew <- matrix(NA, nrow=nx, ncol=ny)
  entropy_post <- matrix(NA, nrow=nx, ncol=ny)
  entropy_diff <- matrix(NA, nrow=nx, ncol=ny)
  prob_reward <-  matrix(0, nrow=nx, ncol=ny)
  prob_no_reward <- matrix(0, nrow=nx, ncol=ny)
  expected_post <-  matrix(NA, nrow=nx, ncol=ny)
  exp_entropy <- matrix(NA, nrow=nx, ncol=ny)
  
  
  #### Probability of Reward ####
  # which is the weight of expected value
  for(a in 1:nx){for(b in 1:ny){
  
    checker <- matrix(NA, nrow=nx, ncol=ny)
  
    for(i in 1:nx){for(j in 1:ny){
      if (sqrt((i-a)^2 + (j-b)^2) < sqrt(sigma[1,1])){ 
        checker[i,j] <- TRUE
      } else {checker[i,j] <- FALSE
      }}}
    
    prob_reward[a,b] <- sum(checker * prior)
    prob_no_reward[a,b] <- 1-prob_reward[a,b]
  }}
  

  #### Entropy of prior ####
  entropy_prior <- sum(prior*log2(1/prior), na.rm=TRUE)
  
  #### Entropy of posterior ####
  for(a in 1:nx){for(b in 1:ny){
 
      inference <- c(a,b)
      likelihood_w_rew <- matrix(data = NA, nrow = nx, ncol = ny, byrow = FALSE, dimnames = NULL)
      posterior_w_rew <- matrix(data = NA, nrow = nx, ncol = ny, byrow = FALSE, dimnames = NULL)
      posterior_wo_rew <- matrix(data = NA, nrow = nx, ncol = ny, byrow = FALSE, dimnames = NULL)
      
      #### Likelihood w/ reward ####
      # defined as gaussian distribution
      for(i in 1:nx){for(j in 1:ny){
          likelihood_w_rew[i,j] <- dmvnorm(c(i,j), c(a,b), sigma)
        }}  
      
      #### Likelihood wo/ reward ####
      # defined as inverted gaussian
      likelihood_wo_rew <- (max(likelihood_w_rew)-likelihood_w_rew)
      
      #### Posteriors in case of reward / non-reward ####
      posterior_w_rew <- likelihood_w_rew * prior
      posterior_wo_rew <- (likelihood_wo_rew) * prior

      posterior_w_rew <- posterior_w_rew  / sum(posterior_w_rew)
      posterior_wo_rew <- posterior_wo_rew  / sum(posterior_wo_rew)
      
      #### Entropy of Posteriors in case of reward / non-reward ####
      entropy_post_w_rew[a,b] <- sum(posterior_w_rew*log2(1/posterior_w_rew), na.rm=TRUE)
      entropy_post_wo_rew[a,b] <- sum(posterior_wo_rew*log2(1/posterior_wo_rew), na.rm=TRUE)
      
      #### Expected Entropy ####
      entropy_post[a,b] <- prob_reward[a,b] * entropy_post_w_rew[a,b] + prob_no_reward[a,b] * entropy_post_wo_rew[a,b]
    
  }}

  
  #### Decide the next Inference
  entropy_diff <- (entropy_prior-entropy_post)
  rank.of.mig <- matrix(rank(-entropy_diff),nrow=nx)
  #inference_mig <- matrix(NA, nrow=10, ncol=2)    

  
  alpha <- which(rank.of.mig == 1, arr.ind = TRUE)
  if (nrow(alpha) > 1){
    alpha <- alpha[sample(1:nrow(alpha), 1),]
    }
  
  inference_mig <- rbind(which(rank.of.mig == 1, arr.ind = TRUE), which(1 < rank.of.mig & rank.of.mig < 11, arr.ind = TRUE))
  
  if (nrow(inference_mig) > 10){
    inference_mig <- inference_mig[1:10,]
    }
  
  if (nrow(inference_mig) < 10){
    namat <- matrix(NA, ncol=2, nrow = 10-nrow(inference_mig))
    inference_mig <- rbind(inference_mig, namat)
    }
  
  #### Return the results
  update_dataset <- list('inference' = inference, 'inference mig' = inference_mig, 'posterior_w_rew' = posterior_w_rew, 'posterior_wo_rew' = posterior_wo_rew, 
                         'entropy_difference' = entropy_diff)
  return(update_dataset)
  
}


###################################
###################################
###################################
######### Bayesian Update #########
###################################
###################################
###################################


bayesian_update <- function(prior, inference, sigma, reward){
  
  likelihood <- matrix(data = NA, nrow = nx, ncol = ny, byrow = FALSE, dimnames = NULL)
  posterior <- matrix(data = NA, nrow = nx, ncol = ny, byrow = FALSE, dimnames = NULL)
  
  for(i in 1:nx){
    for(j in 1:ny){
      mu <- matrix(c(i,j), nrow=2)
      likelihood[i,j] <- dmvnorm(t(inference), mu, sigma)}}
  maxl <- max(likelihood)
  if (reward == 0){
    likelihood <- maxl-likelihood}
  for(i in 1:nx){
    for(j in 1:ny){
      posterior[i,j] <- likelihood[i,j] * prior[i,j]}}
  return(posterior)
  
}

inference_bs <- function(prior){ ##
  inference <- which(prior == max(prior), arr.ind=TRUE) #
  inference <- inference[sample(1:nrow(inference), 1),]
  return(inference)
}

reward_bs <- function(inference, radius, treasure_location){ ##
  if (sqrt((treasure_location[1]-inference[1])^2 + (treasure_location[2]-inference[2])^2) < radius) return(1)
  else return(0)
}



###################################
###################################
###################################
######### Hitmap Plotting #########
###################################
###################################
###################################


gg_heatmap_mixmodel <- function(mat, i){
  dat <- mat[(nx*(i-1)+1):(nx*i),]
  dat <- data.matrix(dat)
  #names(dat) <- paste("X", 1:24) 
  colnames(dat) <- c(1:24)
  rownames(dat) <- c(1:32)
  dat2 <- melt(dat)
  #dat2 <- melt(dat), id.var = "X1")
  dat2$Var1 <- rep(1:nx, times = ny, length.out = NA, each = 1)
  shapes = c(10, 8)
  
  map <- data$`map history`[i,]
  migs <- data$migs.history[(10*(i-1)+1) : (10*i),]
  infer <- data$inference[i,]
  infer_model <- data$`inference of model`[i,]
  
  #dat2[dat2[,1]==19 & dat2[,2]==4,3]
  #dat[19,4]
  ggplot(dat2, aes(Var1, Var2, group=Var2)) +
    geom_tile(aes(fill = value)) + 
    labs(title = sprintf("Trial %s",i),
         x = "",
         y = "") +
    scale_fill_gradient2(low="white", high = "red",
                         oob = scales::squish) +#, midpoint=0.0001, limits=c(0, 0.005), 
    geom_point(size=1, color = "black", aes(x = c(infer[1]), y = c(infer[2]))) +
    
    geom_point(size=1, color = "yellow", aes(x = c(migs[1,1]), y = c(migs[1,2])))  +
    geom_point(size=1, color = "yellow", aes(x = c(migs[2,1]), y = c(migs[2,2])))  +
    geom_point(size=1, color = "yellow", aes(x = c(migs[3,1]), y = c(migs[3,2])))  +
    geom_point(size=1, color = "yellow", aes(x = c(migs[4,1]), y = c(migs[4,2])))  +
    geom_point(size=1, color = "yellow", aes(x = c(migs[5,1]), y = c(migs[5,2])))  +
    geom_point(size=1, color = "yellow", aes(x = c(migs[6,1]), y = c(migs[6,2])))  +
    geom_point(size=1, color = "yellow", aes(x = c(migs[7,1]), y = c(migs[7,2])))  +
    geom_point(size=1, color = "yellow", aes(x = c(migs[8,1]), y = c(migs[8,2])))  +
    geom_point(size=1, color = "yellow", aes(x = c(migs[9,1]), y = c(migs[9,2])))  +
    geom_point(size=1, color = "yellow", aes(x = c(migs[10,1]), y = c(migs[10,2])))  +
    
    
    geom_point(size=1, color = "green", aes(x = c(map[1]), y = c(map[2])))  
  
}



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

gg_multiheatmap_10_mixmodel <- function(mat){
  p1 <- gg_heatmap_mixmodel(mat, 1)
  p2 <- gg_heatmap_mixmodel(mat, 2)
  p3 <- gg_heatmap_mixmodel(mat, 3)
  p4 <- gg_heatmap_mixmodel(mat, 4)
  p5 <- gg_heatmap_mixmodel(mat, 5)
  p6 <- gg_heatmap_mixmodel(mat, 6)
  p7 <- gg_heatmap_mixmodel(mat, 7)
  p8 <- gg_heatmap_mixmodel(mat, 8)
  p9 <- gg_heatmap_mixmodel(mat, 9)
  p10 <- gg_heatmap_mixmodel(mat, 10)
  multiplot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, cols=2)}


#### TD HEATMAP ####
gg_heatmap <- function(matrix, i, col){
  dat <- matrix[(nx*(i-1)+1):(nx*i),]
  #names(dat) <- paste("X", 1:24) 
  colnames(dat) <- c(1:24)
  rownames(dat) <- c(1:32)
  dat2 <- melt(dat)
  #dat2 <- melt(dat, id.var = "X1")
  dat2$Var1 <- rep(1:nx, times = ny, length.out = NA, each = 1)
  #dat2[dat2[,1]==19 & dat2[,2]==4,3]
  #dat[19,4]
  ggplot(dat2, aes(Var1, Var2, group=Var2)) +
    geom_tile(aes(fill = value)) + 
    labs(title = sprintf("Trial %s",i),
         x = "",
         y = "") +
    scale_fill_gradient(low = "white", high = col) #,oob = scales::squish)
  
}

gg_multiheatmap_10 <- function(matrix, col){
  p1 <- gg_heatmap(matrix, 1, col)
  p2 <- gg_heatmap(matrix, 2, col)
  p3 <- gg_heatmap(matrix, 3, col)
  p4 <- gg_heatmap(matrix, 4, col)
  p5 <- gg_heatmap(matrix, 5, col)
  p6 <- gg_heatmap(matrix, 6, col)
  p7 <- gg_heatmap(matrix, 7, col)
  p8 <- gg_heatmap(matrix, 8, col)
  p9 <- gg_heatmap(matrix, 9, col)
  p10 <- gg_heatmap(matrix, 10, col)
  multiplot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, cols=2)}




##############################
##############################
##############################
######### Iterations #########
##############################
##############################
##############################



subjectpool <- c(6)

for (subjectID in subjectpool){
    if (subjectID < 10){
      filename <- paste0("HT0",subjectID)
    } else {filename <- paste0("HT",subjectID)}
    behavdata <- read.csv(paste0("/Volumes/clmnlab/HT/behavior_data/raw_data/",filename,"/",filename,".csv"), header=FALSE)
    names(behavdata) <- c("subjectID", "Trial", "Smiley", "treasureX", "treasureY", "inferX", "inferY", "Radius", "Distance", "RT", "Certainty")
    behavdata$subjectID <- subjectID
    behavdata$Round <- c(rep(1:36, each = 10))
    behavdata$number <- c(1:(length(behavdata$subjectID)))
    
    dir.create(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename))
    
    for (sigma1 in c(0.5, 0.6, 0.7, 0.8)){
      
      
    fitname <- paste0("sigma_",sigma1)
    
    dir.create(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename,"/",fitname))
    
    dir.create(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename,"/",fitname,"/posterior"))
    dir.create(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename,"/",fitname,"/entropy difference"))
    datasetname <- paste0(filename,"_mixed.fitted")
    pilot_mixed.fitted <- matrix(NA,nrow=1,ncol=15)
    prior.history <- matrix(NA,nrow=1,ncol=24)
    ent.history <- matrix(NA,nrow=1,ncol=24)
    
    for (r in 1:36){
      
      data <- modelfitting_MIX((ny/2)*sigma1, 0.78, r, subjectID)

      fitted.data <- cbind(data$inference[,1],data$inference[,2],data$`mig history`[,1],data$`mig history`[,2],data$`map history`[,1],data$`map history`[,2],
                           data$`prob reward history`,data$`reward history`, data$`lkh ratio`, data$`map.likelihood`, data$`mig.likelihood`, 
                           data$`z.map`, data$`z.mig`, data$`f.map`, data$`f.mig`)
      
      prior.history <- rbind(prior.history, data$`prior_history`)
      ent.history <- rbind(ent.history, data$`entropy_diff_history`[-1,])
      pilot_mixed.fitted <- rbind(pilot_mixed.fitted, fitted.data)
      
      dir.create(file.path(dirname(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename,"/",fitname,"/newfit_posterior/Round",r,".pdf"))))
      pdf(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename,"/",fitname,"/newfit_posterior/Round",r,".pdf"), width = 13, height = 20)
      gg_multiheatmap_10_mixmodel(data$`prior history`)
      dev.off()
      
      ### missing value errors are for missing mig points
      dir.create(file.path(dirname(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename,"/",fitname,"/newfit_entropy difference/Round",r,".pdf"))))
      pdf(paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename,"/",fitname,"/newfit_entropy difference/Round",r,".pdf"), width = 13, height = 20)
      gg_multiheatmap_10(data$`entropy diff`,"blue")
      dev.off()
    } 
    pilot_mixed.fitted <- pilot_mixed.fitted[-1,]
    prior.history <- prior.history[-1,]
    ent.history <- ent.history[-1,]
    
    label <- ifelse(pilot_mixed.fitted[,15] > pilot_mixed.fitted[,16], "MAP",
                    ifelse(pilot_mixed.fitted[,15] < pilot_mixed.fitted[,16], "MIG",
                           "NA"))
    pilot_mixed.fitted <- cbind(pilot_mixed.fitted, label)

    csvname <- paste0(filename,"_newfitted")
    write.csv(prior.history, paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename,"/",fitname,"/",csvname,"_prior_history.csv"))
    write.csv(ent.history, paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename,"/",fitname,"/",csvname,"_ent_history.csv"))
    
    csvname <- paste0(filename,"_newfitted")
    write.csv(pilot_mixed.fitted, paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename,"/",fitname,"/",csvname,".csv"))
    pilot_mixed.fitted <-data.frame(pilot_mixed.fitted)
    write_xlsx(pilot_mixed.fitted, paste0("/Users/jisu/Documents/Hidden Target/Model Fit/",filename,"/",fitname,"/",csvname,".xlsx"))
}}




post <- read.csv("/Users/jisu/Documents/Hidden Target/Model Fit/HT06/HT06_newfitted_prior_history.csv")
ent <- read.csv("/Users/jisu/Documents/Hidden Target/Model Fit/HT06/HT06_newfitted_ent_history.csv")


i=1
hist(melt(post[32*(i-1):32*i,])[,2])

dim(melt(post[1:32,]))
