
prior <- matrix(0, nrow=1024, ncol=768)

for (a in 1:1024){
  for (b in 1:768){
    
    for (i in 1:1024){
      for (j in 1:768){
        
        mat <- rbind(c(i,j), c(a,b))      
        if (dist(mat) < (768/2)){
          prior[a,b] <- prior[a,b] + 1
        }
        else {prior[a,b] <- prior[a,b]}
      }
    }
  }
}

gg_heatmap <- function(matrix, i, col){
  dat <- matrix[(nx*(i-1)+1):(nx*i),]
  names(dat) <- paste("X", 1:32) 
  dat2 <- melt(dat, id.var = "X1")
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
pdf("prior_190520", width=10, height=7)
gg_heatmap(prior,1,"red")
dev.off

prior.reduced <- matrix(NA, nrow=64, ncol=32)
for (a in 1:64){
  for (b in 1:32){
    prior.reduced[a,b] <- prior[a*16, b*16]
  }
}

prior.reduced.2 <- matrix(NA, nrow=32, ncol=24)
for (a in 1:32){
  for (b in 1:24){
    prior.reduced.2[a,b] <- prior[a*32, b*24]
  }
}


write.csv(prior.reduced,"prior.reduced.csv")
write.csv(prior.reduced.2,"prior.reduced.2.csv")


a <- read.csv("/Volumes/clmnlab/HT/behavior_data/raw_data/priorreduced.csv")
a <- a[,-c(1)]
dim(a)
which(a == max(a), arr.ind = TRUE)
which(a == min(a), arr.ind = TRUE)





prior <- matrix(NA, ncol=48, nrow=64)
for (a in 1:64){
  for (b in 1:48){
    
    for (i in 1:64){
      for (j in 1:48){
        
        mat <- rbind(c(i,j), c(a,b))      
        if (dist(mat) < (48/2)){
          prior[a,b] <- prior[a,b] + 1
        }
        else {prior[a,b] <- prior[a,b]}
      }
    }
  }
}


prior <- matrix(0, ncol=24, nrow=32)
for (a in 1:32){
  for (b in 1:24){
    
    for (i in 1:32){
      for (j in 1:24){
        
        mat <- rbind(c(i,j), c(a,b))      
        if (sqrt((i-a)^2 + (j-b)^2) < (24/2)){
          prior[a,b] <- prior[a,b] + 1
        } else {prior[a,b] <- prior[a,b]}
      }
    }
  }
}
prior <- prior/sum(prior)
sum(prior)
ent <- sum(prior * log(1/prior))

write.csv(prior, "/Volumes/clmnlab/HT/behavior_data/raw_data/priorreduced2.csv")
