

S=10000
P=100


n1=400 # sample size

# biomarker randomness
bio.sigma = 0.0


out.dir <- c("C:\\Users\\zhaob\\dataset\\")




#########
## setup
#########
library(MASS)
library(Matrix)




#########
## generate
#########

for(s in 1:S)
{
  
  cat("s:",s,"\n")
  
  datnum <- c(s)

  # initialize 4 biomarkers
  biom1.x1.t=runif(n1/2,0,1) + rnorm(n1/2,0,bio.sigma^2) 
  biom1.x1.c=runif(n1/2,0,1) + rnorm(n1/2,0,bio.sigma^2) 
  biom1.x2.t=runif(n1/2,0,1) + rnorm(n1/2,0,bio.sigma^2) 
  biom1.x2.c=runif(n1/2,0,1) + rnorm(n1/2,0,bio.sigma^2) 
  biom1.x3.t=runif(n1/2,0,1) + rnorm(n1/2,0,bio.sigma^2) 
  biom1.x3.c=runif(n1/2,0,1) + rnorm(n1/2,0,bio.sigma^2) 
  biom1.x4.t=runif(n1/2,0,1) + rnorm(n1/2,0,bio.sigma^2) 
  biom1.x4.c=runif(n1/2,0,1) + rnorm(n1/2,0,bio.sigma^2) 
  
  # generate random intercepts for multiple outcomes
  rn.intercept <- matrix(0, nrow = n1/2, ncol = 6)                                  # all 0,  no correlation
  

  cv1.t.index <- sample(n1/2, n1/4) 
  cv2.t.index <- (1:(n1/2))[-cv1.t.index]
  
  cv1.c.index <- sample(n1/2, n1/4) 
  cv2.c.index <- (1:(n1/2))[-cv1.c.index]   
  
  
  out <- list(biom1.x1.t,biom1.x1.c,biom1.x2.t,biom1.x2.c,biom1.x3.t,biom1.x3.c,biom1.x4.t, biom1.x4.c,
              rn.intercept, 
              cv1.t.index, cv2.t.index, cv1.c.index, cv2.c.index)
  
  saveRDS(out,  file = paste(out.dir, "data_", datnum, ".Rds", sep=""), version = 2)

}







