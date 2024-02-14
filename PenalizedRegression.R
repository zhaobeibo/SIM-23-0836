


# type.full <- c(0,1)
# corr.type.full <- c(0,1)
# sce <- c(1,2)  # 4 type and corr.type scenarios 
# 
# outcome.full <- list(c(1,2,3))
# m.out.full <- list(c("1N1B1S"))
# 
# model.full <- list(c("Model1a"),c("Model2b"),c("Model3c"),c("Model4d"),c("Model5"))

# type.full <- c(0,1)
# corr.type.full <- c(0,1)
sce <- c(1,2,3,4)  # 4 type and corr.type scenarios 

# outcome.full <- list(c(1,1,1),c(1,2,3))
# m.out.full <- list(c("3N"),c("1N1B1S"))

outcome.full <- list(c(1,2,3))
m.out.full <- list(c("1N1B1S"))

model.full <- list(c("Model0"),c("Model1a"),c("Model2b"),c("Model3c"),c("Model4d"),c("Model5"))





n1=400 # sample size
S=5 # simulation runs
P<-2# permutation runs

t.rateC <- 0.005
c.rateC <- 0.003



in.dir <- c("/work/users/b/e/beibo/Research/paper1/")
key.dir <- c("/work/users/b/e/beibo/Research/paper1/Jan2024/LP_CV_power_new/")









##############################################################
###### Linear Model Methods with LASSO - 4 biomarkers
##############################################################

library(MASS)
library(Matrix)
library(glmnet)
library(grpreg)
library(grpregOverlap)
library(vennLasso)

library(survival)
library(prodlim)

library(Hmisc)
library(Formula)
library(ggplot2)
library(rms)
library(dplyr)
library(erer)

library(copula)
# set.seed(1)


####################
### Functions
####################


### Binary outcome
simulbin <- function(N, cov, trt, rn.intercept) {
  z = rn.intercept + p.beta0 + p.theta*trt + p.delta*trt*cov       # linear combination with a bias
  pr = 1/(1+exp(-z))         # pass through an inv-logit function
  y = rbinom(N,1,pr)      # bernoulli response variable
}


# Survival function
# baseline hazard: Weibull
# N = sample size
# lambda = scale parameter in h0()
# rho = shape parameter in h0()
# beta = fixed effect parameter
# rateC = rate parameter of the exponential distribution of C
simulWeib <- function(N, lambda, cov, trt, rho, rateC, rn.intercept)
{
  # Weibull latent event times
  v <- runif(n=N)
  Tlat <- (- log(v) / (lambda * exp(cov * trt * s.delta + trt*s.theta + s.beta0 + rn.intercept)))^(1 / rho)
  
  # censoring times, randomly drawn from exponential dist
  C <- rexp(n=N, rate=rateC)
  
  # follow-up times and event indicators
  time <- pmin(Tlat, C)
  status <- as.numeric(Tlat <= C)
  
  # data set
  data.frame(id=1:N,
             time=time,
             status=status,
             x=cov,
             trt=trt)
}

# function to calculate z score from proportion test
z.prop = function(x1,x2,n1,n2){
  numerator = (x1/n1) - (x2/n2)
  p.common = (x1+x2) / (n1+n2)
  denominator = sqrt(p.common * (1-p.common) * (1/n1 + 1/n2))
  z.prop.ris = numerator / denominator
  return(z.prop.ris)
}

# metrics of biomarker selection
# indicate proportion of tests that are subsets of truth
# for example, if truth is x1, x2, what is considered good is:
# only x1
# only x2
# x1 x2
biom1.pred <- function(l1.coeff.1,l1.coeff.2,l1.coeff.3,l1.coeff.4,c.true) {
  l.coeff <- cbind(l1.coeff.1,l1.coeff.2,l1.coeff.3,l1.coeff.4)
  
  # # old standard: miss biomarker by 1
  # l <- apply(l.coeff, 1, function(x) (sum(which(x != 0) %in% c.true)))
  # c.noise <- c.full[!c.full %in% c.true]
  # l.noise <- apply(l.coeff, 1, function(x) (sum(which(x != 0) %in% c.noise)))
  # pred.x <- sum((abs(l - length(c.true)) + l.noise) <= 1)/S
  
  # # exact selection of truth
  l <- apply(l.coeff, 1, function(x) (sum(setequal(which(x != 0),c.true))))
  pred.x.exact <- sum(l)
  
  # subset of truth
  # l <- apply(l.coeff, 1, function(x) (sum(all(which(x != 0) %in% c.true))))
  l <- apply(l.coeff, 1, function(x) (sum(all(ifelse(length(which(x != 0)) == 0, 99, which(x != 0)) %in% c.true))))
  pred.x.sub <- sum(l)
  out <- list(pred.x.exact,pred.x.sub)
  return(out)
}

# function to calculate test statistics, for all 
perc.U.all <- function(n.perc.U.all, ind.t, ind.c) {
  if ((length(ind.t)>=2) & (length(ind.c)>=2) & ((length(ind.t) + length(ind.c))/n.perc.U.all >= prev.cutoff)) {
    # outcome 1
    if (outcome[1] == 2) {
      p.1.a <- z.prop(sum(y1.t[ind.t]), sum(y1.c[ind.c]), length(ind.t), length(ind.c))
    } else if (outcome[1] == 3) {
      y1 <- rbind(y1.t, y1.c)
      ind <- c(ind.t, ind.c + n.perc.U.all)
      model <- coxph(formula = Surv(time, status) ~ trt, data = y1[ind, ])
      if ('try-error' %in% class(model)) next  
      p.1.a <- coef(model)/sqrt(diag(vcov(model)))
    } else {
      p.1.a <- t.test(y1.t[ind.t],y1.c[ind.c])$statistic   
    }  
    # outcome 2
    if (outcome[2] == 2) {
      p.2.a <- z.prop(sum(y2.t[ind.t]), sum(y2.c[ind.c]), length(ind.t), length(ind.c))
    } else if (outcome[2] == 3) {
      y2 <- rbind(y2.t, y2.c)
      ind <- c(ind.t, ind.c + n.perc.U.all)
      model <- coxph(formula = Surv(time, status) ~ trt, data = y2[ind, ])
      if ('try-error' %in% class(model)) next  
      p.2.a <- coef(model)/sqrt(diag(vcov(model)))
    } else {
      p.2.a <- t.test(y2.t[ind.t],y2.c[ind.c])$statistic     
    }    
    
    # outcome 3
    if (outcome[3] == 2) {
      p.3.a <- z.prop(sum(y3.t[ind.t]), sum(y3.c[ind.c]), length(ind.t), length(ind.c))
    } else if (outcome[3] == 3) {
      y3 <- rbind(y3.t, y3.c)
      ind <- c(ind.t, ind.c + n.perc.U.all)
      model <- coxph(formula = Surv(time, status) ~ trt, data = y3[ind, ])
      if ('try-error' %in% class(model)) next  
      p.3.a <- coef(model)/sqrt(diag(vcov(model)))
    } else {
      p.3.a <- t.test(y3.t[ind.t],y3.c[ind.c])$statistic    
    }
  } else {
    p.1.a <-NA
    p.2.a <-NA
    p.3.a <-NA
  }
  out<-list(p.1.a,p.2.a,p.3.a)
  return(out)
}

# function to calculate test statistics, for subgroup
perc.U <- function (n.perc.U, y1.t.perc.U, y1.c.perc.U, y2.t.perc.U, y2.c.perc.U, y3.t.perc.U, y3.c.perc.U, 
                    ab.test, cutoff, sign, type.perc.U,
                    cv.t.index,cv.c.index,
                    weight.vec) {
  # rounding
  ab.test = signif(ab.test, digits = 10)
  cutoff =  signif(cutoff, digits = 10)
  if (sign == "GE") {
    ind.sub<-which(ifelse(ab.test>= cutoff,1,0) == 1) 
  } else {
    ind.sub<-which(ifelse(ab.test< cutoff,1,0) == 1)
  }
  ind.t<-ind.sub[ind.sub<=(n.perc.U/2)]
  ind.c<-ind.sub[ind.sub>(n.perc.U/2)]-(n.perc.U/2)
  if ((length(ind.t)>=2) & (length(ind.c)>=2) 
      # & ((length(ind.t) + length(ind.c))/n.perc.U >= prev.cutoff)
  ) {
    if (outcome[1] == 1) {
      p.1 <- t.test(y1.t.perc.U[ind.t],y1.c.perc.U[ind.c])$p.v  
      z.1 <- t.test(y1.t.perc.U[ind.t],y1.c.perc.U[ind.c])$statistic
    } else if (outcome[1] == 2) {
      p.1 <- prop.test(x=c(sum(y1.t.perc.U[ind.t]), sum(y1.c.perc.U[ind.c])), n=c(length(y1.t.perc.U[ind.t]), length(y1.c.perc.U[ind.c])))$p.v 
      z.1 <- sqrt( prop.test(x=c(sum(y1.t.perc.U[ind.t]), sum(y1.c.perc.U[ind.c])), n=c(length(y1.t.perc.U[ind.t]), length(y1.c.perc.U[ind.c])))$statistic)
    } else {
      y1 <- rbind(y1.t.perc.U, y1.c.perc.U)
      p.1 <- pchisq(survdiff(formula = Surv(time, status) ~ trt, data = y1[ind.sub, ])$chisq, 1, lower.tail=FALSE)    
      z.1 <- sqrt(survdiff(formula = Surv(time, status) ~ trt, data = y1[ind.sub, ])$chisq)
    }
    if (outcome[2] == 1) {
      p.2 <- t.test(y2.t.perc.U[ind.t],y2.c.perc.U[ind.c])$p.v     
      z.2 <- t.test(y2.t.perc.U[ind.t],y2.c.perc.U[ind.c])$statistic
    } else if (outcome[2] == 2) {
      p.2 <- prop.test(x=c(sum(y2.t.perc.U[ind.t]), sum(y2.c.perc.U[ind.c])), n=c(length(y2.t.perc.U[ind.t]), length(y2.c.perc.U[ind.c])))$p.v 
      z.2 <- sqrt(prop.test(x=c(sum(y2.t.perc.U[ind.t]), sum(y2.c.perc.U[ind.c])), n=c(length(y2.t.perc.U[ind.t]), length(y2.c.perc.U[ind.c])))$statistic)
    } else {
      y2 <- rbind(y2.t.perc.U, y2.c.perc.U)
      p.2 <- pchisq(survdiff(formula = Surv(time, status) ~ trt, data = y2[ind.sub, ])$chisq, 1, lower.tail=FALSE)  
      z.2 <- sqrt(survdiff(formula = Surv(time, status) ~ trt, data = y2[ind.sub, ])$chisq)
    }
    if (outcome[3] == 1) {
      p.3 <- t.test(y3.t.perc.U[ind.t],y3.c.perc.U[ind.c])$p.v     
      z.3 <- t.test(y3.t.perc.U[ind.t],y3.c.perc.U[ind.c])$statistic
    } else if (outcome[3] == 2) {
      p.3 <- prop.test(x=c(sum(y3.t.perc.U[ind.t]), sum(y3.c.perc.U[ind.c])), n=c(length(y3.t.perc.U[ind.t]), length(y3.c.perc.U[ind.c])))$p.v   
      z.3 <- sqrt(prop.test(x=c(sum(y3.t.perc.U[ind.t]), sum(y3.c.perc.U[ind.c])), n=c(length(y3.t.perc.U[ind.t]), length(y3.c.perc.U[ind.c])))$statistic)
    } else {
      y3 <- rbind(y3.t.perc.U, y3.c.perc.U)
      p.3 <- pchisq(survdiff(formula = Surv(time, status) ~ trt, data = y3[ind.sub, ])$chisq, 1, lower.tail=FALSE)   
      z.3 <- sqrt(survdiff(formula = Surv(time, status) ~ trt, data = y3[ind.sub, ])$chisq)
    }
    
    p.vec <- c(p.1,p.2,p.3)
    # p.vec.weight <- p.vec / weight.vec
    
    p.vec.weight <- p.vec 
    
    # p.hoch <- min(p.adjust(p.vec, method = "hochberg"), na.rm = TRUE) # hochberg-adjusted p value
    p.hoch <- min(p.adjust(p.vec.weight, method = "hochberg"), na.rm = TRUE) # hochberg-adjusted p value
    
    z.weight <- sum(c(z.1,z.2,z.3)*weight.vec,na.rm = TRUE) 
    
    # p value is two sided, and we keep the effect direction 
    if (type.perc.U == 1) {
      if (outcome[1] == 1) {
        effect = sign(t.test(y1.t.perc.U[ind.t],y1.c.perc.U[ind.c])$statistic)    
      } else if (outcome[1] == 2) {
        effect = sign(z.prop(sum(y1.t.perc.U[ind.t]), sum(y1.c.perc.U[ind.c]), length(y1.t.perc.U[ind.t]), length(y1.c.perc.U[ind.c])))   
      } else {
        model <- coxph(formula = Surv(time, status) ~ trt, data = y1[ind.sub,])
        effect = sign(coef(model)/sqrt(diag(vcov(model))))
        if ('try-error' %in% class(model)) next  
      }
    } else if (type.perc.U == 2) {
      if (outcome[2] == 1) {
        effect = sign(t.test(y2.t.perc.U[ind.t],y2.c.perc.U[ind.c])$statistic)
      } else if (outcome[2] == 2) {
        effect = sign(z.prop(sum(y2.t.perc.U[ind.t]), sum(y2.c.perc.U[ind.c]), length(y2.t.perc.U[ind.t]), length(y2.c.perc.U[ind.c])))
      } else {
        model <- coxph(formula = Surv(time, status) ~ trt, data = y2[ind.sub,])
        effect = sign(coef(model)/sqrt(diag(vcov(model))))
        if ('try-error' %in% class(model)) next  
      }
    } else {
      if (outcome[3] == 1) {
        effect = sign(t.test(y3.t.perc.U[ind.t],y3.c.perc.U[ind.c])$statistic)
      } else if (outcome[3] == 2) {
        effect = sign(z.prop(sum(y3.t.perc.U[ind.t]), sum(y3.c.perc.U[ind.c]), length(y3.t.perc.U[ind.t]), length(y3.c.perc.U[ind.c])))
      } else {
        model <- coxph(formula = Surv(time, status) ~ trt, data = y3[ind.sub,])
        effect = sign(coef(model)/sqrt(diag(vcov(model))))
        if ('try-error' %in% class(model)) next  
      }     
    }
    if (is.na(effect)) {
      z.hoch <- NA
    } else {
      if (effect >= 0) {
        z.hoch <- qnorm(1-p.hoch/2)
      } else {
        z.hoch <- qnorm(p.hoch/2)      
      }
    }
  } else {
    p.1 <- p.2 <- p.3 <- p.hoch <- z.hoch <- z.weight <- NA
  }
  
  # lm.biom1.x1 <- biom1.x1[c(cv.t.index,cv.c.index+n.perc.U)]
  # lm.biom1.x2 <- biom1.x2[c(cv.t.index,cv.c.index+n.perc.U)]
  # lm.biom1.x3 <- biom1.x3[c(cv.t.index,cv.c.index+n.perc.U)]
  # lm.biom1.x4 <- biom1.x4[c(cv.t.index,cv.c.index+n.perc.U)]
  # 
  # plot(lm.biom1.x1[ind.sub],lm.biom1.x2[ind.sub],xlim=c(0,1),ylim=c(0,1)) # plot to see if makes sense
  # plot(lm.biom1.x3[ind.sub],lm.biom1.x4[ind.sub],xlim=c(0,1),ylim=c(0,1)) # plot to see if makes sense
  
  
  
  out<-list(z.hoch, z.weight, 
            p.1, p.2, p.3, p.hoch,
            ind.sub)
  return(out)
}


# function to calculate %U, for subgroup
perc.U.test <- function (input.coef, cutoff, sign, type.U.test, weight.vec) {
  if (dim(input.coef)[1] < dim(x.mat.final.t)[2]) {
    ab.test <- as.matrix(x.mat.final.t)[,-1] %*% input.coef  
  } 
  else {
    ab.test <- as.matrix(x.mat.final.t) %*% input.coef
  }
  # rounding
  ab.test = signif(ab.test, digits = 10)
  cutoff =  signif(cutoff, digits = 10)
  if (sign == "GE") {
    ind.sub<-which(ifelse(ab.test >= cutoff,1,0) == 1) 
  } else {
    ind.sub<-which(ifelse(ab.test < cutoff,1,0) == 1)
  }
  ind.t<-ind.sub[ind.sub<=n.obs]
  ind.c<-ind.sub[ind.sub>n.obs]-n.obs
  ind.true <- which(X1.final.t>c01 & X2.final.t>c02 & X3.final.t + X4.final.t > c03 + c04)
  if ((length(ind.t)>=2) & (length(ind.c)>=2) 
      # & ((length(ind.t) + length(ind.c))/n.obs >= prev.cutoff)
  ) {
    
    # if (type.U.test == 2) {
    #   pi.sub<-(length(ind.sub))/(2*n.obs)
    #   pi.overlap <- length(intersect(ind.t, ind.true))/length(ind.true)
    #   perc.U<-z.prop(sum(b.y.final.t[ind.t]), sum(b.y.final.c[ind.c]), length(ind.t), length(ind.c))/(sqrt(2*n.obs/2)*b.U.best[s])*100 
    # } else if (type.U.test == 3) {
    #   y.final <- rbind(s.y.final.t, s.y.final.c)
    #   pi.sub<-(length(ind.sub))/(2*n.obs)
    #   pi.overlap <- length(intersect(ind.t, ind.true))/length(ind.true)
    #   model <- coxph(formula = Surv(time, status) ~ trt, data = y.final[ind.sub, ])
    #   perc.U<- (coef(model)/sqrt(diag(vcov(model))))/(sqrt(2*n.obs/2)*s.U.best[s])*100 
    # } else {
    #   pi.sub<-(length(ind.sub))/(2*n.obs)
    #   pi.overlap <- length(intersect(ind.t, ind.true))/length(ind.true)
    #   perc.U<-t.test(n.y.final.t[ind.t],n.y.final.c[ind.c])$statistic/(sqrt(2*n.obs/2)*n.U.best[s])*100 
    # }
    
    pi.sub<-(length(ind.sub))/(2*n.obs)
    pi.overlap <- length(intersect(ind.t, ind.true))/length(ind.true)
    
    perc.U.1<-t.test(n.y.final.t[ind.t],n.y.final.c[ind.c])$statistic
    perc.U.2<-z.prop(sum(b.y.final.t[ind.t]), sum(b.y.final.c[ind.c]), length(ind.t), length(ind.c)) 
    
    y.final <- rbind(s.y.final.t, s.y.final.c)
    model <- coxph(formula = Surv(time, status) ~ trt, data = y.final[ind.sub, ])
    perc.U.3<- (coef(model)/sqrt(diag(vcov(model)))) 
    
    perc.U.weight <- sum(c(
      (sqrt(2*n.obs/2)*n.U.best[s]),  
      (sqrt(2*n.obs/2)*b.U.best[s]), 
      (sqrt(2*n.obs/2)*s.U.best[s]))*weight.vec, na.rm = TRUE) 
    
    perc.U <- sum(c(perc.U.1, perc.U.2, perc.U.3)*weight.vec,na.rm = TRUE) / perc.U.weight * 100
    
  } else {
    pi.sub<- NA 
    pi.overlap <- NA
    perc.U<- NA
  }
  # plot(X1.final[ind.sub],X2.final[ind.sub],xlim=c(0,1),ylim=c(0,1)) # plot to see if makes sense
  # plot(X3.final[ind.sub],X4.final[ind.sub],xlim=c(0,1),ylim=c(0,1)) # plot to see if makes sense      
  out<-list(pi.sub,
            perc.U,
            ind.sub,
            pi.overlap,
            
            perc.U.1,
            perc.U.2,
            perc.U.3,
            sum(c(perc.U.1, perc.U.2, perc.U.3)*weight.vec,na.rm = TRUE),
            perc.U.weight
            )
  return(out)
}





# cross-validation for testing
lm.cv <- function(n1.lm,
                  cv.t.index,cv.c.index,
                  lm.y1.t,lm.y1.c,
                  lm.y2.t,lm.y2.c,
                  lm.y3.t,lm.y3.c) {
  
  lm.x.mat <- x.mat[c(cv.t.index,cv.c.index+n1.lm),]
  lm.x.mat.fit.t <- x.mat.fit.t[c(cv.t.index,cv.c.index+n1.lm),]
  
  # initialize
  ab.1 <- pred.1 <-  ab.2 <- pred.2 <- ab.3 <- pred.3 <-  rep(NA,n1.lm)
  uv.1.a <- uv.2.a <- uv.3.a <- rep(NA,n1.lm)  
  uv.1.p <- uv.2.p <- uv.3.p <- rep(NA,n1.lm) 
  uv.1.pl <- uv.2.pl <- uv.3.pl <- rep(NA,n1.lm)  
  p.1.a <- p.2.a <- p.3.a <- rep(NA,n1.lm)  
  p.1.p <- p.2.p <- p.3.p <- rep(NA,n1.lm) 
  p.1.pl <- p.2.pl <- p.3.pl <- rep(NA,n1.lm) 
  p.raw <- p.raw.index <- rep(NA,n1.lm)
  p.weight <- p.weight.index <- rep(NA,n1.lm)
  
  # LMs - normal outcomes
  
  if (outcome[1] == 1) {
    y1.fit <- c(lm.y1.t, lm.y1.c)
    # fit LMs with GLASSO - cross-validated penalty term
    lambda.1 = cv.grpregOverlap(X = as.matrix(lm.x.mat), y = y1.fit, group = grp.mat,  nfolds = fold.num, penalty = penalty.type)$lambda.min
    # lambda.1 = cv.glmnet(x= as.matrix(cbind(1,lm.x.mat)), y = y1.fit, nfolds = fold.num)$lambda.1se
    fit.1=grpregOverlap(X = as.matrix(lm.x.mat), y = y1.fit, group = grp.mat, penalty = penalty.type,lambda = lambda.1)
    # get linear predictors from LMs w GLASSO
    ab.1 <- as.matrix(lm.x.mat.fit.t) %*% as.matrix(predict(fit.1, type = "coefficients"))
    # type of outcome
    y1.type = fit.1$family
  }
  
  if (outcome[2] == 1) {
    y2.fit <- c(lm.y2.t, lm.y2.c)
    # fit LMs with GLASSO - cross-validated penalty term
    lambda.2 = cv.grpregOverlap(X = as.matrix(lm.x.mat), y = y2.fit, group = grp.mat,  nfolds = fold.num, penalty = penalty.type)$lambda.min
    # lambda.2 = cv.glmnet(x= as.matrix(cbind(1,lm.x.mat)), y = y2.fit, nfolds = fold.num)$lambda.1se
    fit.2=grpregOverlap(X = as.matrix(lm.x.mat), y = y2.fit, group = grp.mat, penalty = penalty.type,lambda = lambda.2)
    # get linear predictors from LMs w GLASSO
    ab.2 <- as.matrix(lm.x.mat.fit.t) %*% as.matrix(predict(fit.2, type = "coefficients"))
    # type of outcome
    y2.type = fit.2$family
  }
  
  if (outcome[3] == 1) {
    y3.fit <- c(lm.y3.t, lm.y3.c)
    # fit LMs with GLASSO - cross-validated penalty term
    lambda.3 = cv.grpregOverlap(X = as.matrix(lm.x.mat), y = y3.fit, group = grp.mat,  nfolds = fold.num, penalty = penalty.type)$lambda.min
    # lambda.3 = cv.glmnet(x= as.matrix(cbind(1,lm.x.mat)), y = y3.fit, nfolds = fold.num)$lambda.1se
    fit.3=grpregOverlap(X = as.matrix(lm.x.mat), y = y3.fit, group = grp.mat, penalty = penalty.type,lambda = lambda.3)
    # get linear predictors from LMs w GLASSO
    ab.3 <- as.matrix(lm.x.mat.fit.t) %*% as.matrix(predict(fit.3, type = "coefficients"))
    # type of outcome
    y3.type = fit.3$family
  }
  
  
  
  
  # GLMs - binary outcomes
  
  if (outcome[1] == 2) {
    y1.fit <- c(lm.y1.t, lm.y1.c)
    # fit GLMs with GLASSO - cross-validated penalty term
    lambda.1 = cv.grpregOverlap(X = as.matrix(lm.x.mat), y = y1.fit, group = grp.mat,  family = "binomial", nfolds = fold.num, penalty = penalty.type)$lambda.min
    # lambda.1 = cv.glmnet(x= as.matrix(cbind(1,lm.x.mat)), y = y1.fit,  family = "binomial", nfolds = fold.num)$lambda.1se
    fit.1=grpregOverlap(X = as.matrix(lm.x.mat), y = y1.fit, group = grp.mat, penalty = penalty.type, family = "binomial", lambda = lambda.1)
    # get linear predictors from GLMs w GLASSO
    ab.1 <- as.matrix(lm.x.mat.fit.t) %*% as.matrix(predict(fit.1, type = "coefficients"))
    y1.type = fit.1$family
  }
  
  if (outcome[2] == 2) {
    y2.fit <- c(lm.y2.t, lm.y2.c)
    # fit GLMs with GLASSO - cross-validated penalty term
    lambda.2 = cv.grpregOverlap(X = as.matrix(lm.x.mat), y = y2.fit, group = grp.mat,  family = "binomial", nfolds = fold.num, penalty = penalty.type)$lambda.min
    # lambda.2 = cv.glmnet(x= as.matrix(cbind(1,lm.x.mat)), y = y2.fit,  family = "binomial", nfolds = fold.num)$lambda.1se
    fit.2=grpregOverlap(X = as.matrix(lm.x.mat), y = y2.fit, group = grp.mat, penalty = penalty.type,family = "binomial",lambda = lambda.2)
    # get predicted trt effects from GLMs w GLASSO
    ab.2 <- as.matrix(lm.x.mat.fit.t) %*% as.matrix(predict(fit.2, type = "coefficients"))
    y2.type = fit.2$family    
  }
  
  if (outcome[3] == 2) {
    y3.fit <- c(lm.y3.t, lm.y3.c)
    # fit GLMs with GLASSO - cross-validated penalty term
    lambda.3 = cv.grpregOverlap(X = as.matrix(lm.x.mat), y = y3.fit, group = grp.mat,  family = "binomial", nfolds = fold.num, penalty = penalty.type)$lambda.min
    # lambda.3 = cv.glmnet(x= as.matrix(cbind(1,lm.x.mat)), y = y3.fit,  family = "binomial", nfolds = fold.num)$lambda.1se
    fit.3=grpregOverlap(X = as.matrix(lm.x.mat), y = y3.fit, group = grp.mat, penalty = penalty.type,family = "binomial",lambda = lambda.3)
    # get predicted trt effects from GLMs w GLASSO
    ab.3 <- as.matrix(lm.x.mat.fit.t) %*% as.matrix(predict(fit.3, type = "coefficients"))
    y3.type = fit.3$family
  }
  
  
  
  # GLMs - survival outcomes
  # adjust survival outcomes b.c. of difference in trt direction
  
  if (outcome[1] == 3) {
    y1 <- rbind(lm.y1.t, lm.y1.c)
    # fit GLMs with GLASSO - cross-validated penalty term
    y1.fit <- with(y1, Surv(time,status)) # create survival object
    lambda.1 = cv.grpsurvOverlap(X = as.matrix(lm.x.mat), y = y1.fit, group = grp.mat, penalty = penalty.type, cv.ind = TRUE,  nfolds = fold.num)$lambda.min
    # lambda.1 = cv.glmnet(x= as.matrix(cbind(1,lm.x.mat)), y = y1.fit,  family = "cox", nfolds = fold.num)$lambda.1se
    fit.1=grpregOverlap(X = as.matrix(lm.x.mat), y = y1.fit, group = grp.mat, penalty = penalty.type, family = "cox", lambda = lambda.1)
    # get linear predictors from GLMs w GLASSO
    ab.1 <- as.matrix(lm.x.mat.fit.t[,-1]) %*% as.matrix(predict(fit.1, type = "coefficients"))
    y1.type = fit.1$family
  }
  
  if (outcome[2] == 3) {
    y2 <- rbind(lm.y2.t, lm.y2.c)
    # fit GLMs with GLASSO - cross-validated penalty term
    y2.fit <- with(y2, Surv(time,status)) # create survival object
    lambda.2 = cv.grpsurvOverlap(X = as.matrix(lm.x.mat), y = y2.fit, group = grp.mat, penalty = penalty.type, cv.ind = TRUE,  nfolds = fold.num)$lambda.min
    # lambda.2 = cv.glmnet(x= as.matrix(cbind(1,lm.x.mat)), y = y2.fit,  family = "cox", nfolds = fold.num)$lambda.1se
    fit.2=grpregOverlap(X = as.matrix(lm.x.mat), y = y2.fit, group = grp.mat, penalty = penalty.type, family = "cox", lambda = lambda.2)
    # get linear predictors from GLMs w GLASSO
    ab.2 <- as.matrix(lm.x.mat.fit.t[,-1]) %*% as.matrix(predict(fit.2, type = "coefficients"))
    y2.type = fit.2$family
  }
  
  if (outcome[3] == 3) {
    y3 <- rbind(lm.y3.t, lm.y3.c)
    # fit GLMs with GLASSO - cross-validated penalty term
    y3.fit <- with(y3, Surv(time,status)) # create survival object
    lambda.3 = cv.grpsurvOverlap(X = as.matrix(lm.x.mat), y = y3.fit, group = grp.mat, penalty = penalty.type, cv.ind = TRUE,  nfolds = fold.num)$lambda.min
    # lambda.3 = cv.glmnet(x= as.matrix(cbind(1,lm.x.mat)), y = y3.fit,  family = "cox", nfolds = fold.num)$lambda.1se
    fit.3=grpregOverlap(X = as.matrix(lm.x.mat), y = y3.fit, group = grp.mat, penalty = penalty.type, family = "cox", lambda = lambda.3)
    # get linear predictors from GLMs w GLASSO
    ab.3 <- as.matrix(lm.x.mat.fit.t[,-1]) %*% as.matrix(predict(fit.3, type = "coefficients"))
    y3.type = fit.3$family    
  }
  
  # rounding 
  ab.1 <- signif(ab.1, digits = 10)
  ab.2 <- signif(ab.2, digits = 10)
  ab.3 <- signif(ab.3, digits = 10)
  
  
  
  # Create a matrix to store the vectors
  vector.matrix <- matrix(ncol = 3, nrow = 0)
  
  # Loop over integers from 2 to 50
  for (i in 2:5) {
    n <- 1/i
    if (2*n <= 1) { # this condition ensures that 1 - 2n is nonnegative
      # First configuration
      vec1 <- c(sqrt(n), sqrt(n), sqrt(1 - 2*n))
      vector.matrix <- rbind(vector.matrix, vec1)
      
      # Second configuration
      vec2 <- c(sqrt(1 - 2*n), sqrt(n), sqrt(n))
      vector.matrix <- rbind(vector.matrix, vec2)
      
      # Third configuration
      vec3 <- c(sqrt(n), sqrt(1 - 2*n), sqrt(n))
      vector.matrix <- rbind(vector.matrix, vec3)
    }
  }
  
  
  # Add the vectors (1,0,0), (0,1,0), (0,0,1)
  vector.matrix <- rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), vector.matrix)
  
  
  # Create a list to store the strings
  vector.string.list <- list()
  
  # Generate the strings
  for (i in 1:dim(vector.matrix)[1]) {
    vec <- vector.matrix[i, ]
    string <- paste("c(sqrt(", 1/(vec[1]^2), "), sqrt(", 1/(vec[2]^2), "), sqrt(", 1/(vec[3]^2), "))", sep = "")
    vector.string.list <- c(vector.string.list, string)
  }
  
  # # Print the first 50 strings
  # for (string in vector.string.list) {
  #   print(string)
  # }
  
  p.weight.index.list <- p.weight.max.list <- w.3lm.star.ind.p.list <- rep(NA,dim(vector.matrix)[1])
  
  w.length <- dim(vector.matrix)[1]
  
  for (w in 1:w.length) {
    
    cat("w:",w,"\n")
    
    weight.vec <- vector.matrix[w,]
    
    # weighted z score with actual response to determine subgroup
    # >= threshold
    for (i in 1:n1.lm)   {
      sub.index.1 <- which(ab.1 >= ab.1[i])
      sub.index.2 <- which(ab.2 >= ab.2[i])
      sub.index.3 <- which(ab.3 >= ab.3[i])
      sub.index <- list(sub.index.1, sub.index.2, sub.index.3)
      p.weight.temp <- rep(NA,3)
      
      for (m in 1:3) {
        # coef.t <- unlist(coef[m])
        sub.index.t <- unlist(sub.index[m])
        g.t.t <- sub.index.t[lm.x.mat$t[sub.index.t] == 1]
        g.c.t <- sub.index.t[lm.x.mat$t[sub.index.t] == 0]
        if((length(g.t.t)>=2) & (length(g.c.t)>=2) & ((length(g.t.t) + length(g.c.t))/n1.lm >= prev.cutoff)) {
          if (y1.type=="binomial") {
            p.1.p[i] <- z.prop(sum(y1.fit[g.t.t]), sum(y1.fit[g.c.t]), length(y1.fit[g.t.t]), length(y1.fit[g.c.t]))
          } else if (y1.type=="cox") {
            model <- coxph(formula = Surv(time, status) ~ trt, data = y1[sub.index.t, ])
            if ('try-error' %in% class(model)) next  
            p.1.p[i] <- coef(model)/sqrt(diag(vcov(model)))
          } else {
            p.1.p[i] <- t.test(y1.fit[g.t.t],y1.fit[g.c.t])$statistic       
          }
          if (y2.type=="binomial") {
            p.2.p[i] <-  z.prop(sum(y2.fit[g.t.t]), sum(y2.fit[g.c.t]), length(y2.fit[g.t.t]), length(y2.fit[g.c.t]))   
          } else if (y2.type=="cox") {
            model <- coxph(formula = Surv(time, status) ~ trt, data = y2[sub.index.t, ])
            if ('try-error' %in% class(model)) next  
            p.2.p[i] <- coef(model)/sqrt(diag(vcov(model)))
          } else {
            p.2.p[i] <- t.test(y2.fit[g.t.t],y2.fit[g.c.t])$statistic        
          }
          if (y3.type=="binomial") {
            p.3.p[i] <- z.prop(sum(y3.fit[g.t.t]), sum(y3.fit[g.c.t]), length(y3.fit[g.t.t]), length(y3.fit[g.c.t]))            
          } else if (y3.type=="cox") {
            model <- coxph(formula = Surv(time, status) ~ trt, data = y3[sub.index.t, ])
            if ('try-error' %in% class(model)) next  
            p.3.p[i] <- coef(model)/sqrt(diag(vcov(model)))
          } else {
            p.3.p[i] <- t.test(y3.fit[g.t.t],y3.fit[g.c.t])$statistic            
          }
        } else {
          p.1.p[i] <-NA
          p.2.p[i] <-NA
          p.3.p[i] <-NA
        }    
        # if ((sum(coef.t %notin% coef.1) != 0)) {p.1.p[i] <-NA}
        # if ((sum(coef.t %notin% coef.2) != 0)) {p.2.p[i] <-NA}
        # if ((sum(coef.t %notin% coef.3) != 0)) {p.3.p[i] <-NA}
        p.weight.temp[m] = sum(c((p.1.p[i]), (p.2.p[i]),(p.3.p[i]))*weight.vec, na.rm = TRUE)
      }
      
      if (sum(is.na(p.weight.temp)) == 3) {
        p.weight.index[i] <- NA
        p.weight[i] <- NA
      } else {
        p.weight.index[i] = which.max(abs(p.weight.temp))       # absolute value 
        p.weight[i] = p.weight.temp[p.weight.index[i]]          # select by absolute value, but keep the original
      }
    
      
      #gc()
      }
    p.weight.index.grt <- p.weight.index
    p.weight.grt <- p.weight
    
    # < threshold
    for (i in 1:n1.lm)   {
      sub.index.1 <- which(ab.1 < ab.1[i])
      sub.index.2 <- which(ab.2 < ab.2[i])
      sub.index.3 <- which(ab.3 < ab.3[i])
      sub.index <- list(sub.index.1, sub.index.2, sub.index.3)
      p.weight.temp <- rep(NA,3)
      
      for (m in 1:3) {
        # coef.t <- unlist(coef[m])
        sub.index.t <- unlist(sub.index[m])
        g.t.t <- sub.index.t[lm.x.mat$t[sub.index.t] == 1]
        g.c.t <- sub.index.t[lm.x.mat$t[sub.index.t] == 0]
        if((length(g.t.t)>=2) & (length(g.c.t)>=2) & ((length(g.t.t) + length(g.c.t))/n1.lm >= prev.cutoff)) {
          if (y1.type=="binomial") {
            p.1.p[i] <- z.prop(sum(y1.fit[g.t.t]), sum(y1.fit[g.c.t]), length(y1.fit[g.t.t]), length(y1.fit[g.c.t]))
          } else if (y1.type=="cox") {
            model <- coxph(formula = Surv(time, status) ~ trt, data = y1[sub.index.t, ])
            if ('try-error' %in% class(model)) next  
            p.1.p[i] <- coef(model)/sqrt(diag(vcov(model)))
          } else {
            p.1.p[i] <- t.test(y1.fit[g.t.t],y1.fit[g.c.t])$statistic       
          }
          if (y2.type=="binomial") {
            p.2.p[i] <-  z.prop(sum(y2.fit[g.t.t]), sum(y2.fit[g.c.t]), length(y2.fit[g.t.t]), length(y2.fit[g.c.t]))     
          } else if (y2.type=="cox") {
            model <- coxph(formula = Surv(time, status) ~ trt, data = y2[sub.index.t, ])
            if ('try-error' %in% class(model)) next  
            p.2.p[i] <- coef(model)/sqrt(diag(vcov(model)))
          } else {
            p.2.p[i] <- t.test(y2.fit[g.t.t],y2.fit[g.c.t])$statistic        
          }
          if (y3.type=="binomial") {
            p.3.p[i] <- z.prop(sum(y3.fit[g.t.t]), sum(y3.fit[g.c.t]), length(y3.fit[g.t.t]), length(y3.fit[g.c.t]))              
          } else if (y3.type=="cox") {
            model <- coxph(formula = Surv(time, status) ~ trt, data = y3[sub.index.t, ])
            if ('try-error' %in% class(model)) next  
            p.3.p[i] <- coef(model)/sqrt(diag(vcov(model)))
          } else {
            p.3.p[i] <- t.test(y3.fit[g.t.t],y3.fit[g.c.t])$statistic            
          }
        } else {
          p.1.p[i] <-NA
          p.2.p[i] <-NA
          p.3.p[i] <-NA
        }    
        # if ((sum(coef.t %notin% coef.1) != 0)) {p.1.p[i] <-NA}
        # if ((sum(coef.t %notin% coef.2) != 0)) {p.2.p[i] <-NA}
        # if ((sum(coef.t %notin% coef.3) != 0)) {p.3.p[i] <-NA}
        p.weight.temp[m] = sum(c((p.1.p[i]), (p.2.p[i]),(p.3.p[i]))*weight.vec, na.rm = TRUE)
      }
      
      if (sum(is.na(p.weight.temp)) == 3) {
        p.weight.index[i] <- NA
        p.weight[i] <- NA
      } else {
        p.weight.index[i] = which.max(abs(p.weight.temp))       # absolute value  
        p.weight[i] = p.weight.temp[p.weight.index[i]]          # select by absolute value, but keep the original
      }
    
      #gc()
      }
    
    p.weight.index.les <- p.weight.index
    p.weight.les <- p.weight
    

    p.weight.full <- c(p.weight.grt, p.weight.les)
    p.weight.max <- max(abs(c(p.weight.grt, p.weight.les)), na.rm = TRUE)
    
    p.weight.index <- c(p.weight.index.grt, p.weight.index.les)
    w.3lm.star.ind.p = which.max(abs(p.weight.full)) # maximize weighted test statistics (absolute value)
    
    w.3lm.star.ind.p.list[w] <-   w.3lm.star.ind.p 
    p.weight.index.list[w] <- p.weight.index[w.3lm.star.ind.p] 
    p.weight.max.list[w] <- p.weight.max
    
    #gc()
    
  }
  
  p.index.best <-  p.weight.index.list[which.max(p.weight.max.list)]
  weight.vec.best <- vector.matrix[which.max(p.weight.max.list),]
  w.3lm.star.ind.p <- w.3lm.star.ind.p.list[which.max(p.weight.max.list)]
  
  if (p.index.best == 1) {
    ab.fit.p <- ab.1
    fit.3lm.p <- fit.1
    out.p <- y1.type
    y.p <- 1
  } else if (p.index.best == 2) {
    ab.fit.p <- ab.2
    fit.3lm.p <- fit.2
    out.p <- y2.type
    y.p <- 2
  } else if (p.index.best == 3) {
    ab.fit.p <- ab.3
    fit.3lm.p <- fit.3
    out.p <- y3.type
    y.p <- 3
  }    
  
  if (w.3lm.star.ind.p <= n1.lm) {
    ind.3lm.p = which(ab.fit.p>=ab.fit.p[w.3lm.star.ind.p])     
    sign.p <- "GE"
  } else if (n1.lm < w.3lm.star.ind.p  & w.3lm.star.ind.p  <= (2*n1.lm)) {
    w.3lm.star.ind.p  <- w.3lm.star.ind.p  - n1.lm
    ind.3lm.p = which(ab.fit.p<ab.fit.p[w.3lm.star.ind.p])   
    sign.p <- "LS"
  }
  
  
  # output
  out <- list(lm.x.mat.fit.t,
              # as.matrix(predict(fit.3lm.a, type = "coefficients")),ab.fit.a[w.3lm.star.ind.a],sign.a,out.a, y.a,
              NA,NA,NA,NA,NA,
              as.matrix(predict(fit.3lm.p, type = "coefficients")),ab.fit.p[w.3lm.star.ind.p],sign.p,out.p, y.p,
              # w.3lm.star.ind.a,
              NA,
              w.3lm.star.ind.p,
              weight.vec.best 
              ) 
  return(out)
  
}


# raw metrics, before flip, 0 to indiciate no effect
raw.metrics <- function(cv.coef) {
  biom1.x.test.full <- c("biom1.x1","biom1.x2","biom1.x3","biom1.x4")
  raw.per <-  rep(NA,length(biom1.x.test.full))
  for (i in 1:length(biom1.x.test.full)) {
    biom1.x.test <- biom1.x.test.full[i]
    x.row <- grepl(biom1.x.test, row.names(cv.coef))
    raw.per[i] <- sum(cv.coef[x.row,1] != 0)
  }
  return(raw.per)
}

# flip metrics to evaluate biomarker selection, do not re-evaluate cutoff to obtain flip % 
flip.metrics <- function(n.flip, cv.mat,
                         cv.coef,ind.org, sign, cutoff) {
  
  biom1.x.test.full <- c("biom1.x1","biom1.x2","biom1.x3","biom1.x4")
  flip.per <- rep(NA,length(biom1.x.test.full))
  
  for (i in 1:length(biom1.x.test.full)) {
    biom1.x.test <- biom1.x.test.full[i]
    x.col <- grepl(biom1.x.test, colnames(cv.mat))
    cv.mat.x <- cv.mat
    cv.mat.x[,x.col] <- 0
    ab.test.x <- signif(cv.mat.x %*% cv.coef, digits = 10)
    # do not re-evaluate cutoff
    if (sign == "GE") {
      ind.sub<-which(ifelse(ab.test.x>= cutoff,1,0) == 1)
    } else {
      ind.sub<-which(ifelse(ab.test.x< cutoff,1,0) == 1)
    }
    # flip percent
    org.metric <- rep(0,n.flip)
    sub.metric <- rep(0,n.flip)
    org.metric[ind.org] <- 1 
    sub.metric[ind.sub] <- 1 
    flip.per[i] <- sum(abs(org.metric - sub.metric))/length(ind.org)      
  }
  return(flip.per)
}


# refined subgroups after flip % refine
flip.cv.result <- function(coef.flip, n.flip, cv.mat, 
                           cv.coef, cv.mat.org, sign, ind, threshold) {
  biom1.x.test <- paste(c("biom1.x1","biom1.x2","biom1.x3","biom1.x4")[(coef.flip < threshold)],
                        collapse = "|")
  biom1.x.test <- ifelse(biom1.x.test == "", c("no.replacement"), biom1.x.test) # if empty string, then biomarker change, no.replacement as a placeholder
  x.col <- grep(pattern = biom1.x.test, x=colnames(cv.mat))
  cv.mat.x <- cv.mat
  cv.mat.x[,x.col] <- 0
  ab.test.x <- cv.mat.x %*% cv.coef
  coef.new <- cv.coef
  coef.new[x.col] <- 0
  # to determine the new cutoff, need to use the biomarkers from the original cv
  if (dim(cv.coef)[1] > dim(cv.mat.org)[2]) {
    cutoff.new <- as.matrix(cv.mat.org %*% coef.new[-1])[ind,]  
  } else if (dim(cv.coef)[1] < dim(cv.mat.org)[2]) {
    cutoff.new <- as.matrix(cv.mat.org[,-1] %*% coef.new)[ind,] 
  }
  else {
    cutoff.new <- as.matrix(cv.mat.org %*% coef.new)[ind,]   
  }
  out <- list(coef.new, cutoff.new)
  return(out)
}









####################
### Simulations
####################


# type of penalty
# penalty.type <- c("grLasso")
penalty.type <- c("gel")

# biomarker randomness
bio.sigma = 0.0
mu0 = 0.0

# prevalence cutoff
prev.cutoff <- 0.1
# flip % threshold cutoff
threshold.cutoff <- 0.2

# fold number for cross-validation to determine tuning parameter
fold.num = 10




#########
## setup
#########


for (c in 1:length(sce)) {
  
  sce.in <- sce[c]  
  
  if (sce.in==1) {
    type <- 0  # 0 or 1, if type 1 then one outcome
    corr.type <- 0 # 0 or 1, if 1 then high cor
  } else if (sce.in==2) {
    type <- 0  # 0 or 1, if type 1 then one outcome
    corr.type <- 1 # 0 or 1, if 1 then high cor
  } else if (sce.in==3) {
    type <- 1  # 0 or 1, if type 1 then one outcome
    corr.type <- 0 # 0 or 1, if 1 then high cor
  } else if (sce.in == 4) {
    type <- 1  # 0 or 1, if type 1 then one outcome
    corr.type <- 1 # 0 or 1, if 1 then high cor
  }
  
  
  if (type == 0) {
    type.name <- c("3")
  } else if (type == 1) {
    type.name <- c("1")
  }
  
  if (corr.type == 0) {
    corr.name <- c("indp")
  } else if (corr.type == 1) {
    corr.name <- c("cor")
  }
  


  
  
  # 3 correlated outcomes: y1t, y1c ; y2t, y2c ; y3t, y3c
  # correlation types: 0 = independent, 1 = high-correlation
  if (corr.type == 0) {
    rho1=0.0  #correlation between y1t and y1c ;  y2t and y2c ;  y3t and y3c 
    rho2=0.0  #correlation between y1t and y2t ;  y1t and y3t ;  y2t and y3t ; y1c and y2c ;  y1c and y3c ;  y2c and y3c 
    rho3 = rho1*rho2 #correlation between y1t and y2c ;  y1t and y3c ;  y2t and y1c ; y2t and y3c ;  y3t and y1c ;  y3t and y2c 
    sig2=sig2t=sig2c=1 # common, known variance  
  } else if (corr.type == 1) {
    rho1=0.0  #correlation between y1t and y1c ;  y2t and y2c ;  y3t and y3c 
    rho2=0.8  #correlation between y1t and y2t ;  y1t and y3t ;  y2t and y3t ; y1c and y2c ;  y1c and y3c ;  y2c and y3c 
    rho3 = rho1*rho2 #correlation between y1t and y2c ;  y1t and y3c ;  y2t and y1c ; y2t and y3c ;  y3t and y1c ;  y3t and y2c 
    sig2=sig2t=sig2c=1 # common, known variance 
  }


  # if (corr.type == 0) {
  #   # data.dir <- c("C:\\School\\Research\\thesis\\Paper 1\\dataset\\")
  #   data.dir <- paste(in.dir, "dataset/",  sep='')
  # } else if (corr.type == 1) {
  #   # data.dir <- c("C:\\School\\Research\\thesis\\Paper 1\\dataset_cor\\")
  #   data.dir <- paste(in.dir, "dataset_cor/",  sep='')
  # }
  
  data.dir <- paste(in.dir, "dataset_biocor/",sep='')
  
  # 3 correlated outcomes: y1t, y1c ; y2t, y2c ; y3t, y3c
  # correlation types: 0 = independent, 1 = high-correlation
  if (corr.type == 0) {
    rho1=0.0  #correlation between y1t and y1c ;  y2t and y2c ;  y3t and y3c 
    rho2=0.0  #correlation between y1t and y2t ;  y1t and y3t ;  y2t and y3t ; y1c and y2c ;  y1c and y3c ;  y2c and y3c 
    rho3 = rho1*rho2 #correlation between y1t and y2c ;  y1t and y3c ;  y2t and y1c ; y2t and y3c ;  y3t and y1c ;  y3t and y2c 
    sig2=sig2t=sig2c=1 # common, known variance  
  } else if (corr.type == 1) {
    rho1=0.0  #correlation between y1t and y1c ;  y2t and y2c ;  y3t and y3c 
    rho2=0.8  #correlation between y1t and y2t ;  y1t and y3t ;  y2t and y3t ; y1c and y2c ;  y1c and y3c ;  y2c and y3c 
    rho3 = rho1*rho2 #correlation between y1t and y2c ;  y1t and y3c ;  y2t and y1c ; y2t and y3c ;  y3t and y1c ;  y3t and y2c 
    sig2=sig2t=sig2c=1 # common, known variance 
  }
  

  
  for (o in 1:length(outcome.full)) {
    outcome <- outcome.full[[o]] 
    m.out <- m.out.full[[o]] 
    
    
    for (m in 1:length(model.full)) {
      model.out <- model.full[[m]]
      
      if (m == 1) {
        # 0
        c01=0.0 # true cutoff for biomarker 1
        c02=0.0 # true cutoff for biomarker 2
        c03=0.0 # true cutoff for biomarker 3
        c04=0.0 # true cutoff for biomarker 4
        ### Change point Model: beta0*rn.intercept + t*theta+ t*delta*I(X1>c01 & X2>c02 & X3 + X4 > c03 + c04)
        # normal
        beta0=0.1 # flat control
        theta=0.0 # flat trt
        delta=0.0 # treatment effect for the best subgroup
        # binary
        p.beta0 <- 0.025 # flat control
        p.theta <- 0.0 # flat trt
        p.delta <- 0.0 # treatment effect for the best subgroup
        # survival
        s.beta0 <- 0.016 # flat control
        s.theta <- 0.0 # flat trt
        s.delta <- 0.0 # treatment effect for the best subgroup
      }
      
      if (m == 2) {
        # 1a
        c01=0.5 # true cutoff for biomarker 1
        c02=0.0 # true cutoff for biomarker 2
        c03=0.0 # true cutoff for biomarker 3
        c04=0.0 # true cutoff for biomarker 4
        
        beta0=0.1 # flat control
        theta=0.0 # flat trt
        delta=0.3*2 # treatment effect for the best subgroup
        # binary
        p.beta0 <- 0.025 # flat control
        p.theta <- 0.0 # flat trt
        p.delta <- 0.675*2  # treatment effect for the best subgroup
        # survival
        s.beta0 <- 0.016 # flat control
        s.theta <- 0.0 # flat trt
        s.delta <- 0.5*2  # treatment effect for the best subgroup
      }
      
      if (m == 3) {
        # 2b
        c01=0.0 # true cutoff for biomarker 1
        c02=0.0 # true cutoff for biomarker 2
        c03=0.5 # true cutoff for biomarker 3
        c04=0.5 # true cutoff for biomarker 4
        
        # normal
        beta0=0.1 # flat control
        theta=0.0 # flat trt
        delta=0.3*2 # treatment effect for the best subgroup
        # binary
        p.beta0 <- 0.025 # flat control
        p.theta <- 0.0 # flat trt
        p.delta <- 0.675*2  # treatment effect for the best subgroup
        # survival
        s.beta0 <- 0.016 # flat control
        s.theta <- 0.0 # flat trt
        s.delta <- 0.5*2  # treatment effect for the best subgroup
        
      }
      
      if (m == 4) {
        # 3c
        c01=0.0 # true cutoff for biomarker 1
        c02=1/3 # true cutoff for biomarker 2
        c03=1/2*sqrt(0.5) # true cutoff for biomarker 3
        c04=1/2*sqrt(0.5) # true cutoff for biomarker 4
        
        # normal
        beta0=0.1 # flat control
        theta=0.0 # flat trt
        delta=0.3*2 # treatment effect for the best subgroup
        # binary
        p.beta0 <- 0.025 # flat control
        p.theta <- 0.0 # flat trt
        p.delta <- 0.675*2  # treatment effect for the best subgroup
        # survival
        s.beta0 <- 0.016 # flat control
        s.theta <- 0.0 # flat trt
        s.delta <- 0.5*2  # treatment effect for the best subgroup
        
        
      }
      
      
      if (m == 5) {
        # 4d
        c01=1/5 # true cutoff for biomarker 1
        c02=1/6 # true cutoff for biomarker 2
        c03=1/2*sqrt(0.5) # true cutoff for biomarker 3
        c04=1/2*sqrt(0.5) # true cutoff for biomarker 4
        
        # normal
        beta0=0.1 # flat control
        theta=0.0 # flat trt
        delta=0.3*2 # treatment effect for the best subgroup
        # binary
        p.beta0 <- 0.025 # flat control
        p.theta <- 0.0 # flat trt
        p.delta <- 0.675*2  # treatment effect for the best subgroup
        # survival
        s.beta0 <- 0.016 # flat control
        s.theta <- 0.0 # flat trt
        s.delta <- 0.5*2  # treatment effect for the best subgroup
        
      }
      
      if (m==6) {
        # full
        c01=0.0 # true cutoff for biomarker 1
        c02=0.0 # true cutoff for biomarker 2
        c03=0.0 # true cutoff for biomarker 3
        c04=0.0 # true cutoff for biomarker 4
        
        # normal
        beta0=0.1 # flat control
        theta=0.0 # flat trt
        delta=0.3 # treatment effect for the best subgroup
        # binary
        p.beta0 <- 0.025 # flat control
        p.theta <- 0.0 # flat trt
        p.delta <- 0.675 # treatment effect for the best subgroup
        # survival
        s.beta0 <- 0.016 # flat control
        s.theta <- 0.0 # flat trt
        s.delta <- 0.5 # treatment effect for the best subgroup
        
        
      }
      
      # parameter selection
      c.full <- c(1,2,3,4)
      c.true <- which(c(c01,c02,c03,c04) != 0)
      true.biom1.count <- length(which(c(c01,c02,c03,c04) != 0))
      
      
      
      # Create a list of variable names
      vars <- c("beta0", "theta", "delta", "p.beta0", "p.theta", "p.delta", "s.beta0", "s.theta", "s.delta")
      
      # Loop through each variable in the list
      for (var in vars) {
        # Create a new variable name by appending '.ext'
        new_var <- paste0(var, ".ext")
        
        # Assign the value of the original variable to the new variable
        assign(new_var, get(var))
      }
      
      
      
      if (type == 1) {
        # # only 1 effect
        delta.index <- sample(1:3,1)
        if (delta.index == 1) {
          p.delta <- 0
          s.delta <- 0
        } else if (delta.index == 2) {
          delta <- 0
          s.delta <- 0
        } else {
          delta <- 0
          p.delta <- 0
        }
      }
  
  
  
      # initialize vectors
      zs.a <- zs.p <- z1s.a <- z2s.a <- z1s.p <- z2s.p <- list()
      z.a.raw.hoch <- z.p.raw.hoch <- z.a.raw.weight <- z.p.raw.weight <- rep(NA,S)
      z.a.raw.hoch.1 <- z.a.raw.hoch.2 <-  z.p.raw.hoch.1 <-  z.p.raw.hoch.2 <- rep(NA,S)
      zall.full.1 <- zall.full.2 <- zall.full.3 <- zall.a.cv.1 <- zall.a.cv.2 <- zall.a.cv.3 <- rep(NA,S)
      zall.a.hoch.1 <- zall.a.hoch.2 <-  rep(NA,S)
      zall.a.hoch <- zall.a.max <-  zall.a.avg <-  zall.a.weight <- rep(NA,S)
      z.p.all <- rep(NA,S)
      z.a.hoch <- rep(NA,S)
      z.p.raw.hoch.nocv <- z.p.hoch <- z.p.raw.weight.1 <- z.p.raw.weight.2 <- rep(NA,S)
      p.a.permute <- p.a.permute.1 <- p.a.permute.2 <- p.a.permute.m1 <- p.a.permute.m2 <- rep(NA,S)
      p.p.permute <- p.p.permute.1 <- p.p.permute.2 <- p.p.permute.m1 <- p.p.permute.m2 <- rep(NA,S)
      
      U.all <- U.best <- U.perc <- rep(NA,S) # estimate true test statistics, generic
      n.U.all <- n.U.best <- n.U.perc <- rep(NA,S) # estimate true test statistics, continuous
      b.U.all <- b.U.best <- b.U.perc <- rep(NA,S) # estimate true test statistics, binary
      s.U.all <- s.U.best <- s.U.perc <- rep(NA,S) # estimate true test statistics, survival
      
      perc.U.cv1.a <- perc.U.cv2.a <- perc.U.cv1.p <- perc.U.cv2.p <- rep(NA,S)
      pi.sub.cv1.a <- pi.sub.cv2.a <- pi.sub.cv1.p <- pi.sub.cv2.p <- rep(NA,S)
      perc.U.cv1.a.flip <- perc.U.cv2.a.flip <- perc.U.cv1.p.flip <- perc.U.cv2.p.flip <- rep(NA,S)
      pi.sub.cv1.a.flip <- pi.sub.cv2.a.flip <- pi.sub.cv1.p.flip <- pi.sub.cv2.p.flip <- rep(NA,S)
      
      coeff.cv1.a.1 <- coeff.cv1.a.2 <- coeff.cv1.a.3 <- coeff.cv1.a.4 <- rep(NA,S) 
      coeff.cv1.p.1 <- coeff.cv1.p.2 <- coeff.cv1.p.3 <- coeff.cv1.p.4 <-rep(NA,S) 
      coeff.cv2.a.1 <- coeff.cv2.a.2 <- coeff.cv2.a.3 <- coeff.cv2.a.4 <- rep(NA,S) 
      coeff.cv2.p.1 <- coeff.cv2.p.2 <- coeff.cv2.p.3 <- coeff.cv2.p.4 <-rep(NA,S) 
      perc.x.cv1.a.ex<- perc.x.cv2.a.ex<- perc.x.cv1.p.ex<- perc.x.cv2.p.ex<- rep(NA,S) 
      perc.x.cv1.a.sub<- perc.x.cv2.a.sub<- perc.x.cv1.p.sub<- perc.x.cv2.p.sub<- rep(NA,S) 
      
      
      coeff.cv1.a.1.flip <- coeff.cv1.a.2.flip <- coeff.cv1.a.3.flip <- coeff.cv1.a.4.flip <- rep(NA,S) 
      coeff.cv1.p.1.flip <- coeff.cv1.p.2.flip <- coeff.cv1.p.3.flip <- coeff.cv1.p.4.flip <-rep(NA,S) 
      coeff.cv2.a.1.flip <- coeff.cv2.a.2.flip <- coeff.cv2.a.3.flip <- coeff.cv2.a.4.flip <- rep(NA,S) 
      coeff.cv2.p.1.flip <- coeff.cv2.p.2.flip <- coeff.cv2.p.3.flip <- coeff.cv2.p.4.flip <-rep(NA,S) 
      perc.x.cv1.a.ex.flip<- perc.x.cv2.a.ex.flip<- perc.x.cv1.p.ex.flip<- perc.x.cv2.p.ex.flip<- rep(NA,S) 
      perc.x.cv1.a.sub.flip<- perc.x.cv2.a.sub.flip<- perc.x.cv1.p.sub.flip<- perc.x.cv2.p.sub.flip<- rep(NA,S) 
      
      
      
      weight.vec.out <- list()
      
      
      #########
      ## Runs
      #########
      
      for(s in 1:S)
      {
        # set.seed((mn-1)*S+s)
        
        cat("s:",s,"\n")
        
        datnum <- c((mn-1)*S+s)
        
        data.in <- readRDS(file = paste(data.dir, "data_", datnum, ".Rds", sep=""))
        
        # 2Kx2K variance-covariance matrix 
        sig = matrix(c(sig2t,rho1,rho2,rho3,rho2,rho3,
                       rho1,sig2c,rho3,rho2,rho3,rho2,
                       rho2,rho3,sig2t,rho1,rho2,rho3,
                       rho3,rho2,rho1,sig2c,rho3,rho2,
                       rho2,rho3,rho2,rho3,sig2t,rho1,
                       rho3,rho2,rho3,rho2,rho1,sig2c),6,6)
        # treatment indicator matrix
        t <- c(rep(1,n1/2),rep(0,n1/2)) 
        # initialize 4 biomarkers
        biom1.x1.t=data.in[[1]]
        biom1.x1.c=data.in[[2]]
        biom1.x2.t=data.in[[3]]
        biom1.x2.c=data.in[[4]]
        biom1.x3.t=data.in[[5]]
        biom1.x3.c=data.in[[6]]
        biom1.x4.t=data.in[[7]]
        biom1.x4.c=data.in[[8]]
        biom1.x1 <- c(biom1.x1.t,  biom1.x1.c)
        biom1.x2 <- c(biom1.x2.t,  biom1.x2.c)
        biom1.x3 <- c(biom1.x3.t,  biom1.x3.c)
        biom1.x4 <- c(biom1.x4.t,  biom1.x4.c)
        t.biom1.x1 <- c(biom1.x1.t,  biom1.x1.c)*t
        t.biom1.x2 <- c(biom1.x2.t,  biom1.x2.c)*t
        t.biom1.x3 <- c(biom1.x3.t,  biom1.x3.c)*t
        t.biom1.x4 <- c(biom1.x4.t,  biom1.x4.c)*t
        t.biom1.x1.biom1.x2 <- c(biom1.x1.t,  biom1.x1.c)*c(biom1.x2.t,  biom1.x2.c)*t
        t.biom1.x1.biom1.x3 <- c(biom1.x1.t,  biom1.x1.c)*c(biom1.x3.t,  biom1.x3.c)*t
        t.biom1.x1.biom1.x4 <- c(biom1.x1.t,  biom1.x1.c)*c(biom1.x4.t,  biom1.x4.c)*t
        t.biom1.x2.biom1.x3 <- c(biom1.x2.t,  biom1.x2.c)*c(biom1.x3.t,  biom1.x3.c)*t
        t.biom1.x2.biom1.x4 <- c(biom1.x2.t,  biom1.x2.c)*c(biom1.x4.t,  biom1.x4.c)*t
        t.biom1.x3.biom1.x4 <- c(biom1.x3.t,  biom1.x3.c)*c(biom1.x4.t,  biom1.x4.c)*t
        biom1=c(biom1.x1.t,biom1.x1.c,biom1.x2.t,biom1.x2.c,biom1.x3.t,biom1.x3.c,biom1.x4.t,biom1.x4.c)
        
        
        for (i in 1:length(biom1)) {
          if (biom1[i]<0) {
            biom1[i] = 0  # negative becomes 0
          } else if (biom1[i]>1) {
            biom1[i] = 1 # greater than 1 becomes 1
          }
        } 
        
        #Initialize 3 outcomes 
        y1.t<-y1.c<-y2.t<-y2.c<-y3.t<-y3.c<- mu.t.1<-mu.c.1<-mu.t.2<-mu.c.2<-mu.t.3<-mu.c.3<-rep(NA,n1/2)
        
       
        
        
        # Desired means
        mean.cont.t <- beta0 + theta + (delta*ifelse(biom1.x1.t>c01 & biom1.x2.t>c02 & biom1.x3.t + biom1.x4.t > c03 + c04, 1,0))
        mean.cont.c <- beta0 + theta
        
        desired_mean_beta.t <-  1/(1+exp(-(p.beta0 + p.theta + p.delta*ifelse(biom1.x1.t>c01 & biom1.x2.t>c02 & biom1.x3.t + biom1.x4.t > c03 + c04, 1,0))))
        desired_mean_beta.c <-  1/(1+exp(-(p.beta0 + p.theta)))
        
        Tlat.t <- (- log(runif(n1/2)) / (0.001 * exp(s.beta0 + s.theta + s.delta*ifelse(biom1.x1.t>c01 & biom1.x2.t>c02 & biom1.x3.t + biom1.x4.t > c03 + c04, 1,0))))^(1 / 1)
        Tlat.c <- (- log(runif(n1/2)) / (0.001 * exp(s.beta0 + s.theta)))^(1 / 1)
        
        desired_mean_weibull.t <-  Tlat.t 
        desired_mean_weibull.c <-   Tlat.c  
        
        # Shape parameters
        shape_weibull <- 1
        shape_beta <- 1
        
        # Scale parameter for Weibull to achieve desired mean
        scale_weibull.t <- desired_mean_weibull.t  / gamma(1 + 1 / shape_weibull)
        scale_weibull.c <- desired_mean_weibull.c  / gamma(1 + 1 / shape_weibull)
        
        # Second shape parameter for Beta to achieve desired mean
        shape2_beta.t <- shape_beta * ((1 / desired_mean_beta.t ) - 1)
        shape2_beta.c <- shape_beta * ((1 / desired_mean_beta.c ) - 1)
        
        
        
        copula <- normalCopula(param = P2p(sig), dim = 6, dispstr = "un")
        
        marginals <- mvdc(copula = copula, 
                          margins = c("norm", "norm", "beta",  "beta","weibull","weibull"),
                          paramMargins = list(list(mean = mean.cont.t, sd = 1),
                                              list(mean = mean.cont.c, sd = 1),
                                              list(shape1 = shape_beta, shape2 = shape2_beta.t), 
                                              list(shape1 = shape_beta, shape2 = shape2_beta.c), 
                                              list(shape = shape_weibull, scale = scale_weibull.t),
                                              list(shape = shape_weibull, scale = scale_weibull.c)
                                              ))
        samples <- rMvdc(n1/2, marginals)
        
        # calculate the threshold that makes the mean equal to desired_mean_beta.t
        threshold <- quantile(samples[c(ifelse(biom1.x1.t>c01 & biom1.x2.t>c02 & biom1.x3.t + biom1.x4.t > c03 + c04, 1,0)) == 1, 3], 1 -
                                max(desired_mean_beta.t))
        # convert to binary
        samples[c(ifelse(biom1.x1.t>c01 & biom1.x2.t>c02 & biom1.x3.t + biom1.x4.t > c03 + c04, 1,0)) == 1, 3] <-
          ifelse(samples[c(ifelse(biom1.x1.t>c01 & biom1.x2.t>c02 & biom1.x3.t + biom1.x4.t > c03 + c04, 1,0)) == 1, 3] > threshold, 1, 0)


        # calculate the threshold that makes the mean equal to desired_mean_beta.t
        threshold <- quantile(samples[c(ifelse(biom1.x1.t>c01 & biom1.x2.t>c02 & biom1.x3.t + biom1.x4.t > c03 + c04, 0,1)) == 1, 3], 1 -
                                desired_mean_beta.c)
        # convert to binary
        samples[c(ifelse(biom1.x1.t>c01 & biom1.x2.t>c02 & biom1.x3.t + biom1.x4.t > c03 + c04, 0,1)) == 1, 3] <-
          ifelse(samples[c(ifelse(biom1.x1.t>c01 & biom1.x2.t>c02 & biom1.x3.t + biom1.x4.t > c03 + c04, 0,1)) == 1, 3] > threshold, 1, 0)

        
        # calculate the threshold that makes the mean equal to desired_mean_beta.t
        threshold <- quantile(samples[, 4], 1 - desired_mean_beta.c)
        # convert to binary
        samples[, 4] <- ifelse(samples[, 4] > threshold, 1, 0)
    
        
        q75.t <- quantile(samples[, 5], 0.75)
        q75.c <- quantile(samples[, 6], 0.75)
        
        t.rateC <- 1/(2* q75.t)
        c.rateC <- 1/(2 * q75.c)

        C_time.t <-rexp(n1/2, rate=t.rateC)
        C_time.c <-rexp(n1/2, rate=c.rateC)

        E.t <- as.numeric( samples[, 5] <= C_time.t)  # event indicator
        E.c <- as.numeric( samples[, 6] <= C_time.c)  # event indicator

        samples[, 5] <- pmin(samples[, 5], C_time.t)  # observed time
        samples[, 6] <- pmin(samples[, 6], C_time.c)  # observed time
        
        
        
        y1.t <- samples[,1]
        y1.c <- samples[,2]
        y2.t <- samples[,3]
        y2.c <- samples[,4]
        y3.t <- data.frame(id=1:(n1/2),
                              time= samples[, 5],
                              status=E.t ,
                              x=ifelse(biom1.x1.t>c01 & biom1.x2.t>c02 & biom1.x3.t + biom1.x4.t > c03 + c04,1,0),
                              trt=1)
        y3.c <- data.frame(id=1:(n1/2),
                           time= samples[, 6],
                           status=E.c ,
                           x=ifelse(biom1.x1.c>c01 & biom1.x2.c>c02 & biom1.x3.c + biom1.x4.c > c03 + c04,1,0),
                           trt=0)
        
          

        
        y1.full <- c(y1.t, y1.c)
        y2.full <- c(y2.t, y2.c)
        y3.full <- rbind(y3.t, y3.c)
        
        summary(y3.full[which(y3.full$x == 1 & y3.full$trt == 1),])
        summary(y3.full[which(y3.full$x == 1 & y3.full$trt == 0),])
        
        # prep for GLASSO
        x.mat <- data.frame(cbind(t.biom1.x1 , t.biom1.x2 , t.biom1.x3 , t.biom1.x4 ,
                                  t.biom1.x1.biom1.x2 , t.biom1.x1.biom1.x3 , t.biom1.x1.biom1.x4 , t.biom1.x2.biom1.x3 , t.biom1.x2.biom1.x4 , t.biom1.x3.biom1.x4 , t
                                  , biom1.x1 , biom1.x2 , biom1.x3 , biom1.x4 , biom1.x1*biom1.x2 , biom1.x1*biom1.x3 , biom1.x1*biom1.x4 , biom1.x2*biom1.x3
                                  , biom1.x2*biom1.x4 , biom1.x3*biom1.x4,
                                  biom1.x1*biom1.x1, biom1.x2*biom1.x2, biom1.x3*biom1.x3, biom1.x4*biom1.x4))
        colnames(x.mat)[16:25] <- c("biom1.x1:biom1.x2" , "biom1.x1:biom1.x3" , "biom1.x1:biom1.x4" , "biom1.x2:biom1.x3"
                                    , "biom1.x2:biom1.x4" , "biom1.x3:biom1.x4",
                                    "biom1.x1:biom1.x1" , "biom1.x2:biom1.x2", "biom1.x3:biom1.x3", "biom1.x4:biom1.x4")
        x.mat.fit.t <- data.frame(cbind(1, 1*biom1.x1 , 1*biom1.x2 , 1*biom1.x3 , 1*biom1.x4 ,
                                        1*biom1.x1*biom1.x2 , 1*biom1.x1*biom1.x3 , 1*biom1.x1*biom1.x4 , 1*biom1.x2*biom1.x3 , 1*biom1.x2*biom1.x4 ,
                                        1*biom1.x3*biom1.x4 , 1, biom1.x1 , biom1.x2 , biom1.x3 , biom1.x4 , biom1.x1*biom1.x2 , biom1.x1*biom1.x3 ,
                                        biom1.x1*biom1.x4 , biom1.x2*biom1.x3, biom1.x2*biom1.x4 , biom1.x3*biom1.x4,
                                        biom1.x1*biom1.x1, biom1.x2*biom1.x2, biom1.x3*biom1.x3, biom1.x4*biom1.x4))
        colnames(x.mat.fit.t) <- c("intercept",  colnames(x.mat))
        
        grp.mat <- list(c(11),                       # main effect for treatment
                        c(12), c(13), c(14), c(15),  # main effects for biomarkers
                        c(12,22),c(13,23),c(14,24),c(15,25),   # quadratic terms for biomarkers
                        c(1,12,11),c(2,13,11),c(3,14,11),c(4,15,11),  # main effects and interaction with t
                        c(12,13,16),c(12,14,17),c(12,15,18),c(13,14,19),c(13,15,20),c(14,15,21), # main effects interactions
                        c(1,12,2,13,5,16,11),c(1,12,3,14,6,17,11),c(1,12,4,15,7,18,11),c(2,13,3,14,8,19,11),c(2,13,4,15,9,20,11),c(3,14,4,15,10,21,11)) # main effects interactions, and with t
        
        
        
        
        ### computing %U to report by generating 20000 obsevations
        n.obs = 10000
        X1.final.t=runif(n.obs,0,1) + rnorm(n.obs,0,bio.sigma^2) 
        X1.final.c=runif(n.obs,0,1) + rnorm(n.obs,0,bio.sigma^2) 
        X2.final.t=runif(n.obs,0,1) + rnorm(n.obs,0,bio.sigma^2) 
        X2.final.c=runif(n.obs,0,1) + rnorm(n.obs,0,bio.sigma^2) 
        X3.final.t=runif(n.obs,0,1) + rnorm(n.obs,0,bio.sigma^2) 
        X3.final.c=runif(n.obs,0,1) + rnorm(n.obs,0,bio.sigma^2) 
        X4.final.t=runif(n.obs,0,1) + rnorm(n.obs,0,bio.sigma^2) 
        X4.final.c=runif(n.obs,0,1) + rnorm(n.obs,0,bio.sigma^2) 
        X.final=c(X1.final.t,X1.final.c,X2.final.t,X2.final.c,X3.final.t,X3.final.c,X4.final.t,X4.final.c)
        # treatment indicator matrix
        t <- c(rep(1,n.obs),rep(0,n.obs)) 
        # modify covariate
        X1.final=c(X1.final.t,X1.final.c)
        X2.final=c(X2.final.t,X2.final.c)
        X3.final=c(X3.final.t,X3.final.c)
        X4.final=c(X4.final.t,X4.final.c)
        t.X1.final <- X1.final*t
        t.X2.final <- X2.final*t
        t.X3.final <- X3.final*t
        t.X4.final <- X4.final*t
        t.X1.X2.final <- X1.final*X2.final*t
        t.X1.X3.final <- X1.final*X3.final*t
        t.X1.X4.final <- X1.final*X4.final*t
        t.X2.X3.final <- X2.final*X3.final*t
        t.X2.X4.final <- X2.final*X4.final*t
        t.X3.X4.final <- X3.final*X4.final*t
        
        x.mat.final <-data.frame(cbind(1, t.X1.final, t.X2.final, t.X3.final, t.X4.final ,
                                       t.X1.X2.final, t.X1.X3.final, t.X1.X4.final, t.X2.X3.final, t.X2.X4.final, t.X3.X4.final, t
                                       , X1.final , X2.final , X3.final , X4.final , X1.final*X2.final , X1.final*X3.final , X1.final*X4.final , X2.final*X3.final
                                       , X2.final*X4.final , X3.final*X4.final,
                                       X1.final*X1.final, X2.final*X2.final, X3.final*X3.final, X4.final*X4.final))
        x.mat.final.t <- data.frame(cbind(1, 1*X1.final , 1*X2.final , 1*X3.final , 1*X4.final ,
                                          1*X1.final*X2.final , 1*X1.final*X3.final , 1*X1.final*X4.final , 1*X2.final*X3.final , 1*X2.final*X4.final , 
                                          1*X3.final*X4.final , 1, X1.final , X2.final , X3.final , X4.final , X1.final*X2.final , X1.final*X3.final , 
                                          X1.final*X4.final , X2.final*X3.final, X2.final*X4.final , X3.final*X4.final,
                                          X1.final*X1.final, X2.final*X2.final, X3.final*X3.final, X4.final*X4.final))
        
        
        for (i in 1:length(X.final)) {
          if (X.final[i]<0) {
            X.final[i] = 0  # negative becomes 0
          } else if (X.final[i]>1) {
            X.final[i] = 1 # greater than 1 becomes 1
          }
        }  
        
        # 
        # if (corr.type == 0) {
        #   rn.intercept.final <- matrix(0, nrow = n.obs, ncol = 6)                                  
        # } else if (corr.type == 1) {
        #   rn.intercept.final <- mvrnorm(n.obs,rep(0,6),sig)                                         
        # }
        # 
        
        
        # Desired means
        mean.cont.t <- beta0.ext+theta.ext+(delta.ext*ifelse(X1.final.t>c01 & X2.final.t>c02 & X3.final.t + X4.final.t > c03 + c04, 1,0))
        mean.cont.c <- beta0.ext+theta.ext
        
        desired_mean_beta.t <-  1/(1+exp(-(p.beta0.ext + p.theta.ext + p.delta.ext*ifelse(X1.final.t>c01 & X2.final.t>c02 & X3.final.t + X4.final.t > c03 + c04,1,0))))
        desired_mean_beta.c <-  1/(1+exp(-(p.beta0.ext + p.theta.ext)))
        
        Tlat.t <- (- log(runif(n.obs)) / (0.001 * exp(s.beta0.ext + s.theta.ext + s.delta.ext*ifelse(X1.final.t>c01 & X2.final.t>c02 & X3.final.t + X4.final.t > c03 + c04,1,0))))^(1 / 1)
        Tlat.c <- (- log(runif(n.obs)) / (0.001 * exp(s.beta0.ext + s.theta.ext)))^(1 / 1)
        
        desired_mean_weibull.t <- Tlat.t 
        desired_mean_weibull.c <- Tlat.c  
        
        # Shape parameters
        shape_weibull <- 1
        shape_beta <- 1
        
        # Scale parameter for Weibull to achieve desired mean
        scale_weibull.t <- desired_mean_weibull.t  / gamma(1 + 1 / shape_weibull)
        scale_weibull.c <- desired_mean_weibull.c  / gamma(1 + 1 / shape_weibull)
        
        # Second shape parameter for Beta to achieve desired mean
        shape2_beta.t <- shape_beta * ((1 / desired_mean_beta.t ) - 1)
        shape2_beta.c <- shape_beta * ((1 / desired_mean_beta.c ) - 1)
        
        
        
        copula <- normalCopula(param = P2p(sig), dim = 6, dispstr = "un")
        
        marginals <- mvdc(copula = copula, 
                          margins = c("norm", "norm", "beta",  "beta","weibull","weibull"),
                          paramMargins = list(list(mean = mean.cont.t, sd = 1),
                                              list(mean = mean.cont.c + theta, sd = 1),
                                              list(shape1 = shape_beta, shape2 = shape2_beta.t), 
                                              list(shape1 = shape_beta, shape2 = shape2_beta.c), 
                                              list(shape = shape_weibull, scale = scale_weibull.t),
                                              list(shape = shape_weibull, scale = scale_weibull.c)
                          ))
        samples <- rMvdc(n.obs, marginals)
        
        # calculate the threshold that makes the mean equal to desired_mean_beta.t
        threshold <- quantile(samples[c(ifelse(X1.final.t>c01 & X2.final.t>c02 & X3.final.t + X4.final.t > c03 + c04, 1,0)) == 1, 3], 1 - 
                                max(desired_mean_beta.t))
        # convert to binary
        samples[c(ifelse(X1.final.t>c01 & X2.final.t>c02 & X3.final.t + X4.final.t > c03 + c04, 1,0))== 1, 3] <- 
          ifelse(samples[c(ifelse(X1.final.t>c01 & X2.final.t>c02 & X3.final.t + X4.final.t > c03 + c04, 1,0))== 1, 3] > threshold, 1, 0)
        
        
        # calculate the threshold that makes the mean equal to desired_mean_beta.t
        threshold <- quantile(samples[c(ifelse(X1.final.t>c01 & X2.final.t>c02 & X3.final.t + X4.final.t > c03 + c04, 0,1)) == 1, 3], 1 - 
                                desired_mean_beta.c)
        # convert to binary
        samples[c(ifelse(X1.final.t>c01 & X2.final.t>c02 & X3.final.t + X4.final.t > c03 + c04, 0,1))== 1, 3] <- 
          ifelse(samples[c(ifelse(X1.final.t>c01 & X2.final.t>c02 & X3.final.t + X4.final.t > c03 + c04, 0,1))== 1, 3] > threshold, 1, 0)
        
        
        
        
        # calculate the threshold that makes the mean equal to desired_mean_beta.t
        threshold <- quantile(samples[, 4], 1 - desired_mean_beta.c)
        # convert to binary
        samples[, 4] <- ifelse(samples[, 4] > threshold, 1, 0)
        
        
        q75.t <- quantile(samples[, 5], 0.75)
        q75.c <- quantile(samples[, 6], 0.75)
        
        t.rateC <- 1/(2* q75.t)
        c.rateC <- 1/(2* q75.c)
        
        C_time.t <-rexp(n1/2, rate=t.rateC)
        C_time.c <-rexp(n1/2, rate=c.rateC)
        
        E.t <- as.numeric( samples[, 5] <= C_time.t)  # event indicator
        E.c <- as.numeric( samples[, 6] <= C_time.c)  # event indicator
        
        samples[, 5] <- pmin(samples[, 5], C_time.t)  # observed time
        samples[, 6] <- pmin(samples[, 6], C_time.c)  # observed time
        
        
        
        n.y.final.t <- samples[,1]
        n.y.final.c <- samples[,2]
        b.y.final.t <- samples[,3]
        b.y.final.c <- samples[,4]
        s.y.final.t <- data.frame(id=1:(n.obs),
                           time= samples[, 5],
                           status=E.t ,
                           x=c(ifelse(X1.final.t>c01 & X2.final.t>c02 & X3.final.t + X4.final.t > c03 + c04, 1,0)),
                           trt=1)
        s.y.final.c <- data.frame(id=1:(n.obs),
                           time= samples[, 6],
                           status=E.c ,
                           x=c(ifelse(X1.final.c>c01 & X2.final.c>c02 & X3.final.c + X4.final.c > c03 + c04, 1,0)),
                           trt=0)
        s.y.final <- rbind(s.y.final.t, s.y.final.c)
 
        # summary(s.y.final[which(s.y.final$x == 1 & s.y.final$trt == 1),])
        # summary(s.y.final[which(s.y.final$x == 1 & s.y.final$trt == 0),])
        
        
        
        
        n.U.all[s] <- t.test(n.y.final.t, n.y.final.c)$statistic/sqrt(2*n.obs/2)
        n.U.best[s] <- t.test(n.y.final.t[which(X1.final.t>c01 & X2.final.t>c02 & X3.final.t + X4.final.t > c03 + c04)],
                              n.y.final.c[which(X1.final.c>c01 & X2.final.c>c02 & X3.final.c + X4.final.c > c03 + c04)])$statistic/sqrt(2*n.obs/2)
        n.U.perc[s] <- n.U.all[s]/n.U.best[s]*100

        
        b.U.all[s] <- z.prop(sum(b.y.final.t), sum(b.y.final.c), length(b.y.final.t), length(b.y.final.c))/sqrt(2*n.obs/2)
        b.U.best[s] <- z.prop(sum(b.y.final.t[which(X1.final.t>c01 & X2.final.t>c02 & X3.final.t + X4.final.t > c03 + c04)]), 
                              sum(b.y.final.c[which(X1.final.c>c01 & X2.final.c>c02 & X3.final.c + X4.final.c > c03 + c04)]), 
                              length(b.y.final.t[which(X1.final.t>c01 & X2.final.t>c02 & X3.final.t + X4.final.t > c03 + c04)]), 
                              length(b.y.final.c[which(X1.final.c>c01 & X2.final.c>c02 & X3.final.c + X4.final.c > c03 + c04)]))/sqrt(2*n.obs/2)
        b.U.perc[s] <- b.U.all[s]/b.U.best[s]*100
        
 
        model <- coxph(formula = Surv(time, status) ~ trt, data = s.y.final)
        s.U.all[s] <- (coef(model)/sqrt(diag(vcov(model))))/sqrt(2*n.obs/2)
         model <- coxph(formula = Surv(time, status) ~ trt, data = s.y.final[which(s.y.final$x==1),])
         s.U.best[s] <- (coef(model)/sqrt(diag(vcov(model))))/sqrt(2*n.obs/2)
         s.U.perc[s] <- s.U.all[s]/s.U.best[s]*100
        

        
        # ind.sub <- c(which(X1.final.t>c01 & X2.final.t>c02 & X3.final.t + X4.final.t > c03 + c04), 
        #              which(X1.final.c>c01 & X2.final.c>c02 & X3.final.c + X4.final.c > c03 + c04) + n.obs)
        # plot(X1.final[ind.sub],X2.final[ind.sub],xlim=c(0,1),ylim=c(0,1)) # plot to see if makes sense
        # plot(X3.final[ind.sub],X4.final[ind.sub],xlim=c(0,1),ylim=c(0,1)) # plot to see if makes sense      
        
        
        
        ### apply test to the whole sample ###
        
        
        cv.result <- lm.cv(n1, (1:(n1-1)),  0, y1.t, y1.c, y2.t, y2.c, y3.t, y3.c)
        
        
     
        
        
        # raw coefficients
        coef.1.p <- coef.2.p <- raw.metrics(as.matrix(cv.result[[7]]))
        if (cv.result[[10]] == "cox") {
          cv.mat.p <- as.matrix(cv.result[[1]])[,-1]
        } else {  cv.mat.p <- as.matrix(cv.result[[1]])  }
        zs.p[[s]] <- perc.U(n1,y1.t, y1.c, y2.t, y2.c, y3.t, y3.c, 
                            cv.mat.p %*% as.matrix(cv.result[[7]]),cv.result[[8]],cv.result[[9]],cv.result[[11]],
                            1:(n1-1),0,
                            cv.result[[14]])
        # flip metrics
        # use flip % to refine subgroups
        # note, cv1 flips are calculated using cv1 biomarkers, but cv2 coefficients, cutoffs, signs; vice versa
        coef.1.p.flip <- coef.2.p.flip <- flip.metrics(n1, cv.mat.p, 
                                                       as.matrix(cv.result[[7]]), zs.p[[s]][[7]], cv.result[[9]], cv.result[[8]])
        
        
        
        ## %U
        # u1.a <- u2.a <- perc.U.test(as.matrix(cv.result[[2]]),cv.result[[3]],cv.result[[4]],cv.result[[6]])
        u1.p <- u2.p <- perc.U.test(as.matrix(cv.result[[7]]),cv.result[[8]],cv.result[[9]],cv.result[[11]],
                                    cv.result[[14]])
        u1.a.flip <- u2.a.flip <- u1.p.flip <- u2.p.flip <- 0
        
        
        
        zs.p.hoch.raw <- zs.p[[s]][[1]]
        
        
        
        ### create test statistics
        ## estimated subgroup
        # cross-validation for test statistics Zs
        ## raw statistics
        cv1.t.index <- data.in[[10]]
        cv2.t.index <- data.in[[11]]
        cv1.c.index <- data.in[[12]]
        cv2.c.index <- data.in[[13]]
        # assign cohort by outcome type
        if (outcome[1] == 3) {
          y1.t.cv1 <- y1.t[cv1.t.index,]
          y1.t.cv2 <- y1.t[cv2.t.index,]
          y1.c.cv1 <- y1.c[cv1.c.index,]
          y1.c.cv2 <- y1.c[cv2.c.index,]
        } else {
          y1.t.cv1 <- y1.t[cv1.t.index]
          y1.t.cv2 <- y1.t[cv2.t.index]
          y1.c.cv1 <- y1.c[cv1.c.index]
          y1.c.cv2 <- y1.c[cv2.c.index]
        }
        if (outcome[2] == 3) {
          y2.t.cv1 <- y2.t[cv1.t.index,]
          y2.t.cv2 <- y2.t[cv2.t.index,]
          y2.c.cv1 <- y2.c[cv1.c.index,]
          y2.c.cv2 <- y2.c[cv2.c.index,]
        } else {
          y2.t.cv1 <- y2.t[cv1.t.index]
          y2.t.cv2 <- y2.t[cv2.t.index]
          y2.c.cv1 <- y2.c[cv1.c.index]
          y2.c.cv2 <- y2.c[cv2.c.index]
        }
        if (outcome[3] == 3) {
          y3.t.cv1 <- y3.t[cv1.t.index,]
          y3.t.cv2 <- y3.t[cv2.t.index,]
          y3.c.cv1 <- y3.c[cv1.c.index,]
          y3.c.cv2 <- y3.c[cv2.c.index,]
        } else {
          y3.t.cv1 <- y3.t[cv1.t.index]
          y3.t.cv2 <- y3.t[cv2.t.index]
          y3.c.cv1 <- y3.c[cv1.c.index]
          y3.c.cv2 <- y3.c[cv2.c.index]
        }

        cv1.result <- lm.cv(n1/2, cv1.t.index, cv1.c.index, y1.t.cv1, y1.c.cv1, y2.t.cv1, y2.c.cv1, y3.t.cv1, y3.c.cv1)
        cv2.result <- lm.cv(n1/2, cv2.t.index, cv2.c.index, y1.t.cv2, y1.c.cv2, y2.t.cv2, y2.c.cv2, y3.t.cv2, y3.c.cv2)

        
        

        ## weighted
        cv1.mat.p <- as.matrix(cv1.result[[1]])
        cv2.mat.p <- as.matrix(cv2.result[[1]])
        if (cv2.result[[11]] == 1) {
          y.cv1.t.p <- y1.t.cv1
          y.cv1.c.p <- y1.c.cv1
        } else if (cv2.result[[11]] == 2) {
          y.cv1.t.p <- y2.t.cv1
          y.cv1.c.p <- y2.c.cv1
        } else {
          y.cv1.t.p <- y3.t.cv1
          y.cv1.c.p <- y3.c.cv1
        }
        if (cv1.result[[11]] == 1) {
          y.cv2.t.p <- y1.t.cv2
          y.cv2.c.p <- y1.c.cv2
        } else if (cv1.result[[11]] == 2) {
          y.cv2.t.p <- y2.t.cv2
          y.cv2.c.p <- y2.c.cv2
        } else {
          y.cv2.t.p <- y3.t.cv2
          y.cv2.c.p <- y3.c.cv2
        }
        if (cv2.result[[10]] == "cox") {
          cv1.mat.p <- as.matrix(cv1.result[[1]])[,-1]
        }
        if (cv1.result[[10]] == "cox") {
          cv2.mat.p <- as.matrix(cv2.result[[1]])[,-1]
        }
        # Zs
        z1s.p[[s]] <- perc.U(n1/2,y1.t.cv1, y1.c.cv1, y2.t.cv1, y2.c.cv1, y3.t.cv1, y3.c.cv1,
                             cv1.mat.p %*% as.matrix(cv2.result[[7]]),cv2.result[[8]],cv2.result[[9]],cv2.result[[11]],
                             cv1.t.index,cv1.c.index,
                             cv2.result[[14]] 
                             )
        z2s.p[[s]] <- perc.U(n1/2,y1.t.cv2, y1.c.cv2, y2.t.cv2, y2.c.cv2, y3.t.cv2, y3.c.cv2,
                             cv2.mat.p %*% as.matrix(cv1.result[[7]]),cv1.result[[8]],cv1.result[[9]],cv1.result[[11]],
                             cv2.t.index,cv2.c.index,
                             cv1.result[[14]] 
                             )

        z.p.hoch.1 <- z1s.p[[s]][[1]]
        z.p.hoch.2 <- z2s.p[[s]][[1]]
        zs.p.hoch <- sqrt(0.5)*z1s.p[[s]][[1]]+sqrt(0.5)*z2s.p[[s]][[1]]
        
        
        z.p.weight.1 <- z1s.p[[s]][[2]]
        z.p.weight.2 <- z2s.p[[s]][[2]]
        zs.p.weight <- sqrt(0.5)*z1s.p[[s]][[2]]+sqrt(0.5)*z2s.p[[s]][[2]]
        
        
        z.raw <- list(
          u1.p,u2.p,
          coef.1.p,coef.2.p,
          coef.1.p.flip,coef.2.p.flip,
          
          zs.p.hoch,
          zs.p.weight,
          z.p.hoch.1,z.p.hoch.2,
          
          zs.p.hoch.raw,
          
          z.p.weight.1,
          z.p.weight.2 
        )

        
        perc.U.cv1.p[s] <- z.raw[[1]][[2]]
        perc.U.cv2.p[s] <- z.raw[[2]][[2]]
        
        pi.sub.cv1.p[s] <- z.raw[[1]][[1]]
        pi.sub.cv2.p[s] <- z.raw[[2]][[1]]

        coeff.cv1.p.1[s] <- z.raw[[3]][1]
        coeff.cv1.p.2[s] <- z.raw[[3]][2]
        coeff.cv1.p.3[s] <- z.raw[[3]][3]
        coeff.cv1.p.4[s] <- z.raw[[3]][4]

        coeff.cv2.p.1[s] <- z.raw[[4]][1]
        coeff.cv2.p.2[s] <- z.raw[[4]][2]
        coeff.cv2.p.3[s] <- z.raw[[4]][3]
        coeff.cv2.p.4[s] <- z.raw[[4]][4]
        
        coeff.cv1.p.1.flip[s] <- z.raw[[5]][1]
        coeff.cv1.p.2.flip[s] <- z.raw[[5]][2]
        coeff.cv1.p.3.flip[s] <- z.raw[[5]][3]
        coeff.cv1.p.4.flip[s] <- z.raw[[5]][4]
        
        coeff.cv2.p.1.flip[s] <- z.raw[[6]][1]
        coeff.cv2.p.2.flip[s] <- z.raw[[6]][2]
        coeff.cv2.p.3.flip[s] <- z.raw[[6]][3]
        coeff.cv2.p.4.flip[s] <- z.raw[[6]][4]
        
        
        z.p.raw.hoch[s] <- z.raw[[7]]
        z.p.raw.weight[s] <- z.raw[[8]]
        z.p.raw.hoch.1[s] <- z.raw[[9]]
        z.p.raw.hoch.2[s] <- z.raw[[10]]
        
        z.p.raw.hoch.nocv[s] <- z.raw[[11]]
        
        z.p.raw.weight.1[s] <- z.raw[[12]]
        z.p.raw.weight.2[s] <- z.raw[[13]]
        
        weight.vec.out[[s]] <- cv.result[[14]]
        
        
        ### permutations
        z.a.permute.1 <- z.a.permute.2 <- z.a.permute <- rep(NA,P)
        z.p.permute.1 <- z.p.permute.2 <- z.p.permute <- rep(NA,P)
        
 
      } 
      
      
      U.results <- rbind(cbind(median(n.U.all, na.rm = TRUE), median(n.U.best, na.rm = TRUE), median(n.U.perc, na.rm = TRUE)),
                         cbind(median(b.U.all, na.rm = TRUE), median(b.U.best, na.rm = TRUE), median(b.U.perc, na.rm = TRUE)),
                         cbind(median(s.U.all, na.rm = TRUE), median(s.U.best, na.rm = TRUE), median(s.U.perc, na.rm = TRUE)) )
      
      U.results
      
      
      
      out.final <- paste(model.out,"_", m.out ,"_",
                         paste(type.name,"_",corr.name,sep=""),
                         sep="")
      
      
  
      output.full <- list(cbind(zall.a.cv.1, zall.a.cv.2, zall.a.cv.3, 
                                zall.a.hoch.1, zall.a.hoch.2, zall.a.hoch,
                                z.a.raw.hoch.1, z.a.raw.hoch.2, z.a.raw.hoch, 
                                
                                z.p.raw.hoch.1, z.p.raw.hoch.2, z.p.raw.hoch,
                                
                                z.p.raw.weight.1, z.p.raw.weight.2, z.p.raw.weight,
                                
                                z.p.raw.hoch.nocv,
                                
                                p.a.permute.1, p.a.permute.2, p.a.permute, 
                                p.p.permute.1, p.p.permute.2, p.p.permute,
                                
                                perc.U.cv1.a, pi.sub.cv1.a, 
                                perc.U.cv1.p, pi.sub.cv1.p, 
                                perc.U.cv2.a, pi.sub.cv2.a, 
                                perc.U.cv2.p, pi.sub.cv2.p
                                # ,
                                # perc.U.cv1.a.flip, pi.sub.cv1.a.flip, 
                                # perc.U.cv1.p.flip, pi.sub.cv1.p.flip, 
                                # perc.U.cv2.a.flip, pi.sub.cv2.a.flip, 
                                # perc.U.cv2.p.flip, pi.sub.cv2.p.flip                      
                                
                                # perc.x.cv1.a.ex, perc.x.cv1.a.sub, 
                                # perc.x.cv1.p.ex, perc.x.cv1.p.sub, 
                                # perc.x.cv2.a.ex, perc.x.cv2.a.sub, 
                                # perc.x.cv2.p.ex, perc.x.cv2.p.sub, 
                                
      ),
      
      cbind(coeff.cv1.a.1.flip , coeff.cv1.a.2.flip , coeff.cv1.a.3.flip , coeff.cv1.a.4.flip , 
            coeff.cv1.p.1.flip , coeff.cv1.p.2.flip , coeff.cv1.p.3.flip , coeff.cv1.p.4.flip ,
            coeff.cv2.a.1.flip , coeff.cv2.a.2.flip , coeff.cv2.a.3.flip , coeff.cv2.a.4.flip , 
            coeff.cv2.p.1.flip , coeff.cv2.p.2.flip , coeff.cv2.p.3.flip , coeff.cv2.p.4.flip ),
      
      cbind(coeff.cv1.a.1 , coeff.cv1.a.2 , coeff.cv1.a.3 , coeff.cv1.a.4 , 
            coeff.cv1.p.1 , coeff.cv1.p.2 , coeff.cv1.p.3 , coeff.cv1.p.4 ,
            coeff.cv2.a.1 , coeff.cv2.a.2 , coeff.cv2.a.3 , coeff.cv2.a.4 , 
            coeff.cv2.p.1 , coeff.cv2.p.2 , coeff.cv2.p.3 , coeff.cv2.p.4 ),
      weight.vec.out
      )
      write.list(output.full,  file = paste(key.dir, out.final, "_" ,mn,  ".csv", sep=""))
      
      
    }
  }
  
  
}
  
 
  
