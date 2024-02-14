# mn = 1


type.full <- c(0,1)
corr.type.full <- c(0,1)
sce <- c(1,2,3,4)  # 4 type and corr.type scenarios 

outcome.full <- list(c(1,1,1),c(1,2,3))
m.out.full <- list(c("3N"),c("1N1B1S"))

model.full <- list(c("Model1a"),c("Model2b"),c("Model3c"),c("Model4d"),c("Model5"))





n1=400 # sample size
S=100 # simulation runs
P<-2 # permutation runs


t.rateC <- 0.005
c.rateC <- 0.003

in.dir <- c("/work/users/b/e/beibo/Research/paper1/")
key.dir <- c("/work/users/b/e/beibo/Research/paper1/Jan2024/CART/")





##############################################################
###### CART - 4 biomarkers
##############################################################


library(MASS)
library(Matrix)
library(glmnet)
library(grpreg)
# library(grpregOverlap)
library(vennLasso)

library(survival)
library(prodlim)

library(Hmisc)
library(Formula)
library(ggplot2)
library(rms)
library(dplyr)
library(erer)

library(Matrix)
library(rpart)
# library(rpart.plot)
# library(rpart.utils)
library(data.table)
# set.seed(1)



rpart.rules.table<-function(object)
{
  rules<-rpart.rules(object)
  ff<-object$frame
  ff$rules<-unlist(rules[as.numeric(row.names(ff))])
  ruleList<-lapply(row.names(ff),function (name) setNames(data.frame(name,
                                                                     (strsplit(ff[name,'rules'],split=',')),
                                                                     ff[name,'var']=="<leaf>"
  ),
  c("Rule","Subrule","Leaf")))
  combinedRules<-Reduce(rbind,ruleList)
  
  return(combinedRules)
  
}


rpart.rules<-function(object)
{
  frame<-object$frame
  ruleNums<-as.numeric(row.names(frame))  ##Convert the row names into a list of rule numbers
  is.leaf <- (frame$var == "<leaf>")
  frame[!is.leaf,"order"]<-seq_along(which(!is.leaf)) ##Number the branches to number them for matching with subrule sets
  rules<-replicate(max(ruleNums),NULL)
  rules[1]<-"NULL"
  
  ##The rule numbering convention contains the information to determine branch lineage. 
  ##Most of the potential rule numbers don't actually exist, but this will result in the creation of a NULL rule.
  for (i in as.numeric(row.names(frame))[-1])
  {
    if(i%%2==0)
    {
      rules[i]<-paste(rules[i/2],paste('L',frame[as.character(i/2),"order"],sep=''),sep=',')
    }
    else
    {
      rules[i]<-paste(rules[(i-1)/2],paste('R',frame[as.character((i-1)/2),"order"],sep=''),sep=',')
    }
  }
  rules<-lapply(rules,function (x) gsub("NULL,",'',x))
  return(rules)
}


rpart.subrules.table<-function(object)  
{
  lists<-rpart.lists(object)
  leftCompares<-lapply(lists$L,function (x) attr(x,"compare"))
  rightCompares<-lapply(lists$R,function (x) attr(x,"compare"))
  leftRules<-lapply(seq_along(lists$L),function (i) setNames(data.frame(paste('L',i,sep=''),names(lists$L)[i],as.character(unlist(lists$L[i],use.names=FALSE)),NA,NA),c("Subrule","Variable","Value","Less","Greater")))
  rightRules<-lapply(seq_along(lists$R),function (i) setNames(data.frame(paste('R',i,sep=''),names(lists$R)[i],as.character(unlist(lists$R[i]),use.names=FALSE),NA,NA),c("Subrule","Variable","Value","Less","Greater")))
  
  reassign.columns<-function(object,compare)
  {
    if(grepl("<",compare))
      object$Less<-object$Value
    if(grepl(">",compare))
      object$Greater<-object$Value
    if(!grepl("=",compare))
      object$Value=NA
    return(object)
  }
  
  leftTable<-Reduce(rbind,Map(reassign.columns, leftRules, leftCompares))
  rightTable<-Reduce(rbind,Map(reassign.columns, rightRules, rightCompares))
  
  
  return(rbind(leftTable,rightTable))
}



rpart.lists <- function(object)
{
  
  ff <- object$frame
  n <- nrow(ff)
  if (n == 1L) return("root")            # special case of no splits
  
  
  ##This section  borrowed from the rpart source to identify the appropriate locations from the splits table.
  is.leaf <- (ff$var == "<leaf>")
  whichrow <- !is.leaf
  vnames <- ff$var[whichrow] # the variable names for the primary splits
  
  index <- cumsum(c(1, ff$ncompete + ff$nsurrogate + !is.leaf))
  irow <- index[c(whichrow, FALSE)] # we only care about the primary split
  ncat <- object$splits[irow, 2L]
  ##
  
  lsplit <- rsplit <- list()  
  
  if (any(ncat < 2L)) 
  {               # any continuous vars ?
    
    jrow <- irow[ncat < 2L]
    cutpoint <- object$splits[jrow, 4L]
    temp1 <- (ifelse(ncat < 0, "<", ">="))[ncat < 2L]
    temp2 <- (ifelse(ncat < 0, ">=", "<"))[ncat < 2L]
    lsplit[ncat<2L] <- cutpoint
    #lsplit[ncat<2L] <- lapply(lsplit[ncat<2L],function (x) structure(x, 'numeric'=TRUE))
    
    rsplit[ncat<2L] <- cutpoint
    #rsplit[ncat<2L] <- lapply(rsplit[ncat<2L],function (x) structure(x, 'numeric'=TRUE))
    
  }
  
  if (any(ncat > 1L)) 
  {               # any categorical variables ?
    xlevels <- attr(object, "xlevels")
    ##
    ## jrow will be the row numbers of factors within lsplit and rsplit
    ## crow the row number in "csplit"
    ## and cindex the index on the "xlevels" list
    ##
    jrow <- seq_along(ncat)[ncat > 1L]
    crow <- object$splits[irow[ncat > 1L], 4L] #row number in csplit
    cindex <- (match(vnames, names(xlevels)))[ncat > 1L]
    
    
    lsplit[jrow]<-lapply(seq_along(jrow),function (i) xlevels[[cindex[i]]][object$csplit[crow[i], ]==1L])
    rsplit[jrow]<-lapply(seq_along(jrow),function (i) xlevels[[cindex[i]]][object$csplit[crow[i], ]==3L])
    
  }
  
  
  lsplit<-lapply(seq_along(lsplit), function (i) structure(lsplit[[i]], "compare"=ifelse(ncat[i]<2L,ifelse(ncat[i]<0,"<",">="),"=")))
  rsplit<-lapply(seq_along(lsplit), function (i) structure(rsplit[[i]], "compare"=ifelse(ncat[i]<2L,ifelse(ncat[i]<0,">=","<"),"=")))
  
  
  names(lsplit)<-vnames
  names(rsplit)<-vnames
  
  results<-list("L"=lsplit,"R"=rsplit)  
  
  return(results)
}



# function to calculate z score from proportion test
z.prop = function(x1,x2,n1,n2){
  numerator = (x1/n1) - (x2/n2)
  p.common = (x1+x2) / (n1+n2)
  denominator = sqrt(p.common * (1-p.common) * (1/n1 + 1/n2))
  z.prop.ris = numerator / denominator
  return(z.prop.ris)
}


### Binary outcome
simulbin <- function(N, cov, trt, rn.intercept) {
  z = p.beta0*rn.intercept + p.theta*trt + p.delta*trt*cov       # linear combination with a bias
  pr = 1/(1+exp(-z))         # pass through an inv-logit function
  y = rbinom(N,1,pr)      # bernoulli response variable
}


### Survival outcome
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
  Tlat <- (- log(v) / (lambda * exp(cov * trt * s.delta + trt*s.theta + s.beta0*rn.intercept)))^(1 / rho)
  
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



# CART method

cart <- function(n1,y.t, y.c, y.type,
                 biom1.x1.t, biom1.x1.c,
                 biom1.x2.t, biom1.x2.c,
                 biom1.x3.t, biom1.x3.c,
                 biom1.x4.t, biom1.x4.c,
                 biom1.x1, biom1.x2, biom1.x3, biom1.x4,
                 y2.t, y2.c, y2.type, 
                 y3.t, y3.c, y3.type
                 ) {
  
  if (y.type == "binomial") {
    y.fit <- c(y.t, y.c)
    # set.seed(s+mn*S)
    fit <- rpart(y.fit~biom1.x1+biom1.x2+biom1.x3+biom1.x4, method="class", control=rpart.control(minsplit=prev.cutoff*n1))
    # set.seed(s+mn*S)
  } else if (y.type == "cox") {
    y.t.biom <-  cbind(y.t,biom1.x1.t,biom1.x2.t,biom1.x3.t,biom1.x4.t)
    colnames(y.t.biom)[6:9] <- c("biom1.x1" , "biom1.x2" , "biom1.x3" , "biom1.x4")
    y.c.biom <-  cbind(y.c,biom1.x1.c,biom1.x2.c,biom1.x3.c,biom1.x4.c)
    colnames(y.c.biom)[6:9] <- c("biom1.x1" , "biom1.x2" , "biom1.x3" , "biom1.x4")
    y <- rbind(y.t.biom, y.c.biom)
    # set.seed(s+mn*S)
    fit <- rpart(formula = Surv(time, status) ~ biom1.x1+biom1.x2+biom1.x3+biom1.x4, data = y, control=rpart.control(minsplit=prev.cutoff*n1))
  }else {
    y.fit <- c(y.t, y.c)
    # set.seed(s+mn*S)
    fit <- rpart(y.fit~biom1.x1+biom1.x2+biom1.x3+biom1.x4, method="anova", control=rpart.control(minsplit=prev.cutoff*n1))
  }

  
  # keep top 2-4 splits which minimize cross-validation error
  cp <- subset(data.frame(fit$cptable), nsplit >= 2 & nsplit <= 8)[which.min(subset(data.frame(fit$cptable), nsplit >= 2 & nsplit <= 8)[,"xerror"]),"CP"]
  # prune the tree 
  pfit<- prune(fit, cp=cp) # from cptable   
  
  # plot(pfit, uniform=TRUE)
  # text(pfit, use.n=TRUE, all=TRUE, cex=.8)
  # printcp(pfit) # display the results 
  
  rule_df.0 <- rpart.rules.table(pfit) %>%
    filter(Leaf==TRUE) %>%
    group_by(Rule) %>%
    summarise(Subrules = paste(Subrule, collapse=" & "))
  
  # if no split, then choose all space as the subgroup space
  subrule_df.0 <- try(rpart.subrules.table(pfit))
  
  if("try-error" %in% class(subrule_df.0)) {
    rule_df <- rule_df.0
    rule_df$Subrules <- "biom1.x1 >= 0 & biom1.x2 >= 0 & biom1.x3 >=0 & biom1.x4 >=0"
    loop <-NA
  } else {
    subrule_df.0$Less <- round(suppressWarnings(as.numeric(paste(subrule_df.0$Less))), digits = 2)
    subrule_df.0$Greater <- round(suppressWarnings(as.numeric(paste(subrule_df.0$Greater))), digits = 2)
    k.0 <- round(as.numeric(na.omit(subrule_df.0$Less)), digits = 2)
    
    k.step <- 0.05
    k.grid <- seq(0,1,k.step)
    
    subrule_df.prev <- subrule_df.0
    k <- k.0
    
    loop <- 1
    
    repeat {
      if (loop > 1) {
        prev.k <- k # from last iteration
      }
      
      for (m in 1:length(k)) {
        
        p.a <- matrix(NA, nrow = length(k.grid), ncol = length(rule_df.0$Rule)+2)
        
        for (n in 1:length(k.grid)) {
          
          rule_df <- rule_df.0
          subrule_df <- subrule_df.prev
          
          if (is.na(subrule_df$Less[m])) {
            subrule_df$Less[m+length(k)] <- k.grid[n]
            subrule_df$Greater[m] <- k.grid[n]
          } else {
            subrule_df$Less[m] <- k.grid[n]
            subrule_df$Greater[m+length(k)] <- k.grid[n]
          }
          
          for (i in 1:length(subrule_df$Subrule)) {
            temp <- subrule_df[i,]
            a <- subrule_df$Subrule[i]
            b <- paste(temp$Variable, ifelse(is.na(temp$Less), paste(">=", temp$Greater), paste("<", temp$Less)))
            rule_df$Subrules <- gsub(a, b, rule_df$Subrules)
          }
          
          for (i in 1:length(rule_df$Subrules)) {
            sub.index <- which(ifelse(eval(parse(text = rule_df$Subrules[i])),1,0) == 1)  
            g.t <- sub.index[sub.index<=n1/2]
            g.c <- sub.index[sub.index>n1/2]-n1/2
            
            if((length(g.t)>=2) & (length(g.c)>=2) & ((length(g.t) + length(g.c))/n1 >= prev.cutoff)) {
              if (y.type=="binomial") {
                p.a[n,i] <- prop.test(x=c(sum(y.t[g.t]), sum(y.c[g.c])), n=c(length(y.t[g.t]), length(y.c[g.c])))$p.v
              } else if (y.type=="cox") {
                p.a[n,i] <- pchisq(survdiff(formula = Surv(time, status) ~ trt, data = y[sub.index, ])$chisq, 1, lower.tail=FALSE)
              } else {
                p.a[n,i] <- t.test(y.t[g.t],y.c[g.c])$p.v       
              }
            } else {p.a[n,i] <-99}
          }
          
          p.a[n,length(rule_df.0$Rule)+1] <- which.min(p.a[n,1:length(rule_df.0$Rule)])
          p.a[n,length(rule_df.0$Rule)+2] <- min(p.a[n,1:length(rule_df.0$Rule)])
          
        }
        
        if (is.na(subrule_df.prev$Less[m])) {
          subrule_df.prev$Less[m+length(k)] <-  k.grid[which.min(p.a[,length(rule_df.0$Rule)+2])]
          subrule_df.prev$Greater[m] <-  k.grid[which.min(p.a[,length(rule_df.0$Rule)+2])]
        } else {
          subrule_df.prev$Less[m] <-  k.grid[which.min(p.a[,length(rule_df.0$Rule)+2])]
          subrule_df.prev$Greater[m+length(k)] <-  k.grid[which.min(p.a[,length(rule_df.0$Rule)+2])]
        }
        
      }
      
      k <- round(as.numeric(na.omit(subrule_df.prev$Less)), digits = 2)
      
      if (loop > 1) {
        if (sum(prev.k == k) == length(k)) {
          break # break if repeat jump
        } else if (loop > 10) {
          break # force break
        }
      }
      
      loop <- loop + 1
      
    }
    
  }
  
  p.we <- NA
  
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
  
  rule_df.list <- list()
  p.we.list <- rep(NA,dim(vector.matrix)[1])
  
  w.length <- dim(vector.matrix)[1]
  
  for (w in 1:w.length) {
    
    cat("w:",w,"\n")
    
    weight.vec <- vector.matrix[w,]
    
    for (i in 1:length(rule_df$Subrules)) {
      sub.index <- which(ifelse(eval(parse(text = rule_df$Subrules[i])),1,0) == 1)  
      g.t <- sub.index[sub.index<=n1/2]
      g.c <- sub.index[sub.index>n1/2]-n1/2
      
      if((length(g.t)>=2) & (length(g.c)>=2) & ((length(g.t) + length(g.c))/n1 >= prev.cutoff)) {
        
        
        if (y.type=="binomial") {
          # p1 <- prop.test(x=c(sum(y.t[g.t]), sum(y.c[g.c])), n=c(length(y.t[g.t]), length(y.c[g.c])))$p.v
          p1 <- z.prop(sum(y.t[g.t]), sum(y.c[g.c]), length(y.t[g.t]), length(y.c[g.c]))
        } else if (y.type=="cox") {
          # p1 <- pchisq(survdiff(formula = Surv(time, status) ~ trt, data = y[sub.index, ])$chisq, 1, lower.tail=FALSE)
          model <- coxph(formula = Surv(time, status) ~ trt, data = y[sub.index, ])
          if ('try-error' %in% class(model)) next  
          p1 <- coef(model)/sqrt(diag(vcov(model)))
        } else {
          # p1 <- t.test(y.t[g.t],y.c[g.c])$p.v    
          p1 <- t.test(y.t[g.t],y.c[g.c])$statistic  
        }
        
        if (y2.type=="binomial") {
          # p2 <- prop.test(x=c(sum(y2.t[g.t]), sum(y2.c[g.c])), n=c(length(y2.t[g.t]), length(y2.c[g.c])))$p.v
          p2 <- z.prop(sum(y2.t[g.t]), sum(y2.c[g.c]), length(y2.t[g.t]), length(y2.c[g.c]))
        } else if (y2.type=="cox") {
          
          y2.t.biom <-  cbind(y2.t,biom1.x1.t,biom1.x2.t,biom1.x3.t,biom1.x4.t)
          colnames(y2.t.biom)[6:9] <- c("biom1.x1" , "biom1.x2" , "biom1.x3" , "biom1.x4")
          y2.c.biom <-  cbind(y2.c,biom1.x1.c,biom1.x2.c,biom1.x3.c,biom1.x4.c)
          colnames(y2.c.biom)[6:9] <- c("biom1.x1" , "biom1.x2" , "biom1.x3" , "biom1.x4")
          y2 <- rbind(y2.t.biom, y2.c.biom)
          
          # p2 <- pchisq(survdiff(formula = Surv(time, status) ~ trt, data = y2[sub.index, ])$chisq, 1, lower.tail=FALSE)
          model <- coxph(formula = Surv(time, status) ~ trt, data = y2[sub.index, ])
          if ('try-error' %in% class(model)) next
          p2 <- coef(model)/sqrt(diag(vcov(model)))
        } else {
          # p2 <- t.test(y2.t[g.t],y2.c[g.c])$p.v  
          p2 <- t.test(y2.t[g.t],y2.c[g.c])$statistic 
        }
        
        if (y3.type == "binomial") {
          # p3 <- prop.test(x = c(sum(y3.t[g.t]), sum(y3.c[g.c])), n = c(length(y3.t[g.t]), length(y3.c[g.c])))$p.v
          p3 <- z.prop(sum(y3.t[g.t]), sum(y3.c[g.c]), length(y3.t[g.t]), length(y3.c[g.c]))
        } else if (y3.type == "cox") {
          
          y3.t.biom <-  cbind(y3.t,biom1.x1.t,biom1.x2.t,biom1.x3.t,biom1.x4.t)
          colnames(y3.t.biom)[6:9] <- c("biom1.x1" , "biom1.x2" , "biom1.x3" , "biom1.x4")
          y3.c.biom <-  cbind(y3.c,biom1.x1.c,biom1.x2.c,biom1.x3.c,biom1.x4.c)
          colnames(y3.c.biom)[6:9] <- c("biom1.x1" , "biom1.x2" , "biom1.x3" , "biom1.x4")
          y3 <- rbind(y3.t.biom, y3.c.biom)
          
          # p3 <- pchisq(survdiff(formula = Surv(time, status) ~ trt, data = y3[sub.index, ])$chisq, 1, lower.tail = FALSE)
          model <- coxph(formula = Surv(time, status) ~ trt, data = y3[sub.index, ])
          if ('try-error' %in% class(model)) next
          p3 <- coef(model)/sqrt(diag(vcov(model)))
        } else {
          # p3 <- t.test(y3.t[g.t], y3.c[g.c])$p.v
          p3 <- t.test(y3.t[g.t],y3.c[g.c])$statistic
        }
        
      } else {
        p1 <-NA
        p2 <-NA
        p3 <-NA
      } 
      
      p.we[i] <-  sum(c(p1,p2,p3)*weight.vec, na.rm = TRUE)
    }    
    
    rule_df.list[[w]] <- rule_df[which.max(abs(p.we)),]
    p.we.list[w] <- max(abs(p.we),na.rm = TRUE)
  }
  
  
  p.we.max <- max(abs(p.we.list), na.rm = TRUE)
  rule_df.max <- rule_df.list[[which.max(abs(p.we.list))]]
  weight.vex.max <- vector.matrix[which.max(abs(p.we.list)),]

  out <- list(rule_df.max, p.we.max, loop, weight.vex.max)
  return(out)
} 


# function to calculate test statistics, for subgroup
perc.U <- function (n.obs, cv.t.index, cv.c.index, y1.t, y1.c, y2.t, y2.c, y3.t, y3.c, rule_df, type) {
  ind.sub<-which(ifelse(eval(parse(text = rule_df$Subrules)) & ((1:(n.obs*2)) %in% c(cv.t.index, cv.c.index + n.obs)),1,0) == 1)
  ind.sub<-ind.sub[ind.sub>0]
  
  ind.t<-ind.sub[ind.sub<=(n.obs)]
  ind.c<-ind.sub[ind.sub>(n.obs)]-(n.obs)
  if ((length(ind.t)>=2) & (length(ind.c)>=2) 
      # & ((length(ind.t) + length(ind.c))/n.obs >= prev.cutoff)
      ) {
    
    if (outcome[1] == 1) {
      p.1 <- t.test(y1.t[ind.t],y1.c[ind.c])$p.v   
    } else if (outcome[1] == 2) {
      p.1 <- prop.test(x=c(sum(y1.t[ind.t]), sum(y1.c[ind.c])), n=c(length(y1.t[ind.t]), length(y1.c[ind.c])))$p.v
    } else {
      y1 <- rbind(y1.t, y1.c)
      p.1 <-  pchisq(survdiff(formula = Surv(time, status) ~ trt, data = y1[ind.sub, ])$chisq, 1, lower.tail=FALSE)
    }
    if (outcome[2] == 1) {
      p.2 <- t.test(y2.t[ind.t],y2.c[ind.c])$p.v   
    } else if (outcome[2] == 2) {
      p.2 <- prop.test(x=c(sum(y2.t[ind.t]), sum(y2.c[ind.c])), n=c(length(y2.t[ind.t]), length(y2.c[ind.c])))$p.v
    } else {
      y2 <- rbind(y2.t, y2.c)
      p.2 <-  pchisq(survdiff(formula = Surv(time, status) ~ trt, data = y2[ind.sub, ])$chisq, 1, lower.tail=FALSE)
    }
    if (outcome[3] == 1) {
      p.3 <- t.test(y3.t[ind.t],y3.c[ind.c])$p.v 
    } else if (outcome[3] == 2) {
      p.3 <- prop.test(x=c(sum(y3.t[ind.t]), sum(y3.c[ind.c])), n=c(length(y3.t[ind.t]), length(y3.c[ind.c])))$p.v
    } else {
      y3 <- rbind(y3.t, y3.c)
      p.3 <- pchisq(survdiff(formula = Surv(time, status) ~ trt, data = y3[ind.sub, ])$chisq, 1, lower.tail=FALSE)
    }
    p.hoch <- min(p.adjust(c(p.1,p.2,p.3), method = "hochberg"), na.rm = TRUE) # hochberg-adjusted p value
    # p value is two sided, and we keep the effect direction 
    if (type == 1) {
      if (outcome[1] == 1) {
        effect = sign(t.test(y1.t[ind.t],y1.c[ind.c])$statistic)
      } else if (outcome[1] == 2) {
        effect = sign(z.prop(sum(y1.t[ind.t]), sum(y1.c[ind.c]), length(y1.t[ind.t]), length(y1.c[ind.c])))
      } else {
        model <- coxph(formula = Surv(time, status) ~ trt, data = y1[ind.sub,])
        effect = sign(coef(model)/sqrt(diag(vcov(model))))    
      }
    } else if (type == 2) {
      if (outcome[2] == 1) {
        effect = sign(t.test(y2.t[ind.t],y2.c[ind.c])$statistic)
      } else if (outcome[2] == 2) {
        effect = sign(z.prop(sum(y2.t[ind.t]), sum(y2.c[ind.c]), length(y2.t[ind.t]), length(y2.c[ind.c])))
      } else {
        model <- coxph(formula = Surv(time, status) ~ trt, data = y2[ind.sub,])
        effect = sign(coef(model)/sqrt(diag(vcov(model))))    
      }
    } else {
      if (outcome[3] == 1) {
        effect = sign(t.test(y3.t[ind.t],y3.c[ind.c])$statistic)
      } else if (outcome[3] == 2) {
        effect = sign(z.prop(sum(y3.t[ind.t]), sum(y3.c[ind.c]), length(y3.t[ind.t]), length(y3.c[ind.c])))
      } else {
        model <- coxph(formula = Surv(time, status) ~ trt, data = y3[ind.sub,])
        effect = sign(coef(model)/sqrt(diag(vcov(model))))    
      }   
    }
    if (effect >=0) {
      z.hoch <- qnorm(1-p.hoch/2)
    } else {
      z.hoch <- qnorm(p.hoch/2)    
    }
    if ('try-error' %in% class(model)) next  
  } else {
    p.1 <- p.2 <- p.3 <- p.hoch <- z.hoch <- NA
  }
  out<-list(z.hoch, p.hoch)
  return(out)
}


# function to calculate %U, for subgroup
perc.U.test <- function (rule_df, type) {
  ind.sub<-which(ifelse(eval(parse(text = rule_df$Subrules)),1,0) == 1)  
  ind.sub<-ind.sub[ind.sub>0]
  ind.t<-ind.sub[ind.sub<=n.obs]
  ind.c<-ind.sub[ind.sub>n.obs]-n.obs
  ind.true <- which(X1.final.t>c01 & X2.final.t>c02 & X3.final.t + X4.final.t > c03 + c04)
  
  if ((length(ind.t)>=2) & (length(ind.c)>=2)) {
    if (type == "binomial") {
      pi.sub<-(length(ind.sub))/(2*n.obs)
      pi.overlap <- length(intersect(ind.t, ind.true))/length(ind.true)
      perc.U<-z.prop(sum(b.y.final.t[ind.t]), sum(b.y.final.c[ind.c]), length(ind.t), length(ind.c))/(sqrt(2*n.obs/2)*b.U.best[s])*100 
    } else if (type == "cox") {
      y.final <- rbind(s.y.final.t, s.y.final.c)
      pi.sub<-(length(ind.sub))/(2*n.obs)
      pi.overlap <- length(intersect(ind.t, ind.true))/length(ind.true)
      model <- coxph(formula = Surv(time, status) ~ trt, data = y.final[ind.sub, ])
      perc.U<- (coef(model)/sqrt(diag(vcov(model))))/(sqrt(2*n.obs/2)*s.U.best[s])*100 
    } else {
      pi.sub<-(length(ind.sub))/(2*n.obs)
      pi.overlap <- length(intersect(ind.t, ind.true))/length(ind.true)
      perc.U<-t.test(n.y.final.t[ind.t],n.y.final.c[ind.c])$statistic/(sqrt(2*n.obs/2)*n.U.best[s])*100 
    }
  } else {
    pi.sub<- NA 
    pi.overlap <- NA
    perc.U<- NA
  }
  out<-list(pi.sub,perc.U,ind.sub,pi.overlap)
  return(out)
}





# cross-validation for testing
lm.cv <- function(n1,cv.t.index,cv.c.index,y1.t,y1.c,y2.t,y2.c,y3.t,y3.c,
                  b1.x1.t, b1.x1.c,
                  b1.x2.t, b1.x2.c,
                  b1.x3.t, b1.x3.c,
                  b1.x4.t, b1.x4.c) {
  
  # initialize 4 biomarkers
  temp.b1.x1.t=b1.x1.t[c(cv.t.index)]
  temp.b1.x2.t=b1.x2.t[c(cv.t.index)]
  temp.b1.x1.c=b1.x1.c[c(cv.c.index)]
  temp.b1.x2.c=b1.x2.c[c(cv.c.index)]
  temp.b1.x3.t=b1.x3.t[c(cv.t.index)]
  temp.b1.x3.c=b1.x3.c[c(cv.c.index)]
  temp.b1.x4.t=b1.x4.t[c(cv.t.index)]
  temp.b1.x4.c=b1.x4.c[c(cv.c.index)]
  temp.b1.x1 <- c(temp.b1.x1.t,  temp.b1.x1.c)
  temp.b1.x2 <- c(temp.b1.x2.t,  temp.b1.x2.c)
  temp.b1.x3 <- c(temp.b1.x3.t,  temp.b1.x3.c)
  temp.b1.x4 <- c(temp.b1.x4.t,  temp.b1.x4.c)
  
  out.1 <- cart(n1, y1.t, y1.c, y1.type,
                temp.b1.x1.t, temp.b1.x1.c,
                temp.b1.x2.t, temp.b1.x2.c,
                temp.b1.x3.t, temp.b1.x3.c,
                temp.b1.x4.t, temp.b1.x4.c,
                temp.b1.x1, temp.b1.x2, temp.b1.x3, temp.b1.x4,
                y2.t, y2.c, y2.type,
                y3.t, y3.c, y3.type
                )
  out.2 <- cart(n1, y2.t, y2.c, y2.type,
                temp.b1.x1.t, temp.b1.x1.c,
                temp.b1.x2.t, temp.b1.x2.c,
                temp.b1.x3.t, temp.b1.x3.c,
                temp.b1.x4.t, temp.b1.x4.c,
                temp.b1.x1, temp.b1.x2, temp.b1.x3, temp.b1.x4,
                y1.t, y1.c, y1.type,
                y3.t, y3.c, y3.type
                )
  out.3 <- cart(n1, y3.t, y3.c, y3.type,
                temp.b1.x1.t, temp.b1.x1.c,
                temp.b1.x2.t, temp.b1.x2.c,
                temp.b1.x3.t, temp.b1.x3.c,
                temp.b1.x4.t, temp.b1.x4.c,
                temp.b1.x1, temp.b1.x2, temp.b1.x3, temp.b1.x4,
                y1.t, y1.c, y1.type,
                y2.t, y2.c, y2.type
                )
  
  
  rule.1 <- rbindlist(out.1[1], fill = TRUE)
  p.1.a <-  unlist(out.1[2])
  loop.1 <- unlist(out.1[3])
  rule.2 <- rbindlist(out.2[1], fill = TRUE)
  p.2.a <- unlist(out.2[2])
  loop.2 <- unlist(out.2[3])
  rule.3 <- rbindlist(out.3[1], fill = TRUE)
  p.3.a <- unlist(out.3[2])
  loop.3 <- unlist(out.3[3])
  
  # maximize utility "weighted"
  w.3lm.star.ind.a = which.max(c(p.1.a, p.2.a, p.3.a)) 
  
  # minimize p value "best-outcome"
  # w.3lm.star.ind.a = which.min(c(p.1.a, p.2.a, p.3.a)) 
  
  
  if (w.3lm.star.ind.a <= length(rule.1$Subrules)) {
    w.3lm.star.ind.a  <- w.3lm.star.ind.a 
    ind.3lm.a <- which(ifelse(eval(parse(text = rule.1$Subrules[w.3lm.star.ind.a])),1,0) == 1)  
    rule.best <- rule.1[w.3lm.star.ind.a,]
    out.a <- y1.type
    y.a <- 1
    weight.best <- out.1[[4]]
  } else if (length(rule.1$Subrules) < w.3lm.star.ind.a  & w.3lm.star.ind.a  <= length(rule.1$Subrules) + length(rule.2$Subrules)) {
    w.3lm.star.ind.a  <- w.3lm.star.ind.a  - length(rule.1$Subrules)
    ind.3lm.a <- which(ifelse(eval(parse(text = rule.2$Subrules[w.3lm.star.ind.a])),1,0) == 1)  
    rule.best <- rule.2[w.3lm.star.ind.a,]
    out.a <- y2.type 
    y.a <- 2
    weight.best <- out.2[[4]]
  } else {
    w.3lm.star.ind.a  <- w.3lm.star.ind.a  - (length(rule.1$Subrules) + length(rule.2$Subrules))
    ind.3lm.a <- which(ifelse(eval(parse(text = rule.3$Subrules[w.3lm.star.ind.a])),1,0) == 1)  
    rule.best <- rule.3[w.3lm.star.ind.a,]
    out.a <- y3.type
    y.a <- 3
    weight.best <- out.3[[4]]
  }
  
  # output
  out <- list(rule.best, out.a, y.a,
              weight.best)
  return(out)
  
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
  l <- apply(l.coeff, 1, function(x) (sum(all(which(x != 0) %in% c.true))))
  pred.x.sub <- sum(l)
  out <- list(pred.x.exact,pred.x.sub)
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
      
      if (m == 0) {
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
      
      if (m == 1) {
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
      
      if (m == 2) {
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
      
      if (m == 3) {
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
      
      
      if (m == 4) {
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
      
      if (m==5) {
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
      z.a.raw <- z.a.raw.1 <- z.a.raw.2 <- rep(NA,S)
      z.a.all <- rep(NA,S)
      z.a.hoch <- rep(NA,S)
      p.a.permute <- rep(NA,S)
      
      U.all <- U.best <- U.perc <- rep(NA,S) # estimate true test statistics, generic
      n.U.all <- n.U.best <- n.U.perc <- rep(NA,S) # estimate true test statistics, continuous
      b.U.all <- b.U.best <- b.U.perc <- rep(NA,S) # estimate true test statistics, binary
      s.U.all <- s.U.best <- s.U.perc <- rep(NA,S) # estimate true test statistics, survival
      
      perc.U.cv1.a <- perc.U.cv2.a <- perc.U.cv1.p <- perc.U.cv2.p <- rep(NA,S)
      pi.sub.cv1.a <- pi.sub.cv2.a <- pi.sub.cv1.p <- pi.sub.cv2.p <- rep(NA,S)
      
      coeff.cv1.a.1 <- coeff.cv1.a.2 <- coeff.cv1.a.3 <- coeff.cv1.a.4 <- rep(NA,S) 
      coeff.cv2.a.1 <- coeff.cv2.a.2 <- coeff.cv2.a.3 <- coeff.cv2.a.4 <- rep(NA,S) 
      
      perc.x.cv1.a.ex<- perc.x.cv2.a.ex<- rep(NA,S) 
      perc.x.cv1.a.sub<- perc.x.cv2.a.sub<- rep(NA,S) 

      weight.vec.out <- list()

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
        t.biom1.x1.x2 <- c(biom1.x1.t,  biom1.x1.c)*c(biom1.x2.t,  biom1.x2.c)*t
        t.biom1.x1.x3 <- c(biom1.x1.t,  biom1.x1.c)*c(biom1.x3.t,  biom1.x3.c)*t
        t.biom1.x1.x4 <- c(biom1.x1.t,  biom1.x1.c)*c(biom1.x4.t,  biom1.x4.c)*t
        t.biom1.x2.x3 <- c(biom1.x2.t,  biom1.x2.c)*c(biom1.x3.t,  biom1.x3.c)*t
        t.biom1.x2.x4 <- c(biom1.x2.t,  biom1.x2.c)*c(biom1.x4.t,  biom1.x4.c)*t
        t.biom1.x3.x4 <- c(biom1.x3.t,  biom1.x3.c)*c(biom1.x4.t,  biom1.x4.c)*t
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
        
        # generate random intercepts for multiple outcomes
        rn.intercept <- data.in[[9]]
        
        mu.t<-mu.c<-rep(NA,n1/2)  
        for(i in 1:(n1/2)){
          mu.t[i]<- beta0 + theta+(delta*ifelse(biom1.x1.t[i]>c01 & biom1.x2.t[i]>c02 & biom1.x3.t[i] + biom1.x4.t[i] > c03 + c04, 1,0))
          mu.c[i]<- beta0 # set means in control group to 0
          mu =c(mu.t[i],mu.c[i],mu.t[i],mu.c[i],mu.t[i],mu.c[i])
          if (type == 1) {
            mu =c(mu.t[i],mu.c[i],mu.c[i],mu.c[i],mu.c[i],mu.c[i])
          }
          resp=mvrnorm(1,mu,sig)
          y1.t[i]<-resp[1]
          y1.c[i]<-resp[2]
          y2.t[i]<-resp[3]
          y2.c[i]<-resp[4]
          y3.t[i]<-resp[5]
          y3.c[i]<-resp[6]
          
          y1.full <- c(y1.t, y1.c)
          y2.full <- c(y2.t, y2.c)
          y3.full <- c(y3.t, y3.c)
        }
        
        
        # generate non-normal outcomes - binary, survival
        
        # binary outcome
        if (outcome[1] == 2) {
          y1.t <- simulbin(N = n1/2, cov = ifelse(biom1.x1.t>c01 & biom1.x2.t>c02 & biom1.x3.t + biom1.x4.t > c03 + c04,1,0), trt = 1, rn.intercept = rn.intercept[,1] )
          y1.c <- simulbin(N = n1/2, cov = ifelse(biom1.x1.c>c01 & biom1.x2.c>c02 & biom1.x3.c + biom1.x4.c > c03 + c04,1,0), trt = 0, rn.intercept = rn.intercept[,2] )
          y1.full <- c(y1.t, y1.c)
        }
        
        if (outcome[2] == 2) {
          y2.t <- simulbin(N = n1/2, cov = ifelse(biom1.x1.t>c01 & biom1.x2.t>c02 & biom1.x3.t + biom1.x4.t > c03 + c04,1,0), trt = 1, rn.intercept = rn.intercept[,3] )
          y2.c <- simulbin(N = n1/2, cov = ifelse(biom1.x1.c>c01 & biom1.x2.c>c02 & biom1.x3.c + biom1.x4.c > c03 + c04,1,0), trt = 0, rn.intercept = rn.intercept[,4] )
          y2.full <- c(y2.t, y2.c)
        }
        
        if (outcome[3] == 2) {
          y3.t <- simulbin(N = n1/2, cov = ifelse(biom1.x1.t>c01 & biom1.x2.t>c02 & biom1.x3.t + biom1.x4.t > c03 + c04,1,0), trt = 1, rn.intercept = rn.intercept[,5] )
          y3.c <- simulbin(N = n1/2, cov = ifelse(biom1.x1.c>c01 & biom1.x2.c>c02 & biom1.x3.c + biom1.x4.c > c03 + c04,1,0), trt = 0, rn.intercept = rn.intercept[,6] )
          y3.full <- c(y3.t, y3.c)
        }
        
        
        
        # survival outcome
        if (outcome[1] == 3) {
          y1.t <- simulWeib(N=n1/2, lambda=0.001, cov = ifelse(biom1.x1.t>c01 & biom1.x2.t>c02 & biom1.x3.t + biom1.x4.t > c03 + c04,1,0), 
                            trt = 1, rho=1, rateC=t.rateC, rn.intercept = rn.intercept[,1])
          y1.c <- simulWeib(N=n1/2, lambda=0.001, cov = ifelse(biom1.x1.c>c01 & biom1.x2.c>c02 & biom1.x3.c + biom1.x4.c > c03 + c04,1,0), 
                            trt = 0, rho=1, rateC=c.rateC, rn.intercept = rn.intercept[,2])
          y1.full <- rbind(y1.t, y1.c)
        }
        
        
        if (outcome[2] == 3) {
          y2.t <- simulWeib(N=n1/2, lambda=0.001, cov = ifelse(biom1.x1.t>c01 & biom1.x2.t>c02 & biom1.x3.t + biom1.x4.t > c03 + c04,1,0), 
                            trt = 1, rho=1, rateC=t.rateC, rn.intercept = rn.intercept[,3])
          y2.c <- simulWeib(N=n1/2, lambda=0.001, cov = ifelse(biom1.x1.c>c01 & biom1.x2.c>c02 & biom1.x3.c + biom1.x4.c > c03 + c04,1,0), 
                            trt = 0, rho=1, rateC=c.rateC, rn.intercept = rn.intercept[,4])
          y2.full <- rbind(y2.t, y2.c)
        }
        
        if (outcome[3] == 3) {
          y3.t <- simulWeib(N=n1/2, lambda=0.001, cov = ifelse(biom1.x1.t>c01 & biom1.x2.t>c02 & biom1.x3.t + biom1.x4.t > c03 + c04,1,0), 
                            trt = 1, rho=1, rateC=t.rateC, rn.intercept = rn.intercept[,5])
          y3.c <- simulWeib(N=n1/2, lambda=0.001, cov = ifelse(biom1.x1.c>c01 & biom1.x2.c>c02 & biom1.x3.c + biom1.x4.c > c03 + c04,1,0), 
                            trt = 0, rho=1, rateC=c.rateC, rn.intercept = rn.intercept[,6])
          y3.full <- rbind(y3.t, y3.c)
        }
        
        
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
        
        
        if (corr.type == 0) {
          rn.intercept.final <- matrix(1, nrow = n.obs, ncol = 6)                                  
        } else if (corr.type == 1) {
          rn.intercept.final <- mvrnorm(n.obs,rep(0,6),sig)                                         
        }
        
        n.y.final.t=rnorm(n.obs,beta0.ext*rn.intercept.final[,1]+theta.ext+(delta.ext*ifelse(X1.final.t>c01 & X2.final.t>c02 & X3.final.t + X4.final.t > c03 + c04, 1,0)),sig2)
        n.y.final.c=rnorm(n.obs,beta0.ext*rn.intercept.final[,2],sig2)
        n.U.all[s] <- t.test(n.y.final.t, n.y.final.c)$statistic/sqrt(2*n.obs/2)
        n.U.best[s] <- t.test(n.y.final.t[which(X1.final.t>c01 & X2.final.t>c02 & X3.final.t + X4.final.t > c03 + c04)],
                              n.y.final.c[which(X1.final.c>c01 & X2.final.c>c02 & X3.final.c + X4.final.c > c03 + c04)])$statistic/sqrt(2*n.obs/2)
        n.U.perc[s] <- n.U.all[s]/n.U.best[s]*100
        
        b.y.final.t <- simulbin(N = n.obs, cov = ifelse(X1.final.t>c01 & X2.final.t>c02 & X3.final.t + X4.final.t > c03 + c04,1,0), trt = 1, rn.intercept = rn.intercept.final[,3])
        b.y.final.c <- simulbin(N = n.obs, cov = ifelse(X1.final.c>c01 & X2.final.c>c02 & X3.final.c + X4.final.c > c03 + c04,1,0), trt = 0, rn.intercept = rn.intercept.final[,4])
        b.U.all[s] <- z.prop(sum(b.y.final.t), sum(b.y.final.c), length(b.y.final.t), length(b.y.final.c))/sqrt(2*n.obs/2)
        b.U.best[s] <- z.prop(sum(b.y.final.t[which(X1.final.t>c01 & X2.final.t>c02 & X3.final.t + X4.final.t > c03 + c04)]), 
                              sum(b.y.final.c[which(X1.final.c>c01 & X2.final.c>c02 & X3.final.c + X4.final.c > c03 + c04)]), 
                              length(b.y.final.t[which(X1.final.t>c01 & X2.final.t>c02 & X3.final.t + X4.final.t > c03 + c04)]), 
                              length(b.y.final.c[which(X1.final.c>c01 & X2.final.c>c02 & X3.final.c + X4.final.c > c03 + c04)]))/sqrt(2*n.obs/2)
        b.U.perc[s] <- b.U.all[s]/b.U.best[s]*100
        
        s.y.final.t <- simulWeib(N=n.obs, lambda=0.001, cov = ifelse(X1.final.t>c01 & X2.final.t>c02 & X3.final.t + X4.final.t > c03 + c04,1,0), 
                                 t = 1, rho=1, rateC=t.rateC, rn.intercept = rn.intercept.final[,5])
        s.y.final.c <- simulWeib(N=n.obs, lambda=0.001, cov = ifelse(X1.final.c>c01 & X2.final.c>c02 & X3.final.c + X4.final.c > c03 + c04,1,0), 
                                 t = 0, rho=1, rateC=c.rateC, rn.intercept = rn.intercept.final[,6])
        s.y.final <- rbind(s.y.final.t, s.y.final.c)
        model <- coxph(formula = Surv(time, status) ~ trt, data = s.y.final)
        s.U.all[s] <- (coef(model)/sqrt(diag(vcov(model))))/sqrt(2*n.obs/2)
        model <- coxph(formula = Surv(time, status) ~ trt, data = s.y.final[which(s.y.final$x==1),])
        s.U.best[s] <- (coef(model)/sqrt(diag(vcov(model))))/sqrt(2*n.obs/2)
        s.U.perc[s] <- s.U.all[s]/s.U.best[s]*100
      

        
        #########
        # CART
        #########
      
        # initialize
        # normal outcome
        if (outcome[1] == 1) {
          y1.type = "gaussian"  
        }
        
        if (outcome[2] == 1) {
          y2.type = "gaussian"  
        }
        
        if (outcome[3] == 1) {
          y3.type = "gaussian" 
        }
        
        # binary outcome
        if (outcome[1] == 2) {
          y1.type = "binomial" 
        }
        
        if (outcome[2] == 2) {
          y2.type = "binomial" 
        }
        
        if (outcome[3] == 2) {
          y3.type = "binomial" 
        }
        
        # survival outcome
        if (outcome[1] == 3) {
          y1.type = "cox" 
        }
        
        if (outcome[2] == 3) {
          y2.type = "cox" 
        }
        
        if (outcome[3] == 3) {
          y3.type = "cox" 
        }
        
        ### apply to whole sample
        cv.result <- lm.cv(n1, 1:(n1/2), 1:(n1/2), y1.t, y1.c, y2.t, y2.c, y3.t, y3.c,
                           biom1.x1.t, biom1.x1.c,
                           biom1.x2.t, biom1.x2.c,
                           biom1.x3.t, biom1.x3.c,
                           biom1.x4.t, biom1.x4.c) 
        
        
        ## %U
        rule.sub <- function(rule.best) {
          rule.best$Subrules <- gsub("biom1.x1", "X1.final", rule.best$Subrules)
          rule.best$Subrules <- gsub("biom1.x2", "X2.final", rule.best$Subrules)
          rule.best$Subrules <- gsub("biom1.x3", "X3.final", rule.best$Subrules)
          rule.best$Subrules <- gsub("biom1.x4", "X4.final", rule.best$Subrules) 
          return(rule.best)
        }
        u1.a <- u2.a <- perc.U.test(rule.sub(cv.result[[1]]), cv.result[[2]])
        
        # main effects
        coeff.sub <- function(rule.best) {
          coeff.3lm.a.1 <- ifelse(grepl("X1.final",rule.best$Subrules),1,0)
          coeff.3lm.a.2 <- ifelse(grepl("X2.final",rule.best$Subrules),1,0)
          coeff.3lm.a.3 <- ifelse(grepl("X3.final",rule.best$Subrules),1,0)
          coeff.3lm.a.4 <- ifelse(grepl("X4.final",rule.best$Subrules),1,0)
          out <- list(coeff.3lm.a.1, coeff.3lm.a.2, coeff.3lm.a.3,  coeff.3lm.a.4)
        }
        coef.1.a <- coef.2.a <- coeff.sub(rule.sub(cv.result[[1]]))
        
        
      
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
          
          cv1.result <- lm.cv(n1/2, cv1.t.index, cv1.c.index, y1.t.cv1, y1.c.cv1, y2.t.cv1, y2.c.cv1, y3.t.cv1, y3.c.cv1,
                              biom1.x1.t, biom1.x1.c,
                              biom1.x2.t, biom1.x2.c,
                              biom1.x3.t, biom1.x3.c,
                              biom1.x4.t, biom1.x4.c)
          cv2.result <- lm.cv(n1/2, cv2.t.index, cv2.c.index, y1.t.cv2, y1.c.cv2, y2.t.cv2, y2.c.cv2, y3.t.cv2, y3.c.cv2,
                              biom1.x1.t, biom1.x1.c,
                              biom1.x2.t, biom1.x2.c,
                              biom1.x3.t, biom1.x3.c,
                              biom1.x4.t, biom1.x4.c)
          
          ## unweighted
          if (cv2.result[[3]] == 1) {
            y.cv1.t.a <- y1.t
            y.cv1.c.a <- y1.c
          } else if (cv2.result[[3]] == 2) {
            y.cv1.t.a <- y2.t
            y.cv1.c.a <- y2.c
          } else {
            y.cv1.t.a <- y3.t
            y.cv1.c.a <- y3.c
          }
          if (cv1.result[[3]] == 1) {
            y.cv2.t.a <- y1.t
            y.cv2.c.a <- y1.c
          } else if (cv1.result[[3]] == 2) {
            y.cv2.t.a <- y2.t
            y.cv2.c.a <- y2.c
          } else {
            y.cv2.t.a <- y3.t
            y.cv2.c.a <- y3.c
          }
          # Zs
          z1s.a <- perc.U(n1/2, cv1.t.index, cv1.c.index, 
                          y1.t, y1.c, y2.t, y2.c, y3.t, y3.c, 
                          cv2.result[[1]], cv2.result[[3]])[[1]]
          z2s.a <- perc.U(n1/2, cv2.t.index, cv2.c.index, 
                          y1.t, y1.c, y2.t, y2.c, y3.t, y3.c, 
                          cv1.result[[1]], cv1.result[[3]])[[1]]
          zs.a <- sqrt(0.5)*z1s.a+sqrt(0.5)*z2s.a 
          
          # ## %U
          # rule.sub <- function(rule.best) {
          #   rule.best$Subrules <- gsub("biom1.x1", "X1.final", rule.best$Subrules)
          #   rule.best$Subrules <- gsub("biom1.x2", "X2.final", rule.best$Subrules)
          #   rule.best$Subrules <- gsub("biom1.x3", "X3.final", rule.best$Subrules)
          #   rule.best$Subrules <- gsub("biom1.x4", "X4.final", rule.best$Subrules) 
          #   return(rule.best)
          # }
          # u2.a <- perc.U.test(rule.sub(cv2.result[[1]]), cv2.result[[2]])
          # u1.a <- perc.U.test(rule.sub(cv1.result[[1]]), cv1.result[[2]])
          # 
          # # main effects
          # coeff.sub <- function(rule.best) {
          #   coeff.3lm.a.1 <- ifelse(grepl("X1.final",rule.best$Subrules),1,0)
          #   coeff.3lm.a.2 <- ifelse(grepl("X2.final",rule.best$Subrules),1,0)
          #   coeff.3lm.a.3 <- ifelse(grepl("X3.final",rule.best$Subrules),1,0)
          #   coeff.3lm.a.4 <- ifelse(grepl("X4.final",rule.best$Subrules),1,0)
          #   out <- list(coeff.3lm.a.1, coeff.3lm.a.2, coeff.3lm.a.3,  coeff.3lm.a.4)
          # }
          # coef.1.a <- coeff.sub(rule.sub(cv1.result[[1]]))
          # coef.2.a <- coeff.sub(rule.sub(cv2.result[[1]]))
          
          z.raw <-  list(zs.a,
                      u1.a,u2.a,
                      coef.1.a, coef.2.a,
                      z1s.a,z2s.a)
      
      
        ### CART
        z.a.raw[s] <- z.raw[[1]]
        z.a.raw.1[s] <- z.raw[[6]]
        z.a.raw.2[s] <- z.raw[[7]]
        
        ## %U
        perc.U.cv1.a[s] <- z.raw[[2]][[2]]
        perc.U.cv2.a[s] <- z.raw[[3]][[2]]
        
        pi.sub.cv1.a[s] <- z.raw[[2]][[1]]
        pi.sub.cv2.a[s] <- z.raw[[3]][[1]]
        
        ## coeffcients
        coeff.cv1.a.1[s] <- z.raw[[4]][1]
        coeff.cv1.a.2[s] <- z.raw[[4]][2]
        coeff.cv1.a.3[s] <- z.raw[[4]][3]
        coeff.cv1.a.4[s] <- z.raw[[4]][4]
        
        coeff.cv2.a.1[s] <- z.raw[[5]][1]
        coeff.cv2.a.2[s] <- z.raw[[5]][2]
        coeff.cv2.a.3[s] <- z.raw[[5]][3]
        coeff.cv2.a.4[s] <- z.raw[[5]][4]
        
        ## variable selection
        # exact selection
        perc.x.cv1.a.ex[s] <- biom1.pred(coeff.cv1.a.1[s], coeff.cv1.a.2[s], coeff.cv1.a.3[s], coeff.cv1.a.4[s], c.true)[[1]]
        perc.x.cv2.a.ex[s] <- biom1.pred(coeff.cv2.a.1[s], coeff.cv2.a.2[s], coeff.cv2.a.3[s], coeff.cv2.a.4[s], c.true)[[1]]
        # subset selection
        perc.x.cv1.a.sub[s] <- biom1.pred(coeff.cv1.a.1[s], coeff.cv1.a.2[s], coeff.cv1.a.3[s], coeff.cv1.a.4[s], c.true)[[2]]
        perc.x.cv2.a.sub[s] <- biom1.pred(coeff.cv2.a.1[s], coeff.cv2.a.2[s], coeff.cv2.a.3[s], coeff.cv2.a.4[s], c.true)[[2]]
        
        
        weight.vec.out[[s]] <-   cv.result[[4]] 
        
        ### permutations
        z.a.permute.1 <- z.a.permute.2 <- z.a.permute <- rep(NA,P)
      

      
        
      }


      U.results <- rbind(cbind(median(n.U.all, na.rm = TRUE), median(n.U.best, na.rm = TRUE), median(n.U.perc, na.rm = TRUE)),
                         cbind(median(b.U.all, na.rm = TRUE), median(b.U.best, na.rm = TRUE), median(b.U.perc, na.rm = TRUE)),
                         cbind(median(s.U.all, na.rm = TRUE), median(s.U.best, na.rm = TRUE), median(s.U.perc, na.rm = TRUE)) )
      
      
      U.results
      
      out.final <- paste(model.out,"_", m.out ,"_",
                         paste(type.name,"_",corr.name,sep=""),
                         sep="")
      
      
      output.full <- list(cbind(z.a.raw.1, z.a.raw.2, z.a.raw,
                                p.a.permute,
                                perc.U.cv1.a, pi.sub.cv1.a, 
                                perc.U.cv2.a, pi.sub.cv2.a),
                          
                          cbind(coeff.cv1.a.1 , coeff.cv1.a.2 , coeff.cv1.a.3 , coeff.cv1.a.4 , 
                                coeff.cv2.a.1 , coeff.cv2.a.2 , coeff.cv2.a.3 , coeff.cv2.a.4),
                          weight.vec.out
      )
      write.list(output.full,  file = paste(key.dir, out.final, "_" ,mn,  ".csv", sep=""))

    }
  }
  
  
}




