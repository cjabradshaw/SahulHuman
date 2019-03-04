#####################################################################################
## global sensitivity analysis using Latin hypercube sampling for stochastic 
## ancient human founding population model for Sahul
##
## Corey Bradshaw
## September 2018 / updated February 2019
##
## testing effect on probability of quasi-extinction (< 25 females)
## of modifying the following variables:
##
## Bentley fertility total
## set in model: 'fert.bentley <- 4.69/2' = 2.345; test random uniform sample between 2.1105 and 2.5795 (10%)
##
## change assumed SD of survival
## set in model: Sx.sd <- 0.05; test random uniform sample between 0.025 and 0.10
##
## change initial infant mortality rate — Siler hazard model parameter a1
## set in model: a1 <- 0.422; test random uniform sample between 0.3 and 0.5
##
## change rate of mortality decline — Siler hazard model parameter a2
## set in model: a2 <- 0.013; test random uniform sample between 0.010 and 0.015
##
## change initial adult mortality rate — Siler hazard model parameter a3
## set in model: a3 <- 1.47e-04; test random uniform sample between 0.0001323 and 0.0001617 (10%)
##
## change age-independent mortality — Siler hazard model parameter b1
## set in model: b1 <- 1.131; test random uniform sample between 1.0 and 2.0
##
## change rate of mortality increase — Siler hazard model parameter b3
## set in model: b3 <- 0.086; test random uniform sample between 0.06 and 0.095
##
## change modifier to survival with > carrying capacity
## set in model: surv.min <- 0.98; test random uniform sample between 0.95 and 0.99
##
## change threshold of extinction probability
## set in model: min.thresh <- 50; test random uniform sample between 25 and 75 (/2 for females)
##
## change probability of catastrophe
## set in model to 'cat.numerator <- 0.14'; test random uniform sample between 0.1 and 0.2
##
## change intensity of catastrophe
## set in model to 'cat.intensity <- 0.5'; test random uniform sample between 0.25 and 0.75
##
## change low-density estimate for calculation of carrying capacity
## set in model 'start.K' <- 47000; test random uniform sample between 50% (23500) and 200% (94000)
#####################################################################################

## Remove everything
rm(list = ls())

## Load libraries
library(iterators)
library(snow)
library(doSNOW)
library(foreach)
library(lhs)
library(data.table)

## source (update when appropriate)
source("~/Documents/Papers/Other/Global human population/ms/PNAS/R1/matrixOperators.r")

## set working directory
setwd("~/Documents/Papers/Palaeo/Sahul/Aus human age-structured model/analysis/gsa")

## Set up parallel processing (nproc is the number of processing cores to use)
nproc <- 6
cl.tmp = makeCluster(rep('localhost',nproc), type='SOCK')
registerDoSNOW(cl.tmp)
getDoParWorkers()

## ancient human founding population simulation function
ahfp_sim <- function(input, dir.nm, rowNum) {
  
  ## assign all parameter values
  for (i in 1:ncol(input)) {assign(names(input)[i], input[,i])}
  
  ####################################################
  ## necessary input calculations for stochastic model
  ####################################################
  
  ## source (update when appropriate)
  source("matrixOperators.r") 
  
  ## define necessary functions
  # stochastic beta sampler for a survival vector (replaces above function for faster processing)
  stoch.surv.func <- function(mu, var) {
    Sx <- rbeta(length(mu), (((1 - mu) / var - 1 / mu) * mu ^ 2), ((((1 - mu) / var - 1 / mu) * mu ^ 2)*(1 / mu - 1)))
    return(params=Sx)
  }
  
  ## import necessary data (place in appropriate directory)
  npp.sah <- read.table("ClimateSahul_Npp.csv", header=T, sep=",") 
  dat.world13 <- read.table("world2013lifetable.csv", header=T, sep=",") #fertility data

  # net primary productivity (NPP; kg C/yr/m2)
  npp.nsah <- subset(npp.sah, Lat..degree. <= 0 & Lat..degree. >= -14)
  npp.nsah.med <- as.numeric(apply(npp.nsah[,-c(1,2)], MARGIN=2, quantile, probs=0.5))
  npp.nsah.lo <- as.numeric(apply(npp.nsah[,-c(1,2)], MARGIN=2, quantile, probs=0.25))
  npp.nsah.up <- as.numeric(apply(npp.nsah[,-c(1,2)], MARGIN=2, quantile, probs=0.75))
  npp.nsah.dat <- as.data.frame(t(data.frame(npp.nsah.lo,npp.nsah.med,npp.nsah.up)))
  rownames(npp.nsah.dat) <- c("lo","med","up")
  colnames(npp.nsah.dat) <- colnames(npp.nsah[,-c(1,2)])
  year.vec <- 120000:0
  yr1000 <- seq(120,0,-1)
  npp.med.interp <- approx(yr1000*1000, npp.nsah.dat[2,], xout=year.vec, method="linear")
  npp.lo.interp <- approx(yr1000*1000, npp.nsah.dat[1,], xout=year.vec, method="linear")
  npp.up.interp <- approx(yr1000*1000, npp.nsah.dat[3,], xout=year.vec, method="linear")
  K.dat <- data.frame(year.vec,npp.lo.interp$y,npp.med.interp$y,npp.up.interp$y)
  colnames(K.dat) <- c("year","npp.lo","npp.med","npp.up")
  K.dat$npp.lo.sc1 <- scale(K.dat$npp.lo + abs(min(K.dat$npp.lo)), scale=F, center=F)
  K.dat$npp.med.sc1 <- scale(K.dat$npp.med + abs(min(K.dat$npp.lo)), scale=F, center=F) 
  K.dat$npp.up.sc1 <- scale(K.dat$npp.up + abs(min(K.dat$npp.lo)), scale=F, center=F)
  K.dat$npp.lo.sc <- K.dat$npp.lo.sc1/max(K.dat$npp.med.sc1)
  K.dat$npp.med.sc <- K.dat$npp.med.sc1/max(K.dat$npp.med.sc1)
  K.dat$npp.up.sc <- K.dat$npp.up.sc1/max(K.dat$npp.med.sc1)
  #start.K <- 47000 # Guatney & Holliday 2015
  K.dat$K.lo <- round(start.K * 1/K.dat$npp.lo.sc, 0)
  K.dat$K.med <- (K.dat$npp.med.sc/K.dat$npp.lo.sc) * K.dat$K.lo
  K.dat$K.up <- (K.dat$npp.up.sc/K.dat$npp.lo.sc) * K.dat$K.lo

  # Siler hazard h(x) (Gurven et al. 2007)
  #a1 <- 0.422 # initial infant mortality rate (also known as αt)
  #b1 <- 1.131 # rate of mortality decline (also known as bt)
  #a2 <- 0.013 # age-independent mortality (exogenous mortality due to environment); also known as ct
  #a3 <- 1.47e-04 # initial adult mortality rate (also known as βt)
  #b3 <- 0.086 # rate of mortality increase
  longev <- 80
  x <- seq(0,longev,1) # age vector
  h.x <- a1 * exp(-b1*x) + a2 + a3 * exp(b3 * x) # Siler's hazard model
  l.x <- exp((-a1/b1) * (1 - exp(-b1*x))) * exp(-a2 * x) * exp(a3/b3 * (1 - exp(b3 * x))) # Siler's survival (proportion surviving) model
  l.inf <- exp(-a1/b1) # survival at infinite time
  T.m <- 1/b1 # time constant at which maturity is approached
  h.m <- a2 # hazard for mature animals
  l.m <- exp(-a2*x) # survival
  h.s <- a3*exp(b3*x) # hazard for senescence
  l.s <- exp((a3/b3)*(1 - exp(b3*x))) # survival for senescence
  f.x <- a3*exp(b3*x)*exp((a3/b3)/(1-exp(b3*x))) # probability density function
  T.s <- (1/b3) # modal survival time
  init.pop <- 10000
  lx <- round(init.pop*l.x,0)
  len.lx <- length(lx)
  dx <- lx[1:(len.lx-1)]-lx[2:len.lx]
  qx <- dx/lx[1:(length(lx)-1)]
  Sx <- 1 - qx
  sx <- lx[2:len.lx]/lx[1:(len.lx-1)]
  mx <- 1 - sx
  Lx <- (lx[1:(len.lx-1)] + lx[2:len.lx])/2
  ex <- rev(cumsum(rev(Lx)))/lx[-len.lx]
  ex.avg <- ex + x[-len.lx]

  # set SD for Sx
  #Sx.sd <- 0.05 # can set to any value
  
  # fertility (Walker et al. 2006)
  primiparity.walker <- c(17.7,18.7,19.5,18.5,18.5,18.7,25.7,19,20.5,18.8,17.8,18.6,22.2,17,16.2,18.4)
  prim.mean <- round(mean(primiparity.walker),0)
  prim.lo <- round(quantile(primiparity.walker,probs=0.025),0)
  prim.hi <- round(quantile(primiparity.walker,probs=0.975),0)
  fert.world13 <- dat.world13$m.f
  fert.trunc <- fert.world13[1:(longev+1)]
  pfert.trunc <- fert.trunc/sum(fert.trunc)
  #fert.bentley <- 4.69/2 # Bentley 1985 for !Kung
  fert.vec <- fert.bentley * pfert.trunc
  
  ## construct base matrices
  stages <- len.lx
  popmat <- matrix(0,nrow=stages,ncol=stages)
  colnames(popmat) <- x
  rownames(popmat) <- x
  popmat[1,] <- fert.vec
  diag(popmat[2:stages,]) <- Sx
  popmat[stages,stages] <- 0 # Sx[stages-1]
  popmat.orig <- popmat ## save original matrix
  
  ## initial population vector
  pop.found <- 700 # founding population size
  init.vec <- stable.stage.dist(popmat) * pop.found
  
  ## use natural-mortality popmat for set-up parameters
  gen.l <- G.val(popmat.orig, stages) # mean generation length
  #cat.numerator <- 0.14
  cat.pr <- cat.numerator/gen.l # probability of catastrophe (Reed et al. 2003)
  #cat.intensity <- 0.5
  
  #initial colonisation window
  col.old <- 60000
  col.yng <- 50000
  
  # generations to project
  gen.proj <- 100
  
  # set quasi-extinction threshold
  #min.thresh <- 50

  ## set time limit for projection in 1-yr increments
  yr.st <- round(runif(iter,col.yng,col.old), 0) 
  #************************
  yr.end <- yr.st - round(gen.proj*gen.l, 0) # set projection end date
  #************************
  t <- (yr.st - yr.end)
  
  ## set population storage matrices
  n.mat <- matrix(0,nrow=stages,ncol=(t[1]+1))
  n.mat[,1] <- init.vec
  popmat <- popmat.orig
  
  ## iterate projection
  #iter <- 1000
  itdiv <- iter/100
  
  # set storage matrices & vectors
  n.sums.mat <- matrix(data = 0, nrow = iter, ncol = (t[1]+1))

  for (e in 1:iter) {
    
    # K
    yr.st.sub <- which(K.dat$year == yr.st[e])
    npp.sc.st <- runif(1, K.dat$npp.lo.sc[yr.st.sub], K.dat$npp.up.sc[yr.st.sub])
    K.run <- round((npp.sc.st / K.dat$npp.med.sc[yr.st.sub]) * K.dat$K.med[yr.st.sub:(yr.st.sub+t[e])], 0) # this iteration's realised K series
    yr.vec.run <- yr.st[e]:yr.end[e]
    
    ## reset popmat
    popmat <- popmat.orig
    
    ## set up projection loop
    for (i in 1:t[e]) {
      
      ## reconstruct popmat with stochastic elements
      popmat[1, ] <- pfert.trunc * rnorm(1, fert.bentley,0.05*fert.bentley) # fertility sampler
      diag(popmat[2:(stages), ]) <- ((ifelse(sum(n.mat[,i], na.rm=T) > K.run[i], surv.min, 1)) * (stoch.surv.func(Sx, Sx.sd^2))[-stages]) # survival sampler
      
      # survival (+ catastrophic mortality at 50%)
      if ((rbinom(1, 1, cat.pr)) == 1) {
        diag(popmat[2:(stages), ]) <- (0.5* (diag(popmat[2:(stages), ])))}
      popmat[stages,stages] <- 0
      
      ## project over interval
      n.mat[,i+1] <- popmat %*% (n.mat[,i])
      
      #print(i)
    }
    
    n.sums.mat[e,] <- as.vector(colSums(n.mat))
    if (e %% itdiv==0) print(e) 
  }
  
  # N confidence limits
  n.mn <- apply(n.sums.mat, MARGIN=2, mean, na.rm=T) # mean over all iterations
  n.up <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  n.lo <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
  n.min <- apply(n.sums.mat, MARGIN=1, min, na.rm=T) # minimum over all projected years
  pr.ext <- pr.ext <- sum(ifelse(n.min < (min.thresh/2), 1, 0))/iter

  # save
  input$pr.ext <- pr.ext
  save.nm <- paste0('res',sprintf("%09.0f", rowNum))
  assign(save.nm, input)
  save(list=save.nm,file=paste(dir.nm,save.nm,sep='/'))

}

## parameter ranges
ranges <- list()
  ranges$fert.bentley <- c(2.1105, 2.5795)
  ranges$Sx.sd <- c(0.025, 0.100)
  ranges$a1 <- c(0.3, 0.5)
  ranges$a2 <- c(0.010, 0.015)
  ranges$a3 <- c(0.0001323, 0.0001617)
  ranges$b1 <- c(1.0, 2.0)
  ranges$b3 <- c(0.060, 0.095)
  ranges$surv.min <- c(0.95, 0.99)
  ranges$min.thresh <- c(25, 75)
  ranges$cat.numerator <- c(0.1, 0.2)
  ranges$cat.intensity <- c(0.25, 0.75)
  ranges$start.K <- c(23500, 94000)

## create hypercube
nSamples <- 1000
lh <- data.frame(randomLHS(n=nSamples, k=length(ranges)))
names(lh) <- names(ranges)

## convert parameters to required scale
for (i in 1:ncol(lh)) {
  par <- names(lh)[i]
  lh[,par] <- qunif(lh[,i], min=ranges[[par]][1], max=ranges[[par]][2]) ## continuous
}

## number of iterations for each parameter set
lh$iter <- 100

## folder for saving the results of each row
## we could just store in memory, but then if something breaks we will lose the lot
dir.nm <- 'test2'
dir.create(dir.nm)

## run in parallel
#res <- foreach(rowNum=1:nrow(lh),.verbose=F) %dopar% {ahfp_sim(input=lh[rowNum,],dir.nm=dir.nm,rowNum=rowNum)}
res <- foreach(rowNum=1:nrow(lh),.verbose=F) %do% {ahfp_sim(input=lh[rowNum,],dir.nm=dir.nm,rowNum=rowNum)}

## retrieve results
res.nms <- list.files(dir.nm)
res.list <- lapply(res.nms, function(x) {load(paste(dir.nm,x,sep='/')) ; print(x) ; return(eval(as.name(x)))})
dat <- rbindlist(res.list)
head(dat)
sum(is.na(dat$pr.ext))

## BRT
dat.nona <- na.omit(dat) # remove NAs
library(dismo)
library(gbm)
brt.fit <- gbm.step(dat.nona, gbm.x = attr(dat.nona, "names")[1:12], gbm.y = attr(dat.nona, "names")[14], family="gaussian", tolerance = 0.0001, learning.rate = 0.01, bag.fraction=0.75, tree.complexity = 2)
summary(brt.fit)
gbm.plot(brt.fit)
gbm.plot.fits(brt.fit)


## post-hoc relationships for quantifying magnitude of change in pr.ext with changes in explanatory variables
# set in model: b1 <- 1.131; test random uniform sample between 1.0 and 2.0
plot(dat$b1, dat$pr.ext, pch=19, xlab="age-independent mortality (b1)", ylab="probability of extinction")
#plot(dat$b1, logit(dat$pr.ext), pch=19, xlab="age-independent mortality (b1)", ylab="logit probability of extinction")
fit.b1 <- lm((dat$pr.ext) ~ dat$b1)
abline(fit.b1, lty=2, col="red")
# 50% higher b1 ...
b1 <- 1.131
b1.prime <- 1.5*b1
abline(v=b1.prime, lty=3)
pr.ext.b1.prime <- dat$pr.ext[which(dat$b1 >= 0.95*b1.prime & dat$b1 <= 1.05*b1.prime)]
mean.b1.prime <- mean(pr.ext.b1.prime)
mean.b1.prime
abline(v=b1, lty=3)
pr.ext.b1 <- dat$pr.ext[which(dat$b1 >= 0.95*b1 & dat$b1 <= 1.05*b1)]
mean.b1 <- mean(pr.ext.b1)
mean.b1
(mean.b1-mean.b1.prime)/mean.b1 * 100 # % reduction
# bootstrap
iter <- 100000
pc.red.b1 <- rep(0,iter)
for (i in 1:iter) {
  b1.boot <- sample(pr.ext.b1, 1, replace=T)
  b1.prime.boot <- sample(pr.ext.b1.prime, 1, replace=T)
  pc.red.b1[i] <- (b1.boot-b1.prime.boot)/b1.boot * 100
}
hist(pc.red.b1)
mean(pc.red.b1[is.finite(pc.red.b1)], na.rm=T)
# conclusion: 50% reduction in age-independent mortality -> 68% reduction in probability of extinction

# set in model: 'fert.bentley <- 4.69/2' = 2.345; test random uniform sample between 2.1105 and 2.5795 (10%)
plot(dat$fert.bentley, dat$pr.ext, pch=19, xlab="total fertility (F)", ylab="probability of extinction")
fit.F <- lm((dat$pr.ext) ~ dat$fert.bentley)
abline(fit.F, lty=2, col="red")
# 10% lower F ...
F <- 2.345
F.prime <- 0.9*F
abline(v=F.prime, lty=3)
pr.ext.F.prime <- dat$pr.ext[which(dat$fert.bentley >= 0.95*F.prime & dat$fert.bentley <= 1.05*F.prime)]
mean.F.prime <- mean(pr.ext.F.prime)
mean.F.prime
abline(v=F, lty=3)
pr.ext.F <- dat$pr.ext[which(dat$fert.bentley >= 0.95*F & dat$fert.bentley <= 1.05*F)]
mean.F <- mean(pr.ext.F)
mean.F
(mean.F.prime-mean.F)/mean.F*100 # % increase
# bootstrap
iter <- 100000
pc.red.F <- rep(0,iter)
for (i in 1:iter) {
  F.boot <- sample(pr.ext.F, 1, replace=T)
  F.prime.boot <- sample(pr.ext.F.prime, 1, replace=T)
  pc.red.F[i] <- (F.prime.boot-F.boot)/F.boot * 100
}
hist(pc.red.F)
mean(pc.red.F[is.finite(pc.red.F)], na.rm=T)
# conclusion: 10% reduction in fertility -> 773% increase in probability of extinction

## set in model: a1 <- 0.422; test random uniform sample between 0.3 and 0.5
plot(dat$a1, dat$pr.ext, pch=19, xlab="initial infant mortality (a1)", ylab="probability of extinction")
fit.a1 <- lm((dat$pr.ext) ~ dat$a1)
abline(fit.a1, lty=2, col="red")

## set in model: a2 <- 0.013; test random uniform sample between 0.010 and 0.015
plot(dat$a2, dat$pr.ext, pch=19, xlab="mortality decline rate (a2)", ylab="probability of extinction")
fit.a2 <- lm((dat$pr.ext) ~ dat$a2)
abline(fit.a2, lty=2, col="red")

## set in model to 'cat.numerator <- 0.14'; test random uniform sample between 0.1 and 0.2
plot(dat$cat.numerator, dat$pr.ext, pch=19, xlab="Pr(catastrophe) (Pc)", ylab="probability of extinction")
fit.Pc <- lm((dat$pr.ext) ~ dat$cat.numerator)
abline(fit.Pc, lty=2, col="red")

## set in model 'start.K' <- 47000; test random uniform sample between 50% (23500) and 200% (94000)
plot(dat$start.K, dat$pr.ext, pch=19, xlab="nadir population density (Dmin)", ylab="probability of extinction")
fit.Dmin <- lm((dat$pr.ext) ~ dat$start.K)
abline(fit.Dmin, lty=2, col="red")
summary(fit.Dmin)
