## Human arrival population model (Bird et al.)
## Corey Bradshaw
## Sep 2018

## Remove everything
rm(list = ls())

## libraries
library(boot)
library(tcltk)
library(plotly)    

## source
source("~/Documents/Papers/Other/Global human population/ms/PNAS/R1/matrixOperators.r")

# Set functions
# beta distribution shape parameter estimator function
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

# stochastic beta sampler for a survival vector (replaces above function for faster processing)
stoch.surv.func <- function(mu, var) {
  Sx <- rbeta(length(mu), (((1 - mu) / var - 1 / mu) * mu ^ 2), ((((1 - mu) / var - 1 / mu) * mu ^ 2)*(1 / mu - 1)))
  return(params=Sx)
}

## set working directory
setwd("~/Documents/Papers/Palaeo/Sahul/Aus human age-structured model/data")

# Siler hazard h(x) (Gurven et al. 2007)
# average hunter-gatherer
a1 <- 0.422 # initial infant mortality rate (also known as αt)
b1 <- 1.131 # rate of mortality decline (also known as bt)
a2 <- 0.013 # age-independent mortality (exogenous mortality due to environment); also known as ct
a3 <- 1.47e-04 # initial adult mortality rate (also known as βt)
b3 <- 0.086 # rate of mortality increase
longev <- 80
x <- seq(0,longev,1) # age vector
h.x <- a1 * exp(-b1*x) + a2 + a3 * exp(b3 * x) # Siler's hazard model
plot(x,h.x,pch=19,type="l")
plot(x,log(h.x),pch=19,type="l")
l.x <- exp((-a1/b1) * (1 - exp(-b1*x))) * exp(-a2 * x) * exp(a3/b3 * (1 - exp(b3 * x))) # Siler's survival (proportion surviving) model
plot(x,l.x,type="l")

l.inf <- exp(-a1/b1) # survival at infinite time
T.m <- 1/b1 # time constant at which maturity is approached
h.m <- a2 # hazard for mature animals
l.m <- exp(-a2*x) # survival
h.s <- a3*exp(b3*x) # hazard for senescence
l.s <- exp((a3/b3)*(1 - exp(b3*x))) # survival for senescence
f.x <- a3*exp(b3*x)*exp((a3/b3)/(1-exp(b3*x))) # probability density function
(log(a3) - log(a1)) / a3
T.s <- (1/b3) # modal survival time

## survival
init.pop <- 10000
lx <- round(init.pop*l.x,0)
len.lx <- length(lx)
dx <- lx[1:(len.lx-1)]-lx[2:len.lx]
dx
qx <- dx/lx[1:(length(lx)-1)]
qx
Sx <- 1 - qx
Sx
sx <- lx[2:len.lx]/lx[1:(len.lx-1)]
mx <- 1 - sx
Lx <- (lx[1:(len.lx-1)] + lx[2:len.lx])/2
ex <- rev(cumsum(rev(Lx)))/lx[-len.lx]
ex
ex.avg <- ex + x[-len.lx]
ex.avg

par(mfrow=c(2,1))
plot(x[-1], Sx, pch=19, type="l", xlab="age (years)", ylab="Sx")
plot(x[-1], ex, pch=19, type="l", xlab="age (years)", ylab="life expectancy")
par(mfrow=c(1,1))
surv.out <- data.frame(x[-1],Sx,ex)
colnames(surv.out)[1] <- "x"

# fertility (Walker et al. 2006)
primiparity.walker <- c(17.7,18.7,19.5,18.5,18.5,18.7,25.7,19,20.5,18.8,17.8,18.6,22.2,17,16.2,18.4)
prim.mean <- round(mean(primiparity.walker),0)
prim.lo <- round(quantile(primiparity.walker,probs=0.025),0)
prim.hi <- round(quantile(primiparity.walker,probs=0.975),0)
print(c(prim.lo, prim.mean, prim.hi))
setwd("~/Documents/Papers/Other/Global human population/data/import data/")
dat.world13 <- read.table("world2013lifetable.csv", header=T, sep=",")
fert.world13 <- dat.world13$m.f
fert.trunc <- fert.world13[1:(longev+1)]
pfert.trunc <- fert.trunc/sum(fert.trunc)
fert.bentley <- 4.69/2 # Bentley 1985 for !Kung
fert.vec <- fert.bentley * pfert.trunc
plot(x,fert.vec, type="l", xlab="age (years)", ylab="fertility")

fert.out <- data.frame(x,fert.vec)
colnames(fert.out)[2] <- "fert"

## construct matrix
stages <- len.lx
popmat <- matrix(0,nrow=stages,ncol=stages)
colnames(popmat) <- x
rownames(popmat) <- x

## populate matrix
popmat[1,] <- fert.vec
diag(popmat[2:stages,]) <- Sx
popmat[stages,stages] <- 0 # Sx[stages-1]
popmat.orig <- popmat ## save original matrix

## matrix properties
max.lambda(popmat) ## 1-yr lambda
max.r(popmat) # rate of population change, 1-yr
stable.stage.dist(popmat) ## stable stage distribution
plot(x,stable.stage.dist(popmat),type="l")
R.val(popmat,stages) # reproductive value
gen.l <- G.val(popmat,stages) # mean generation length
cat.pr <- 0.14/gen.l # probability of catastrophe (Reed et al. 2003)

## initial population vector
pop.found <- 125 # founding population size
init.vec <- stable.stage.dist(popmat) * pop.found


#################
## project
## set time limit for projection in 1-yr increments
yr.now <- 1 # update if more data available post-2010
#************************
yr.end <- 50000 # set projection end date
#************************
t <- (yr.end - yr.now)

tot.F <- sum(popmat.orig[1,])
popmat <- popmat.orig
yr.vec <- seq(yr.now,yr.end)

## set population storage matrices
n.mat <- matrix(0,nrow=stages,ncol=(t+1))
n.mat[,1] <- init.vec

## set up projection loop
for (i in 1:t) {
  n.mat[,i+1] <- popmat %*% n.mat[,i] 
}

## year projection vector
yrs <- seq(yr.now,yr.end,1)

# plot
plot(yrs,as.vector(colSums(n.mat)),type="l",xlab="year",ylab="N",xlim=c(yr.end,yr.now), xaxt="n")
axis(1, at = seq(0, yr.end, 1000), las=1)


####################################################
## relative density, carrying capacity & feedbacks
setwd("~/Documents/Papers/Palaeo/Sahul/Aus human age-structured model/data/updated K data/")

## NORTH OF SAHUL (between 10 and 18 degrees S Latitude)
# net primary productivity (NPP; kg C/yr/m2)
# ALL SAHUL FIRST
npp.sah <- read.table("ClimateSahul_Npp.csv", header=T, sep=",") 

## JUST NORTH (0-10 DEGREES)
npp.nsah <- subset(npp.sah, Lat..degree. <= 0 & Lat..degree. >= -10)
npp.sah60 <- npp.sah[,c(1,2,63)]
head(npp.sah60)

# plot raster
library(sp)
library(rgdal)
library(raster)
coordinates(npp.sah60) = ~ Lon..degree. + Lat..degree.
proj4string(npp.sah60)=CRS("+proj=longlat +datum=WGS84") # set it to lat-long
gridded(npp.sah60) = TRUE
npp60 = raster(npp.sah60)
plot(npp60, xlab="longitude", ylab="latitude",col=)

npp.nsah.med <- as.numeric(apply(npp.nsah[,-c(1,2)], MARGIN=2, quantile, probs=0.5))
npp.nsah.lo <- as.numeric(apply(npp.nsah[,-c(1,2)], MARGIN=2, quantile, probs=0.25))
npp.nsah.up <- as.numeric(apply(npp.nsah[,-c(1,2)], MARGIN=2, quantile, probs=0.75))

npp.nsah.dat <- as.data.frame(t(data.frame(npp.nsah.lo,npp.nsah.med,npp.nsah.up)))
rownames(npp.nsah.dat) <- c("lo","med","up")
colnames(npp.nsah.dat) <- colnames(npp.nsah[,-c(1,2)])
head(npp.nsah.dat)

yr1000 <- seq(120,0,-1)
plot(yr1000,npp.nsah.dat[2,],type="l",lwd=2,ylim=c(min(npp.nsah.dat[1,]),max(npp.nsah.dat[3,])), xlab="ka", ylab="NPP")
lines(yr1000,npp.nsah.dat[1,],type="l",lty=2,lwd=1)
lines(yr1000,npp.nsah.dat[3,],type="l",lty=2,lwd=1)
abline(v=55,lty=3,lwd=3)
abline(v=65,lty=3,lwd=3)


# interpolate values between successive points
year.vec <- 120000:0
npp.med.interp <- approx(yr1000*1000, npp.nsah.dat[2,], xout=year.vec, method="linear")
npp.lo.interp <- approx(yr1000*1000, npp.nsah.dat[1,], xout=year.vec, method="linear")
npp.up.interp <- approx(yr1000*1000, npp.nsah.dat[3,], xout=year.vec, method="linear")
plot(year.vec,npp.med.interp$y,type="l",lwd=2,ylim=c(min(npp.lo.interp$y),max(npp.up.interp$y)), xlab="ka", ylab="NPP")
lines(year.vec,npp.lo.interp$y,type="l",lty=2,lwd=1,col="red")
lines(year.vec,npp.up.interp$y,type="l",lty=2,lwd=1,col="red")
abline(v=60000,lty=3,lwd=3)
abline(v=70000,lty=3,lwd=3)
K.dat <- data.frame(year.vec,npp.lo.interp$y,npp.med.interp$y,npp.up.interp$y)
colnames(K.dat) <- c("year","npp.lo","npp.med","npp.up")

## Calculation of K (n people)
K.dat$npp.lo.sc1 <- scale(K.dat$npp.lo + abs(min(K.dat$npp.lo)), scale=F, center=F)
K.dat$npp.med.sc1 <- scale(K.dat$npp.med + abs(min(K.dat$npp.lo)), scale=F, center=F) 
K.dat$npp.up.sc1 <- scale(K.dat$npp.up + abs(min(K.dat$npp.lo)), scale=F, center=F)

K.dat$npp.lo.sc <- K.dat$npp.lo.sc1/max(K.dat$npp.med.sc1)
K.dat$npp.med.sc <- K.dat$npp.med.sc1/max(K.dat$npp.med.sc1)
K.dat$npp.up.sc <- K.dat$npp.up.sc1/max(K.dat$npp.med.sc1)

plot(K.dat$year, K.dat$npp.med.sc, type="l", xlab="year", ylab="scaled NPP", lty=1, lwd=2, ylim=c(min(K.dat$npp.lo.sc), max(K.dat$npp.up.sc)))
lines(K.dat$year, K.dat$npp.lo.sc, lty=2, col="red")
lines(K.dat$year, K.dat$npp.up.sc, lty=2, col="red")
abline(h=0,lty=2)

start.K <- 47000 # Guatney & Holliday 2015
K.dat$K.med <- round(start.K * K.dat$npp.med.sc, 0)
K.dat$K.lo <- round(start.K * K.dat$npp.lo.sc, 0)
K.dat$K.up <- round(start.K * K.dat$npp.up.sc, 0)

K.6070 <- subset(K.dat, year <= 70000 & year >= 65000)

plot(K.dat$year, K.dat$K.med, type="l", xlab="year", ylab="relative K", lty=1, lwd=2, ylim=c(min(K.dat$K.lo), max(K.dat$K.up)))
lines(K.dat$year, K.dat$K.lo, lty=2, col="red")
lines(K.dat$year, K.dat$K.up, lty=2, col="red")
abline(v=60000,lty=3,lwd=3)
abline(v=70000,lty=3,lwd=3)


## plot during sampled interval only
plot(K.dat$year, K.dat$K.med, type="l", xlab="year", ylab="relative K", lty=1, lwd=2, xlim=c(0.95*55000, 1.05*65000), ylim=c(min(K.dat$K.lo), max(K.dat$K.up)))
lines(K.dat$year, K.dat$K.lo, lty=2, col="red")
lines(K.dat$year, K.dat$K.up, lty=2, col="red")
abline(v=55000,lty=3,lwd=3)
abline(v=65000,lty=3,lwd=3)

## density feedback function on survival
popmat <- popmat.orig
Sx.mod <- 0.99615*Sx
diag(popmat[2:stages,]) <- Sx.mod
popmat[stages,stages] <- 0 # Sx[stages-1]
max.r(popmat) # rate of population change, 1-yr

surv.min <- 0.98




#####################################################################################
## NORTHERN ROUTE
## Mangoli—Buru—Ceram—Sahul
## Pr(random) = 0.054, 0.187, 0.251
## areas = 8966—9302—21495—10000000 (km^2)
## coords = (125.8,-1.8)—(126.6-3.5)—(129.4,-3)—(133,-3.5)
#####################################################################################

# initial colonisation window
col.old <- 70000
col.yng <- 60000

# iterations
iter <- 1000
itdiv <- iter/100

# relative island sizes (to Sahul)
mangoli.A <- 8966
buru.A <- 9302
ceram.A <- 21495
sahul.A <- 786000 # kmˆ2 (only New Guinea)

mangoli.rA <- mangoli.A/sahul.A
buru.rA <- buru.A/sahul.A
ceram.rA <- ceram.A/sahul.A

# population densities (from Tallavaara et al. 2018 PNAS)
mangoli.D <- 0.1750582
buru.D <- 0.1852028
ceram.D <- 0.2796982
sahul.D <- 0.7452437

# generations to project
#gen.proj <- 200
gen.proj <- 100

# set probabilities of randomly drifting across divide between islands (from M. Bird)
mangoli.buru.drift.pr <- 0.054
buru.ceram.drift.pr <- 0.187
ceram.sahul.drift.pr <- 0.251

# set quasi-extinction threshold (min number of females)
q.ext <- 10

# set SD for Sx
Sx.sd <- 0.05 # can set to any value

## set colonist sample size (washed off) vector
col.size.vec <- 2*seq(10,50,5)

## set set probability of a random colonisation event occurring vector
pr.col.vec <- seq(0.01, 0.25, 0.02)

# set storage matrix (col.size x pr.col)
pr.sahul.arrive.stor <- pr.sahul.ext.stor <- matrix(0, nrow=length(pr.col.vec), ncol=length(col.size.vec))

for (c in 1:length(col.size.vec)) {
  for (p in 1:length(pr.col.vec)) {
    
    # storage vectors
    sahul.arrive <- sahul.ext <- rep(0,iter)
    
    # progress counter
    pb <- txtProgressBar(min=1, max=iter, style=3)
    
    # iterate over attempts
    for (e in 1:iter) {
    
      ## set time limit for projection in 1-yr increments
      yr.st <- round(runif(iter,col.yng,col.old), 0) 
      #************************
      yr.end <- yr.st - round(gen.proj*gen.l, 0) # set projection end date
      #************************
      t <- (yr.st - yr.end)
      
      # reset all popmats
      popmat.mangoli <- popmat.buru <- popmat.ceram <- popmat.sahul <- popmat.orig
      
      ## set K run series for each island
      yr.st.sub <- which(K.dat$year == yr.st[e])
      npp.sc.st <- runif(1, K.dat$npp.lo.sc[yr.st.sub], K.dat$npp.up.sc[yr.st.sub])
      mangoli.K1.run <- mangoli.rA * round((npp.sc.st / K.dat$npp.med.sc[yr.st.sub]) * K.dat$K.med[yr.st.sub:(yr.st.sub+t[e]+4)], 0) # this iteration's realised K series
      mangoli.K.run <- mangoli.K1.run * 1/(mean(mangoli.K1.run)/(mangoli.D*mangoli.A))
      buru.K.run <- (buru.D/mangoli.D) * (buru.A/mangoli.A) * mangoli.K.run
      ceram.K.run <- (ceram.D/buru.D) * (ceram.A/buru.A) * buru.K.run
      sahul.K.run <- (sahul.D/ceram.D) * (sahul.A/ceram.A) * ceram.K.run
      
      ## initial population vector
      ## assume Mangoli Island is at K
      mangoli.init.vec <- stable.stage.dist(popmat.mangoli) * mangoli.K.run[e] # stable stage distribution x founding population size
      
      ## set population storage matrices
      n.mat.mangoli <- matrix(0,nrow=stages,ncol=(t[e]+1))
      n.mat.buru <- matrix(0,nrow=stages,ncol=(t[e]+2)) # must increment by one to make i loop work
      n.mat.ceram <- matrix(0,nrow=stages,ncol=(t[e]+3)) # must increment by one to make i loop work
      n.mat.sahul <- matrix(0,nrow=stages,ncol=(t[e]+4)) # must increment by one to make i loop work
      
      # start everything off in Mangoli
      n.mat.mangoli[,1] <- mangoli.init.vec
        
      ## set up projection loop
      for (i in 1:t[1]) { # i projection loop
        
        #############
        ## MANGOLI ##
        #############
        ## reconstruct popmat with stochastic elements
        popmat.mangoli[1, ] <- pfert.trunc * rnorm(1, fert.bentley,0.05*fert.bentley) # fertility sampler
        diag(popmat.mangoli[2:(stages), ]) <- ((ifelse(sum(n.mat.mangoli[,i], na.rm=T) > mangoli.K.run[i], surv.min, 1)) * (stoch.surv.func(Sx, Sx.sd^2))[-stages]) # survival sampler
        
        # survival (+ catastrophic mortality at 50%)
        if ((rbinom(1, 1, cat.pr)) == 1) {
          diag(popmat.mangoli[2:(stages), ]) <- (0.5* (diag(popmat.mangoli[2:(stages), ])))}
        popmat.mangoli[stages,stages] <- 0
        
        ## project over interval
        n.mat.mangoli[,i+1] <- popmat.mangoli %*% n.mat.mangoli[,i]
        
        ## did the Mangoli population go extinct?
        if (sum(n.mat.mangoli[,i+1]) < q.ext) {
          n.mat.mangoli[,i+1] <- 0
        }
        
        ##########################################
        ## colonisation event from Mangoli-Buru #
        ##########################################
        if ((rbinom(1, 1, pr.col.vec[p]*mangoli.buru.drift.pr)) == 1) { # is it a colonisation event, and did they make it?
          col.n <- round(rnorm(1, mean=col.size.vec[c]/2, sd=0.1*col.size.vec[c]/2),0) ## size of Mangoli-Buru colonisation event
          if (round(sum(n.mat.mangoli[(16:stages),i+1]), 0) >= col.n) { # are they enough adults in Mangoli to wash off to Buru?
            n.mat.buru[16:stages,i+1] <- n.mat.buru[16:stages,i+1] + (stable.stage.dist(popmat.buru)[16:stages] * col.n) # add new adult colonists to Buru
            n.mat.mangoli[16:stages,i+1] <- n.mat.mangoli[16:stages,i+1] - n.mat.buru[16:stages,i+1] # take away washed-off colonists from Mangoli 
            if (length(which(n.mat.mangoli[,i+1] < 0)) > 0) {
              n.mat.mangoli[which(n.mat.mangoli[,i+1] < 0), i+1] <- 0   }}} # make sure no negative values in Mangoli
        
        ##########
        ## BURU ##
        ##########
        if (sum(n.mat.buru[,i+1]) > 0) { # only project if there are people in this iteration on Buru
          
          ## reconstruct popmat with stochastic elements
          popmat.buru[1, ] <- pfert.trunc * rnorm(1, fert.bentley,0.05*fert.bentley) # fertility sampler
          diag(popmat.buru[2:(stages), ]) <- ((ifelse(sum(n.mat.buru[,i+1], na.rm=T) > buru.K.run[i+1], surv.min, 1)) * (stoch.surv.func(Sx, Sx.sd^2))[-stages]) # survival sampler
          
          # survival (+ catastrophic mortality at 50%)
          if ((rbinom(1, 1, cat.pr)) == 1) {
            diag(popmat.buru[2:(stages), ]) <- (0.5* (diag(popmat.buru[2:(stages), ])))}
          popmat.buru[stages,stages] <- 0
          
          ## project over interval
          n.mat.buru[,i+2] <- popmat.buru %*% n.mat.buru[,i+1]
          
          ## did the Buru population go extinct?
          if (sum(n.mat.buru[,i+2]) < q.ext) {
            n.mat.buru[,i+2] <- 0}}
        
        ########################################
        ## colonisation event from Buru-Ceram ##
        ########################################
        if ((rbinom(1, 1, pr.col.vec[p]*buru.ceram.drift.pr)) == 1) { # is it a colonisation event, and did they make it?
          col2.n <- round(rnorm(1, mean=col.size.vec[c]/2, sd=0.1*col.size.vec[c]/2),0) ## size of Buru-Ceram colonisation event
          if (round(sum(n.mat.buru[(16:stages),i+1]), 0) >= col2.n) { # are they enough adults in Buru to wash off to Ceram?
            n.mat.ceram[16:stages,i+2] <- n.mat.ceram[16:stages,i+2] + (stable.stage.dist(popmat.ceram)[16:stages] * col2.n) # add new adult colonists to Ceram
            n.mat.buru[16:stages,i+2] <- n.mat.buru[16:stages,i+2] - n.mat.ceram[16:stages,i+2] # take away washed-off colonists from Buru 
            if (length(which(n.mat.buru[,i+2] < 0)) > 0) {
              n.mat.buru[which(n.mat.buru[,i+2] < 0), i+2] <- 0   }}} # make sure no negative values in Buru
    
        ###########
        ## CERAM ##
        ###########
        if (sum(n.mat.ceram[,i+2]) > 0) { # only project if there are people in this iteration on Ceram
          
          ## reconstruct popmat with stochastic elements
          popmat.ceram[1, ] <- pfert.trunc * rnorm(1, fert.bentley,0.05*fert.bentley) # fertility sampler
          diag(popmat.ceram[2:(stages), ]) <- ((ifelse(sum(n.mat.ceram[,i+2], na.rm=T) > ceram.K.run[i+2], surv.min, 1)) * (stoch.surv.func(Sx, Sx.sd^2))[-stages]) # survival sampler
          
          # survival (+ catastrophic mortality at 50%)
          if ((rbinom(1, 1, cat.pr)) == 1) {
            diag(popmat.ceram[2:(stages), ]) <- (0.5* (diag(popmat.ceram[2:(stages), ])))}
          popmat.ceram[stages,stages] <- 0
          
          ## project over interval
          n.mat.ceram[,i+3] <- popmat.ceram %*% n.mat.ceram[,i+2]
          
          ## did the Ceram population go extinct?
          if (sum(n.mat.ceram[,i+3]) < q.ext) {
            n.mat.ceram[,i+3] <- 0}}
      
        #########################################
        ## colonisation event from Ceram-Sahul ##
        #########################################
        if ((rbinom(1, 1, pr.col.vec[p]*ceram.sahul.drift.pr)) == 1) { # is it a colonisation event, and did they make it?
          col3.n <- round(rnorm(1, mean=col.size.vec[c]/2, sd=0.1*col.size.vec[c]/2),0) ## size of Ceram-Sahul colonisation event
          if (round(sum(n.mat.ceram[(16:stages),i+3]), 0) >= col3.n) { # are they enough adults in Ceram to wash off to Sahul?
            n.mat.sahul[16:stages,i+3] <- n.mat.sahul[16:stages,i+3] + (stable.stage.dist(popmat.sahul)[16:stages] * col3.n) # add new adult colonists to Sahul
            n.mat.ceram[16:stages,i+3] <- n.mat.ceram[16:stages,i+3] - n.mat.sahul[16:stages,i+3] # take away washed-off colonists from Ceram 
            if (length(which(n.mat.ceram[,i+3] < 0)) > 0) {
              n.mat.ceram[which(n.mat.ceram[,i+3] < 0), i+3] <- 0   }}} # make sure no negative values in Ceram
    
        ###########
        ## SAHUL ##
        ###########
        if (sum(n.mat.sahul[,i+3]) > 0) { # only project if there are people in this iteration on Sahul
          
          ## reconstruct popmat with stochastic elements
          popmat.sahul[1, ] <- pfert.trunc * rnorm(1, fert.bentley,0.05*fert.bentley) # fertility sampler
          diag(popmat.sahul[2:(stages), ]) <- ((ifelse(sum(n.mat.sahul[,i+3], na.rm=T) > sahul.K.run[i+3], surv.min, 1)) * (stoch.surv.func(Sx, Sx.sd^2))[-stages]) # survival sampler
          
          # survival (+ catastrophic mortality at 50%)
          if ((rbinom(1, 1, cat.pr)) == 1) {
            diag(popmat.sahul[2:(stages), ]) <- (0.5* (diag(popmat.sahul[2:(stages), ])))}
          popmat.sahul[stages,stages] <- 0
          
          ## project over interval
          n.mat.sahul[,i+4] <- popmat.sahul %*% n.mat.sahul[,i+3]
          
          ## did the Sahul population go extinct?
          if (sum(n.mat.sahul[,i+4]) < q.ext) {
            n.mat.sahul[,i+4] <- 0}}
        
          #print(i)
      } # end i loop  
      
      ## did people make it to Sahul?
      sahul.arrive[e] <- ifelse(length(which((round(colSums(n.mat.sahul),0)) > 0)) > 0, 1, 0)
      
      # if arrived in Sahul, what was the probability that the population went extinct?
      sahul.ext[e] <- ifelse(sahul.arrive[e] == 1 & (ifelse(sahul.arrive[e] == 1, (round(colSums(n.mat.sahul),0))[length(round(colSums(n.mat.sahul),0))], NA)) == 0, 1, NA)
      
      setTxtProgressBar(pb, e)
      # if (e %% itdiv==0) print(e)
    } # end e loop
    
    # probability of arrival in Sahul
    pr.sahul.arrive.stor[p, c] <- sum(sahul.arrive)/iter
    
    # probability of arrival & extinction in Sahul
    pr.sahul.ext.stor[p, c] <- sum(sahul.ext,na.rm=T)/sum(sahul.arrive,na.rm=T)
  
    print("")
    print(paste("probability of a colonisation event = ", pr.col.vec[p], sep=""))
  } # end p loop (pr.vec)

  print("##################################################")
  print(paste("colonisation event size (male & female) = ", col.size.vec[c], sep=""))
  print("##################################################")
  
} # end c loop (col.size)

# plot 3D surfaces (pr.col x col.size)
# probability of arrival to Sahul
pr.sahul.arrive.stor
persp3D(x=pr.col.vec, y=col.size.vec, z=pr.sahul.arrive.stor, border="white", bty="f", resfac=2, theta = 230, phi= 20, col=grey.colors(20), image=F, contour=F, cex.axis=0.6, ticktype="detailed", xlab="Pr(colonisation event)", ylab="event size (n)", zlab="Pr(reach Sahul)", expand=0.7)

b <- list(title="Pr(colonisation event)", tickmode="array", nticks=length(pr.col.vec), tickvals=pr.col.vec)
a <- list(title="event size (n)", tickmode="array", nticks=length(col.size.vec), tickvals=col.size.vec)
p1 <- plot_ly(y=~pr.col.vec, x=~col.size.vec, z=pr.sahul.arrive.stor, type = "contour", colorscale='Viridis', reversescale=T) %>%
  colorbar(title = "Pr(reach Sahul)") %>%
  layout(xaxis = a, yaxis=b)
p1

pr.sahul.ext.stor
persp3D(x=pr.col.vec, y=col.size.vec, z=pr.sahul.ext.stor, border="white", bty="f", resfac=2, theta = 120, phi= 20, col=rev(grey.colors(20)), image=F, contour=F, cex.axis=0.6, ticktype="detailed", xlab="Pr(colonisation event)", ylab="event size (n)", zlab="Pr(Sahul extinction)", expand=0.7)

b <- list(title="Pr(colonisation event)", tickmode="array", nticks=length(pr.col.vec), tickvals=pr.col.vec)
a <- list(title="event size (n)", tickmode="array", nticks=length(col.size.vec), tickvals=col.size.vec)
p <- plot_ly(y=~pr.col.vec, x=~col.size.vec, z=pr.sahul.ext.stor, type = "contour", colorscale='Viridis', reversescale=F) %>%
  colorbar(title = "Pr(Sahul extinction)") %>%
  layout(xaxis = a, yaxis=b)
p



#####################################################################################
## SOUTHERN ROUTE
## Alor—Timor—Sahul
## Pr(random) = 0.013, 0.013
## areas = 3803—35527—10000000 (km^2)
## coords = (-8.3,124.7)—(126.6-3.5)—(129.4,-3)—(133,-3.5)
#####################################################################################

# initial colonisation window
col.old <- 70000
col.yng <- 60000

# iterations
iter <- 1000
itdiv <- iter/100

# relative island sizes (to Sahul)
alors.A <- 3803
timor.A <- 35527
sahul.A <- 1000000 # 'northern' Aus

alor.rA <- alor.A/sahul.A
timor.rA <- timor.A/sahul.A

# population densities (from Tallavaara et al. 2018 PNAS)
alor.D <- 0.1016968
timor.D <- 0.112876
sahul.D <- 0.1427313

# generations to project
#gen.proj <- 200
gen.proj <- 100

# set probabilities of randomly drifting across divide between islands (from M. Bird)
alor.timor.drift.pr <- 0.013
timor.sahul.drift.pr <- 0.013

# set quasi-extinction threshold (min number of females)
q.ext <- 10

# set SD for Sx
Sx.sd <- 0.05 # can set to any value

## set colonist sample size (washed off) vector
col.size.vec <- 2*seq(10,50,5)

## set set probability of a random colonisation event occurring vector
pr.col.vec <- seq(0.01, 0.25, 0.02)

# set storage matrix (col.size x pr.col)
pr.sahul.arrive.stor <- pr.sahul.ext.stor <- matrix(0, nrow=length(pr.col.vec), ncol=length(col.size.vec))

for (c in 1:length(col.size.vec)) {
  for (p in 1:length(pr.col.vec)) {
    
    # storage vectors
    sahul.arrive <- sahul.ext <- rep(0,iter)
    
    # progress counter
    pb <- txtProgressBar(min=1, max=iter, style=3)
    
    # iterate over attempts
    for (e in 1:iter) {
      
      ## set time limit for projection in 1-yr increments
      yr.st <- round(runif(iter,col.yng,col.old), 0) 
      #************************
      yr.end <- yr.st - round(gen.proj*gen.l, 0) # set projection end date
      #************************
      t <- (yr.st - yr.end)
      
      # reset all popmats
      popmat.alor <- popmat.timor <- popmat.sahul <- popmat.orig
      
      ## set K run series for each island
      yr.st.sub <- which(K.dat$year == yr.st[e])
      npp.sc.st <- runif(1, K.dat$npp.lo.sc[yr.st.sub], K.dat$npp.up.sc[yr.st.sub])
      alor.K1.run <- alor.rA * round((npp.sc.st / K.dat$npp.med.sc[yr.st.sub]) * K.dat$K.med[yr.st.sub:(yr.st.sub+t[e]+3)], 0) # this iteration's realised K series
      alor.K.run <- alor.K1.run * 1/(mean(alor.K1.run)/(alor.D*alor.A))
      timor.K.run <- (timor.D/alor.D) * (timor.A/alor.A) * alor.K.run
      sahul.K.run <- (sahul.D/timor.D) * (sahul.A/timor.A) * timor.K.run
      
      ## initial population vector
      ## assume Alor Island is at K
      alor.init.vec <- stable.stage.dist(popmat.alor) * alor.K.run[e] # stable stage distribution x founding population size
      
      ## set population storage matrices
      n.mat.alor <- matrix(0,nrow=stages,ncol=(t[e]+1))
      n.mat.timor <- matrix(0,nrow=stages,ncol=(t[e]+2)) # must increment by one to make i loop work
      n.mat.sahul <- matrix(0,nrow=stages,ncol=(t[e]+3)) # must increment by one to make i loop work
      
      # start everything off in Mangoli
      n.mat.alor[,1] <- alor.init.vec
      
      ## set up projection loop
      for (i in 1:t[1]) { # i projection loop
        
        ##########
        ## ALOR ##
        ##########
        ## reconstruct popmat with stochastic elements
        popmat.alor[1, ] <- pfert.trunc * rnorm(1, fert.bentley,0.05*fert.bentley) # fertility sampler
        diag(popmat.alor[2:(stages), ]) <- ((ifelse(sum(n.mat.alor[,i], na.rm=T) > alor.K.run[i], surv.min, 1)) * (stoch.surv.func(Sx, Sx.sd^2))[-stages]) # survival sampler
        
        # survival (+ catastrophic mortality at 50%)
        if ((rbinom(1, 1, cat.pr)) == 1) {
          diag(popmat.alor[2:(stages), ]) <- (0.5* (diag(popmat.alor[2:(stages), ])))}
        popmat.alor[stages,stages] <- 0
        
        ## project over interval
        n.mat.alor[,i+1] <- popmat.alor %*% n.mat.alor[,i]
        
        ## did the Alor population go extinct?
        if (sum(n.mat.alor[,i+1]) < q.ext) {
          n.mat.alor[,i+1] <- 0
        }
        
        ########################################
        ## colonisation event from Alor-Timor ##
        ########################################
        if ((rbinom(1, 1, pr.col.vec[p]*alor.timor.drift.pr)) == 1) { # is it a colonisation event, and did they make it?
          col.n <- round(rnorm(1, mean=col.size.vec[c]/2, sd=0.1*col.size.vec[c]/2),0) ## size of Alor-Timor colonisation event
          if (round(sum(n.mat.alor[(16:stages),i+1]), 0) >= col.n) { # are they enough adults in Alor to wash off to Timor?
            n.mat.timor[16:stages,i+1] <- n.mat.timor[16:stages,i+1] + (stable.stage.dist(popmat.timor)[16:stages] * col.n) # add new adult colonists to Timor
            n.mat.alor[16:stages,i+1] <- n.mat.alor[16:stages,i+1] - n.mat.timor[16:stages,i+1] # take away washed-off colonists from Alor 
            if (length(which(n.mat.alor[,i+1] < 0)) > 0) {
              n.mat.alor[which(n.mat.alor[,i+1] < 0), i+1] <- 0   }}} # make sure no negative values in Alor
        
        ###########
        ## TIMOR ##
        ###########
        if (sum(n.mat.timor[,i+1]) > 0) { # only project if there are people in this iteration on Timor
          
          ## reconstruct popmat with stochastic elements
          popmat.timor[1, ] <- pfert.trunc * rnorm(1, fert.bentley,0.05*fert.bentley) # fertility sampler
          diag(popmat.timor[2:(stages), ]) <- ((ifelse(sum(n.mat.timor[,i+1], na.rm=T) > buru.K.run[i+1], surv.min, 1)) * (stoch.surv.func(Sx, Sx.sd^2))[-stages]) # survival sampler
          
          # survival (+ catastrophic mortality at 50%)
          if ((rbinom(1, 1, cat.pr)) == 1) {
            diag(popmat.timor[2:(stages), ]) <- (0.5* (diag(popmat.timor[2:(stages), ])))}
          popmat.timor[stages,stages] <- 0
          
          ## project over interval
          n.mat.timor[,i+2] <- popmat.timor %*% n.mat.timor[,i+1]
          
          ## did the Timor population go extinct?
          if (sum(n.mat.timor[,i+2]) < q.ext) {
            n.mat.timor[,i+2] <- 0}}
        
        #########################################
        ## colonisation event from Timor-Sahul ##
        #########################################
        if ((rbinom(1, 1, pr.col.vec[p]*timor.sahul.drift.pr)) == 1) { # is it a colonisation event, and did they make it?
          col2.n <- round(rnorm(1, mean=col.size.vec[c]/2, sd=0.1*col.size.vec[c]/2),0) ## size of Timor-Sahul colonisation event
          if (round(sum(n.mat.timor[(16:stages),i+2]), 0) >= col2.n) { # are they enough adults in Timor to wash off to Sahul?
            n.mat.sahul[16:stages,i+2] <- n.mat.sahul[16:stages,i+2] + (stable.stage.dist(popmat.sahul)[16:stages] * col2.n) # add new adult colonists to Sahul
            n.mat.timor[16:stages,i+2] <- n.mat.timor[16:stages,i+2] - n.mat.sahul[16:stages,i+2] # take away washed-off colonists from Timor 
            if (length(which(n.mat.timor[,i+2] < 0)) > 0) {
              n.mat.timor[which(n.mat.timor[,i+2] < 0), i+2] <- 0   }}} # make sure no negative values in Timor
        
        ###########
        ## SAHUL ##
        ###########
        if (sum(n.mat.sahul[,i+2]) > 0) { # only project if there are people in this iteration on Sahul
          
          ## reconstruct popmat with stochastic elements
          popmat.sahul[1, ] <- pfert.trunc * rnorm(1, fert.bentley,0.05*fert.bentley) # fertility sampler
          diag(popmat.sahul[2:(stages), ]) <- ((ifelse(sum(n.mat.sahul[,i+2], na.rm=T) > sahul.K.run[i+2], surv.min, 1)) * (stoch.surv.func(Sx, Sx.sd^2))[-stages]) # survival sampler
          
          # survival (+ catastrophic mortality at 50%)
          if ((rbinom(1, 1, cat.pr)) == 1) {
            diag(popmat.sahul[2:(stages), ]) <- (0.5* (diag(popmat.sahul[2:(stages), ])))}
          popmat.sahul[stages,stages] <- 0
          
          ## project over interval
          n.mat.sahul[,i+3] <- popmat.sahul %*% n.mat.sahul[,i+2]
          
          ## did the Sahul population go extinct?
          if (sum(n.mat.sahul[,i+3]) < q.ext) {
            n.mat.sahul[,i+3] <- 0}}
        
      } # end i loop  

      ## did people make it to Sahul?
      sahul.arrive[e] <- ifelse(length(which((round(colSums(n.mat.sahul),0)) > 0)) > 0, 1, 0)
      
      # if arrived in Sahul, what was the probability that the population went extinct?
      sahul.ext[e] <- ifelse(sahul.arrive[e] == 1 & (ifelse(sahul.arrive[e] == 1, (round(colSums(n.mat.sahul),0))[length(round(colSums(n.mat.sahul),0))], NA)) == 0, 1, NA)
      
      setTxtProgressBar(pb, e)
      # if (e %% itdiv==0) print(e)
    } # end e loop
    
    # probability of arrival in Sahul
    pr.sahul.arrive.stor[p, c] <- sum(sahul.arrive)/iter
    
    # probability of arrival & extinction in Sahul
    pr.sahul.ext.stor[p, c] <- sum(sahul.ext,na.rm=T)/sum(sahul.arrive,na.rm=T)
    
    print("")
    print(paste("probability of a colonisation event = ", pr.col.vec[p], sep=""))
  } # end p loop (pr.vec)
  
  print("##################################################")
  print(paste("colonisation event size (male & female) = ", col.size.vec[c], sep=""))
  print("##################################################")
  
} # end c loop (col.size)

# plot 3D surfaces (pr.col x col.size)
# probability of arrival to Sahul
pr.sahul.arrive.stor
persp3D(x=pr.col.vec, y=col.size.vec, z=pr.sahul.arrive.stor, border="white", bty="f", resfac=2, theta = 30, phi= 20, col=grey.colors(20), image=F, contour=F, cex.axis=0.6, ticktype="detailed", xlab="Pr(colonisation event)", ylab="event size (n)", zlab="Pr(reach Sahul)", expand=0.7)

m = list(l=100, r=20, b=60, t=0)
b <- list(title="Pr(colonisation event)", tickmode="array", nticks=length(pr.col.vec), tickvals=pr.col.vec, tickfont=(list(size=18)), titlefont=list(size=23))
a <- list(title="event size (n)", tickmode="array", nticks=length(col.size.vec), tickvals=col.size.vec, tickfont=(list(size=18)), titlefont=list(size=23))
p1 <- plot_ly(y=~pr.col.vec, x=~col.size.vec, z=pr.sahul.arrive.stor, type = "contour", contours=list(start=0.05, end=0.95, size=0.05, showlabels=T, labelfont=list(size=18, color="white")), colorscale='Viridis', reversescale=T) %>%
  colorbar(title = "Pr(reach Sahul)") %>%
  layout(xaxis = a, yaxis=b, margin=m)
p1

pr.sahul.ext.stor
persp3D(x=pr.col.vec, y=col.size.vec, z=pr.sahul.ext.stor, border="white", bty="f", resfac=2, theta = 230, phi= 20, col=rev(grey.colors(20)), image=F, contour=F, cex.axis=0.6, ticktype="detailed", xlab="Pr(colonisation event)", ylab="event size (n)", zlab="Pr(Sahul extinction)", expand=0.7)

m = list(l=100, r=20, b=60, t=0)
b <- list(title="Pr(colonisation event)", tickmode="array", nticks=length(pr.col.vec), tickvals=pr.col.vec, tickfont=(list(size=18)), titlefont=list(size=23))
a <- list(title="event size (n)", tickmode="array", nticks=length(col.size.vec), tickvals=col.size.vec, tickfont=(list(size=18)), titlefont=list(size=23))
p <- plot_ly(y=~pr.col.vec, x=~col.size.vec, z=pr.sahul.ext.stor, type = "contour", contours=list(start=0.05, end=0.95, size=0.05, showlabels=T, labelfont=list(size=18, color="white")), colorscale='Viridis', reversescale=F) %>%
  colorbar(title = "Pr(Sahul extinction)") %>%
  layout(xaxis = a, yaxis=b, margin=m)
p

