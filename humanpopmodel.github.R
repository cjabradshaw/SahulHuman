## Ancient human demographic matrix projections in Australia
## Corey Bradshaw
## February 2019

## Remove everything
rm(list = ls())

## libraries
library(boot)
library(tcltk)

## source
source("matrixOperators.r")

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

# average forager-horticuluralist
a1.fh <- 0.418; b1.fh <- 1.657; a2.fh <- 0.012; a3.fh <- 3.65e-04; b3.fh <- 0.074
h.x.fh <- a1.fh * exp(-b1.fh*x) + a2.fh + a3.fh * exp(b3.fh * x)
plot(x,h.x.fh,pch=19,type="l")
plot(x,log(h.x.fh),pch=19,type="l")
l.x.fh <- exp((-a1.fh/b1.fh) * (1 - exp(-b1.fh*x))) * exp(-a2.fh * x) * exp(a3.fh/b3.fh * (1 - exp(b3.fh * x)))
plot(x,l.x.fh,type="l")

# NT Aboriginal
a1.nta <- 0.242; b1.nta <- 1.031; a2.nta <- 0.000; a3.nta <- 7.13e-04; b3.nta <- 0.063
h.xnta <- a1.nta * exp(-b1.nta*x) + a2.nta + a3.nta * exp(b3.nta * x)
plot(x,h.x.nta,pch=19,type="l")
plot(x,log(h.x.nta),pch=19,type="l")
l.x.nta <- exp((-a1.nta/b1.nta) * (1 - exp(-b1.nta*x))) * exp(-a2.nta * x) * exp(a3.nta/b3.nta * (1 - exp(b3.nta * x)))
plot(x,l.x.nta,type="l")

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

# average forager-horticulturalist
lx.fh <- round(init.pop*l.x.fh,0)
len.lx.fh <- length(lx.fh)
dx.fh <- lx.fh[1:(len.lx.fh-1)]-lx.fh[2:len.lx.fh]
qx.fh <- dx.fh/lx.fh[1:(length(lx.fh)-1)]
Sx.fh <- 1 - qx.fh
sx.fh <- lx.fh[2:len.lx.fh]/lx.fh[1:(len.lx.fh-1)]
mx.fh <- 1 - sx.fh
Lx.fh <- (lx.fh[1:(len.lx.fh-1)] + lx.fh[2:len.lx.fh])/2
ex.fh <- rev(cumsum(rev(Lx.fh)))/lx[-len.lx.fh]
ex.avg.fh <- ex.fh + x[-len.lx.fh]

# average NT aboriginal
lx.nta <- round(init.pop*l.x.nta,0)
len.lx.nta <- length(lx.nta)
dx.nta <- lx.nta[1:(len.lx.nta-1)]-lx.nta[2:len.lx.nta]
qx.nta <- dx.nta/lx.nta[1:(length(lx.nta)-1)]
Sx.nta <- 1 - qx.nta
sx.nta <- lx.nta[2:len.lx.nta]/lx.nta[1:(len.lx.nta-1)]
mx.nta <- 1 - sx.nta
Lx.nta <- (lx.nta[1:(len.lx.nta-1)] + lx.nta[2:len.lx.nta])/2
ex.nta <- rev(cumsum(rev(Lx.nta)))/lx[-len.lx.nta]
ex.avg.nta <- ex.nta + x[-len.lx.nta]

# set SD for Sx
Sx.sd <- 0.05 # can set to any value

par(mfrow=c(2,1))
plot(x[-1], Sx, pch=19, type="l", xlab="age (years)", ylab="Sx")
plot(x[-1], ex, pch=19, type="l", xlab="age (years)", ylab="life expectancy")
par(mfrow=c(1,1))
surv.out <- data.frame(x[-1],Sx,ex)
colnames(surv.out)[1] <- "x"

# plot average hunter-gatherer, average forager-horticulturalist, and NT aboriginal together
plot(x[-1], Sx, pch=19, type="l", xlab="age (years)", ylab="Sx") # average hunter-gatherer
lines(x[-1], Sx.fh, lty=2, lwd=3, col="red") # average forager-horticulturalist
lines(x[-1], Sx.fh, lty=3, lwd=3, col="green") # NT aboriginal


# fertility (Walker et al. 2006)
primiparity.walker <- c(17.7,18.7,19.5,18.5,18.5,18.7,25.7,19,20.5,18.8,17.8,18.6,22.2,17,16.2,18.4)
prim.mean <- round(mean(primiparity.walker),0)
prim.lo <- round(quantile(primiparity.walker,probs=0.025),0)
prim.hi <- round(quantile(primiparity.walker,probs=0.975),0)
print(c(prim.lo, prim.mean, prim.hi))
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

## different Siler model parameters' influence on base matrix
popmat.fh <- popmat.orig # forager-horticulturalist
diag(popmat.fh[2:stages,]) <- Sx.fh
popmat.fh[stages,stages] <- 0 # Sx[stages-1]
max.lambda(popmat.fh)

popmat.nta <- popmat.orig # NT Aboriginal
diag(popmat.nta[2:stages,]) <- Sx.nta
popmat.nta[stages,stages] <- 0 # Sx[stages-1]
max.lambda(popmat.nta)

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

## NORTH OF SAHUL (between 0 and 14 degrees S Latitude)
# net primary productivity (NPP; kg C/yr/m2)
# ALL SAHUL FIRST
npp.sah <- read.table("ClimateSahul_Npp.csv", header=T, sep=",") 
## JUST NORTH (0-14 DEGREES)
npp.nsah <- subset(npp.sah, Lat..degree. <= 0 & Lat..degree. >= -14)
dim(npp.nsah)

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
abline(v=55000,lty=3,lwd=3)
abline(v=65000,lty=3,lwd=3)
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
SahulN.D <- 0.6572452 # density according to Tallavaara et al. 2018 from 0-14 S Lat / 129.2-151 Lon
Aus.A <- 7692000  # today's terrestrial size of Australia (km^2)
NG.A <- 786000 # today's new guinea area (km^2)
NAus.prop <- 27.75/695.5 # proportion of Australia north of 14 degrees S
SahulN.A <- 1.25 * ((Aus.A*NAus.prop) + NG.A) # assuming 25% bigger than today
SahulN.K <- SahulN.A * SahulN.D
SahulN.K

K.dat$K.lo <- round(start.K * 1/K.dat$npp.lo.sc, 0)
K.dat$K.med <- (K.dat$npp.med.sc/K.dat$npp.lo.sc) * K.dat$K.lo
K.dat$K.up <- (K.dat$npp.up.sc/K.dat$npp.lo.sc) * K.dat$K.lo

K.5565 <- subset(K.dat, year <= 65000 & year >= 55000)
min.med5565 <- min(K.5565$K.med)
K.5565[which(K.5565$K.med == min.med5565),]
max.med5565 <- max(K.5565$K.med)
maxK.dat <- K.5565[which(K.5565$K.med == max.med5565),]
range(maxK.dat$year)



plot(K.dat$year, K.dat$K.med, type="l", xlab="year", ylab="relative K", lty=1, lwd=2, ylim=c(min(K.dat$K.lo), max(K.dat$K.up)))
lines(K.dat$year, K.dat$K.lo, lty=2, col="red")
lines(K.dat$year, K.dat$K.up, lty=2, col="red")
abline(v=55000,lty=3,lwd=3)
abline(v=65000,lty=3,lwd=3)


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
## INITIAL FOUNDING POPULATION SIZE REQUIRED
## stochastic projection 
## resampling proportion breeding & fertility
## set x% SD on Sx values (beta distribution)
## catastrophe sampler
#####################################################################################

# initial colonisation window
col.old <- 60000
col.yng <- 50000

# generations to project
gen.proj <- 100

# iterations per founding population size
iter <- 10000

# set quasi-extinction threshold
min.thresh <- 50

# set founding population vector
found.min <- 875
found.max <- 1000
found.vec <- seq(found.min, found.max, by=25)

# Pr(ext) storage vector
pr.ext.vec <- rep(0,length(found.vec))

for (f in 1:length(found.vec)) { # founding population size loop
  
  #################
  ## project
  ## set time limit for projection in 1-yr increments
  yr.st <- round(runif(iter,col.yng,col.old), 0) 
  #************************
  yr.end <- yr.st - round(gen.proj*gen.l, 0) # set projection end date
  #************************
  t <- (yr.st - yr.end)
  
  ## initial population vector
  popmat <- popmat.orig
  init.vec <- stable.stage.dist(popmat) * found.vec[f] # stable stage distribution x founding population size
  
  ## set population storage matrices
  n.mat <- matrix(0,nrow=stages,ncol=(t[1]+1))
  n.mat[,1] <- init.vec
  
  ## iterate projection
  itdiv <- iter/100
  
  # set storage matrices & vectors
  n.sums.mat <- matrix(data = 0, nrow = iter, ncol = (t[1]+1))
  #r.mat <- matrix(data = 0, nrow = iter, ncol = t)
  
  # progress counter
  pb <- txtProgressBar(min=1, max=iter, style=3)
  
  for (e in 1:iter) { # e iterations loop
    
    yr.st.sub <- which(K.dat$year == yr.st[e])
    npp.sc.st <- runif(1, K.dat$npp.lo.sc[yr.st.sub], K.dat$npp.up.sc[yr.st.sub])
    K.run <- round((npp.sc.st / K.dat$npp.med.sc[yr.st.sub]) * K.dat$K.med[yr.st.sub:(yr.st.sub+t[e])], 0) # this iteration's realised K series
    yr.vec.run <- yr.st[e]:yr.end[e]
    #plot(yr.vec.run, K.run, type="l", xlab="year",ylab="K (N)", xlim=c(55000-t[1],65000))
    
    ## reset popmat to original values
    popmat <- popmat.orig
    
    ## set up projection loop
    for (i in 1:t[e]) { # i projection loop
      
      ## reconstruct popmat with stochastic elements
      popmat[1, ] <- pfert.trunc * rnorm(1, fert.bentley,0.05*fert.bentley) # fertility sampler
      diag(popmat[2:(stages), ]) <- ((ifelse(sum(n.mat[,i], na.rm=T) > K.run[i], surv.min, 1)) * (stoch.surv.func(Sx, Sx.sd^2))[-stages]) # survival sampler
      
      # survival (+ catastrophic mortality at 50%)
      if ((rbinom(1, 1, cat.pr)) == 1) {
        diag(popmat[2:(stages), ]) <- (0.5* (diag(popmat[2:(stages), ])))}
      popmat[stages,stages] <- 0
      
      ## project over interval
      n.mat[,i+1] <- popmat %*% n.mat[,i] 
      
      #print(i)
    } # end i projection-length loop
    
    #r.mat[e,] <- r.stoch
    n.sums.mat[e,] <- as.vector(colSums(n.mat))
    
    setTxtProgressBar(pb, e)
    #if (e %% itdiv==0) print(e) 
  } # end e iterations loop
  
  # N confidence limits
  n.mn <- apply(n.sums.mat, MARGIN=2, mean, na.rm=T) # mean over all iterations
  n.up <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  n.lo <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
  n.min <- apply(n.sums.mat, MARGIN=1, min, na.rm=T) # minimum over all projected years
  q.ext.vec <- ifelse(n.min < (min.thresh/2), 1, 0)
  pr.ext.vec[f] <- pr.ext <- sum(q.ext.vec)/iter
  
  # plot
  #yrs.vec <- 0:t[1]
  #plot(yrs.vec, log10(n.mn), type="l", xlab="years into future",ylab="log10(N)", ylim=c(log10(min(n.lo)),log10(max(n.up))))
  #lines(yrs.vec, log10(n.lo), lty=2, col="red")
  #lines(yrs.vec, log10(n.up), lty=2, col="red")
  
  print("##################################################")
  print(paste("founding population size = ", found.vec[f], sep=""))
  print(paste("Pr(quasi-extinction) = ", pr.ext.vec[f], sep=""))
  print("##################################################")
  
} # end founding population size (f) loop

plot(found.vec, pr.ext.vec, type="l", xlab="founding population minimum size (N)", ylab="Pr(min N quasi-extinct)") # females only
plot(2*found.vec, pr.ext.vec, type="l", xlab="founding population minimum size (N)", ylab="Pr(min N quasi-extinct)") # all only
abline(h=0.2,lty=2)
abline(h=0.1,lty=2)

prminext.out <- data.frame(2*found.vec,pr.ext.vec)
colnames(prminext.out) <- c("founderN","PrExt")



#####################################################################################
## FREQUENCY OF PEOPLE ENTERING AUSTRALIA (CONSTANT AMOUNT)
## stochastic projection 
## resampling proportion breeding & fertility
## set x% SD on Sx values (beta distribution)
## catastrophe sampler
#####################################################################################

# initial colonisation window
col.old <- 65000
col.yng <- 55000

# generations to project
gen.proj <- 100

# iterations per founding population size
iter <- 10000

# set quasi-extinction threshold
min.thresh <- 50

## set interval vector
int.vec <- seq(110,120,10)

## set female sample-size window for minimum viable population size (from previous step)
mvp.lo <- round(1300/2, 0)
mvp.hi <- round(1550/2, 0)

# storage vector
pr.ext.vec <- rep(0,length(int.vec))

for (s in 1:length(int.vec)) {
  
  #################
  ## project
  ## set time limit for projection in 1-yr increments
  yr.st <- round(runif(iter,col.yng,col.old), 0) 
  #************************
  yr.end <- yr.st - round(gen.proj*gen.l, 0) # set projection end date
  #************************
  t <- (yr.st - yr.end)
  
  ## initial population vector
  popmat <- popmat.orig
  found.st <- round(runif(1, min=mvp.lo/10, max=mvp.hi/10),0)
  
  init.vec <- stable.stage.dist(popmat) * min.thresh/2 # if not greater
  if (found.st > min.thresh/2) {
    init.vec <- stable.stage.dist(popmat) * found.st} # stable stage distribution x founding population size
  
  ## set population storage matrices
  n.mat <- matrix(0,nrow=stages,ncol=(t[1]+1))
  n.mat[,1] <- init.vec
  
  ## iterate projection
  itdiv <- iter/100
  
  # set storage matrices & vectors
  n.sums.mat <- matrix(data = 0, nrow = iter, ncol = (t[1]+1))
  
  # progress counter
  pb <- txtProgressBar(min=1, max=iter, style=3)
  
  for (e in 1:iter) { # e iterations loop
    
    yr.st.sub <- which(K.dat$year == yr.st[e])
    npp.sc.st <- runif(1, K.dat$npp.lo.sc[yr.st.sub], K.dat$npp.up.sc[yr.st.sub])
    K.run <- round((npp.sc.st / K.dat$npp.med.sc[yr.st.sub]) * K.dat$K.med[yr.st.sub:(yr.st.sub+t[e])], 0) # this iteration's realised K series
    yr.vec.run <- yr.st[e]:yr.end[e]
    #plot(yr.vec.run, K.run, type="l", xlab="year",ylab="K (N)", xlim=c(55000-t[1],65000))
    
    ## reset popmat to original values
    popmat <- popmat.orig
    
    ## add to n.mat for first years of invasion
    add.n <- round(runif(int.vec[s]*10, min=mvp.lo/10, max=mvp.hi/10),0)
    use.sub <- rep(0,int.vec[s]*10)
    use.sub[seq(int.vec[s], (int.vec[s]*10), (int.vec[s]+1))] <- 1
    add.n <- ifelse(use.sub==1, add.n, 0)
    
    ## set up projection loop
    for (i in 1:t[e]) { # i projection loop
      
      ## reconstruct popmat with stochastic elements
      popmat[1, ] <- pfert.trunc * rnorm(1, fert.bentley,0.05*fert.bentley) # fertility sampler
      diag(popmat[2:(stages), ]) <- ((ifelse(sum(n.mat[,i], na.rm=T) > K.run[i], surv.min, 1)) * (stoch.surv.func(Sx, Sx.sd^2))[-stages]) # survival sampler
      
      # survival (+ catastrophic mortality at 50%)
      if ((rbinom(1, 1, cat.pr)) == 1) {
        diag(popmat[2:(stages), ]) <- (0.5* (diag(popmat[2:(stages), ])))}
      popmat[stages,stages] <- 0
      
      ## project over interval
      n.mat[,i+1] <- popmat %*% n.mat[,i] 
      if (i <= length(add.n)) {
        n.mat[,i+1] <- n.mat[,i+1] + as.vector((stable.stage.dist(popmat) * add.n[i]))} # add new founders for length of add.n vector
      
      #print(i)
    } # end i projection-length loop
    
    #r.mat[e,] <- r.stoch
    n.sums.mat[e,] <- as.vector(colSums(n.mat))
    
    setTxtProgressBar(pb, e)
    #if (e %% itdiv==0) print(e) 
  } # end e iterations loop
  
  # probability of extinction
  pr.ext.vec[s] <- sum(ifelse((apply(n.sums.mat, MARGIN=1, min, na.rm=T)) < (min.thresh/2), 1, 0))/iter
  
  # N confidence limits
  n.mn <- apply(n.sums.mat, MARGIN=2, mean, na.rm=T) # mean over all iterations
  n.up <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  n.lo <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  

  print("##################################################")
  print(paste("interval of arrival (years) = ", int.vec[s], sep=""))
  print(paste("Pr(quasi-extinction) = ", pr.ext.vec[s], sep=""))
  print("##################################################")
  
} # end s loop  

plot(int.vec, pr.ext.vec, type="l", xlab="interval of arrival", ylab="Pr(q-extinction)")
abline(h=0.1,lty=2)




#####################################################################################
## FREQUENCY OF PEOPLE ENTERING AUSTRALIA (RANDOM SERIES)
## stochastic projection 
## resampling proportion breeding & fertility
## set x% SD on Sx values (beta distribution)
## catastrophe sampler
#####################################################################################

# initial colonisation window
col.old <- 65000
col.yng <- 55000

# generations to project
gen.proj <- 100

# iterations per founding population size
iter <- 10000

# set quasi-extinction threshold
min.thresh <- 50

## set female sample-size window for minimum viable population size (from previous step)
mvp.lo <- round(1300/2, 0)
mvp.hi <- round(1550/2, 0)

#################
## project
## set time limit for projection in 1-yr increments
yr.st <- round(runif(iter,col.yng,col.old), 0) 
#************************
yr.end <- yr.st - round(gen.proj*gen.l, 0) # set projection end date
#************************
t <- (yr.st - yr.end)

## initial population vector
popmat <- popmat.orig

## iterate projection
itdiv <- iter/100

# interval vector (years)
interval.vec <- seq(10, 300, 10)
linterv.vec <- length(interval.vec)

# store pr.ext for each interval value
pr.ext.int.vec <- rep(0,linterv.vec)

for (v in 1:linterv.vec) {
  
  # set storage matrices & vectors
  n.sums.mat <- matrix(data = 0, nrow = iter, ncol = (t[1]+1))
  
  # progress counter
  pb <- txtProgressBar(min=1, max=iter, style=3)
  
  for (e in 1:iter) { # e iterations loop
    
    yr.st.sub <- which(K.dat$year == yr.st[e])
    npp.sc.st <- runif(1, K.dat$npp.lo.sc[yr.st.sub], K.dat$npp.up.sc[yr.st.sub])
    K.run <- round((npp.sc.st / K.dat$npp.med.sc[yr.st.sub]) * K.dat$K.med[yr.st.sub:(yr.st.sub+t[e])], 0) # this iteration's realised K series
    yr.vec.run <- yr.st[e]:yr.end[e]
    
    ## reset popmat to original values
    popmat <- popmat.orig
    
    ## add to n.mat for first years of invasion
    tot.init.pop <- round(runif(1,mvp.lo,mvp.hi),0)
    pop.n <- 0 # series of additions
    pop.n[1] <- round(runif(1,min.thresh/2,tot.init.pop/2))
    j <- 1
    while (sum(pop.n) < tot.init.pop-(min.thresh/2)) {
      j = j+1
      pop.n[j] <- round(runif(1,min.thresh/2,(tot.init.pop - sum(pop.n))))
    }
    if (sum(pop.n) > tot.init.pop-(min.thresh/2)) {
      fin.sub <- length(pop.n) + 1
      pop.n[fin.sub] <- tot.init.pop - sum(pop.n)
    }
    if (pop.n[length(pop.n)] == 0) {
      pop.n <- pop.n[-length(pop.n)]
    }
    
    # create interval between additions
    interv.samp <- c(0, (round(rnorm(length(pop.n)-1, mean=interval.vec[v], sd=(0.1*interval.vec[v])), 0))) + 1
    add.n <- rep(0,sum(interv.samp))
    add.n[cumsum(interv.samp)] <- pop.n
    
    ## initial population vector
    init.vec <- stable.stage.dist(popmat) * add.n[1] # stable stage distribution x founding population size
    
    ## remove first entry of add.n (already accounted for in init.vec)
    add.n <- add.n[-1]
    
    ## set population storage matrices
    n.mat <- matrix(0,nrow=stages,ncol=(t[1]+1))
    n.mat[,1] <- init.vec
    
    ## set up projection loop
    for (i in 1:t[e]) { # i projection loop
      
      ## reconstruct popmat with stochastic elements
      popmat[1, ] <- pfert.trunc * rnorm(1, fert.bentley,0.05*fert.bentley) # fertility sampler
      diag(popmat[2:(stages), ]) <- ((ifelse(sum(n.mat[,i], na.rm=T) > K.run[i], surv.min, 1)) * (stoch.surv.func(Sx, Sx.sd^2))[-stages]) # survival sampler
      
      # survival (+ catastrophic mortality at 50%)
      if ((rbinom(1, 1, cat.pr)) == 1) {
        diag(popmat[2:(stages), ]) <- (0.5* (diag(popmat[2:(stages), ])))}
      popmat[stages,stages] <- 0
      
      ## project over interval
      n.mat[,i+1] <- popmat %*% n.mat[,i] 
      if (i <= length(add.n)) {
        n.mat[,i+1] <- n.mat[,i+1] + as.vector((stable.stage.dist(popmat) * add.n[i]))} # add new founders for length of add.n vector
      
    } # end i projection-length loop
    
    n.sums.mat[e,] <- as.vector(colSums(n.mat))
    
    setTxtProgressBar(pb, e)
    #if (e %% itdiv==0) print(e) 
  } # end e iterations loop
  
  # N confidence limits
  n.mn <- apply(n.sums.mat, MARGIN=2, mean, na.rm=T) # mean over all iterations
  n.up <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  n.lo <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
  n.min <- apply(n.sums.mat, MARGIN=1, min, na.rm=T) # minimum over all projected years
  q.ext.vec <- ifelse(n.min < (min.thresh/2), 1, 0)
  pr.ext.int.vec[v] <- sum(q.ext.vec)/iter
  

  print("--------------")
  print("--------------")
  print(paste("interval = ", interval.vec[v], " years (complete)", sep=""))
  print(paste("Pr(quasi-extinction) = ", pr.ext.int.vec[v], sep=""))
  print("--------------")
  print("--------------")
  
} # end v interval loop

plot(interval.vec, pr.ext.int.vec,type="l", xlab="average interval of introduction (years)", ylab="Pr(quasi-extinction)")
abline(y=0.1, lty="2")

