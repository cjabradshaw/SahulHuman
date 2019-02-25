# Matrix operators for population models
# place in R 'Resources/R/' folder

## maximum lambda function
max.lambda <- function(x) Re((eigen(x)$values)[1]) ## where 'x' is a Leslie matrix

## Maximum r function
max.r <- function(x) log(Re((eigen(x)$values)[1])) ## where 'x' is a Leslie matrix

## Stable stage distribution
stable.stage.dist <- function(x) ((x %*% (Re((eigen(x)$vectors)[,1])))/(sum((x %*% (Re((eigen(x)$vectors)[,1]))))))[,1]

## Generation length function
R.val <- function(X,age.max) ## reproductive value (R0) where X = Leslie matrix; age.max = maximum age of females
{		
		## define the transition matrix
		T <- X[1:age.max,1:age.max]
		T[1,1:(age.max)] <- 0

		## define the fertility matrix
		F <- X[1:age.max,1:age.max]
		diag(F[2:age.max,1:(age.max-1)]) <- 0

		## define the identity matrix
		I <- matrix(data<-0,nrow<-age.max,ncol<-age.max)
		diag(I) <- 1

		## define the fundamental matrix
		library(MASS)
		N.fund <- ginv(I-T)

		## define the reproductive matrix
		R <- F %*% N.fund

		## define R0 (number of female offspring produced per female during lifetime)
		R0 <- Re((eigen(R)$values)[1])
		
		## output
		print("number of female offspring produced per female during its lifetime")
		print("_________________________________________________________________")
		print(R0)

}

## Mean generation time function
G.val <- function (X,age.max) ## where X is a Leslie Matrix
{	
		G <- (log(R.val(X,age.max)))/(log(Re((eigen(X)$values)[1])))
		print("mean generation time")
		print("____________________")
		print(G)
}
