
## Fitting functions:

## Information bound 
## Information bound with lapse
## Weber
## Weber with lapse
## Weber nonlinear
## Weber nonlienar with laspe

# setwd("/home/piantado/Desktop/Science/Projects/Number-Stan/c++")

library(Rcpp)

sourceCpp("R-Test.cpp")

WEBER_RATE <- 1.0 # rate on the prior for Weber

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
estimation.likelihood <- function(shown, response, param, model="weber") {
	if(model=="weber") {
		stopifnot(param > 0); stopifnot(length(param)==1)
		## This weber likelihood will integrate from [n-0.5,n+0.5] instead of using the density, 
		## because we really need the true probability that it assigns to each discrete outcome
		## (not the relative probability)
		sum(mapply(logminusexp, pnorm(response+0.5,shown, shown*param, log=T),
								pnorm(response-0.5,shown, shown*param, log=T)))
	}
	else if(model == "cheyette") {
		estimation_likelihood(shown,response, param) 
	}
	else {
		print("Invalid model", model)
	}
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
estimation.fit <- function(shown, response, model="weber") {
	if(model == "weber") {
		params <- optimize(interval=c(0,5),  
					function(w) { 
							-estimation.likelihood(shown,response,w,model=model) - dexp(w, WEBER_RATE, log=T) 
					 })
						
		
		# recompute the likelihood only (score w/o prior)
		list(param=params$minimum, ll=estimation.likelihood(shown,response,params$minimum,model=model), model=model)	
	}
	else if (model == "cheyette") {
		ib <- estimation_fit_ib_bound(shown,response)
		list(param=ib, ll=estimation.likelihood(shown,response,ib,model=model), model=model)
	}
	else {
		print("Invalid model", model)
	}
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
estimation.full.fit <- function(shown, response, model="weber", plot=TRUE) {

	mode <- estimation.fit(shown, response, model=model)$param
	
	# now evaluate at a bunch of places to compute the fit 
	parms <- seq(0, 3*mode, length.out=100) 
	lls <- sapply(parms, function(ib){ estimation.likelihood(shown,response,ib, model=model)})
	z <- logsumexp(lls)
	lls.m <- exp(lls - z) # and normalize
	
	# a density we will fit to this thing -- in this case, renormalized on x
	lg <- function(x, m, s) {
		v <- dnorm(x, m, s, log=T) # normalize only on ibs
		(v-logsumexp(v))
	}
	
	# Now fit a normal to these points. We are going to fit in terms of MSE
	# which may not be ideal but works pretty well. 
	# we will start with the mean on the right value 
	params <- optimize(interval=c(0,5),  function(x) { 
		sum((lls.m - exp(lg(parms, mode, x)))^2)
	})
		
	fitsd <- params$minimum 
		
	if(plot) {
		plot(parms, lls.m, type="l")
		lines(parms, exp(lg(parms, mode, fitsd)), col=2  )
		points(mode, 0, col=3)
	}
	
	list(param=mode, sd=fitsd, ll=estimation.likelihood(response, shown, mode, model=model), model=model )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
estimation.sample <- function(shown, param, model="weber") {
	if(model=="weber") {
		round(rnorm(length(shown), shown, shown*param))
	}
	else if(model == "cheyette") {
		estimation_sample(shown, param)
	}
	else {
		print("Invalid model", model)
	}
	
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# =========================================================
discrimination.likelihood <- function(n1, n2, accuracy, param, model="weber") {
	if(model=="weber") {
		stopifnot(param > 0); stopifnot(length(param)==1)
		pcorrect <- 1.0-pnorm(0, abs(n1-n2), param * sqrt(n1^2+n2^2), log=F)
		sum(ifelse(accuracy,log(pcorrect),log(1.0-pcorrect)))
	}
	else if(model == "cheyette") {
		discrimination_likelihood(n1,n2,accuracy,param)
	}
	else {
		print("Invalid model", model)
	}
}

# =========================================================
discrimination.fit <- function(n1, n2, accuracy, model="weber") {
	if(model=="weber") {
		best <- optimize(interval=c(0.001,5.0),  
			function(w) { 
				-discrimination.likelihood(n1,n2,accuracy,w,model=model) - dexp(w, WEBER_RATE, log=T) 
		})	
	
		# recompute the likelihood only (score w/o prior)
		list(param=best$minimum, ll=discrimination.likelihood(n1,n2,accuracy,best$minimum), model=model)
	}
	else if(model == "cheyette") {
		ib <- discrimination_fit_ib_bound(n1,n2,accuracy)
		list(param=ib, ll=discrimination.likelihood(n1,n2,accuracy,ib,model=model), model=model)	
	}
	else {
		print("Invalid model", model)
	}
}

# =========================================================
discrimination.full.fit <- function(n1, n2, accuracy, model="weber", plot=TRUE){

	mode <- discrimination.fit(n1, n2, accuracy, model=model)$param
	
	# now evaluate at a bunch of places to compute the fit 
	parms <- seq(1e-3, 3*mode, length.out=100) 
	lls <- sapply(parms, function(x){ discrimination.likelihood(n1, n2, accuracy, x, model=model)})
	z <- logsumexp(lls)
	lls.m <- exp(lls - z) # and normalize
	
	# a density we will fit to this thing -- in this case, renormalized on x
	lg <- function(x, m, s) {
		v <- dnorm(x, m, s, log=T) # normalize only on ibs
		(v-logsumexp(v))
	}
	
	# Now fit a normal to these points. We are going to fit in terms of MSE
	# which may not be ideal but works pretty well. 
	# we will start with the mean on the right value 
	params <- optimize(interval=c(0,5),  function(x) { 
		sum((lls.m - exp(lg(parms, mode, x)))^2)
	})
		
	fitsd <- params$minimum 
		
	if(plot) {
		plot(parms, lls.m, type="l")
		lines(parms, exp(lg(parms, mode, fitsd)), col=2  )
		points(mode, 0, col=3)
	}
	
	list(param=mode, sd=fitsd, ll=discrimination.likelihood(n1, n2, accuracy, mode, model=model), model=model )

}

# =========================================================
discrimination.sample <- function(n1, n2, param, model="weber") {
	if(model=="weber") {
		stopifnot(param > 0); stopifnot(length(param)==1)
		pcorrect <- 1.0-pnorm(0, abs(n1-n2), w * sqrt(n1^2+n2^2), log=F)
		rbinom(length(n1), 1, pcorrect)
	}
	else if(model == "cheyette") {
		discrimination_sample(n1,n2, param)
	}
	else {
		print("Invalid model", model)
	}	
}

# # Find a gamma which is close to this likelihood
# estimation.get.likelihood.params <- function(shown,response, plot=FALSE) {
# 	
# 	# We will fit hte mode and center the normal on it
# 	mode <- estimation_fit_ib_bound(shown,response) 
# 
# 	ibs <- seq(0.0, 10.0, length.out=100) # compute on these points to fit a normal 
# 	lls <- sapply(ibs, function(ib){ estimation_likelihood(shown,response,ib)})
# 	z <- logsumexp(lls)
# 	lls.m <- exp(lls - z) # and normalize	
# 	
# 	# a log gamma density renormalized on ibs
# 	lg <- function(x, m, s) {
# 		v <- dnorm(x, m, s, log=T) # normalize only on ibs
# 		(v-logsumexp(v))
# 	}
# 	
# 	# Now fit a normal to these points. We are going to fit in terms of MSE
# 	# which may not be ideal but works pretty well. 
# 	# we will start with the mean on the right value 
# 	params <- optimize(interval=c(0,5),  function(x) { 
# 			ly <- lg(ibs, mode, x)
# 			sum((lls.m - exp(ly-logsumexp(ly)))^2)
# 	})
# 		
# 	fitsd <- params$minimum 
# 		
# 	if(plot) {
# 		plot(ibs, lls.m, type="l")
# 		lines(ibs, exp(lg(ibs, mode, fitsd)), col=2  )
# 		points(mode,0, color=3)
# 	}
# 	
# 	list(mode=mode, sd=fitsd, ll=estimation_likelihood(response, shown, mode) )
# }

# 
# # Find a gamma which is close to this likelihood
# discrimination.get.likelihood.params <- function(shown,response,accuracy,plot=FALSE) {
# 	
# 	# We will fit hte mode and center the normal on it
# 	mode <- discrimination_fit_ib_bound(shown,response,accuracy) 
# 
# 	ibs <- seq(0, mode+2, length.out=100) # compute on these points to fit a normal 
# 	lls <- sapply(ibs, function(ib){ discrimination_likelihood(shown,response,accuracy,ib)})
# 	z <- logsumexp(lls)
# 	lls.m <- exp(lls - z) # and normalize	
# 	
# 	if(is.nan(max(lls))) {
# 		return(list(mode=mode, sd=NA, ll=discrimination_likelihood(response, shown, accuracy, mode) ))
# 	}
# 	
# 	# a log gamma density renormalized on ibs
# 	lg <- function(x, m, s) {
# 		v <- dnorm(x, m, s, log=T) # normalize only on ibs
# 		(v-logsumexp(v))
# 	}
# 	
# 	# Now fit a normal to these points. We are going to fit in terms of MSE
# 	# which may not be ideal but works pretty well. 
# 	# we will start with the mean on the right value 
# 	params <- optimize(interval=c(0,5),  function(x) { 
# 			ly <- lg(ibs, mode, x)
# 			sum((lls.m - exp(ly-logsumexp(ly)))^2)
# 	})
# 		
# 	fitsd <- params$minimum 
# 		
# 	if(plot) {
# 		plot(ibs, lls.m, type="l")
# 		lines(ibs, exp(lg(ibs, mode, fitsd)), col=2  )
# 		points(mode,0, color=3)
# 	}
# 	
# 	list(mode=mode, sd=fitsd, ll=discrimination_likelihood(response, shown, accuracy, mode) )
# }



# cauchy.noise.likelihood <- function(shown, response, s=10) {	
# 	sum(mapply(logminusexp, pcauchy(response+0.5,0,s, log=T),
# 							pcauchy(response-0.5,0,s, log=T)))
# }
