source("ColalaNumberPsychophysics.R")

## Simulate and test recovery for estimation
D <- NULL
model <- "cheyette" ## NOTE: if you change, you must change the seq
for(x in rep(seq(0.1, 3.0, 0.1), 25)) {

	shown <- rep(2:25,3) 
	sim.response <- estimation.sample(shown, x, model=model)
	
	#recovery 
	fit <- estimation.full.fit(shown, sim.response, model=model, plot=TRUE)

	D <- rbind(D, data.frame(true=x, recovered=fit$param))
	print(D)
}
