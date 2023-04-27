library(tidyverse)
library(ggplot2)
source("ColalaNumberPsychophysics.R")

d <- read.csv("data/Cheyette-NHB.csv", header=T)
d <- subset(d, cond == "DURATION" & control=="NONE" & d$prez_num == 0)
d <- subset(d, n_guess < 30 & d$rt < 5000) # Toss crazy stuff 

# We are going to remove outliers in responses since those are likely 
# from some other process. Also the Cheyette model is not robust
# to these kinds of outliers
d <- d %>% group_by(n_shown, exp_dur) %>%
  filter(between(n_guess, quantile(n_guess,0.05), quantile(n_guess,0.95)))
ggplot(d, aes(x=n_shown, y=n_guess)) + 
	geom_jitter() + 
	facet_wrap(~exp_dur) +
	theme_bw()
	

# Loop through the data by subject and duration, and estimate psychophysical
# models for each subject. 
D <- NULL
for(q in split(d, paste0(d$uniqueid,".", d$exp_dur))) {
	if(nrow(q) < 20) next # skip conditions with too few data 
	
	# Fit information bound with certainty 
	ib <- estimation.full.fit(q$n_shown, q$n_guess, model="cheyette", plot=TRUE)
	
	# Fit simple Weber model
	weber <- estimation.fit(q$n_shown, q$n_guess, model="weber")
	
	D <- rbind(D, data.frame(n=nrow(q), 
						     uniqueid=q$uniqueid[1], 
							 exp_dur=q$exp_dur[1], 
							 ib=ib$param, 
							 ib.sd=ib$sd, 
							 ib.ll=ib$ll,
							 w=weber$param,
							 w.ll=weber$ll,
							 mad=median(abs(q$n_shown-q$n_guess)),
							 mean=mean(abs(q$n_shown-q$n_guess))))
		print(D)
}

# estimation_likelihood(q$n_shown, q$n_guess, 2.82)
# 
# L <- NULL
# for(i in 1:nrow(q)) {
# 	L <- rbind(L, data.frame(n_shown=q$n_shown[i], cheyette=estimation.likelihood(q$n_shown[i], q$n_guess[i], 2.82, model="cheyette"),weber=estimation.likelihood(q$n_shown[i], q$n_guess[i], 0.18, model="weber")))
# }


# plot the relative scores of these models as a function of
# 
plt <- ggplot(D, aes(x=ib.ll, y=w.ll)) + 
	geom_point() + 
	geom_abline(intercept=0, slope=1) +
	xlim(-200,0) + ylim(-200,0) + ## NOTE: this removes some outliers for which Cheyette does not fit very well
	theme_bw() + 
	facet_wrap(~exp_dur)

# Look at the relationship between information and W
plt <- ggplot(D, aes(x=w, y=ib)) + 
	geom_point() + 
	facet_wrap(~exp_dur)
	
# NOTE: In Cheyette & Piantadosi, we parametrically fit exp_dur as a function
# of the display time. Here, we will just visualize it
plt <- ggplot(D, aes(x=exp_dur, y=ib)) + 
	stat_summary(fun.data=mean_cl_boot, geom="errorbar", color="red") + 
	geom_point()
