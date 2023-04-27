library(tidyverse)
library(ggplot2)
source("ColalaNumberPsychophysics.R")

d <- read.csv("data/WData.csv", header=T)

# Loop through the data by subject and duration, and estimate psychophysical
# models for each subject. 
D <- NULL
for(q in split(d, paste0(d$species,".", d$study, ".", d$task))) {
	
	# Fit information bound with certainty 
	ib <- discrimination.full.fit(q$n1, q$n2, q$ntrials, q$ncorrect, model="cheyette", plot=TRUE)
	
	# Fit simple Weber model
	weber <- discrimination.fit(q$n1, q$n2, q$ntrials, q$ncorrect, model="weber")
	
	D <- rbind(D, data.frame(n=sum(q$ntrials), 
						     species=q$species[1], 
						     study=q$study[1],
						     task=q$task[1],
							 ib=ib$param, 
							 ib.sd=ib$sd, 
							 ib.ll=ib$ll,
							 w=weber$param,
							 w.ll=weber$ll,
							 mean=mean(q$ntrials/q$ncorrect)
							 ))
		print(D)
}
