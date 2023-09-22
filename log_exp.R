
# Darrah et al. 2018
#Supporting Information Appendix S3

# R code for analysis of real data set with Bayesian indicator variable model selection

##read in nest data, which contains one row per observation per #nest and the following columns:##
#Nest.ID : factor with unique codes for each nest
#Interval: length of nest check (exposure) interval
#Status: nest status at each nest check. 1 = active, 2 = #depredated, 3 = flooded, 4 = abandoned, 5 = unknown cause #of failure
#Veg.Cover: numeric column with cover categories 0-3
#Habitat: factor column with habitat types
#Exclosed: numeric column with 0 = no exclosure and 1 = #exclosure
#Nest.Initiation: numeric column with Julian day of nest #initiation

# nestdata <- read.csv("allTheData.csv")

#remove date of discovery from each nest
#IMPORTANT: this method requires that the data be sorted
#by nest and then by date
#equals vector: is the previous observation of the same nest?
equals <- rep(NA, length(nestdata$Nest.ID))
for(i in 2:length(nestdata$Nest.ID)){
  equals[i] <- nestdata$Nest.ID[i] == nestdata$Nest.ID[i-1]
}
equals[1] <- FALSE
nestdata$equals <- equals
nestdata <- subset(nestdata, equals=="TRUE")

#define fate matrix
n <- length(nestdata$Interval)
Surv <- rep(0, n)
Aban <- rep(0, n)
Flood <- rep(0, n)
Pred <- rep(0, n)

for (i in 1:n){
  Surv[i][nestdata$Status[i] == 1 | nestdata$Status[i] == 6] <- 1
  Aban[i][nestdata$Status[i] == 4] <- 1
  Pred[i][nestdata$Status[i] == 2] <- 1
  Flood[i][nestdata$Status[i] == 3] <- 1
}

Fate <- cbind(Surv, Aban, Pred, Flood)

#change other/unknown fates from row of 0s(which causes JAGS to #crash) to NAs
for (i in 1:length(Fate[,1])){
  if (sum(Fate[i,]) == 0) {Fate[i,] <- NA}
}

#reduce vegetation cover categories from 4 levels to 3
veg.cover.low <- rep(0, n) #1-20% vegetation cover
veg.cover.mod <- rep(0, n) #> 20% vegetation cover
for (i in 1:n){
  veg.cover.low[i][nestdata$Veg.Cover[i] == 1] <- 1
  veg.cover.mod[i][nestdata$ Veg.Cover[i] == 2 | nestdata$Veg.Cover[i] == 3] <- 1
}

#define habitat covariate codes
#lump rare interior habitats with Dune
nestdata$Habitat[nestdata$Habitat == "Blowout" | nestdata$Habitat == "Dry bayside flats"] <- "Dune"
dune <- rep(0, n)
overwash <- rep(0, n)
dune[][nestdata$Habitat[] == "Dune"] <- 1
overwash[][nestdata$Habitat[] == "Overwash"] <- 1

#create BUGS text file

sink("BUGS_model_subset_select.txt ") 
cat("
model {
  #HyperPriors
   for (i in 1:nSite){
    eta.p.s[i] ~ dnorm(0, tau.p)   #site effect for predation
    eta.a.s[i] ~ dnorm(0, tau.a)   #site effect for abandonment
    eta.f.s[i] ~ dnorm(0, tau.f)   #site effect flooding   
    }
  #Precision hyperparameters
    tau.p <- 1/(pow(sigma.p,2)) #predation
    tau.a <- 1/(pow(sigma.a,2)) #abandonment
    tau.f <- 1/(pow(sigma.f,2)) #flooding

  #Standard Deviation hyperparameters
    sigma.p ~ dunif(0,50)  
    sigma.a ~ dunif(0,50)  
    sigma.f ~ dunif(0,50)

  #Model Priors for total model variance
    tau.V ~ dgamma(3.29, 7.8) 
#number of parameters entering model (sum of inclusion   #variables)
   K <- (w.p.ex 
   	+ w.p.veg
	   + w.a.hab
	   + w.a.ex
	   + w.a.init
    	)
#total model variance
  tau.model <- tau.V*K  
    
 # Priors for predation sources of mortality
    alpha.p ~ dnorm(0, 0.001) #intercept
    beta.p.ex ~ dnorm(0, tau.model)  #exclosure slope
    beta.p.lowveg ~ dnorm(0, tau.model) #veg slope 1
    beta.p.modveg ~ dnorm(0, tau.model) #veg slope 2

 #Priors for predation inclusion variables
    w.p.ex ~ dbern(0.5)  #exclosure slope
    w.p.veg ~ dbern(0.5) #vegetation slopes

 #priors for abandonment
    alpha.a ~ dnorm(0, 0.001)  #intercept
    beta.a.ex ~ dnorm(0, tau.model)  #exclosure slope
    beta.a.init ~ dnorm(0, tau.model) #init. date
    beta.a.dune ~ dnorm(0, tau.model) #habitat - dune
    beta.a.over ~ dnorm(0, tau.model) #habitat - overwash

  #priors for abandonment inclusion variables
    w.a.ex ~ dbern(0.5) #exclosure
    w.a.init ~ dbern(0.5) #nest initiation date
    w.a.hab ~ dbern(0.5) #habitat

   #priors for flooding
    alpha.f ~ dnorm(0,0.001)  #intercept
   
# Likelihood
  for (i in 1:n) {
   #linear predictors (Equation 2)
   #flooding: 
    ctf[i] <- exp(alpha.f + eta.f.s[site[i]])

   #predation:
    ctp[i] <- exp(alpha.p + eta.p.s[site[i]] 
		                + w.p.ex*beta.p.ex*ex[i]
		                + w.p.veg*beta.p.lowveg*veg.cov.low[i] 
		                + w.p.veg*beta.p.modveg*veg.cov.mod[i]
	                	)

   #abandonment:
    cta[i] <- exp(alpha.a + eta.a.s[site[i]] 
		                + w.a.ex*beta.a.ex*ex[i] 
		                + w.a.init*beta.a.init*init[i]
		                + w.a.hab*beta.a.dune*dune[i]
	                	+ w.a.hab*beta.a.over*over[i]
	                	)
    cts[i] <- 1

   #Equation 5
    den[i] <- ctf[i] + ctp[i] + cta[i] + cts[i]  
    survp[i] <- cts[i]/den[i]
#interval nest loss probabilities (Equation 4)
#flooding: 
p[i,4] <- ((ctf[i]/(den[i]))/(1 - survp[i]))*(1 - pow(survp[i], interval[i])) 

#abandonment
p[i,2] <- ((cta[i]/(den[i]))/(1 - survp[i]))*(1 - pow(survp[i], interval[i])) 
 
#predation
p[i,3] <- ((ctp[i]/(den[i]))/(1 - survp[i]))*(1 - pow(survp[i], interval[i])) 

#interval survival probability 
p[i,1] <- pow(survp[i],interval[i])

#Equation 1    
Fate[i,1:4] ~ dmulti(p[i,] , 1 )
    
}

#Derived quantities 
  #Intermediate calculations for exclosed nests:
    ctp.ex <- exp(alpha.p + beta.p.ex)
    cta.ex <- exp(alpha.a + beta.a.ex)
    den.ex <- (ctp.ex + cta.ex + ctf.mean + 1)

  #Intermediate calculations for unexclosed nests:  
    ctp.un <- exp(alpha.p)
    cta.un <- exp(alpha.a)
    den.un <-(ctp.un + cta.un + ctf.mean + 1)

  #Intermediate calculations for flooding
    ctf.mean <- exp(alpha.f)

#Daily fate probabilities for exclosed nests
    Survp.ex<-1/den.ex
    Abanp.ex<-cta.ex/den.ex
    Predp.ex<-ctp.ex/den.ex
    
#Daily fate probabilties for unexclosed nests
    Survp.un<-1/den.un
    Abanp.un<-cta.un/den.un
    Predp.un<-ctp.un/den.un

#Mean daily flooding probability   
    Floodp<-ctf.mean/den.un
    
    }
    
", fill=T)
sink()

#package data for analysis in JAGS
win.data.sc <- list(ex = scale(nestdata$Exclosed), n = n,    
                    site = as.numeric(droplevels(nestdata$Site.x)), nSite = 46, interval = as.numeric(nestdata$Interval),Fate = Fate, init = scale(nestdata$Nest.Initiation), dune = scale(dune),over = scale(overwash),veg.cov.low = scale(veg.cover.low), veg.cov.mod = scale(veg.cover.mod))

#define function to draw initial values for MCMC chains
inits <- function() {list(alpha.p = rnorm(1, 0, 1), alpha.a =   rnorm(1, 0, 1), alpha.f = rnorm(1,0,1), beta.a.ex = rnorm(1, 0, 1), sigma.a = rlnorm(1), sigma.p = rlnorm(1), sigma.f = rlnorm(1), beta.p.ex = rnorm(1,0,1), beta.a.init = rnorm(1,0,1), beta.a.dune = rnorm(1,0,1), beta.a.over = rnorm(1,0,1), beta.p.lowveg = rnorm(1,0,1), beta.p.modveg = rnorm(1,0,1),	w.p.ex = 1, w.p.veg  = 1, w.a.ex = 1, w.a.init = 1,w.a.hab =1)}

#list parameters to monitor
params <- c("sigma.a", "sigma.p", "sigma.f", "alpha.p", "alpha.f", "alpha.a", "beta.a.ex", "beta.p.ex", "beta.p.lowveg", "beta.p.modveg", "beta.a.init", "beta.a.dune", "beta.a.over", "Survp.ex", "Survp.un", "Abanp.ex", "Abanp.un",  
            "Predp.ex", "Predp.un","Floodp", "w.p.ex", "w.p.veg", "w.a.ex",   "w.a.init", "w.a.hab")

ni <- 70000
nb <- 40000
nt <- 1
nc <- 3

out.sub3.sc <- jags(win.data.sc, inits, params, "BUGS_model_subset_select.txt", n.iter=ni, n.thin=nt, n.burnin=nb, n.chains=nc, paral
