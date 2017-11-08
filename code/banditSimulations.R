# ----------------------------------------------------------------------
# Description
# ----------------------------------------------------------------------

# Script for running the Bernoulli bandit simulations from 
# Chapelle, Li 2011 article


# ----------------------------------------------------------------------
# Loading the librarires and functions
# ----------------------------------------------------------------------

# set the directory if using the script interactively, 
# below you see example of my path to the folder
# setwd("/home/hstojic/Research/replications/gap_Chapelle_Li_2011/code")

# house keeping
rm(list = ls())

# load libraries and functions
source("utils.R")
source("banditAlgorithms.R")


# ----------------------------------------------------------------------
# Setting the parameters
# ----------------------------------------------------------------------

# parameters for simulation
noSim <- 100
subsample <- 100  # we dont save all the data at the end

# defining the bandit problem characteristics
noTrials <- 10^7
noArms <- c(10, 100)
epsilon <- c(0.1, 0.02)
prob <- 0.5

# algorithm info
algoInfo <- list(
    "Thompson" = list(Thompson, pars = list(a = 1, b = 1)),
    "UCB" = list(UCB, pars = NA)
)
parUpdate <- 1


# to speed up simulations, we run it in parallel 
seed <- 1234
noCores <-  ifelse(detectCores() == 1, 1, detectCores() - 1)
registerDoParallel(noCores)
registerDoRNG(seed)  


# ----------------------------------------------------------------------
# Computing asymptotic lower bounds
# ----------------------------------------------------------------------

# we can compute them explicitly, no need for simulations

# setting conditions info
algo <- "Asymptotic Lower Bound"
constant <- 0
condInfo <- expand.grid(noArms, epsilon, stringsAsFactors = FALSE)
colnames(condInfo) <- c("noArms", "epsilon")

# data frame for all results
boundResults <- data.frame()

# iterating over conditions
for (cond in 1:nrow(condInfo)) {
    # cond = 1
    print(cond)

    # extract number of arms and epsilon from condInfo
    probs <- c(
        prob, 
        rep(prob - condInfo$epsilon[cond], condInfo$noArms[cond] - 1)
    )
    
    # iterating over simulation iterations in parallel
    algoRes <- 
        log(1:noTrials) * 
        (sum((prob - probs) / 
        sapply( probs, function(x) klDivergence(prob, x)), na.rm = TRUE) + 
        constant)

    # computing regret
    results <- data.frame(
        trial = 1:noTrials,
        reg_mean = algoRes,
        reg_sd = NA,
        algo = algo,
        noArms = condInfo$noArms[cond],
        epsilon = condInfo$epsilon[cond]
    ) %>%
    filter(trial %in% c(1, seq(subsample, noTrials, subsample))) 
    boundResults <- rbind(boundResults, results)
}


# ----------------------------------------------------------------------
# Simulation
# ----------------------------------------------------------------------

# setting condition info for simulation iterations
condInfo <- expand.grid(noArms, epsilon, names(algoInfo), 
    stringsAsFactors = FALSE)
colnames(condInfo) <- c("noArms", "epsilon", "algo")

# data frame for all results
simResults <- data.frame()

# iterating over conditions
for (cond in 1:nrow(condInfo)) {
    # cond = 1
    cat("-----------------------------------------------------------------\n
         -----------------------------------------------------------------\n
         Condition: ", cond, " out of ", nrow(condInfo), "\n 
         -----------------------------------------------------------------\n
         -----------------------------------------------------------------\n")

    # extract number of arms and epsilon from condInfo
    probs <- c(
        prob, 
        rep(prob - condInfo$epsilon[cond], condInfo$noArms[cond] - 1)
    )
    algo <- algoInfo[[condInfo$algo[cond]]][[1]]
    algoPars <- algoInfo[[condInfo$algo[cond]]][["pars"]]
    noArms <- condInfo$noArms[cond]
    env <- list(
        reward = genBernoulliBandit(noArms, probs),
        noTrials = noTrials,
        noArms = noArms
    )

    # iterating over simulation iterations in parallel
    results <- foreach(sim = 1:noSim, .combine = rbind) %dorng% {

        # single simulation of an algorithm
        startTime <-  Sys.time()
        algoRes <- algo(env, algoPars, parUpdate)
        endTime <-  Sys.time()
        runTime <- endTime - startTime
        if (attr(runTime, "units") == "secs") runTime <- runTime/60
        cat("Simulation: ", sim, " | Execution time: ", 
            as.numeric(runTime, "mins"), "min\n")

        # computing regret
        algoCumRegret <- data.frame(
            iter = sim, 
            trial = 1:noTrials,
            cumreg = cumsum(0.5 - probs[algoRes$choices])
        ) %>%
        filter(trial %in% c(1, seq(subsample, noTrials, subsample))) 

        # provide some estimate of remaining runs duration
        cat("Estimated remaining time: ", 
            as.numeric(runTime, "mins")*(noSim-sim)/noCores*nrow(condInfo), 
            "min\n")

        return(algoCumRegret)
    }

    # computing averages across simulation runs
    results <- results %>%
        group_by(trial) %>% 
        summarise(reg_mean = mean(cumreg), reg_sd = sd(cumreg)) %>% 
        mutate(
            algo = condInfo$algo[cond],
            noArms = condInfo$noArms[cond],
            epsilon = condInfo$epsilon[cond]
        )
    simResults <- rbind(simResults, results)
        
    # saving interim results
    save(results, file = paste0("../data/bernoulliBandit_", cond, ".RData"))
}


# ----------------------------------------------------------------------
# Saving
# ----------------------------------------------------------------------

# combining first bound and simulation results
results <- rbind(boundResults, simResults)

# if everything run succesfully, saving the data
save(results, file = "../data/bernoulliBandit.RData")
