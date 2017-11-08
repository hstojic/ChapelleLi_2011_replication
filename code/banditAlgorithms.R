# ----------------------------------------------------------------------
# Information
# ----------------------------------------------------------------------

# This script defines the algorithms from the paper for the 
# Bernoulli bandit problem


# ----------------------------------------------------------------------
# Auxiliary functions
# ----------------------------------------------------------------------

# finding a max value in a vector with real numbers, randomly breaking
# the ties
idxMax <- function (x, tol = 1e-16) {
    stopifnot(is.numeric(x), is.numeric(tol), all(!is.na(x)))
    if (all(is.infinite(x)) || abs(max(x) - min(x)) < tol) {
        idx <- sample(1:length(x), 1)
    }
    else {
        idx <- which(x == max(x))
        if (length(idx) > 1) {
            idx <- sample(idx, 1)
        }
    }
    return(idx)
}

# help function for KL divergence
klDivergence <- function(p, q) {
    sum(p*log(p/q))
}

# a function for generating the Bernoulli bandit observations 
genBernoulliBandit <- function(noArms, prob) {
    # basic checks
    stopifnot(noArms == length(prob))
    # make draws 
    draws <- function() rbinom(noArms, 1, prob)
    return(draws)
}


# ----------------------------------------------------------------------
# Algorithms
# ----------------------------------------------------------------------


Thompson <- function(env, pars, parUpdate = 100) {

    # basic vars
    aPrior <- pars[["a"]]
    bPrior <- pars[["b"]]
    noTrials <- env$noTrials
    noArms <- env$noArms
    rewardFnc <- env$reward

    # objects for results
    choices <- rep(NA, noTrials)
    success <- rep(0, noArms)
    failure <- rep(0, noArms)

    for (trial in 1:noTrials) {
        
        # draw prob of reward based on the model of the environment 
        # i.e. from the posterior distribution P(theta|data)
        thetat <- rbeta(noArms, success + aPrior, failure + bPrior)

        # choice according to the draws from the posterior
        choice <- which.max(thetat)
        choices[trial] <- choice

        # environment returns a reward for a given choice
        reward <- rewardFnc()[choice]

        # update the model of the environment, e.g. only every 100 trials  
        if (trial %% parUpdate == 0) {
            if(reward == 0) {
                failure[choice] <- failure[choice] + 1
            } else {
                success[choice] <- success[choice] + 1
            }
        }
    }

    return(list(choices = choices))
}


UCB <- function(env, pars, parUpdate = 100) {

    # basic vars
    noTrials <- env$noTrials
    noArms <- env$noArms
    rewardFnc <- env$reward

    # objects for results
    choices <- rep(NA, noTrials)
    k <- rep(0, noArms)  
    m <- rep(0, noArms)

    for (trial in 1:noTrials) {
        # trial = 11
        
        idxNontried <- m == 0
        if (any(idxNontried)) {
            ucbEst <- idxNontried*1
        } else {
            # update delta 
            delta <- sqrt(1 / trial)

            # computing the upper confidence bound
            ucbEst <- k / m + 
                      sqrt((2*(k/m)*log(1/delta)) / m) + 
                      (2*log(1/delta)) / m
        }

        # choice according to the best bound
        choice <- idxMax(ucbEst)
        choices[trial] <- choice

        # environment returns a reward for a given choice
        reward <- rewardFnc()[choice]

        # update the model of the environment, e.g. only every 100 trials 
        if (trial %% parUpdate == 0) {
            k[choice] <- k[choice] + reward
            m[choice] <- m[choice] + 1
        }
    }

    return(list(choices = choices))
}

