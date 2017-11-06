# ----------------------------------------------------------------------
# Information
# ----------------------------------------------------------------------

# This script defines various models for our Bernoulli MAB task 
# used in simulations


# ----------------------------------------------------------------------
# Auxiliary functions
# ----------------------------------------------------------------------


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

klDivergence <- function(p, q) {
    sum(p*log(p/q))
}


# a function for generating a Bernoulli bandit observations 
genBernoulliBandit <- function(noTrials = 10, noArms = 2, 
                               a = 1, b = 1, probs) {

    if (missing(probs)) probs <- rbeta(noArms, a, b)

    bandit <- matrix(NA, noTrials, noArms) 
    for (trial in 1:noTrials) {
        bandit[trial,] <- rbinom(noArms, 1, prob = probs)
    }

    return(bandit)
}


# ----------------------------------------------------------------------
# Algorithms
# ----------------------------------------------------------------------


Thompson <- function(rewards, pars, parUpdate = 100) {

    # basic vars
    aPrior <- pars[["a"]]
    bPrior <- pars[["b"]]
    noTrials <- nrow(rewards)
    noArms <- ncol(rewards)

    # frames and counters 
    choices <- rep(NA, noTrials)
    success <- rep(0, noArms)
    failure <- rep(0, noArms)

    for (trial in 1:noTrials) {
        
        # draw prob of reward based on the model of the environment 
        # i.e. from P(theta|data)
        # sapply(1:noArms, function(arm) rbeta())
        thetat <- rbeta(noArms, success + aPrior, failure + bPrior)

        # choice according to the best estimated probability
        choice <- which.max(thetat)
        choices[trial] <- choice
        reward <- rewards[trial, choice]

        # update the model of the environment 
        if (trial %% parUpdate == 0) {
            if(reward == 0) {
                failure[choice] <- failure[choice] + 1
            } else {
                success[choice] <- success[choice] + 1
            }
        }
    }

    payoff <- rewards[cbind(1:noTrials, choices)]
    return(list(choices = choices, payoff = payoff))
}


UCB <- function(rewards, pars, parUpdate = 100) {

    # basic vars
    noTrials <- nrow(rewards)
    noArms <- ncol(rewards)

    # frames 
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

            # choosing the arm with highest upper confidence bound
            ucbEst <- k / m + 
                      sqrt((2*(k/m)*log(1/delta)) / m) + 
                      (2*log(1/delta)) / m
        }

        choice <- idxMax(ucbEst)
        choices[trial] <- choice

        # choice according to the best estimated probability
        reward <- rewards[trial, choice]

        # update the model of the environment, e.g. only every 100 trials 
        if (trial %% parUpdate == 0) {
            k[choice] <- k[choice] + reward
            m[choice] <- m[choice] + 1
        }
    }

    payoff <- rewards[cbind(1:noTrials, choices)]
    return(list(choices = choices, payoff = payoff))
}

