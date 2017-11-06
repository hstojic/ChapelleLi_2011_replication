# ----------------------------------------------------------------------
# Description
# ----------------------------------------------------------------------

# Script for Chapelle Li 2011 article
# It produces the figure 1 from the article. Note that I use the same
# scales as in the article, to faciliate the comparison, even though 
# it makes the current figures a bit less readable.


# ----------------------------------------------------------------------
# Loading the data
# ----------------------------------------------------------------------

# you will need to set the directory if using the script interactively, 
# below you see example of my path to the folder
# setwd("/home/hstojic/Research/replications/gap_Chapelle_Li_2011/code")

# house keeping
rm(list = ls())

# load libraries and auxiliary functions
source("utils.R")

# load the dataset
load(file = "../data/results.RData")

# outDir
outDir <- ""


# ----------------------------------------------------------------------
# Reshaping the data
# ----------------------------------------------------------------------

# select models to display
whichAlgo <- c("ALB", "Thompson", "UCB")

# subsample the data
noTrials <- max(results$trial)
subsample <- 100
resultsSub <- results %>%
    filter(
        trial %in% c(1, seq(subsample, noTrials, subsample)),
        algo %in% whichAlgo
    ) %>%
    mutate(lineLabel = ifelse(trial == 0.59*noTrials, as.character(algo), NA))


# ----------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------

# plotting function
genPlot <- function(plotData, yMax, yStep) {
    ggplot(data = plotData, aes(x = trial, y = reg_mean)) + 
    geom_line(aes(group = algo, color = algo, linetype = algo)) + 
    facet_wrap(epsilon ~ noArms, 
        labeller = label_bquote(
            cols = K ~ "=" ~ .(noArms) ~ "," ~ epsilon ~ "=" ~ .(epsilon))
        ) +
    geom_label_repel( 
        aes(x = trial, y = reg_mean, label = lineLabel),
        colour = "black", segment.alpha = 1,
        segment.colour = "black", segment.size = 0.5, force = 3,
        min.segment.length = unit(0.1, "lines"),
        size = fontSize,
        inherit.aes = FALSE) + 
    scale_x_log10("\nTime (in log units)",
        limits = c(10^2, 10^7),
        breaks = 10^c(2,3,4,5,6,7)) +
    scale_y_continuous("Cumulative regret\n", 
        limits = c(0,yMax),
        breaks = seq(0, yMax, yStep)) +
    scale_color_manual("", values = cbbPalette[c(1,6,7)]) +
    pdftheme + 
    theme(legend.position = "none") 
}

# generating plots for each experiment condition
epsilon01_n10 <- genPlot(
    plotData = filter(resultsSub, epsilon == "0.1", noArms == "10"),
    yMax = 900, yStep = 100
    )

epsilon002_n10 <- genPlot(
    plotData = filter(resultsSub, epsilon == "0.02", noArms == "10"),
    yMax = 4000, yStep = 500
    )

epsilon01_n100 <- genPlot(
    plotData = filter(resultsSub, epsilon == "0.1", noArms == "100"),
    yMax = 10000, yStep = 2000
    )

epsilon002_n100 <- genPlot(
    plotData = filter(resultsSub, epsilon == "0.02", noArms == "100"),
    yMax = 50000, yStep = 10000
    )


# ----------------------------------------------------------------------
# Saving figures
# ----------------------------------------------------------------------

### PDFs

filename <- paste0(outDir, "epsilon01_n10.pdf")
cairo_pdf(filename, height = 4, width = 7, onefile = TRUE)
print(epsilon01_n10)
dev.off()


filename <- paste0(outDir, "epsilon002_n10.pdf")
cairo_pdf(filename, height = 4, width = 7, onefile = TRUE)
print(epsilon002_n10)
dev.off()

filename <- paste0(outDir, "epsilon01_n100.pdf")
cairo_pdf(filename, height = 4, width = 7, onefile = TRUE)
print(epsilon01_n100)
dev.off()

filename <- paste0(outDir, "epsilon002_n100.pdf")
cairo_pdf(filename, height = 4, width = 7, onefile = TRUE)
print(epsilon002_n100)
dev.off()



### SVGs

svg("../article/Figure2.svg", height = 4, width = 7, onefile = TRUE)
print(relFrequencies50)
dev.off()

svg("../article/Figure2_1000.svg", height = 4, width = 7, onefile = TRUE)
print(relFrequencies1000)
dev.off()