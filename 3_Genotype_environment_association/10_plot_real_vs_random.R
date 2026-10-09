################################################################################
# Real vs. randomized LFMM genotype-environment associations
#
# Produces the two-panel figure:
#   "Data randomization demonstrates statistically significant but non-causal
#    signal in genotype-environment association analysis using LFMM"
#
#   (a) Boxplots of -log10(p) for each of the 19 bioclim variables, for the
#       real data and three randomized datasets (random1-3). In the random
#       datasets, population coordinates (and therefore bioclim values) were
#       randomly reassigned among populations
#       (see 01_randomize_env_run_lfmm.R). Only associations with
#       FDR-corrected p-values (q-values) < 0.01 are shown.
#   (b) Ratio of random to real signal (number of associations with random
#       data / number with real data) for cumulative -log10(p) bins, starting
#       from the most significant associations (> 42, > 41, ..., > 3).
#       Points are means across the three random runs, error bars are the
#       standard error of the mean, the blue line is a fitted curve through
#       the means.
#
# The significance thresholds drawn in (a) (-log10(p) = 21 and 35) correspond
# to the cutoffs at which the random/real ratio in (b) reaches ~50% and ~5%
# (1/20), respectively. They are set by hand from (b), not computed here.
#
# Input (one column per bioclim variable, one row per SNP, header line):
#   all.populations.all.bio.txt            - LFMM p-values, real data
#   all.populations.all.bio.random.{1,2,3}.txt
#                                          - LFMM p-values, randomized data
#   These are written by 01_randomize_env_run_lfmm.R.
#
# Output:
#   boxplot.all.bio.real.vs.random+ratio_new.png
################################################################################


## ---- Packages ---------------------------------------------------------------

library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(splines)   # ns() for the fitted curve in panel b
library(ggpubr)    # ggarrange() to combine panels
library(qvalue)    # Bioconductor: BiocManager::install("qvalue")


## ---- Settings ---------------------------------------------------------------

real_file    <- "all.populations.all.bio.txt"
random_files <- sprintf("all.populations.all.bio.random.%d.txt", 1:3)
out_file     <- "boxplot.all.bio.real.vs.random+ratio_new.png"

q_cutoff     <- 0.01       # keep associations with q-value below this
thresholds   <- c(35, 21)  # -log10(p) thresholds read off panel b (5% / 50%)
cutoff_range <- 42:3       # -log10(p) cutoffs evaluated in panel b

bioclim <- paste0("bio", 1:19)


## ---- Read p-values ----------------------------------------------------------

# Real data: columns bio1..bio19
p_real <- read.table(real_file, header = TRUE)
names(p_real) <- bioclim

# Random data: columns bio1.random.1 .. bio19.random.3
p_random <- lapply(seq_along(random_files), function(r) {
  p_r <- read.table(random_files[r], header = TRUE)
  names(p_r) <- paste0(bioclim, ".random.", r)
  p_r
})

p_all <- do.call(cbind, c(list(p_real), p_random))


## ---- Long format and FDR filter ----------------------------------------------

# One row per SNP x bioclim variable x run
data <- melt(p_all, measure.vars = names(p_all))

# q-values are computed jointly over all p-values (real and random runs pooled),
# then only associations with q < q_cutoff are kept.
q <- qvalue(data$value)
data <- data[which(q$qvalues < q_cutoff), ]

# Annotate run (real / random1-3) and bioclim variable from the column name
data$run <- "real"
data$run[grep("random.1", data$variable)] <- "random1"
data$run[grep("random.2", data$variable)] <- "random2"
data$run[grep("random.3", data$variable)] <- "random3"
data$run <- factor(data$run, levels = c("real", "random1", "random2", "random3"))

data$bio <- factor(gsub("\\..*", "", data$variable), levels = bioclim)


## ---- Panel a: p-value distributions -----------------------------------------

# 19 colours interpolated from the "Paired" palette
bio_colours <- colorRampPalette(brewer.pal(12, "Paired"))(length(bioclim))

panel_a <- ggplot(data, aes(x = run, y = -log10(value),
                            group = interaction(run, bio), colour = bio)) +
  geom_boxplot() +
  scale_colour_manual(values = bio_colours) +
  geom_hline(yintercept = thresholds, linetype = c(2, 1)) +  # dashed 5%, solid 50%
  xlab("Geographic coordinates") +
  ylab("-log10(p-value)") +
  labs(colour = "Bioclim variable")


## ---- Panel b: random/real signal ratio --------------------------------------

# For each cutoff i, count associations with -log10(p) > i in the real data and
# in each random run, and express each random count as a fraction of the real
# count. Because the cutoff is a lower bound only, the bins are cumulative
# (e.g. 42-43, 41-43, 40-43, ...).
neglog_p <- -log10(data$value)

ratio_stats <- t(sapply(cutoff_range, function(i) {
  n_real   <- sum(neglog_p[data$run == "real"] > i)
  n_random <- sapply(c("random1", "random2", "random3"),
                     function(r) sum(neglog_p[data$run == r] > i))
  ratios <- n_random / n_real
  c(mean = mean(ratios), se = sd(ratios) / sqrt(length(ratios)))
}))

ratio_df <- data.frame(pvalue = cutoff_range,
                       ratio  = ratio_stats[, "mean"] * 100,   # in %
                       se     = ratio_stats[, "se"]   * 100)

panel_b <- ggplot(ratio_df, aes(x = pvalue, y = ratio)) +
  geom_point() +
  geom_errorbar(aes(ymin = ratio - se, ymax = ratio + se)) +
  # Smooth curve through the means: quasi-Poisson GLM on a natural spline
  # with 2 degrees of freedom
  stat_smooth(method = "glm", method.args = list(family = "quasipoisson"),
              formula = y ~ ns(x, 2), fullrange = FALSE, se = FALSE) +
  geom_hline(yintercept = 50) +               # 50% line (solid)
  geom_hline(yintercept = 5, linetype = 2) +  # 5% line (dashed)
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  xlab("-log10(p-value) cutoff") +
  ylab("Random signal/real signal (%)")


## ---- Combine panels and save -------------------------------------------------

png(out_file, width = 12.5, height = 5.5, res = 300, units = "in")
print(ggarrange(panel_a, panel_b, ncol = 2, labels = c("a", "b"),
                widths = c(1.8, 1)))
dev.off()


## ---- Optional: ANOVA on retained p-values ------------------------------------
# Tests whether p-values of retained associations differ between real and
# random runs and among bioclim variables. Not part of the figure.

data$real_or_random <- factor(ifelse(data$run == "real", "real", "random"))

summary.lm(aov(value ~ bio * real_or_random, data = data))
summary.lm(aov(value ~ run, data = data))
