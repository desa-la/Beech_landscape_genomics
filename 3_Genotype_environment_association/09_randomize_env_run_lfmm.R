################################################################################
# LFMM genotype-environment association analysis with real and randomized
# environmental data
#
# Runs LFMM (LEA package; latent factor mixed model, MCMC version) for all 19
# bioclim variables
#   1. with the real environmental data, and
#   2. with three randomized datasets in which the bioclim values (i.e. the
#      population coordinates) were randomly reassigned among populations.
#
# Randomization is done at the population level: every population receives the
# complete set of 19 bioclim values of another, randomly chosen population, so
# all individuals of a population keep sharing the same (now wrong)
# environment. Genotypes and population structure stay unchanged; only the
# genotype-environment link is broken. Any associations found with random data
# are therefore statistically significant but non-causal.
#
# Input
#   geno_file : genotype matrix in LFMM format (individuals x SNPs, missing
#               data imputed)
#   env_file  : environmental matrix (individuals x 19 bioclim variables, no
#               header), same individual order as geno_file
#   pop_file  : sample table with header; column "prov" = population ID,
#               same individual order as geno_file
#
# Output
#   all.populations.all.bio.txt               LFMM p-values, real data
#   all.populations.all.bio.random.{1,2,3}.txt LFMM p-values, random data
#   env.random.{1,2,3}.env                    randomized environment files
#   random.env.{1,2,3}.unique.txt             randomized environment per population
#   (one row per SNP, one column per bioclim variable V1..V19)
#
# The p-value files are the input for 02_plot_real_vs_random.R.
################################################################################


## ---- Packages ---------------------------------------------------------------

library(LEA)   # Bioconductor: BiocManager::install("LEA")


## ---- Settings ---------------------------------------------------------------

geno_file <- "final.filtered.bim.indep.lfmm_imputed.lfmm"
env_file  <- "bio.present.env"
pop_file  <- "populations.final.UAcorrected.txt"

n_random  <- 3          # number of randomized datasets
seeds     <- c(1, 2, 3) # one seed per randomization, for reproducibility

# LFMM parameters
lfmm_K           <- 3     # number of latent factors
lfmm_repetitions <- 5     # independent MCMC runs, combined by lfmm.pvalues()
lfmm_iterations  <- 10000  # MCMC iterations
lfmm_burnin      <- 5000  # burn-in iterations
#   the number of repetitions, iterations and burnin can be reduced to speed up
#   calculations (e.g., to 3, 5000 and 2500)
lfmm_CPU         <- 48    # threads

n_bio <- 19


## ---- Functions --------------------------------------------------------------

# Reassign environments among populations at random.
#   env  : data frame, individuals x bioclim variables
#   prov : population ID of each individual (same order as env)
# Returns a data frame of the same shape in which each population carries the
# bioclim values of another, randomly drawn population.
randomize_env_by_population <- function(env, prov, seed) {
  set.seed(seed)

  # One row of bioclim values per population
  env_pop <- unique(cbind(prov = prov, env))
  if (any(duplicated(env_pop$prov))) {
    stop("Individuals of the same population have different bioclim values.")
  }

  # Shuffle population labels: the bioclim row of population Y is now
  # labelled with a random population X
  env_pop$prov_random <- sample(env_pop$prov)

  # Give every individual of population X the bioclim row labelled X
  rows <- match(prov, env_pop$prov_random)
  env_random <- env_pop[rows, names(env)]
  rownames(env_random) <- NULL
  env_random
}

# Run LFMM for all bioclim variables and return a SNP x variable matrix of
# p-values (calibrated and combined across repetitions by lfmm.pvalues()).
run_lfmm_all_bio <- function(geno_file, env_file) {
  project <- lfmm(input.file       = geno_file,
                  environment.file = env_file,
                  K                = lfmm_K,
                  d                = 1:n_bio,
                  repetitions      = lfmm_repetitions,
                  iterations       = lfmm_iterations,
                  burnin           = lfmm_burnin,
                  CPU              = lfmm_CPU,
                  project          = "new")

  pvalues <- sapply(1:n_bio, function(d) {
    lfmm.pvalues(project, K = lfmm_K, d = d)$pvalues
  })
  colnames(pvalues) <- paste0("V", 1:n_bio)
  pvalues
}


## ---- Read data --------------------------------------------------------------

pops <- read.table(pop_file, header = TRUE)
env  <- read.table(env_file)

stopifnot(nrow(pops) == nrow(env), ncol(env) == n_bio)


## ---- LFMM with real data ----------------------------------------------------

p_real <- run_lfmm_all_bio(geno_file, env_file)
write.table(p_real, "all.populations.all.bio.txt", row.names = FALSE)


## ---- LFMM with randomized data ----------------------------------------------

for (r in 1:n_random) {

  # Randomize environment among populations and save it
  env_random <- randomize_env_by_population(env, pops$prov, seed = seeds[r])
  env_random_file <- sprintf("env.random.%d.env", r)
  write.env(as.matrix(env_random), env_random_file)

  # Keep a per-population table of the randomized environment for the record
  write.table(unique(cbind(pops[, -1], env_random)),
              sprintf("random.env.%d.unique.txt", r), row.names = FALSE)

  # Run LFMM on the randomized environment
  p_random <- run_lfmm_all_bio(geno_file, env_random_file)
  write.table(p_random, sprintf("all.populations.all.bio.random.%d.txt", r),
              row.names = FALSE)
}
