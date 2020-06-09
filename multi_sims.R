
# Title: Simple models to understand the expected effects of biodiversity on multifunctionality

# load relevant libraries
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(broom)
library(RColorBrewer)
library(viridis)
library(here)

# define number of species and environments
specnum <- 2
funnum <- 2
envnum <- 2

# set up species matrix
spec <- paste("species_", c(1:specnum), sep = "")
spec_comb <- lapply(c(1:specnum), function(x) combn(spec, x)) 
names(spec_comb) <- paste("n_spp", c(1:specnum))

spec_comb

# set up number of species combinations
n_combs <- 
  lapply(spec_comb, function(x) apply(x, 2, function(y) paste(y, collapse = "_"))) %>%
  unlist(., use.names = FALSE) %>%
  length()

# set up environmental matrix
env <- paste("env_", c(1:envnum), sep = "")
env_comb <- lapply(c(1:envnum), function(x) combn(env, x)) 
names(env_comb) <- paste("heterogeneity", c(1:envnum))

env_comb

# set up number of environment combinations
n_envs <- 
  lapply(env_comb, function(x) apply(x, 2, function(y) paste(y, collapse = "_"))) %>%
        unlist(., use.names = FALSE) %>%
        length()

spec_comb

env_comb

# create a species matrix

spp_env_mat <- 
  tibble(
  environment = rep(lapply(env_comb, function(x) apply(x, 2, function(y) paste(y, collapse = "_"))) %>%
                      unlist(., use.names = FALSE), 
                    each = specnum*n_combs),
  richness = rep( rep( rep(1:specnum, unlist(lapply(spec_comb, ncol))), each = specnum ), n_envs),
  comp = rep( lapply(spec_comb, function(x) apply(x, 2, function(y) paste(y, collapse = "_"))) %>%
    unlist(., use.names = FALSE) %>%
    rep(., each = specnum), n_envs) ,
  spp = rep(rep(spec, times = n_combs), n_envs) 
)

# set up functional matrix
spp_abun <- 
  matrix(
    nrow = length(env),
    ncol = length(spec),
    dimnames = list(env, spec))

# populate this with values from a multivariate normal i.e. species must trade-off in abundance in different habitats
spp_abun[1:nrow(spp_abun), 1:ncol(spp_abun)] <- 
  runif( n = prod( dim( spp_abun)), 0, 1)

data.frame(spp_abun, row.names = "env")

spp_env_mat




### run original code

# define number of species and environments
specnum <- 2
envnum <- 2

# set up species matrix
spec <- paste("species_", c(1:specnum), sep = "")
spec_comb <- lapply(c(1:specnum), function(x) combn(spec, x)) 
names(spec_comb) <- paste("richness", c(1:specnum))

spec_comb

# set up environmental matrix
env <- paste("env_", c(1:envnum), sep = "")
env_comb <- lapply(c(1:envnum), function(x) combn(env, x)) 
names(env_comb) <- paste("heterogeneity", c(1:envnum))

env_comb

spp_traits <- matrix(
  nrow = length(env),
  ncol = length(spec),
  dimnames = list(env, spec))

# make species trade-off here
spp_traits[1:nrow(spp_traits), 1:ncol(spp_traits)] <-
  runif( n = prod( dim( spp_traits)), 0, 1)

spp_traits[1:nrow(spp_traits), 1:ncol(spp_traits)] <-
  c(0,1,1,0)


spp_traits


# 1 = Dominance (# pure selection effect and perfect environmental sorting)
# 2 = Weighted mean (# partial selection effect)
# 3 = Mean abundance

Scenario <- 3

plot_values <- sapply(env, #for each environment
                      function(x) {
                        lapply(spec_comb, #for each richness level
                               function(y) {
                                 apply(y, 2, function(z){ #for each species combination
                                   if(Scenario == 1){
                                     max(spp_traits[x,z])
                                   } else if(Scenario == 2){
                                     mean(sum(spp_traits[x,z]^2) / abs(sum(spp_traits[x,z])))
                                   } else {
                                     mean(spp_traits[x,z])
                                   } 
                                 }
                                 )
                               }
                        )
                      }
)



plot_values_df <- 
  tibble(
    richness = rep( rep(1:specnum, unlist(lapply(spec_comb, ncol))), envnum),
    spec_comb = rep(lapply(spec_comb, 
                           function(x) apply(x, 2, function(y) 
                             paste(y, collapse = " "))) %>% unlist(),
                    envnum),
    environment = rep(env, each = sum(choose(specnum,1:specnum))),
    functioning = unlist(plot_values)
  )

plot_values_df

plot_values_df <- 
  tibble(
    richness = rep( rep( rep(1:specnum, unlist(lapply(spec_comb, ncol))), each = specnum), envnum),
    spec_comb = rep( rep(lapply(spec_comb, 
                           function(x) apply(x, 2, function(y) 
                             paste(y, collapse = " "))) %>% unlist(),
                    each = specnum), envnum ),
    spec_abun = rep(rep(spec, times = n_combs), times = specnum),
    environment = rep(rep(env, each = sum(choose(specnum,1:specnum))), each = specnum) ,
    functioning = rep(unlist(plot_values), each = specnum)
  )

plot_values_df %>%
  mutate(functioning = if_else(richness == 1 & spec_comb != spec_abun, 0, functioning)) %>%
  group_by(environment, spec_comb) %>%
  mutate(sum_functioning = sum(functioning))

  
  functioning = if_else(richness > 1, functioning/richness, functioning))










