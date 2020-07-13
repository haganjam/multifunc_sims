
# Title: Many dimensions of multifunctionality

# Re-analysing the Jena data to examine how different components of multifunctionality respond to different aspects of diversity

# Next steps: 
# (1) hypothesise about which aspects of diversity affect which dimension of multifunctionality
# (2) quantify species composition and multifunctional composition (and correlate)

# load relevant libraries
library(readr)
library(tibble)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(broom)
library(RColorBrewer)
library(viridis)
library(here)
library(vegan)
library(betapart)


# load the Jena community data
jena_comm <- read_delim(here("data/Jena_Community_02-08.csv"), delim = ",")
head(jena_comm)

# remove the first plots that were not sown with any species
jena_comm <- 
  jena_comm %>%
  filter(!sowndiv %in% c(0))

# create a vector of species names
sp_names <- names(jena_comm)[51:110]

# the multifunctionality measurements are based on the final year of the data (2007)
# if functions were measured in multiple years, they were averaged
# also, they exclude plotcode = B4A03
unique(jena_comm$year)

plot_codes <- 
  jena_comm %>%
  select(-sp_names) %>%
  filter(plotcode != "B4A03") %>%
  filter(year == 2007) %>%
  filter(sowndiv < 60, sowndiv > 1) %>%
  pull(plotcode) %>%
  unique()

# separate out only the species data
sp_dat <- 
  jena_comm %>% 
  select(plotcode, year, month, sp_names)

# replace NAs in sp_dat with zeros to reflect absence
sp_dat <- 
  sp_dat %>% 
  mutate_at(vars(sp_names), ~replace(., is.na(.), 0))

# check if there are any -9999's in the species data
sp_dat %>% 
  filter_at(vars(sp_names), any_vars(. < 0))

# subset out the relevant data from the sp_dat data frame
sp_dat <- 
  sp_dat %>%
  filter(plotcode %in% plot_codes) %>%
  filter(year == 2007) %>%
  select(-month) %>%
  group_by(plotcode, year) %>%
  summarise(across(.cols = everything(), ~mean(x = .x, na.rm = TRUE)),
            .groups = "drop") %>%
  select(-year)

sp_dat
  
  
# load the multifunctionality data (Meyers et al. 2017)
multi_dat <- read_delim(here("data/Jena_data_functions_final_year_available.csv"), delim = ";")
multi_dat

# check the plotcodes
multi_dat$Plot %>%
  unique() %>%
  length()

# plot B4A03 (a monoculture which was removed from the analysis Meyer et al. 2017)

# ecosystem functions to remove from Google drive (n = 14) decided during the workshop
r_funcs <- c("Diam_root", "Height_targ", "Height_weeds", "S_hymenopt", "S_parastie_hymenopt",
             "S_pollin", "S_seedb_targ", "S_seedl_targ", "S_weeds", "S_seedb_weed",
             "S_seedl_weed", "Soil_Dens", "Soil_ph_CaCl", "Soil_ph_H2O")

length(r_funcs) #yes there are 14

# remove the r_func functions from the func_dat dataset
multi_dat <- 
  multi_dat %>%
  select(-r_funcs)

# rename the Plot column to plotcode
multi_dat <- 
  multi_dat %>%
  rename(plotcode = Plot)

# which functions should be inverted?
func_inv <- c("BM_weed_DW", "Cov_bare", "Cov_weed", "Mortality_hymenopt", "N_weeds",
              "N_seedb_weed", "N_seedl_weed", "Soil_NH4_AfterGrowth",
              "Soil_NO3_AfterGrowth", "Soil_P_AfterGrowth")

# invert functions in the spp_func_dat dataset
multi_dat <- 
  multi_dat %>% 
  mutate_at(func_inv, ~(.*-1))

# examine the ecosystem function measured
names(multi_dat)

# choose a subset of functions to remove (in my opinion)
rem_func <- 
  c("BM_weed_DW",
    "Cov_weed",
    "N_weeds",
    "N_seedb_weed",
    "N_seedl_weed")

multi_dat <- 
  multi_dat %>%
  select(-all_of(rem_func) )


# remove plots to match with community data
multi_dat <- 
  multi_dat %>%
  filter(plotcode %in% plot_codes)

multi_dat

# standardise the functions by the standard deviation

# write a function to do the standardisation
scale_this <- function(x) {
  (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
}

multi_dat <- 
  multi_dat %>%
  mutate(across(where(is.numeric), ~scale_this(.x))) %>%
  mutate(across(where(is.numeric), ~(.x + abs(min(.x)) ) ) )

# get function names
f_names <- 
  multi_dat %>%
  select(-plotcode) %>%
  names()

# draw landscapes of n = 4 plots and write them into a list
f_names
sp_names

samp_dat <- full_join(sp_dat, multi_dat, by = c("plotcode"))

# set number of landscapes
n <- 1000

# set plots in landscapes
p <- 8

# set up output lists for: diversity
l_scapes_div <- vector("list", length = n)
names(l_scapes_div) <- paste0(c("l"), seq_along(1:n))

# set up output lists for: multifunctionality
l_scapes_mf <- vector("list", length = n)
names(l_scapes_mf) <- paste0(c("l"), seq_along(1:n))

for (i in seq_along(1:n)) {
  
    x <- samp_dat[ sample(x = seq_along(1:nrow(samp_dat)), size = p), replace = TRUE, ]
    
    l_scapes_div[[i]] <- 
      x %>%
      select(all_of(sp_names) )
    
    l_scapes_mf[[i]] <- 
      x %>%
      select(all_of(f_names) )
  
}

# calculate diversity and multifunctionality metrics

# diversity metrics
# alpha, beta, gamma diversity
div <- 
  l_scapes_div %>%
  lapply(., function(x) {
    
    # alpha diversity (a)
    v <- rowSums(decostand(x = x, method = "pa"))
    
    alpha <- mean(v, na.rm = TRUE)
    
    alpha_cv <- sd(v, na.rm = TRUE)/alpha
    
    w <- exp(diversity(x = x, index = "shannon"))
    
    alpha_ens <- mean(w, na.rm = TRUE)
    
    alpha_ens_cv <- sd(w, na.rm = TRUE)/alpha_ens
    
    a <- tibble(alpha, alpha_cv, alpha_ens, alpha_ens_cv)
      
    # beta diversity (b)
    b <- 
      beta.multi.abund(x = x, index.family = "bray") %>% 
      unlist(.) %>%
      enframe(., name = "beta", value = "bc_dis") %>%
      spread(key = "beta", value = "bc_dis")
    
    names(b) <- c("beta", "beta_bal", "beta_gra")
    
    # gamma diversity (c)
    z <- summarise(x, across(.cols = everything(), sum)) 
    
    gamma <- 
      rowSums(decostand(x = z, method = "pa"))
    
    gamma_ens <- exp(diversity(x = z, index = "shannon"))
    
    c <- tibble(gamma , gamma_ens)
    
    # bind these diversity metrics
    bind_cols(a, b, c)
      
  } )

div <- bind_rows(div, .id = "landscape")


# multifunctionality metrics
# alpha, beta, gamma multifunctionality
mf <- 
  l_scapes_mf %>%
  lapply(., function(x) {
    
    # alpha multifunctionality
    v <- (rowSums(x)/ncol(x) ) 
    
    alpha_mf <- mean(v, na.rm = TRUE)
    
    alpha_cv_mf <- sd(v, na.rm = TRUE)/alpha_mf
    
    alpha_s_mf <- 
      ( v/apply(x, 1, sd) ) %>%
      mean(., na.rm = TRUE)
    
    a <- tibble(alpha_mf, alpha_cv_mf, alpha_s_mf)
    
    # beta multifunctionality
    b <- 
      beta.multi.abund(x = x, index.family = "bray") %>% 
      unlist(.) %>%
      enframe(., name = "beta", value = "bc_dis") %>%
      spread(key = "beta", value = "bc_dis")
    
    names(b) <- c("beta_mf", "beta_bal_mf", "beta_gra_mf")
    
    # gamma multifunctionality
    z <- summarise(x, across(.cols = everything(), sum)) 
    
    gamma_mf <- 
      rowSums(z)/ncol(z)
    
    gamma_s_mf <- 
      (rowSums(z)/ncol(z))/apply(z, MARGIN = 1, sd)
    
    c <- tibble(gamma_mf  , gamma_s_mf)
    
    bind_cols(a, b, c)
    
  } )

mf <- bind_rows(mf, .id = "landscape")  


### join the diversity metrics and multifunctionality metrics and examine relationships

mf_sims <- full_join(div, mf, by = "landscape")
names(mf_sims)

# check alpha relationships
mf_sims %>% 
  select(contains("alpha")) %>%
  pairs()

# check beta relationships
mf_sims %>%
  select(contains("beta")) %>%
  mutate(across(.cols = everything(), ~ car::logit(.x)) ) %>%
  pairs()

mf_sims %>%
  select(contains("beta")) %>%
  mutate(across(.cols = everything(), ~ car::logit(.x)) ) %>%
  cor(method = "spearman")

# check gamma relationships
mf_sims %>%
  select(contains("gamma")) %>%
  pairs()


# check specific relationships
mf_sims %>% names()

ggplot(data = mf_sims, 
       mapping = aes(x = alpha_cv_mf, y = car::logit(beta_gra_mf) )) +
  geom_point() +
  geom_smooth()

ggplot(data = mf_sims, 
       mapping = aes(x = alpha_cv_mf, y = car::logit(beta_bal_mf) )) +
  geom_point() +
  geom_smooth()

ggplot(data = mf_sims, 
       mapping = aes(x = car::logit(beta_bal), y = car::logit(beta_bal_mf), colour = alpha)) +
  geom_point() +
  scale_colour_viridis_c() +
  geom_smooth()

ggplot(data = mf_sims, 
       mapping = aes(x = car::logit(beta_gra), y = car::logit(beta_gra_mf), colour = alpha_cv_mf)) +
  geom_point() +
  geom_smooth(method = "lm")



### quantify species composition and multifunctionality composition (identity)

# does turnover in species composition drive turnover in multifunctional composition?

# use a mantel test
multi_dat

f_pca <- 
  princomp(formula = reformulate(paste(f_names[1:30], sep = "+")), data = multi_dat, cor = FALSE)

mantel(vegdist(x = select(sp_dat, -plotcode), method = "bray"),
       vegdist(x = select(multi_dat, -plotcode), method = "euclidean"),
       method = "pearson")

plot(x = vegdist(x = select(sp_dat, -plotcode), method = "bray"),
     y = vegdist(x = select(multi_dat, -plotcode), method = "euclidean"))


### do different species contribute to different functions?

sp_dat

func_list <- 
  multi_dat %>%
  pivot_longer(., cols = where(is.numeric), names_to = c("eco_function"), values_to = c("value") ) %>%
  arrange(eco_function, plotcode) %>%
  split(., .$eco_function)

length(func_list)

# let's try this with the first 10 functions
func_list_test <- func_list[1:10]

# set the number of replications for the null expectation
reps <- 10

ses_spp <- 
  lapply(func_list_test, function(x) {
  
  ses <- vector("list", length = length(sp_names))
  
  for (i in seq_along(1:length(sp_names))) {
    
    u <- 
      sp_dat %>%
      select(plotcode, all_of(sp_names[i])) %>%
      rename(group = sp_names[i]) %>%
      mutate(group = as.numeric(decostand(x = group, method = "pa")) )
    
    u
    
    # get the observed SES value for each species  
    w <- 
      full_join(u, x, by = "plotcode") %>%
      group_by(group) %>%
      summarise(mean_value = mean(value, na.rm = TRUE),
                n = n(), .groups = "drop") %>%
      summarise(SES = diff(mean_value),
                n = abs(diff(n)) ) %>%
      mutate(species = sp_names[i],
             eco_function = f_names[i]) %>%
      select(eco_function, species, SES, n)
    
    w
    
    # get the null expectation of observed values among species
    
    null_func <- 
      
      function(x) {
        
        s <- sample(x = seq_along(1:nrow(x)), 
                    size = nrow(x), replace = FALSE)
        
        bind_cols(u, select(x[s, ], value) ) %>%
          group_by(group) %>%
          summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
          summarise(SES_null = diff(mean_value)) %>%
          mutate(species = sp_names[i]) %>%
          select(species, SES_null) %>%
          pull(SES_null)
      }
    
    ses[[i]] <- bind_cols(w, null = replicate(n = reps, expr = null_func(x)))
    
  }
  
  bind_rows(ses) %>%
    group_by(species) %>%
    summarise(SES = first(SES),
              low_thres = quantile(x = null, probs = 0.025),
              upp_thres = quantile(x = null, probs = 0.975), .groups = "keep") %>%
    mutate(sig = if_else(SES > upp_thres | SES < low_thres, 1, 0)) %>%
    ungroup()
  
})

ses_spp

spp <- 
  lapply(ses_spp, function(x) {
  
  x %>%
    select(species, sig) %>%
    pivot_wider(., names_from = "species", values_from = "sig")
  
})

bind_rows(spp, .id = "func") %>%
  select(-func) %>%
  vegdist(x = ., method = "bray")

# if they do, then we can make some predictions about beta diversity

# otherwise, if similar species contribute to the different functions, then this is useless
















