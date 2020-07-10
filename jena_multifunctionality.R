
# Title: Many dimensions of multifunctionality

# Re-analysing the Jena data to examine how different components of multifunctionality respond to different aspects of diversity

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

# seperate the data into spp and site characteristics matrix
site_dat <- 
  jena_comm %>%
  select(-sp_names)

sp_dat <- 
  jena_comm %>% 
  select(plotcode, year, month, sp_names)

# replace NAs in sp_dat with zeros to reflect absence
sp_dat <- 
  sp_dat %>% 
  mutate_at(vars(sp_names), ~replace(., is.na(.), 0))

# check if there are any -9999's in the species data
sp_dat %>% 
  filter_at(vars(sp_names), any_vars(. < 0)) # there are none!

# calculate the number of species in the surveyed 3 x 3 m plots
rowSums( decostand(select(sp_dat, sp_names), method = "pa") ) # number of species

# add these observed species richness values to the site_dat dataset
site_dat <- 
  site_dat %>% 
  mutate(observed_species = rowSums(decostand(select(sp_dat, sp_names), method = "pa")),
         ens = exp(diversity(x = select(sp_dat, sp_names), index = "shannon")) )


# remove unnecessary site_dat variables
names(site_dat)
site_dat <- 
  site_dat %>%
  select(plotcode:leg.ef, observed_species, ens)

# check the plotcodes to see how many there are
site_dat$plotcode %>% 
  unique() %>%
  length()

# the multifunctionality measurements are based on the final year of the data (2007)
# if functions were measured in multiple years, they were averaged
# also, they exclude plotcode = B4A03
unique(site_dat$year)

site_dat <- 
  site_dat %>%
  filter(plotcode != "B4A03") %>%
  filter(year == 2007) %>%
  filter(sowndiv < 60) %>%
  group_by(block, plot, plotcode, year) %>%
  summarise(across(.cols = c("sowndiv", "numfg", "numgrass", "numsherb",
                             "numtherb", "numleg", "gr.ef", "sh.ef", 
                             "th.ef", "leg.ef"), ~first(x = .x)),
            across(.cols = c("observed_species", "ens"), ~mean(x = .x, na.rm = TRUE)),
            .groups = "drop")

site_dat$plotcode

sp_dat <- 
  sp_dat %>%
  filter(plotcode %in% site_dat$plotcode) %>%
  filter(year == 2007) %>%
  select(-month) %>%
  group_by(plotcode, year) %>%
  summarise(across(.cols = everything(), ~mean(x = .x, na.rm = TRUE)),
            .groups = "drop") %>%
  select(-year) %>%
  separate(col = plotcode, into = c("block", "plot"), sep = 2)

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

# remove plots with sowndiv = 60
multi_dat <- 
  multi_dat %>%
  filter(plotcode %in% site_dat$plotcode)

# split plotcode into plot and code
multi_dat <- 
  multi_dat %>%
  separate(col = plotcode, into = c("block", "plot"), sep = 2)

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
  select(-block, -plot) %>%
  names()


# draw landscapes of n = 4 plots and write them into a list
f_names
sp_names

samp_dat <- full_join(sp_dat, multi_dat, by = c("block", "plot"))

# set number of landscapes
n <- 100

# set plots in landscapes
p <- 4

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

l_scapes_div
l_scapes_mf

# calculate diversity and multifunctionality metrics

# diversity metrics
# alpha, beta, gamma diversity
div <- 
  l_scapes_div %>%
  lapply(., function(x) {
    
    # alpha diversity (a)
    alpha <- 
      rowSums(decostand(x = x, method = "pa")) %>%
      mean(., na.rm = TRUE)
    
    alpha_ens <- exp(diversity(x = x, index = "shannon")) %>%
      mean(., na.rm = TRUE)
    
    a <- tibble(alpha, alpha_ens)
      
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
    alpha_mf <- 
      (rowSums(x)/ncol(x) ) %>%
      mean(., na.rm = TRUE)
    
    alpha_s_mf <- 
      ( (rowSums(x)/ncol(x))/apply(x, 1, sd) ) %>%
      mean(., na.rm = TRUE)
    
    a <- tibble(alpha_mf, alpha_s_mf)
    
    # beta multifunctionality
    b <- 
      mean(vegdist(x = x, method = "euclidean"))
    
    b <- tibble(beta_euc = b)
    
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
  pairs()

# check gamma relationships
mf_sims %>%
  select(contains("gamma")) %>%
  pairs()

ggplot(data = mf_sims, 
       mapping = aes(x = gamma, y = alpha)) +
  geom_point()






