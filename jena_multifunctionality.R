
# Title: Many dimensions of multifunctionality

# Re-analysing the Jena data to examine how different components of multifunctionality respond to different aspects of diversity

# load relevant libraries
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(broom)
library(RColorBrewer)
library(viridis)
library(here)
library(vegan)


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
  group_by(block, plot, plotcode, year) %>%
  summarise(across(.cols = c("sowndiv", "numfg", "numgrass", "numsherb",
                             "numtherb", "numleg", "gr.ef", "sh.ef", 
                             "th.ef", "leg.ef"), ~first(x = .x)),
            across(.cols = c("observed_species", "ens"), ~mean(x = .x, na.rm = TRUE)),
            .groups = "drop")

sp_dat <- 
  sp_dat %>%
  filter(plotcode != "B4A03") %>%
  filter(year == 2007) %>%
  select(-month) %>%
  group_by(plotcode, year) %>%
  summarise(across(.cols = everything(), ~mean(x = .x, na.rm = TRUE)),
            .groups = "drop")
  
  


# load the multifunctionality data (Meyers et al. 2017)
multi_dat <- read_delim(here("data/Jena_data_functions_final_year_available.csv"), delim = ";")
multi_dat

# check the plotcodes
multi_dat$Plot %>%
  unique() %>%
  length()

# which plot is missing from the multifunctionality data?
unique(site_dat$plotcode)[!( unique(site_dat$plotcode) %in% (multi_dat$Plot) )]

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

# Invert functions in the spp_func_dat dataset
multi_dat <- 
  multi_dat %>% 
  mutate_at(func_inv, ~(.*-1))





# remove the B4A03 plot 


# join the two datasets by their plotcode
spp_multi <- full_join(spp_dat, func_dat, by = c("plotcode"))

# Check the join
# filter(spp_func_dat, plotcode == "B1A05") %>% select(10:15)
# filter(spp_dat, plotcode == "B1A05")
# filter(func_dat, plotcode == "B1A05")

# Check colnames
colnames(spp_func_dat)

# Visualise some of the data: without the 60 diversity plot
ggplot(data = spp_func_dat %>% filter(sowndiv < 60), 
       mapping = aes(x = sowndiv, y = (-BM_weed_DW))) +
  geom_point() +
  geom_smooth(method = "lm")



# Did I invert the correct functions?
spp_func_dat %>% select(func_inv) #yes

# Remove the sowndiv = 60 plots because there are only four of them in the data set
spp_func_dat16 <- spp_func_dat %>% filter(sowndiv < 60)

# Remove the plot variables besides plotcode and sowndiv
colnames(spp_func_dat16)
plot_data <- c("block", "plot", "numfg", "numgrass", "numsherb", "numtherb",
               "numleg", "gr.ef", "sh.ef", "th.ef", "leg.ef") 

spp_func_dat16 <- spp_func_dat16 %>% select(-plot_data)
spp_func_dat16



