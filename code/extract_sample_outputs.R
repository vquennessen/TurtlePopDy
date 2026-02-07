# get vlines for plotting 10 yr median lambdas over time + IQR

# empty environment
rm(list = ls())

# set working directory
# setwd('~/Projects/iliketurtles3/code')

# load libraries
library(matrixStats)
library(readr)
library(ggpattern)
library(dplyr)
library(tidyr)
library(magrittr)

# source functions
source('mating function/OSRs_to_betas.R')

# which computer we using
computer <- 'cluster'

# path based on computer being used
user <- ifelse(computer == 'cluster', '/home/quennessenv/iliketurtles3/output/',
               ifelse(computer == 'desktop',
                      'C:/Users/Vic/Box/Quennessen_Thesis/PhD Thesis/model output/iliketurtles3/',
                      'C:/Users/vique/Box/Quennessen_Thesis/PhD Thesis/model output/iliketurtles3/'))

# number of sims
nsims <- 10000

# evolution
traits <- c('emergence_success_t0', 'T_piv')
# rates <- c('effective', 'high')
rates <- c('effective')
TRTs <- c('narrow')
input_folder <- c('2025_12_21_evolution_n10_b5000')
paths <- expand.grid(input_folder, TRTs, traits, rates) %>%
    mutate(foo = paste(Var2, Var3, Var4, sep = ' - ')) %>%
  mutate(path = paste(Var1, Var3, Var4, Var2, sep = '/'))
# scenarios <- c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)
# osrs <- c(0.49, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05)
scenarios <- c(1, 3.5)
osrs <- c(0.4, 0.2)

# # conservation
# TRTs <- c('narrow')
# intensities <- paste('i', c(0.2, 0.4, 0.6, 0.8, 1), sep = '')
# frequencies <- paste('F', c(2, 4), sep = '')
# input_folder <- c('2025_12_21_conservation_n10_b5000')
# paths <- expand.grid(input_folder, TRTs, intensities, frequencies) %>%
#   mutate(foo = paste(Var2, Var3, Var4, sep = ' - ')) %>%
#   mutate(path = paste(Var1, Var3, Var4, Var2, sep = '/'))

# # base
# TRTs <- c('narrow', 'wide')
# # input_folder <- c('2025_12_21_base_n10_b5000')
# input_folder <- c('2026_01_22_base_red_noise_n10_b5000')
# paths <- expand.grid(input_folder, TRTs) %>%
#   mutate(foo = paste(Var1, Var2, sep = ' - ')) %>%
#   mutate(path = paste(Var1, Var2, sep = '/'))
# scenarios <- c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)
# osrs <- c(0.49, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05)

# # conservation AND evolution
# traits <- c('emergence_success_t0')
# rates <- c('effective')
# intensities <- paste('i', c(0.2, 0.4, 0.6, 0.8, 1), sep = '')
# frequencies <- paste('F', c(2, 4), sep = '')
# TRTs <- c('narrow')
# input_folder <- c('2025_12_21_both_n10_b5000')
# paths <- expand.grid(input_folder, TRTs, traits, rates, intensities, frequencies) %>%
#   mutate(foo = paste(Var2, Var3, Var4, Var5, Var6, sep = ' - ')) %>%
#   mutate(path = paste(Var1, Var3, Var4, Var5, Var6, Var2, sep = '/'))

# scenarios <- c(1, 2, 3, 4.5)
# osrs <- c(0.4, 0.3, 0.2, 0.05)

betas <- OSRs_to_betas(osrs)

# years to average over for lambdas
average_over <- 10

# cutoff point - at what persistence value should we not track lambda?
cutoff <- 0.1

# total years run per sim
years <- 1:100

# abundances / sexes
abundances <- c('Immature Females', 'Immature Males', 'Mature Females',
                'Mature Males', 'Hatchlings', 'Mature', 'Total')
indices <- list(c(1), c(2), c(3), c(4), c(1:2), c(3:4), c(1:4))

# dimensions
# P <- length(paths)
P <- nrow(paths)
S <- length(scenarios)
B <- length(betas)
Y <- length(years)
A <- length(abundances)

numrows <- P*S*B*A*Y

# maturity ogive
max_age <- 85
age_maturity_mu <- 25
age_maturity_sd <- 2.5
M <- pnorm(q = 1:max_age, mean = age_maturity_mu, sd = age_maturity_sd)

# other demographic parameters
temp_mu <- 31.8

# initialize super data frame
SDF <- data.frame(
  Input_folder = NULL,
  TRT = NULL,
  Trait = NULL,
  Rate = NULL,
  Intensity = NULL,
  Frequency = NULL,
  Scenario = NULL,
  Beta = NULL,
  OSR = NULL,
  Year = NULL,
  Abundance = NULL,
  Persist_mean = NULL,
  Temperature = NULL,
  Emergence_success_median = NULL, 
  Emergence_success_Q25 = NULL, 
  Emergence_success_Q75 = NULL, 
  Operational_sex_ratio_median = NULL,
  Operational_sex_ratio_Q25 = NULL, 
  Operational_sex_ratio_Q75 = NULL, 
  Abundance_median = NULL,
  Abundance_Q25 = NULL,
  Abundance_Q75 = NULL,
  Effective_median = NULL, 
  Effective_Q25 = NULL, 
  Effective_Q75 = NULL,
  Sex_Ratio_median = NULL,
  Sex_Ratio_Q25 = NULL,
  Sex_Ratio_Q75 = NULL,
  Breeding_success_median = NULL,
  Breeding_success_Q25 = NULL,
  Breeding_success_Q75 = NULL,
  Lambda_mean = NULL,
  Lambda_median = NULL,
  Lambda_Q25 = NULL,
  Lambda_Q75 = NULL,
  Lambda_10yr_mean = NULL,
  Lambda_10yr_median = NULL,
  Lambda_10yr_Q25 = NULL,
  Lambda_10yr_Q75 = NULL
)

# for each path (folder and model/TRT)
for (p in 1:P) {
  
  # for each temp increase scenario
  for (s in 1:S) {
    
    # temperature increases
    temps <- seq(from = temp_mu,
                 to = temp_mu + scenarios[s],
                 length = Y)
    
    # for each mating function
    for (b in 1:B) {
      
      # load in appropriate output files
      N <- paste(user, paths[p, ]$path, '/', scenarios[s], 'C/beta', betas[b],
                 '/', nsims, '_N.Rda', sep = '')
      
      Osr <- paste(user, paths[p, ]$path, '/', scenarios[s], 'C/beta', betas[b],
                   '/', nsims, '_OSR.Rda', sep = '')
      
      ES <- paste(user, paths[p, ]$path, '/', scenarios[s], 'C/beta', betas[b],
                  '/', nsims, '_ES.Rda', sep = '')
      
      # if the file exists
      if ( file.exists(N) & file.exists(Osr) & file.exists(ES) ) {

        load(N)
        
        load(Osr)
        
        load(ES)
        
        # iterate through and get output
        for (a in 1:A) {        
          
          # extract abundances
          
          if (a < 5) {
            
            N_abundances <- colSums(sims_N[indices[[a]], , , ], dim = 1)
            
          } else if (a == 5) {
            
            N_abundances <- sims_N[1, 1, , ] + sims_N[2, 1, , ]
            
          } else if (a == 6) {
            
            N_abundances <- colSums(sims_N[indices[[a]], , , ], dim = c(2))
            
            }
          
          # initial and final abundances to calculate mean persistence to year y
          # for abundance/sex a
          inits <- N_abundances[1, ]
          finals <- N_abundances
          
          # no sex ratios for most abundances / sexes
          sex_ratio_medians <- rep(NA, Y)
          sex_ratio_Q25s <- rep(NA, Y)
          sex_ratio_Q75s <- rep(NA, Y)
          
          # hatchlings
          if (a == 5) {
            
            F_hatchlings <- sims_N[1, 1, , ]
            M_hatchlings <- sims_N[2, 1, , ]
            sex_ratios <- M_hatchlings / (F_hatchlings + M_hatchlings)
            sex_ratio_medians <- rowMedians(sex_ratios, na.rm = TRUE)
            sex_ratio_Q25s <- rowQuantiles(sex_ratios, na.rm = TRUE, probs = c(0.25))
            sex_ratio_Q75s <- rowQuantiles(sex_ratios, na.rm = TRUE, probs = c(0.75))
            
            emergence_success_medians = rowMedians(sims_ES, na.rm = TRUE)
            emergence_success_Q25s = rowQuantiles(sims_ES, probs = 0.25, na.rm = TRUE)
            emergence_success_Q75s = rowQuantiles(sims_ES, probs = 0.75, na.rm = TRUE)
            
            
          } else {
            
            emergence_success_medians = rep(NA, Y)
            emergence_success_Q25s = rep(NA, Y)
            emergence_success_Q75s = rep(NA, Y)
            
          }
          
          # mature
          if (a == 6) {
            
            F_mature <- colSums(sims_N[3, , , ], dim = 1)
            M_mature <- colSums(sims_N[4, , , ], dim = 1)
            
            # adult sex ratios
            sex_ratios <- M_mature / (F_mature + M_mature)
            sex_ratio_medians <- rowMedians(sex_ratios, na.rm = TRUE)
            sex_ratio_Q25s <- rowQuantiles(sex_ratios, na.rm = TRUE, probs = c(0.25))
            sex_ratio_Q75s <- rowQuantiles(sex_ratios, na.rm = TRUE, probs = c(0.75))
            
            # operational sex ratios
            operational_sex_ratio_medians = rowMedians(sims_OSR, na.rm = TRUE)
            operational_sex_ratio_Q25s = rowQuantiles(sims_OSR, probs = 0.25, na.rm = TRUE)
            operational_sex_ratio_Q75s = rowQuantiles(sims_OSR, probs = 0.75, na.rm = TRUE)
            
            # breeding success
            breeding_success_medians <- pbeta(q = 2*sex_ratio_medians,
                                              shape1 = 1,
                                              shape2 = as.numeric(betas[b]))
            breeding_success_Q25s <- pbeta(q = 2*sex_ratio_Q25s,
                                           shape1 = 1,
                                           shape2 = as.numeric(betas[b]))
            breeding_success_Q75s <- pbeta(q = 2*sex_ratio_Q75s,
                                           shape1 = 1,
                                           shape2 = as.numeric(betas[b]))
            
            # sex ratio effective population size
            effective <- (4 * F_mature * M_mature) / (F_mature + M_mature)
            effective_medians <- rowMedians(effective, na.rm = TRUE)
            effective_Q25s <- rowQuantiles(effective, na.rm = TRUE, probs = 0.25)
            effective_Q75s <- rowQuantiles(effective, na.rm = TRUE, probs = 0.75)
            
          }  else {
            
            operational_sex_ratio_medians = rep(NA, Y)
            operational_sex_ratio_Q25s = rep(NA, Y)
            operational_sex_ratio_Q75s = rep(NA, Y)
            
            breeding_success_medians <- rep(NA, Y)
            breeding_success_Q25s <- rep(NA, Y)
            breeding_success_Q75s <- rep(NA, Y)
            
            effective_medians <- rep(NA, Y)
            effective_Q25s <- rep(NA, Y)
            effective_Q75s <- rep(NA, Y)
            
          }
          
          # abundances
          abundance_medians <- rowMedians(N_abundances, na.rm = TRUE)
          abundance_Q25s <- rowQuantiles(N_abundances, na.rm = TRUE, probs = c(0.25))
          abundance_Q75s <- rowQuantiles(N_abundances, na.rm = TRUE, probs = c(0.75))
          
          # calculate lambdas across all years
          lambdas <- N_abundances[2:Y, ] / N_abundances[1:(Y - 1), ]
          
          # replace infinite / NaN with NA
          lambdas[!is.finite(lambdas)] <- NA
          
          # lambdas with NA in front for year 1
          lambda_means <- c(NA, rowMeans(lambdas, na.rm = TRUE))
          lambda_medians <- c(NA, rowMedians(lambdas, na.rm = TRUE))
          lambda_Q25s <- c(NA, rowQuantiles(lambdas, prob = c(0.25),
                                            na.rm = TRUE))
          lambda_Q75s <- c(NA, rowQuantiles(lambdas, prob = c(0.75),
                                            na.rm = TRUE))
          
          # subset to add to super data frame
          sub_SDF <- data.frame(
            Folder = rep(paths[p, ]$Var1, Y),
            TRT = rep(paths[p, ]$Var2, Y),
            
            # # conservation
            # Trait = rep(NA, Y),
            # Rate = rep(NA, Y),
            # Intensity = rep(paths[p, ]$Var3, Y),
            # Frequency = rep(paths[p, ]$Var4, Y),
            
            # # evolution
            # Trait = rep(paths[p, ]$Var3, Y),
            # Rate = rep(paths[p, ]$Var4, Y),
            # Intensity = rep(NA, Y),
            # Frequency = rep(NA, Y),
            
            # base
            Trait = rep(NA, Y),
            Rate = rep(NA, Y),
            Intensity = rep(NA, Y),
            Frequency = rep(NA, Y),
            
            # # both
            # Trait = rep(paths[p, ]$Var3, Y),
            # Rate = rep(paths[p, ]$Var4, Y),
            # Intensity = rep(paths[p, ]$Var5, Y),
            # Frequency = rep(paths[p, ]$Var6, Y),
            
            Scenario = rep(scenarios[s], Y),
            Beta = rep(betas[b], Y),
            OSR = rep(osrs[b], Y),
            Year = years,
            Abundance = rep(abundances[a], Y),
            Persist_mean = rowMeans(finals >= (0.10 * inits), na.rm = TRUE),
            Temperature = temps,
            Emergence_success_median = emergence_success_medians,
            Emergence_success_Q25 = emergence_success_Q25s, 
            Emergence_success_Q75 = emergence_success_Q75s, 
            Operational_sex_ratio_median = operational_sex_ratio_medians,
            Operational_sex_ratio_Q25 = operational_sex_ratio_Q25s, 
            Operational_sex_ratio_Q75 = operational_sex_ratio_Q75s, 
            Abundance_median = abundance_medians,
            Abundance_Q25 = abundance_Q25s,
            Abundance_Q75 = abundance_Q75s,
            Effective_median = effective_medians, 
            Effective_Q25 = effective_Q25s, 
            Effective_Q75 = effective_Q75s, 
            Sex_ratio_median = sex_ratio_medians,
            Sex_ratio_Q25 = sex_ratio_Q25s,
            Sex_ratio_Q75 = sex_ratio_Q75s,
            Breeding_success_median = breeding_success_medians,
            Breeding_success_Q25 = breeding_success_Q25s,
            Breeding_success_Q75 = breeding_success_Q75s,
            Lambda_mean = lambda_means,
            Lambda_median = lambda_medians,
            Lambda_Q25 = lambda_Q25s,
            Lambda_Q75 = lambda_Q75s,
            Lambda_10yr_mean = rep(NA, Y),
            Lambda_10yr_median = rep(NA, Y),
            Lambda_10yr_Q25 = rep(NA, Y),
            Lambda_10yr_Q75 = rep(NA, Y)
          )
          
          # calculate lambdas averaged over previous 10 years
          lambda_10yr_means <- rep(NA, Y)
          lambda_10yr_medians <- rep(NA, Y)
          lambda_10yr_Q25s <- rep(NA, Y)
          lambda_10yr_Q75s <- rep(NA, Y)
          
          for (y in (average_over + 1):Y) {
            
            lambda_10yr_means[y] <- mean(lambdas[(y - average_over):(y - 1), ],
                                         na.rm = TRUE)
            
            lambda_10yr_medians[y] <- median(lambdas[(y - average_over):(y - 1), ],
                                             na.rm = TRUE)
            
            lambda_10yr_Q25s[y] <- quantile(lambdas[(y - average_over):(y - 1), ],
                                            na.rm = TRUE, probs = c(0.25))
            
            lambda_10yr_Q75s[y] <- quantile(lambdas[(y - average_over):(y - 1), ],
                                            na.rm = TRUE, probs = c(0.75))
            
          }
          
          # add 10 year average lambdas to subset df
          sub_SDF$Lambda_10yr_mean <- lambda_10yr_means
          sub_SDF$Lambda_10yr_median <- lambda_10yr_medians
          sub_SDF$Lambda_10yr_Q25 <- lambda_10yr_Q25s
          sub_SDF$Lambda_10yr_Q75 <- lambda_10yr_Q75s
          
          SDF <- rbind(SDF, sub_SDF)
          
          # print progress update
          prop <- round(nrow(SDF) / (numrows) * 100, 2)
          boop <- format(lubridate::now())
          
          print(paste(boop, ' - ', paths$foo[p], ' - ', scenarios[s], 
                      'C - beta ',  betas[b], ' - ', abundances[a], 
                      ' - done - ', prop, '% of total done!', sep = ''))
          
        }
        
      } else { print(paste(paths[p, ]$Var2, paths[p, ]$Var3, paths[p, ]$Var4, 
                           'file does not exist???', '\n', N, sep = ' - '))}
      
    }
    
  }
  
}

# save(SDF,
#      file = paste(user, folder, '_SDF_test.Rdata', sep = ''))

# adjust some stuff
all_outputs <- SDF %>%
  mutate(Emergence_Success = 0.86 / (1 + exp(1.7 * (Temperature - 32.7)))) %>%
  mutate(OSR = as.numeric(as.character(OSR))) %>%
  mutate(Mating_Function = if_else(OSR < 0.26, 'Steep', 'Shallow')) %>%
  mutate(Facet_label = paste(Abundance, ' abundance', sep = '')) %>%
  mutate(Beta = factor(Beta)) %>%
  mutate(OSR = factor(OSR)) %>%
  mutate(Scenario = factor(Scenario)) %>%
  mutate(Operational_sex_ratio_median = case_when(Persist_mean < cutoff ~ NA, TRUE ~ Operational_sex_ratio_median)) %>%
  mutate(Operational_sex_ratio_Q25 = case_when(Persist_mean < cutoff ~ NA, TRUE ~ Operational_sex_ratio_Q25)) %>%
  mutate(Operational_sex_ratio_Q75 = case_when(Persist_mean < cutoff ~ NA, TRUE ~ Operational_sex_ratio_Q75)) %>%
  mutate(Abundance_median = case_when(Persist_mean < cutoff ~ NA, TRUE ~ Abundance_median)) %>%
  mutate(Abundance_Q25 = case_when(Persist_mean < cutoff ~ NA, TRUE ~ Abundance_Q25)) %>%
  mutate(Abundance_Q75 = case_when(Persist_mean < cutoff ~ NA, TRUE ~ Abundance_Q75)) %>%
  mutate(Sex_ratio_median = case_when(Persist_mean < cutoff ~ NA, TRUE ~ Sex_ratio_median)) %>%
  mutate(Sex_ratio_Q25 = case_when(Persist_mean < cutoff ~ NA, TRUE ~ Sex_ratio_Q25)) %>%
  mutate(Sex_ratio_Q75 = case_when(Persist_mean < cutoff ~ NA, TRUE ~ Sex_ratio_Q75)) %>%
  mutate(Breeding_success_median = case_when(Persist_mean < cutoff ~ NA, TRUE ~ Breeding_success_median)) %>%
  mutate(Breeding_success_Q25 = case_when(Persist_mean < cutoff ~ NA, TRUE ~ Breeding_success_Q25)) %>%
  mutate(Breeding_success_Q75 = case_when(Persist_mean < cutoff ~ NA, TRUE ~ Breeding_success_Q75)) %>%
  mutate(Lambda_median = case_when(Persist_mean < cutoff ~ NA, TRUE ~ Lambda_median)) %>%
  mutate(Lambda_Q25 = case_when(Persist_mean < cutoff ~ NA, TRUE ~ Lambda_Q25)) %>%
  mutate(Lambda_Q75 = case_when(Persist_mean < cutoff ~ NA, TRUE ~ Lambda_Q75)) %>%
  mutate(Lambda_10yr_mean = case_when(Persist_mean < cutoff ~ NA, TRUE ~ Lambda_10yr_mean)) %>%
  mutate(Lambda_10yr_median = case_when(Persist_mean < cutoff ~ NA, TRUE ~ Lambda_10yr_median)) %>%
  mutate(Lambda_10yr_Q25 = case_when(Persist_mean < cutoff ~ NA, TRUE ~ Lambda_10yr_Q25)) %>%
  mutate(Lambda_10yr_Q75 = case_when(Persist_mean < cutoff ~ NA, TRUE ~ Lambda_10yr_Q75))

# name to save to
name <- paste(input_folder, '_n', nsims, sep = '')

# save to file
save(all_outputs,
     file = paste('../output/', name, '_all_outputs.Rdata', 
                  sep = ''))
