# extract hatchling phenotype values

rm(list = ls())

# set working directory
# setwd('~/Projects/iliketurtles3/code')

# source functions
source('mating function/OSRs_to_betas.R')

# load libraries
library(ggplot2)
library(matrixStats)
library(dplyr)
library(tidyr)
library(abind)
library(zoo)

##### to modify ################################################################

# which computer we using
computer <- 'cluster'

# path based on computer being used
user <- ifelse(computer == 'cluster', '/home/quennessenv/iliketurtles3/output/',
               ifelse(computer == 'desktop',
                      'C:/Users/Vic/Box/Quennessen_Thesis/PhD Thesis/model output/iliketurtles3/',
                      'C:/Users/vique/Box/Quennessen_Thesis/PhD Thesis/model output/iliketurtles3/'))

# name of folder for current runs
input_folders <- c('2025_12_21_evolution_n10_b5000')

# number of sims
nsims <- 1000

# name to save to
name <- paste(input_folders, '_n', nsims, sep = '')

# model(s)
traits <- c('emergence_success_t0', 'T_piv')
rates <- c('effective', 'high')
TRTs <- c('narrow')
intensities <- 'i0.2'
frequencies <- 'F1'

paths <- expand.grid(input_folders, TRTs, traits, rates) %>%
  mutate(path = paste(Var1, Var3, Var4, Var2, sep = '/'))

# # conservation AND evolution
# traits <- c('emergence_success_t0')
# rates <- c('effective', 'high')
# intensities <- paste('i', c(0.2, 0.4, 0.6, 0.8, 1), sep = '')
# frequencies <- paste('F', c(2, 4, 6, 8, 10), sep = '')
# TRTs <- c('narrow')
# input_folder <- c('2025_12_11_evolution_conservation_n1_b5000')
# paths <- expand.grid(input_folder, TRTs, traits, rates, intensities, frequencies) %>%
#   mutate(path = paste(Var1, Var3, Var4, Var5, Var6, Var2, sep = '/'))

# years to average over
average_over <- 10
years <- 1:100

# temperature increase scenarios
scenarios <- c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)
# scenarios <- c(0.5, 4.5)

# operational sex ratios / betas
osrs <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.49)
# osrs <- c(0.1, 0.35)
betas <- OSRs_to_betas(osrs)

# dimensions
P <- nrow(paths)
S <- length(scenarios)
B <- length(betas)
Y <- length(years)
numrows <- P*S*B*Y

# initialize super data frame
SDF <- data.frame(Folder = NULL,
                  TRT = NULL,
                  Trait = NULL, 
                  Rate = NULL, 
                  Scenario = NULL, 
                  OSR = NULL, 
                  Year = NULL,
                  G_median = NULL,
                  G_Q25 = NULL,
                  G_Q75 = NULL,
                  G_mean = NULL,
                  G_var = NULL, 
                  P_median = NULL,
                  P_Q25 = NULL,
                  P_Q75 = NULL,
                  P_mean = NULL,
                  P_var = NULL)

for (p in 1:P) {
  
  for (s in 1:S) {
    
    for (b in 1:B) {
      
      # initialize empty dataframe, one for each filepath
      sub_DF <- data.frame(Folder = rep(paths[p, ]$Var1, Y),
                           TRT = rep(paths[p, ]$Var2, Y),
                           Trait = rep(paths[p, ]$Var3, Y),
                           Rate = rep(paths[p, ]$Var4, Y),
                           Scenario = rep(scenarios[s], Y),
                           Beta = rep(betas[b], Y),
                           Year = years, 
                           G_median = rep(NA, Y),
                           G_Q25 = rep(NA, Y),
                           G_Q75 = rep(NA, Y),                           
                           G_mean = rep(NA, Y),
                           G_var = rep(NA, Y), 
                           P_median = rep(NA, Y),
                           P_Q25 = rep(NA, Y),
                           P_Q75 = rep(NA, Y),
                           P_mean = rep(NA, Y),
                           P_var = rep(NA, Y))
      
      # load in appropriate output file
      
      G_stats <- paste(user, paths[p, ]$path, '/', scenarios[s], 'C/beta', 
                            betas[b], '/', nsims, '_G_stats.Rda', sep = '')
      
      P_stats <- paste(user, paths[p, ]$path, '/', scenarios[s], 'C/beta', 
                            betas[b], '/', nsims, '_P_stats.Rda', sep = '')
      
      # if the file exists - cluster
      if (file.exists(G_stats) & file.exists(P_stats))  {
        
        # load in genetics stats
        load(G_stats)
        
        # load in phenotypes stats
        load(P_stats)
        
      }
      
      # extract mean g stats for hatchlings by year, across simulations
      # GF_medians <- rowMedians(sims_G_stats[1, 1, , 1, ], na.rm = TRUE)
      # GM_medians <- rowMedians(sims_G_stats[2, 1, , 1, ], na.rm = TRUE)
      # sub_DF$G_median <- rowMeans(data.frame(X = GF_medians, Y = GM_medians), 
      #                               na.rm = TRUE) 
      sub_DF$G_median <- rowMedians(cbind(sims_G_stats[1, 1, , 1, ], 
                                          sims_G_stats[2, 1, , 1, ]), 
                                    na.rm = TRUE) 
      
      # GF_Q25s <- rowMedians(sims_G_stats[1, 1, , 2, ], na.rm = TRUE)
      # GM_Q25s <- rowMedians(sims_G_stats[2, 1, , 2, ], na.rm = TRUE)
      # sub_DF$G_Q25 <- rowMeans(data.frame(X = GF_Q25s, Y = GM_Q25s), 
      #                             na.rm = TRUE) 
      sub_DF$G_Q25 <- rowMedians(cbind(sims_G_stats[1, 1, , 2, ], 
                                       sims_G_stats[2, 1, , 2, ]), na.rm = TRUE) 
      
      # GF_Q75s <- rowMedians(sims_G_stats[1, 1, , 3, ], na.rm = TRUE)
      # GM_Q75s <- rowMedians(sims_G_stats[2, 1, , 3, ], na.rm = TRUE)
      # sub_DF$G_Q75 <- rowMeans(data.frame(X = GF_Q75s, Y = GM_Q75s), 
      #                             na.rm = TRUE) 
      sub_DF$G_Q75 <- rowMedians(cbind(sims_G_stats[1, 1, , 3, ], 
                                       sims_G_stats[2, 1, , 3, ]), na.rm = TRUE) 
      
      GF_means <- rowMeans(sims_G_stats[1, 1, , 4, ], na.rm = TRUE)
      GM_means <- rowMeans(sims_G_stats[2, 1, , 4, ], na.rm = TRUE)
      sub_DF$G_mean <- rowMeans(data.frame(X = GF_means, Y = GM_means), 
                                na.rm = TRUE)
      
      GF_vars <- rowMeans(sims_G_stats[1, 1, , 5, ], na.rm = TRUE)
      GM_vars <- rowMeans(sims_G_stats[2, 1, , 5, ], na.rm = TRUE)
      sub_DF$G_var <- rowMeans(data.frame(X = GF_vars, GY = GM_vars), 
                               na.rm = TRUE)
      
      # extract mean p stats for hatchlings by year, across simulations
      # PF_medians <- rowMedians(sims_P_stats[1, 1, , 1, ], na.rm = TRUE)
      # PM_medians <- rowMedians(sims_P_stats[2, 1, , 1, ], na.rm = TRUE)
      # sub_DF$P_median <- rowMeans(data.frame(X = PF_medians, Y = PM_medians), 
      #                             na.rm = TRUE) 
      sub_DF$P_median <- rowMedians(cbind(sims_P_stats[1, 1, , 1, ], 
                                          sims_P_stats[2, 1, , 1, ]), 
                                    na.rm = TRUE) 
      
      # PF_Q25s <- rowMedians(sims_P_stats[1, 1, , 2, ], na.rm = TRUE)
      # PM_Q25s <- rowMedians(sims_P_stats[2, 1, , 2, ], na.rm = TRUE)
      # sub_DF$P_Q25 <- rowMeans(data.frame(X = PF_Q25s, Y = PM_Q25s), 
      #                          na.rm = TRUE) 
      sub_DF$P_Q25 <- rowMedians(cbind(sims_P_stats[1, 1, , 2, ], 
                                       sims_P_stats[2, 1, , 2, ]), na.rm = TRUE) 
      
      # PF_Q75s <- rowMedians(sims_P_stats[1, 1, , 3, ], na.rm = TRUE)
      # PM_Q75s <- rowMedians(sims_P_stats[2, 1, , 3, ], na.rm = TRUE)
      # sub_DF$P_Q75 <- rowMeans(data.frame(X = PF_Q75s, Y = PM_Q75s), 
      #                          na.rm = TRUE) 
      sub_DF$P_Q75 <- rowMedians(cbind(sims_P_stats[1, 1, , 3, ], 
                                       sims_P_stats[2, 1, , 3, ]), na.rm = TRUE) 
      
      PF_means <- rowMeans(sims_P_stats[1, 1, , 4, ], na.rm = TRUE)
      PM_means <- rowMeans(sims_P_stats[2, 1, , 4, ], na.rm = TRUE)
      sub_DF$P_mean <- rowMeans(data.frame(X = PF_means, Y = PM_means), 
                                na.rm = TRUE)
      
      PF_vars <- rowMeans(sims_P_stats[1, 1, , 5, ], na.rm = TRUE)
      PM_vars <- rowMeans(sims_P_stats[2, 1, , 5, ], na.rm = TRUE)
      sub_DF$G_var <- rowMeans(data.frame(X = PF_vars, Y = PM_vars), 
                               na.rm = TRUE)
      
      # add DF to SDF
      SDF <- rbind(SDF, sub_DF)
      
      prop <- round(nrow(SDF) / (numrows) * 100, 2)
      boop <- format(lubridate::now())
      print(paste(boop, ' - ', paths[p, ]$Var2, ' - ', paths[p, ]$Var3, 
                  ' - ', paths[p, ]$Var4, ' - ', scenarios[s],
                  'C - beta ', betas[b], ' - done - ', prop, 
                  '% of total done!', sep = ''))
      
    }
    
  }
  
}

traits <- SDF %>%
  mutate(G_median_10yr = rollapply(G_median, width = average_over, median, 
                                   na.rm = TRUE, fill = NA, align = 'right')) %>% 
  mutate(G_Q25_10yr = rollapply(G_Q25, width = average_over, median, 
                                   na.rm = TRUE, fill = NA, align = 'right')) %>% 
  mutate(G_Q75_10yr = rollapply(G_Q75, width = average_over, median, 
                                na.rm = TRUE, fill = NA, align = 'right')) %>% 
  mutate(G_mean_10yr = rollapply(G_mean, width = average_over, mean, 
                                 na.rm = TRUE, fill = NA, align = 'right')) %>%
  mutate(G_var_10yr = rollapply(G_var, width = average_over, mean, 
                                 na.rm = TRUE, fill = NA, align = 'right')) %>%
  
  mutate(P_median_10yr = rollapply(P_median, width = average_over, median, 
                                   na.rm = TRUE, fill = NA, align = 'right')) %>% 
  mutate(P_Q25_10yr = rollapply(P_Q25, width = average_over, median, 
                                na.rm = TRUE, fill = NA, align = 'right')) %>% 
  mutate(P_Q75_10yr = rollapply(P_Q75, width = average_over, median, 
                                na.rm = TRUE, fill = NA, align = 'right')) %>% 
  mutate(P_mean_10yr = rollapply(P_mean, width = average_over, mean, 
                                 na.rm = TRUE, fill = NA, align = 'right')) %>%
  mutate(P_var_10yr = rollapply(P_var, width = average_over, mean, 
                                na.rm = TRUE, fill = NA, align = 'right'))

# save dataframe as R object
save(traits, 
     file = '../output/evolution_trait_values.Rdata')
