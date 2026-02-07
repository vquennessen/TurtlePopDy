# reproduction

reproduction <- function(N, M, y, beta, max_age, 
                         F_remigration_int, M_remigration_int,
                         clutches_mu, clutches_sd, eggs_mu, eggs_sd, 
                         emergence_success_A, emergence_success_k, 
                         emergence_success_t0, 
                         season_temp_mus, clutch_temp_sd,
                         k_piv, T_piv, T_threshold, 
                         evolve, trait, male_probs, contributions,
                         h2, varGenetic, varPhenotypic, 
                         G, P, G_stats, P_stats,
                         conserve, conservation_years, intensity, effect_size) {
  
  # breeding females this year
  n_available_F <- sum(rbinom(n = max_age, 
                              size = N[3, , y], 
                              prob = 1 / F_remigration_int), na.rm = TRUE)
  
  # breeding males this year
  n_available_M <- sum(rbinom(n = max_age, 
                              size = N[4, , y], 
                              prob = 1 / M_remigration_int), na.rm = TRUE)
  
  # check that there is at least one available female and one available male
  if (n_available_F < 1 | n_available_M < 1) {
    
    OSR <- NA
    female_hatchlings <- 0
    male_hatchlings <- 0 
    hatchlings <- 0
    PSR <- NA
    ES <- NA
    
    if (evolve == TRUE) {
      G[1, 1] <- list(NA)
      P[1, 1] <- list(NA)
      G_stats[1, 1, y, ] <- rep(NA, 5)
      P_stats[1, 1, y, ] <- rep(NA, 5)
    } 
    
  } else {
    
    # operational sex ratio - proportion of males
    # multiply by 2 to transform to beta function with x from 0 to 0.5 instead 
    # of 0 to 1
    OSR <- n_available_M / (n_available_M + n_available_F)
    
    # calculate reproductive success
    # use beta function to calculate breeding success multiply OSR by 2 to 
    # transform to beta function with x from 0 to 0.5 instead of 0 to 1
    breeding_success <- pbeta(2 * OSR, 
                              shape1 = 1, 
                              shape2 = beta) 
    
    # how many females actually find a male to mate with and then nest
    # set.seed(seed)
    n_breeding_F <- sum(rbinom(n = n_available_F, 
                               size = 1, 
                               prob = breeding_success), na.rm = TRUE)
    
    # n_breeding_F <- sum(breeding_F)
    
    if (n_breeding_F < 1) {
      
      female_hatchlings <- 0
      male_hatchlings <- 0 
      hatchlings <- 0
      PSR <- NA
      ES <- NA
      
      if (evolve == TRUE) {
        G[1, 1] <- list(NA)
        P[1, 1] <- list(NA)
        G_stats[1, 1, y, ] <- rep(NA, 5)
        P_stats[1, 1, y, ] <- rep(NA, 5)
      } 
      
    } else {
      
      # vector of number of clutches per female (round to nearest integer)
      # set.seed(seed)
      clutches <- round(rnorm(n = n_breeding_F, 
                              mean = clutches_mu, 
                              sd = clutches_sd)) 
      
      # replace any number < 1 with +1
      clutches[which(clutches < 1)] <- 1
      
      # eggs list, one number for each clutch
      # set.seed(seed)
      eggs <- lapply(clutches,
                     rnorm,
                     mean = eggs_mu,
                     sd = eggs_sd) %>%
        lapply(pmax, 0) %>%
        lapply(round)
      
      # list of clutch temperatures, one number for each clutch 
      # set.seed(seed)
      clutch_temps <- lapply(clutches, 
                             rnorm, 
                             mean = season_temp_mus[y], 
                             sd = clutch_temp_sd) %>%
        lapply(pmax, 0)
      
      # adjust clutch temperatures based on conservation measures
      if (conserve == TRUE & y %in% conservation_years) {
        
        clutch_temps <- conservation(initial_temps = clutch_temps, 
                                     intensity = intensity, 
                                     effect_size = effect_size)
        
      }
      
      # if evolution
      if (evolve == TRUE) {
        
        evo_output <- evolution(max_age, G, P, 
                                n_breeding_F, n_available_M, 
                                trait, male_probs, contributions,
                                h2, varGenetic, varPhenotypic,
                                clutches, eggs, clutch_temps, 
                                emergence_success_A, emergence_success_k, 
                                emergence_success_t0, 
                                T_threshold, k_piv, T_piv)
        
        G[1, 1] <- evo_output[[1]]
        G[2, 1] <- evo_output[[2]]
        
        P[1, 1] <- evo_output[[3]]
        P[2, 1] <- evo_output[[4]]
        
        G_stats[1:2, 1, y, ] <- evo_output[[5]]
        P_stats[1:2, 1, y, ] <- evo_output[[6]]
        
        # numbers of hatchlings produced
        female_hatchlings <- length(unlist(G[1, 1]))
        male_hatchlings <- length(unlist(G[2, 1]))
        
        # emergence success
        ES <- (female_hatchlings + male_hatchlings) / sum(unlist(eggs))
        
      } else {
        
        # list of probability of emergence, one number for each clutch 
        # set.seed(seed)
        probs_emerged <- lapply(
          clutch_temps, 
          function(x) {
            a <- emergence_success_A / (
              1 + exp(-emergence_success_k * (x - emergence_success_t0)))
            a[x > T_threshold] <- 0
            return(a)
          }
        ) %>% lapply(pmax, 0)
        
        # list of probabilities of developing as male, one for each clutch
        # set.seed(seed)
        probs_male <- lapply(clutch_temps, 
                             function(x) {
                               1 / (1 + exp(-k_piv * (x - T_piv)))
                             }) %>%
          lapply(pmax, 0)
        
        # list of numbers of emerged hatchlings, one for each clutch
        # set.seed(seed)
        hatchlings <- map2(eggs, probs_emerged,
                           ~ unlist(map2(.x, .y,
                                         ~ sum(rbinom(n = .x, size = 1, prob = .y),
                                               na.rm = TRUE))))
        
        # list of number of males, one for each clutch
        males <- map2(hatchlings, probs_male,
                      ~ unlist(map2(.x, .y,
                                    ~ sum(rbinom(n = .x, size = 1, prob = .y),
                                          na.rm = TRUE))))
        
        # list of number of females, one for each clutch
        females <- Map(`-`, hatchlings, males) 
        
        # total number of male hatchlings
        male_hatchlings <- sum(unlist(males))
        
        # total number of female hatchlings
        female_hatchlings <- sum(unlist(females))
        
        # emergence success
        ES <- (female_hatchlings + male_hatchlings) / sum(unlist(eggs))
        
      }
      
    }
    
  }
  
  # output
  output <- list(OSR, ES, female_hatchlings, male_hatchlings, 
                 G, P, G_stats, P_stats)
  
  return(output)
  
}
