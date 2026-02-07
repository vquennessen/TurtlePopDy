# base model

base_model <- function(scenario, beta, yrs, max_age,
                       IF_survival, IM_survival, MF_survival, MM_survival,
                       IF_init, IM_init, MF_init, MM_init,
                       M, F_remigration_int, M_remigration_int,
                       clutches_mu, clutches_sd, eggs_mu, eggs_sd,
                       emergence_success_A, emergence_success_k,
                       emergence_success_t0, T_piv, k_piv, T_threshold, 
                       evolve, trait, value, male_probs, contributions,
                       h2, varGenetic, varPhenotypic, 
                       temp_mu, climate_stochasticity,
                       season_temp_sd, clutch_temp_sd, noise, AC,
                       conserve, frequency, intensity, effect_size) {
  
  ##### source initialized arrays ##############################################
  # seed <- 1325633
  # set.seed(seed)
  init_output <- initialize_arrays(scenario, yrs, max_age, 
                                   IF_init, IM_init, MF_init, MM_init,
                                   M, F_remigration_int, M_remigration_int, 
                                   T_piv, k_piv, T_threshold, 
                                   temp_mu, climate_stochasticity, 
                                   season_temp_sd, clutch_temp_sd, noise, AC, 
                                   evolve, trait, value, 
                                   varGenetic, varPhenotypic,
                                   conserve, frequency)
  
  N                  <- init_output[[1]]   # population size array
  season_temp_mus    <- init_output[[2]]   # mean temps at the season level
  OSRs               <- init_output[[3]]   # operational sex ratio
  ESs                <- init_output[[4]]
  G                  <- init_output[[5]]
  P                  <- init_output[[6]]
  G_stats            <- init_output[[7]]   # genotype stats, to keep
  P_stats            <- init_output[[8]]   # phenotype stats, to keep  
  conservation_years <- init_output[[9]]   # years for conservation action
  
  ##### model ##################################################################
  # set.seed(seed)
  for (y in 2:yrs) {
    
    # population dynamics
    # survival for each age 
    # set.seed(seed)
    popdy_output <- pop_dynamics(N, max_age, y, M,
                                 IF_survival, IM_survival, 
                                 MF_survival, MM_survival, 
                                 evolve, G, P, G_stats, P_stats)
    
    N       <- popdy_output[[1]]
    G       <- popdy_output[[2]]
    P       <- popdy_output[[3]]
    G_stats <- popdy_output[[4]]
    P_stats <- popdy_output[[5]]
    
    # reproduction
    # set.seed(seed)
    rep_output <- reproduction(N, M, y, beta, max_age, 
                               F_remigration_int, M_remigration_int,
                               clutches_mu, clutches_sd, eggs_mu, eggs_sd, 
                               emergence_success_A, emergence_success_k, 
                               emergence_success_t0, 
                               season_temp_mus, clutch_temp_sd,
                               k_piv, T_piv, T_threshold, 
                               evolve, trait, male_probs, contributions,
                               h2, varGenetic, varPhenotypic, 
                               G, P, G_stats, P_stats,
                               conserve, conservation_years, 
                               intensity, effect_size)
    
    # add recruits to population size array
    OSRs[y]    <- rep_output[[1]]
    ESs[y]     <- rep_output[[2]]
    N[1, 1, y] <- rep_output[[3]]
    N[2, 1, y] <- rep_output[[4]]
    G          <- rep_output[[5]]
    P          <- rep_output[[6]]
    G_stats    <- rep_output[[7]]
    P_stats    <- rep_output[[8]]
    
    # break out of loop if there are zero males or females at any age
    if (sum(N[c(1, 3), , y], na.rm = TRUE) < 1 || 
        sum(N[c(2, 4), , y], na.rm = TRUE) < 1) {
      
      # print('ran out of females or males!')
      
      break 
    }
    
    # print(y)
    
  }
  
  ##### output #################################################################
  
  if (evolve == TRUE) {
    
    output <- list(N, 
                   OSRs,
                   ESs,
                   G_stats, 
                   P_stats)
  } else {  
    output <- list(N, 
                   OSRs, 
                   ESs) 
  }
  
  return(output)
  
}
