evolution <- function(max_age, G, P, 
                      n_breeding_F, n_available_M, 
                      trait, male_probs, contributions, 
                      h2, varGenetic, varPhenotypic,
                      clutches, eggs, clutch_temps, 
                      emergence_success_A, emergence_success_k, 
                      emergence_success_t0, 
                      T_threshold, k_piv, T_piv) {
  
  # extract maternal genotypes
  GM <- as.list(resample(unlist(G[3, 2:max_age]), size = n_breeding_F))
  
  # extract potential paternal genotypes
  potential_GP <- resample(unlist(G[4, 2:max_age]), size = n_available_M)       
  
  # how many males does each female mate with
  nMales <- as.list(resample(1:length(male_probs), 
                             size = n_breeding_F, 
                             prob = male_probs, 
                             replace = TRUE))
  
  # if there are more males assigned to a female than there are available, 
  # reduce it with the maximum number of males available
  nMales[nMales > n_available_M] <- n_available_M
  
  # assign male genotypes to each female
  GP <- map(nMales, ~ resample(potential_GP, size = .x))
  
  # assign male genotypes to eggs
  GP_eggs <- map2(GP, eggs, 
                  ~ lapply(.y, function(x) {
                    resample(.x, 
                             size = x, 
                             prob = contributions[[length(.x)]], 
                             replace = TRUE) } ))
  
  # assign female genotypes to eggs
  GM_eggs <- map2(GM, eggs, 
                  ~ lapply(.y, function(x) {
                    rep(.x, times = x)
                  }))
  
  # egg genotypes
  G_eggs <- map2(GP_eggs, GM_eggs, 
                 ~ map2(.x, .y, 
                        ~ (.x + .y) / 2 + rnorm(n = length(.x), 
                                                mean = 0, 
                                                sd = sqrt(varGenetic / 2))))
  
  # egg phenotypes
  P_eggs <- lapply(G_eggs, function(x) {
    lapply(x, function(x) {
      rnorm(n = length(x), 
            mean = x, 
            sd = sqrt(varGenetic * (1 - h2) / h2)) } ) } )
  
  if (trait == 'emergence_success_t0') {
    
    # list of probability of emergence, one for each egg 
    probs_emerged <- map2(P_eggs, clutch_temps, 
                          ~ map2(.x, .y, 
                                 ~ if (.y < T_threshold) {
                                   emergence_success_A / (
                                     1 + exp(-emergence_success_k * (.y - .x)))
                                   } else { 0 } )) 
    
  } else {
    
    probs_emerged <- lapply(
      clutch_temps, 
      function(x) {
        unlist(lapply(x, function(x) {
        if (x < T_threshold) {
          emergence_success_A / (
            1 + exp(-emergence_success_k * (x - emergence_success_t0)))
        } else { 0 } } ) ) } ) %>% 
      lapply(pmax, 0)
    
  }
  
  # which eggs emerge as hatchlings?
  indices_hatchlings <- map2(eggs, probs_emerged, 
                             ~ map2(.x, .y, 
                                    ~ as.logical(rbinom(n = .x, 
                                                        size = 1, 
                                                        prob = .y))))
  
  # how many hatchlings are there?
  hatchlings <- lapply(indices_hatchlings, 
                       function(x) { unlist(lapply(x, sum, na.rm = TRUE)) } )
  
  # hatchling genotypes and phenotypes
  G_hatchlings <- map2(G_eggs, indices_hatchlings, 
                       ~ map2(.x, .y, ~ .x[as.logical(.y)]))
  P_hatchlings <- map2(P_eggs, indices_hatchlings, 
                       ~ map2(.x, .y, ~ .x[as.logical(.y)]))
  
  if (trait == 'T_piv') {
    
    # probability of developing as male, one for each egg
    probs_male <- map2(P_hatchlings, clutch_temps, 
                       ~ map2(.x, .y, 
                              ~ 1 / (1 + exp(-k_piv * (.y - .x)))))
  } else {
        
    # list of probability of developing as male, one for each clutch 
    probs_male <- lapply(
      clutch_temps, function(x) {
        unlist(lapply(x, function(x) {
          1 / (1 + exp(-k_piv * (x - T_piv)))
          } ) ) } ) %>% 
      lapply(pmax, 0)
    
  }
  
  # which hatchlings developed as male?
  indices_males <- map2(hatchlings, probs_male, 
                        ~ map2(.x, .y, 
                               ~ as.logical(rbinom(n = .x, 
                                                   size = 1, 
                                                   prob = .y))))
  
  # which hatchlings developed as female?
  indices_females <- lapply(indices_males, function(x) {
    lapply(x, function(x) { as.logical(Map(`-`, 1, x)) } ) } )
  
  # female genotypes
  G_females <- map2(G_hatchlings, indices_females, 
                    ~ map2(.x, .y, ~ .x[as.logical(.y)]))
  
  # male genotypes
  G_males <- map2(G_hatchlings, indices_males, 
                    ~ map2(.x, .y, ~ .x[as.logical(.y)]))
  
  # female phenotypes
  P_females <- map2(P_hatchlings, indices_females, 
                    ~ map2(.x, .y, ~ .x[as.logical(.y)]))
  
  # male phenotypes
  P_males <- map2(P_hatchlings, indices_males, 
                  ~ map2(.x, .y, ~ .x[as.logical(.y)]))
  
  # assign hatchling genotypes
  EGF <- list(unlist(lapply(G_females, unlist)))
  EGM <- list(unlist(lapply(G_males, unlist)))
  
  # assign hatchling phenotypes
  EPF <- list(unlist(lapply(P_females, unlist)))
  EPM <- list(unlist(lapply(P_males, unlist)))
  
  # G stats values
  G_medians <- lapply(c(EGF, EGM), function(x) median(unlist(x), na.rm = TRUE))
  G_Q25s <- lapply(c(EGF, EGM), function(x) quantile(unlist(x), prob = 0.25, na.rm = TRUE))
  G_Q75s <- lapply(c(EGF, EGM), function(x) quantile(unlist(x), prob = 0.75, na.rm = TRUE))
  G_means <- lapply(c(EGF, EGM), function(x) mean(unlist(x), na.rm = TRUE))    
  G_variances <- lapply(c(EGF, EGM), function(x) var(unlist(x), na.rm = TRUE))
  
  # P stats values
  P_medians <- lapply(c(EPF, EPM), function(x) median(unlist(x), na.rm = TRUE))
  P_Q25s <- lapply(c(EPF, EPM), function(x) quantile(unlist(x), prob = 0.25, na.rm = TRUE))
  P_Q75s <- lapply(c(EPF, EPM), function(x) quantile(unlist(x), prob = 0.75, na.rm = TRUE))
  P_means <- lapply(c(EPF, EPM), function(x) mean(unlist(x), na.rm = TRUE))
  P_variances <- lapply(c(EPF, EPM), function(x) var(unlist(x), na.rm = TRUE))
  
  EG_stats <- unlist(c(G_medians, G_Q25s, G_Q75s, G_means, G_variances))
  
  EP_stats <- unlist(c(P_medians, P_Q25s, P_Q75s, P_means, P_variances))
  
  ##### output #################################################################
  
  output <- list(EGF, EGM, EPF, EPM, EG_stats, EP_stats)
  
  return(output)
  
}
