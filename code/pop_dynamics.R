# population dynamics with demographic stochasticity

pop_dynamics <- function(N, max_age, y, M,
                         IF_survival, IM_survival, MF_survival, MM_survival,
                         evolve, G, P, G_stats, P_stats) {
  
  ##############################################################################
  ##### population dynamics ####################################################
  ##############################################################################
  
  #########################
  ##### females ###########
  #########################
  
  # which immature females survived
  survived_IF <- map2(N[1, 1:(max_age - 1), y - 1], 
                      IF_survival[1:(max_age - 1)], 
                      ~ rbinom(n = .x, size = 1, prob = .y))
  
  # number of immature females that survived
  n_sIF <- unlist(lapply(survived_IF, sum))
  
  # which females matured
  matured_IF <- map2(n_sIF, M[1:(max_age - 1)], 
                     ~ as.logical(rbinom(n = .x, size = 1, prob = .y)))
  
  # number of females that matured
  n_mIF <- unlist(lapply(matured_IF, sum))   
  
  # which females did not mature
  not_matured_IF <- map(matured_IF, ~ as.logical(Map(`-`, 1, .x)))
  
  # number of females that didn't mature
  n_nmIF <- unlist(lapply(not_matured_IF, sum))  
  
  # new immature females: number that survived - number that matured
  N[1, 2:max_age, y] <- n_nmIF
  
  # mature females that survived
  survived_MF <- map(N[3, 1:(max_age - 1), y - 1], 
                     ~ rbinom(n = .x, size = 1, prob = MF_survival))
  
  # number of mature females that survived
  n_sMF <- unlist(lapply(survived_MF, sum))
  
  # new mature females: 
  # number of immature females that matured + mature females that survived
  N[3, 2:max_age, y] <- n_mIF + n_sMF
  
  #########################
  ##### males #############
  #########################
  
  # which immature males survived
  survived_IM <- map2(N[2, 1:(max_age - 1), y - 1], 
                      IM_survival[1:(max_age - 1)], 
                      ~ rbinom(n = .x, size = 1, prob = .y))

  
  # number of immature males that survived
  n_sIM <- unlist(lapply(survived_IM, sum, na.rm = TRUE))
  
  # which males matured
  matured_IM <- map2(n_sIM, M[1:(max_age - 1)], 
                     ~ as.logical(rbinom(n = .x, size = 1, prob = .y)))
  
  # number of males that matured
  n_mIM <- unlist(lapply(matured_IM, sum, na.rm = TRUE))   
  
  # which males did not mature
  not_matured_IM <- map(matured_IM, ~ as.logical(Map(`-`, 1, .x)))
  
  # number of males that did not mature
  n_nmIM <- unlist(lapply(not_matured_IM, sum, na.rm = TRUE))
  
  # new immature males
  N[2, 2:max_age, y] <- n_nmIM
  
  # mature males that survived
  survived_MM <- map(N[4, 1:(max_age - 1), y - 1], 
                     ~ rbinom(n = .x, size = 1, prob = MM_survival)
                     )
  
  # number of mature males that survived
  n_sMM <- unlist(lapply(survived_MM, sum))
  
  # new mature males
  # number of immature males that matured + mature males that survived
  N[4, 2:max_age, y] <- n_mIM + n_sMM
  
  ##############################################################################
  ##### evolution ##############################################################
  ##############################################################################
  
  if (evolve == TRUE) {
    
    # initialize new G and P
    G_new <- asplit(array(rep(NA, 4 * max_age * 1), dim = c(4, max_age, 1)),
                    c(1, 2))
    
    P_new <- asplit(array(rep(NA, 4 * max_age * 1), dim = c(4, max_age, 1)),
                    c(1, 2))
    
    #########################
    ##### females ###########
    #########################
    
    # genotypes and phenotypes of immature females that survived
    G_sIF <- map2(G[1, 1:(max_age - 1)], survived_IF,
                  ~ .x[which(as.logical(.y))])
    P_sIF <- map2(P[1, 1:(max_age - 1)], survived_IF,
                  ~ .x[which(as.logical(.y))])
    
    ##### check 1 #####
    if (sum(n_sIF - unlist(lapply(G_sIF, length))) > 0) {
      stop('number of genotypes more than number of IF survived')
    }
    if (sum(n_sIF - unlist(lapply(P_sIF, length))) > 0) {
      stop('number of phenotypes more than number of IF survived')
    }
    
    # genotypes of immature females that survived and matured
    G_mIF <- map2(G_sIF, matured_IF, ~ .x[which(.y)])
    P_mIF <- map2(P_sIF, matured_IF, ~ .x[which(.y)])
    
    ##### check 2 #####
    if (sum(n_mIF - unlist(lapply(G_mIF, length))) > 0) {
      stop('number of genotypes more than number of IF matured')
    }
    if (sum(n_mIF - unlist(lapply(P_mIF, length))) > 0) {
      stop('number of phenotypes more than number of IF matured')
    }
    
    # immature females that survived but didn't mature
    G_nmIF <- map2(G_sIF, not_matured_IF, ~ .x[which(.y)])
    P_nmIF <- map2(P_sIF, not_matured_IF, ~ .x[which(.y)])
    
    ##### check that the matured and didn't mature add up to survived
    
    # new immature females
    G_new[1, 2:max_age] <- G_nmIF
    P_new[1, 2:max_age] <- P_nmIF
    
    ##### check 3 #####
    if (sum(N[1, 2:max_age, y] - unlist(
      lapply(G_new[1, 2:max_age], function(x) length(unlist(x))))) > 0) {
      stop('number of new IF genotypes more than new IF abundance')
    }
    if (sum(N[1, 2:max_age, y] - unlist(
      lapply(P_new[1, 2:max_age], function(x) length(unlist(x))))) > 0) {
      stop('number of new IF phenotypes more than new IF abundance')
    }
    
    # genotypes of mature females that survived
    G_sMF <- map2(G[3, 1:(max_age - 1)], n_sMF, ~ resample(.x, size = .y))
    P_sMF <- map2(P[3, 1:(max_age - 1)], n_sMF, ~ resample(.x, size = .y))
    
    # new mature females
    G_new[3, 2:max_age] <- Map(`c`, G_mIF, G_sMF)
    P_new[3, 2:max_age] <- Map(`c`, P_mIF, P_sMF)
    
    ##### check 4 #####
    if (sum(N[3, 2:max_age, y] - unlist(lapply(G_new[3, 2:max_age], length))) > 0) {
      stop('number of MF genotypes more than MF abundance')
    }
    if (sum(N[3, 2:max_age, y] - unlist(lapply(P_new[3, 2:max_age], length))) > 0) {
      stop('number of MF phenotypes more than MF abundance')
    }
    
    #########################
    ##### males #############
    #########################
    
    # genotypes and phenotypes of immature males that survived
    G_sIM <- map2(G[2, 1:(max_age - 1)], survived_IM,
                  ~ .x[which(as.logical(.y))])
    P_sIM <- map2(P[2, 1:(max_age - 1)], survived_IM,
                  ~ .x[which(as.logical(.y))])
    
    ##### check 1 #####
    if (sum(n_sIM - unlist(lapply(G_sIM, length))) > 0) {
      stop('number of genotypes more than number of IM survived')
    }
    if (sum(n_sIM - unlist(lapply(P_sIM, length))) > 0) {
      stop('number of phenotypes more than number of IM survived')
    }
    
    # genotypes of immature males that survived and matured and didn't mature
    G_mIM <- map2(G_sIM, matured_IM, ~ .x[which(.y)])
    P_mIM <- map2(P_sIM, matured_IM, ~ .x[which(.y)])
    
    ##### check 2 #####
    if (sum(n_mIM - unlist(lapply(G_mIM, length))) > 0) {
      stop('number of genotypes more than number of IM matured')
    }
    if (sum(n_mIM - unlist(lapply(P_mIM, length))) > 0) {
      stop('number of phenotypes more than number of IM matured')
    }
    
    # immature males that survived but didn't mature
    G_nmIM <- map2(G_sIM, not_matured_IM, ~ .x[which(.y)])
    P_nmIM <- map2(P_sIM, not_matured_IM, ~ .x[which(.y)])
    
    # new immature males
    G_new[2, 2:max_age] <- G_nmIM
    P_new[2, 2:max_age] <- P_nmIM
    
    ##### check 3 #####
    if (sum(N[2, 2:max_age, y] - unlist(
      lapply(G_new[2, 2:max_age], function(x) length(unlist(x))))) > 0) {
      stop('number of new IM genotypes more than new IM abundance')
    }
    if (sum(N[2, 2:max_age, y] - unlist(
      lapply(P_new[2, 2:max_age], function(x) length(unlist(x))))) > 0) {
      stop('number of new IM phenotypes more than new IM abundance')
    }
    
    # genotypes of mature males that survived
    G_sMM <- map2(G[4, 1:(max_age - 1)], n_sMM, ~ resample(.x, size = .y))
    P_sMM <- map2(P[4, 1:(max_age - 1)], n_sMM, ~ resample(.x, size = .y))
    
    # new mature males
    G_new[4, 2:max_age] <- Map(`c`, G_mIM, G_sMM)
    P_new[4, 2:max_age] <- Map(`c`, P_mIM, P_sMM)
    
    ##### check 4 #####
    if (sum(N[4, 2:max_age, y] - unlist(
      lapply(G_new[4, 2:max_age], length))) > 0) {
      stop('number of MM genotypes more than MM abundance')
    }
    if (sum(N[4, 2:max_age, y] - unlist(
      lapply(P_new[4, 2:max_age], length))) > 0) {
      stop('number of MM phenotypes more than MM abundance')
    }
    
    #########################
    ##### stats #############
    #########################
    
    # genotype stats
    
    G_stats[, , y, 1] <- apply(G_new, c(1, 2), 
                               function(x)
                                 ifelse(is.numeric(unlist(x)), 
                                        median(unlist(x), na.rm = TRUE), 
                                        NA))  
    
    G_stats[, , y, 2] <- apply(G_new, c(1, 2), 
                               function(x)
                                 ifelse(is.numeric(unlist(x)), 
                                          quantile(unlist(x), prob = 0.25, 
                                                   na.rm = TRUE), 
                                        NA)) 
    
    G_stats[, , y, 3] <- apply(G_new, c(1, 2), 
                               function(x)
                                 ifelse(is.numeric(unlist(x)), 
                                        quantile(unlist(x), prob = 0.75, 
                                                 na.rm = TRUE), 
                                        NA))  
    
    G_stats[, , y, 4] <- apply(G_new, c(1, 2), 
                               function(x) 
                                 ifelse(is.numeric(unlist(x)), 
                                        mean(unlist(x), na.rm = TRUE), 
                                        NA))    

    G_stats[, , y, 5] <- apply(G_new, c(1, 2), 
                               function(x)
                                 ifelse(is.numeric(unlist(x)), 
                                        var(unlist(x), na.rm = TRUE), 
                                        NA))
    
    # phenotype stats
    P_stats[, , y, 1] <- apply(P_new, c(1, 2), 
                               function(x)
                                 ifelse(is.numeric(unlist(x)), 
                                        median(unlist(x), na.rm = TRUE), 
                                        NA)) 
    
    P_stats[, , y, 2] <- apply(P_new, c(1, 2), 
                               function(x)
                                 ifelse(is.numeric(unlist(x)), 
                                        quantile(unlist(x), prob = 0.25, 
                                                 na.rm = TRUE), 
                                        NA)) 
    
    P_stats[, , y, 3] <- apply(P_new, c(1, 2), 
                               function(x)
                                 ifelse(is.numeric(unlist(x)), 
                                        quantile(unlist(x), prob = 0.75, 
                                                 na.rm = TRUE), 
                                        NA))  
    P_stats[, , y, 4] <- apply(P_new, c(1, 2), 
                               function(x)
                                 ifelse(is.numeric(unlist(x)), 
                                        mean(unlist(x), na.rm = TRUE), 
                                        NA))
 
    P_stats[, , y, 5] <- apply(P_new, c(1, 2), 
                               function(x) 
                                 ifelse(is.numeric(unlist(x)), 
                                        var(unlist(x), na.rm = TRUE), 
                                        NA))    
    
  } else {
    
    G_new   <- NULL
    P_new   <- NULL
    G_stats <- NULL
    P_stats <- NULL
    
  }
  
  # what objects to return
  output <- list(N, G_new, P_new, G_stats, P_stats)
  
  # output
  return(output)
  
}
