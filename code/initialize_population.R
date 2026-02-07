### test version of function where we don't record previous years, only last year

initialize_population <- function(arguments) {
  
  # pull model parameters from arguments list (varies)
  trt         <- as.character(arguments$Var1)
  beta        <- arguments$Var2
  evolve      <- arguments$Var3
  trait       <- as.character(arguments$Var4)
  rate        <- as.character(arguments$Var5)
  nsims       <- arguments$Var6
  burn_in     <- arguments$Var7
  date_to_use <- arguments$Var8

  # rm(list = ls())
  # # setwd('~/Projects/iliketurtles3/code')
  # library(purrr)
  # library(magrittr)
  # library(abind)
  # library(slider)
  # 
  # resample <- function(x, ...) x[sample.int(length(x), ...)]
  # 
  # trt        <- 'narrow'
  # beta       <- 1.17
  # evolve     <- FALSE
  # trait      <- 'T_piv'
  # rate       <- 'high'
  # burn_in    <- 50
  # nsims      <- 2
  # date_to_use <- gsub('-', '_', Sys.Date())
  
  # model parameters (always the same)
  sexes   <- c('IF', 'IM', 'MF', 'MM')
  S <- length(sexes)
  
    # turtle demographics
  max_age <- 85
  F_survival_values <- c(0.35, 0.8, 0.8, 0.85, 0.799)
  M_survival_values <- c(0.35, 0.8, 0.8, 0.85, 0.799)
  
  F_years_in_stage <- c(1, 2, 7, 12, 63)# lifespan
  M_years_in_stage <- c(1, 2, 7, 12, 63)# lifespan
  
  # max_age - (1 + 2 + 7 + 12)              # years for last ageclass
  IF_survival <- c(rep(F_survival_values[1], F_years_in_stage[1]), 
                   rep(F_survival_values[2], F_years_in_stage[2]), 
                   rep(F_survival_values[3], F_years_in_stage[3]), 
                   rep(F_survival_values[4], 
                       (F_years_in_stage[4] + F_years_in_stage[5])))  
  IM_survival <- c(rep(M_survival_values[1], M_years_in_stage[1]), 
                   rep(M_survival_values[2], M_years_in_stage[2]), 
                   rep(M_survival_values[3], M_years_in_stage[3]), 
                   rep(M_survival_values[4], 
                       (M_years_in_stage[4] + M_years_in_stage[5])))  
  MF_survival <- 0.799
  MM_survival <- 0.799
  
  age_maturity_mu <- 25                     # age at first reproduction, mean
  age_maturity_sd <- 2.5                    # age at first reproduction, SD
  F_remigration_int <- 3.87                 # remigration interval - females
  M_remigration_int <- 1.47                 # remigration interval - males
  clutches_mu <- 4.95                       # mean # of clutches/F/season
  eggs_mu <- 100.58                         # mean # of eggs/clutch - 100.58
  emergence_success_A <- 0.86               # logistic by temp - A
  emergence_success_k <- -1.7               # logistic by temp - beta
  emergence_success_t0 <- 32.7              # logistic by temp - t0
  T_piv <- 29.4                             # thermal reaction norm midpoint
  T_threshold <- 35                         # lethal temperature threshold

  temp_mu <- 31.8
  seasonal_temp_mus <- seq(from = 29.4, 
                           to = temp_mu,
                           length = burn_in)
  clutch_temp_sd <- 0.790                 # variance in temp at clutch level
  
  # thermal reaction norm values
  k_piv <- ifelse(trt == 'narrow', 
                  -1.54, 
                  -0.77)
  
  ##### maturity ogive
  M <- round(pnorm(q = 1:max_age, 
                   mean = age_maturity_mu, 
                   sd = age_maturity_sd), 
             3)  
  
  if (evolve == TRUE) {
    
    # probabilities of females mating with 1-10 males
    male_probs <- c(0.188, 0.280, 0.236, 0.150, 0.080, 
                    0.038, 0.017, 0.007, 0.003, 0.001)
    
    # male fertilization contributions
    raw_contributions <- 0.687 * (c(1:10))^(-1.710)
    
    contributions <- list(1)
    
    for (i in 2:length(raw_contributions)) {
      
      contributions[i] <- list(c(
        raw_contributions[1:i]/sum(raw_contributions[1:i])))
      
    }
    
    if (trait == 'T_piv') {
      
      value <- T_piv
      
      if (rate == 'effective') { 
        h2 <- 0.221
        varGenetic <- 0.926
        
      } else { 
        h2 <- 0.576
        varGenetic <- 2.41 }
      
    } 
    
    # or, if the evolvable trait is the emergence success midpoint (t0)
    if (trait == 'emergence_success_t0') {
      
      value <- emergence_success_t0
      
      if (rate == 'effective') { 
        h2 <- 0.75
        varGenetic <- 1.19
        
      } else { 
        h2 <- 0.88
        varGenetic <- 1.39 }
      
    }
    
    # phenotypic variance, error term for offspring phenotype, one for each year
    varPhenotypic <- varGenetic / h2  
    
  } else {
    
    value         <- NULL
    h2            <- NULL
    varGenetic    <- NULL
    varPhenotypic <- NULL
    G             <- NULL
    P             <- NULL
    male_probs    <- NULL
    contributions <- NULL
    
  }
  
  # immature ogive
  Mi <- 1 - M
  
  ##### first update ###########################################################
  
  # write to progress text file
  time1 <- lubridate::now()
  TIME1 <- format(time1)
  
  if (evolve == TRUE) {
    update1 <- paste(TIME1, ' - evolution - ', trait, ' - ', rate, ' - ', 
                     trt, ' - beta ', beta, ' - ', nsims, 
                     ' sims - ', burn_in, ' burn-in years', sep = '')
  } else {
    update1 <- paste(TIME1, ' - ', trt, ' - beta ', beta, ' - ', nsims, 
                     ' sims - ', burn_in, ' burn-in years', sep = '')
  }
  
  write(update1, file = '../output/SAD_progress.txt', append = TRUE)
  
  ##############################################################################
  
  # initialize N dataframe
  N <- array(rep(0, times = 4 * max_age), 
             dim = c(4, max_age))
  
  # initial pop size
  N[1, 1] <- 1000
  N[2, 1] <- 100
  N[1, 2:max_age] <- round(100 * Mi[2:max_age])
  N[2, 2:max_age] <- round(5 * Mi[2:max_age])
  N[3, 2:max_age] <- round(10 * M[2:max_age])
  N[4, 2:max_age] <- round(1 * M[2:max_age])
  
  # initialize results dataframe
  SAD <- data.frame(Nsims = NULL, 
                    Burn_in = NULL, 
                    TRT = NULL, 
                    Beta = NULL,
                    Evolution = NULL,
                    Trait = NULL,
                    Rate = NULL,
                    Year = NULL,
                    Sex = NULL,
                    Age = NULL,
                    Abundance = NULL,
                    Proportion = NULL, 
                    PSR = NULL,
                    OSR = NULL, 
                    G_mean = NULL, 
                    G_var = NULL,
                    P_mean = NULL, 
                    P_var = NULL
  )
  
  for (i in 1:nsims) {
    
    # initialize results dataframe
    sub1_SAD <- data.frame(Nsims = rep(nsims, S * max_age), 
                           Burn_in = rep(burn_in, S * max_age), 
                           TRT = rep(trt, S * max_age), 
                           Beta = rep(beta, S * max_age),
                           Evolution = rep(evolve, S * max_age),
                           Trait = rep(trait, S * max_age),
                           Rate = rep(rate, S * max_age),
                           Year = rep(1, S * max_age),
                           Sex = rep(sexes, each = max_age),
                           Age = rep(1:max_age, times = S),
                           Abundance = rep(NA, S * max_age),
                           Proportion = rep(NA, S * max_age), 
                           PSR = rep(NA, S * max_age),
                           OSR = rep(NA, S * max_age), 
                           G_mean = rep(NA, S * max_age),
                           G_var = rep(NA, S * max_age),
                           P_mean = rep(NA, S * max_age), 
                           P_var = rep(NA, S * max_age)
    )
    
    # abundance by age
    sub1_SAD$Abundance <- c(N[1, ], N[2, ], N[3, ], N[4, ])
    
    # total population size in year 1
    total <- sum(N, na.rm = TRUE)
    
    # add proportions to SAD
    sub1_SAD$Proportion <- c(N[1, ]/total, 
                             N[2, ]/total, 
                             N[3, ]/total, 
                             N[4, ]/total)
    
    sub1_SAD$PSR <- rep(sum(N[2, 1], na.rm = TRUE) / 
                          sum(N[1:2, 1], na.rm = TRUE), 
                        times = S * max_age)
    
    ##### initialize population ################################################
    
    # if evolution is true
    
    if (evolve == TRUE) {
      
      G <- lapply(t(N), rnorm, mean = value, sd = sqrt(varGenetic))
      
      sub1_SAD$G_mean <- unlist(lapply(G, mean, na.rm = TRUE))
      
      sub1_SAD$G_var <- unlist(lapply(G, var, na.rm = TRUE))
      
      P <- lapply(G, function(x) 
        rnorm(n = length(unlist(x)), mean = x, sd = sqrt(varPhenotypic)))
      
      sub1_SAD$P_mean <- unlist(lapply(P, mean, na.rm = TRUE))
      
      sub1_SAD$P_var <- unlist(lapply(P, var, na.rm = TRUE))
      
    }
    
    # move population forward in time burn_in years
    for (y in 2:burn_in) {  
      
      # initialize results dataframe
      sub2_SAD <- data.frame(Nsims = rep(nsims, S * max_age), 
                             Burn_in = rep(burn_in, S * max_age), 
                             TRT = rep(trt, S * max_age), 
                             Beta = rep(beta, S * max_age),
                             Evolution = rep(evolve, S * max_age),
                             Trait = rep(trait, S * max_age),
                             Rate = rep(rate, S * max_age),
                             Year = rep(y, S * max_age),
                             Sex = rep(sexes, each = max_age),
                             Age = rep(1:max_age, times = S),
                             Abundance = NA,
                             Proportion = NA, 
                             PSR = NA,
                             OSR = NA, 
                             G_mean = NA, 
                             G_var = NA,
                             P_mean = NA, 
                             P_var = NA
      )
      
      ##### population dynamics ##############################################
      
      # initialize N dataframe
      N_new <- array(rep(0, times = 4 * max_age), 
                     dim = c(4, max_age))
      
      # how many immature females survived in each age class
      survived_IF <- round(
        N[1, 1:(max_age - 1)] * IF_survival[1:(max_age - 1)])
      
      # how many immature females matured in each age class
      matured_IF <- round(survived_IF * M[1:(max_age - 1)])
      
      # how many didn't mature
      not_matured_IF <- survived_IF - matured_IF
      
      # new immature females: number that survived - number that matured
      N_new[1, 2:max_age] <- not_matured_IF
      
      # mature females that survived
      survived_MF <- round(N[3, 1:(max_age - 1)] * MF_survival)
      
      # new mature females: number of immature females that matured + mature 
      # females that survived
      N_new[3, 2:max_age] <- matured_IF + survived_MF
      
      # how many immature males survived in each age class
      survived_IM <- round(
        N[2, 1:(max_age - 1)] * IM_survival[1:(max_age - 1)])
      
      # how many immature males matured in each age class
      matured_IM <- round(survived_IM * M[1:(max_age - 1)])
      
      # how many didn't mature
      not_matured_IM <- survived_IM - matured_IM
      
      # new immature males: number that survived - number that matured
      N_new[2, 2:max_age] <- not_matured_IM
      
      # mature males that survived
      survived_MM <- round(N[4, 1:(max_age - 1)] * MM_survival)
      
      # new mature females: number of immature males that matured + mature 
      # males that survived
      N_new[4, 2:max_age] <- matured_IM + survived_MM
      
      # if we need to keep track of genotypes
      if (evolve == TRUE) {
        
        # initialize new G and P
        G_new <- list()
        P_new <- list()
        
        # indices, genotypes, and phenotypes of immature females that survived
        i_sIF <- map2(G[1:(max_age - 1)], survived_IF, 
                      ~ resample(1:length(.x), size = .y))
        
        G_sIF <- map2(G[1:(max_age - 1)], i_sIF, ~ .x[.y])
        P_sIF <- map2(P[1:(max_age - 1)], i_sIF, ~ .x[.y])
        
        # indices, genotypes, and phenotypes of immature females that survived 
        # and matured
        i_mIF <- map2(G_sIF, matured_IF, ~ resample(1:length(.x), size = .y))
        
        G_mIF <- map2(G_sIF, i_mIF, ~ .x[.y])
        P_mIF <- map2(P_sIF, i_mIF, ~ .x[.y])
        
        # genotypes and phenotypes of immature females that survived but 
        # didn't mature
        G_nmIF <- map2(G_sIF, i_mIF, 
                       ~ if(sum(.y > 0)) {.x[-.y]} else {.x})
        P_nmIF <- map2(P_sIF, i_mIF, 
                       ~ if(sum(.y > 0)) {.x[-.y]} else {.x})  
        
        # new immature females
        G_new[2:max_age] <- G_nmIF
        P_new[2:max_age] <- P_nmIF
        
        # genotypes of mature females that survived
        i_sMF <- map2(G[(2*max_age + 1):(3*max_age - 1)], survived_MF, 
                      ~ resample(1:length(.x), size = .y))
        
        G_sMF <- map2(G[(2*max_age + 1):(3*max_age - 1)], i_sMF, ~ .x[.y])
        P_sMF <- map2(P[(2*max_age + 1):(3*max_age - 1)], i_sMF, ~ .x[.y])
        
        # new mature females
        G_new[(2*max_age + 1):(3*max_age)] <- append(NA, Map(`c`, G_mIF, G_sMF))
        P_new[(2*max_age + 1):(3*max_age)] <- append(NA, Map(`c`, P_mIF, P_sMF))
        
        # indices, genotypes, and phenotypes of immature males that survived
        i_sIM <- map2(G[(max_age + 1):(2 * max_age - 1)], survived_IM, 
                      ~ resample(1:length(.x), size = .y))
        
        G_sIM <- map2(G[(max_age + 1):(2 * max_age - 1)], i_sIM, ~ .x[.y])
        P_sIM <- map2(P[(max_age + 1):(2 * max_age - 1)], i_sIM, ~ .x[.y])
        
        # indices, genotypes, and phenotypes of immature males that survived 
        # and matured
        i_mIM <- map2(G_sIM, matured_IM, ~ resample(1:length(.x), size = .y))
        
        G_mIM <- map2(G_sIM, i_mIM, ~ .x[.y])
        P_mIM <- map2(P_sIM, i_mIM, ~ .x[.y])
        
        # genotypes and phenotypes of immature males that survived but 
        # didn't mature
        G_nmIM <- map2(G_sIM, i_mIM, 
                       ~ if(sum(.y > 0)) {.x[-.y]} else {.x})
        P_nmIM <- map2(P_sIM, i_mIM, 
                       ~ if(sum(.y > 0)) {.x[-.y]} else {.x})  
        
        # new immature males
        G_new[(max_age + 2):(2*max_age)] <- G_nmIM
        P_new[(max_age + 2):(2*max_age)] <- P_nmIM
        
        # genotypes of mature males that survived
        i_sMM <- map2(G[(3*max_age + 1):(4*max_age - 1)], survived_MM, 
                      ~ resample(1:length(.x), size = .y))
        
        G_sMM <- map2(G[(3*max_age + 1):(4*max_age - 1)], i_sMM, ~ .x[.y])
        P_sMM <- map2(P[(3*max_age + 1):(4*max_age - 1)], i_sMM, ~ .x[.y])
        
        # new mature males
        G_new[(3*max_age + 1):(4*max_age)] <- append(NA, Map(`c`, G_mIM, G_sMM))
        P_new[(3*max_age + 1):(4*max_age)] <- append(NA, Map(`c`, P_mIM, P_sMM))
        
      }
      
      # breeding females this year
      n_available_F <- round(
        sum(survived_MF, na.rm = TRUE) / F_remigration_int)
      
      # breeding males this year
      n_available_M <- round(
        sum(survived_MM, na.rm = TRUE) / M_remigration_int)
      
      # check that there are at least 1 available male and female for breeding
      if (n_available_F < 1 | n_available_M < 1) {
        
        female_hatchlings <- 0
        male_hatchlings <- 0  
        hatchlings <- 0
        PSR <- NA
        genotype_means <- rep(NA, S * max_age)
        phenotype_means <- rep(NA, S * max_age)
        genotype_vars <- rep(NA, S * max_age)
        phenotype_vars <- rep(NA, S * max_age)
        
      } else {
        
        # operational sex ratio - proportion of males
        OSR <- n_available_M / (n_available_M + n_available_F)
        
        # calculate reproductive success
        breeding_success <- pbeta(2 * OSR, shape1 = 1, shape2 = beta) 
        
        # how many females actually find a male to mate with and then nest
        n_breeding_F <- round(n_available_F * breeding_success)
        
        # check that there is at least one female that will successfully breed
        if (n_breeding_F < 1) {
          
          female_hatchlings <- 0
          male_hatchlings <- 0  
          hatchlings <- 0
          PSR <- NA
          genotype_means <- rep(NA, S * max_age)
          phenotype_means <- rep(NA, S * max_age)
          genotype_vars <- rep(NA, S * max_age)
          phenotype_vars <- rep(NA, S * max_age)
          
        } else {
          
          # number of clutches total
          clutches <- n_breeding_F * round(clutches_mu)
          
          # vector of clutch temperatures, one number for each clutch
          clutch_temps <- rnorm(n = clutches, 
                                # mean = seasonal_temp_mus[y],
                                mean = temp_mu,
                                sd = clutch_temp_sd)
          
          # eggs vector, one number for each clutch
          eggs <- round(rep(eggs_mu, clutches))
          
          if (evolve == TRUE) {
            
            GM <- resample(unlist(G[(2*max_age + 2):(3*max_age)]), 
                           size = n_breeding_F)
            
            # extract potential paternal genotypes
            potential_GP <- resample(unlist(G[(3*max_age + 2):(4*max_age)]), 
                                     size = n_available_M)       
            
            # how many males does each female mate with
            nMales <- resample(1:length(male_probs), 
                               size = n_breeding_F, 
                               prob = male_probs, 
                               replace = TRUE)
            
            # if there are more males assigned to a female than there are available, 
            # reduce it with the maximum number of males available
            nMales[nMales > n_available_M] <- n_available_M
            
            # assign male genotypes to each female
            GP <- map(nMales, ~ resample(potential_GP, size = .x))
            
            # assign male genotypes to each egg
            GP_eggs <- lapply(GP, function(x) {
              replicate(n = round(clutches_mu), 
                        resample(x, 
                                 size = round(eggs_mu), 
                                 prob = contributions[[length(x)]], 
                                 replace = TRUE))
            })
            
            GM_eggs <- lapply(GM, function(x) {
              replicate(n = round(clutches_mu), 
                        rep(x, times = round(eggs_mu))) } )
            
            # egg genotypes
            G_eggs <- map2(GP_eggs, GM_eggs, 
                           ~ (.x + .y)/2 + rnorm(n = length(.x), 
                                                 mean = 0, 
                                                 sd = sqrt(varGenetic/2))) %>%
              abind(along = 2) %>%
              asplit(MARGIN = 2)
            
            # egg phenotypes
            P_eggs <- lapply(G_eggs, function(x) {
              rnorm(n = length(x), 
                    mean = x, 
                    sd = sqrt(varGenetic*(1 - h2)/h2)) } ) %>%
              abind(along = 2) %>%
              asplit(MARGIN = 2)
            
            # calculate how many emerged
            if (trait == 'emergence_success_t0') {
              
              # list of probability of emergence, one for each egg 
              probs_emerged <- map2(as.list(clutch_temps), P_eggs, 
                                    ~ if (.x < T_threshold) {
                                      emergence_success_A / (
                                        1 + exp(-emergence_success_k * (.x - .y)))
                                    } else { 0 } ) %>%
                lapply(pmax, 0)
              
            } else {
              
              probs_emerged <- lapply(
                clutch_temps, 
                function(x) {
                  if (x < T_threshold) {
                    emergence_success_A / (
                      1 + exp(-emergence_success_k * (x - emergence_success_t0)))
                  } else { 0 } } ) %>% 
                lapply(pmax, 0)
              
            }
            
            # which eggs emerge as hatchlings?
            indices_hatchlings <- map2(eggs, probs_emerged, 
                                       ~ as.logical(
                                         rbinom(n = .x, size = 1, prob = .y)))
            
            # how many hatchlings are there?
            hatchlings <- unlist(lapply(indices_hatchlings, sum, na.rm = TRUE))
            
            # hatchling genotypes and phenotypes
            G_hatchlings <- map2(G_eggs, indices_hatchlings, 
                                 ~ .x[as.logical(.y)])
            P_hatchlings <- map2(P_eggs, indices_hatchlings, 
                                 ~ .x[as.logical(.y)])
            
            if (trait == 'T_piv') {
              
              # probability of developing as male, one for each egg
              probs_male <- map2(as.list(clutch_temps), 
                                 P_hatchlings, 
                                 ~ 1 / (1 + exp(-k_piv * (.x - .y)))) %>%
                lapply(pmax, 0)
              
            } else {
              
              # list of probability of developing as male, one for each clutch 
              probs_male <- lapply(clutch_temps, 
                                   function(x) {
                                     1 / (1 + exp(-k_piv * (x - T_piv)))
                                   }) %>%
                lapply(pmax, 0)
              
            }
            
            # which hatchlings developed as male?
            indices_males <- map2(hatchlings, probs_male, 
                                  ~ as.logical(rbinom(n = .x, size = 1, prob = .y)))
            indices_females <- map(indices_males, ~ as.logical(Map(`-`, 1, .x)))
            
            # genotypes of females and males
            G_new[1] <- list(unlist(map2(G_hatchlings, 
                                         indices_females, 
                                         ~ .x[as.logical(.y)])))
            G_new[max_age + 1] <- list(unlist(map2(G_hatchlings, 
                                                   indices_males, 
                                                   ~ .x[as.logical(.y)])))            
            # phenotypes of females and males
            P_new[1] <- list(unlist(map2(P_hatchlings, 
                                         indices_females, 
                                         ~ .x[as.logical(.y)])))
            P_new[max_age + 1] <- list(unlist(map2(P_hatchlings, 
                                                   indices_males, 
                                                   ~ .x[as.logical(.y)])))   
            
            # number of female hatchlings
            female_hatchlings <- length(unlist(G_new[1]))
            male_hatchlings <- length(unlist(G_new[max_age + 1]))
            
            # genetics summary stats
            genotype_means <- c(unlist(lapply(G_new, mean, na.rm = TRUE)))
            genotype_vars <- c(unlist(lapply(G_new, var, na.rm = TRUE)))
            
            phenotype_means <- c(unlist(lapply(P_new, mean, na.rm = TRUE)))
            phenotype_vars <- c(unlist(lapply(P_new, var, na.rm = TRUE)))
            
            G <- G_new
            P <- P_new
            
            # genetic summary stats for dataframe
            genotype_means[1] <- mean(unlist(G_new[1]), na.rm = TRUE)
            genotype_means[max_age + 1] <- mean(unlist(G_new[max_age + 1]), 
                                                na.rm = TRUE)
            genotype_vars[1] <- var(unlist(G_new[1]), na.rm = TRUE)
            genotype_vars[max_age + 1] <- var(unlist(G_new[max_age + 1]), 
                                              na.rm = TRUE)            
            
            phenotype_means[1] <- mean(unlist(P_new[1]), na.rm = TRUE)
            phenotype_means[max_age + 1] <- mean(unlist(P_new[max_age + 1]), 
                                                 na.rm = TRUE)
            phenotype_vars[1] <- var(unlist(P_new[1]), na.rm = TRUE)
            phenotype_vars[max_age + 1] <- var(unlist(P_new[max_age + 1]), 
                                               na.rm = TRUE)  
            
            # otherwise, if there's no evolution:
          } else {
            
            # calculate expected emergence success
            emergence_success <- emergence_success_A / 
              (1 + exp(-emergence_success_k * (
                clutch_temps - emergence_success_t0)))
            
            # vector of probabilities of developing as male, one for each clutch
            probs_male <- 1 / (1 + exp(-k_piv * (clutch_temps - (T_piv))))          
            
            # hatchlings vector, one for each clutch
            hatchlings <- round(eggs * emergence_success)
            
            # number of males
            male_hatchlings <- sum(round(hatchlings * probs_male), na.rm = TRUE)
            
            # number of females
            female_hatchlings <- sum(hatchlings, na.rm = TRUE) - male_hatchlings
            
            # no genetics info without evolution
            genotype_means <- rep(NA, S * max_age)
            genotype_vars <- rep(NA, S * max_age)
            phenotype_means <- rep(NA, S * max_age)
            phenotype_vars <- rep(NA, S * max_age)
            
          }
          
        }
        
      }
      
      # add hatchlings to N
      N_new[1, 1] <- female_hatchlings
      N_new[2, 1] <- male_hatchlings
      
      # add population abundances
      sub2_SAD$Abundance <- c(N_new[1, ], N_new[2, ], N_new[3, ], N_new[4, ])
      
      # total population size 
      total <- sum(N_new, na.rm = TRUE)      
      PSR <- male_hatchlings / sum(hatchlings, na.rm = TRUE)
      
      # add proportions to SAD
      sub2_SAD$Proportion <- c(N[1, ]/total, 
                               N[2, ]/total, 
                               N[3, ]/total, 
                               N[4, ]/total)
      
      sub2_SAD$PSR <- rep(PSR, times = S * max_age) 
      sub2_SAD$OSR <- rep(OSR, times = S * max_age) 
      
      sub2_SAD$G_mean <- genotype_means
      sub2_SAD$G_var <- genotype_vars
      sub2_SAD$P_mean <- phenotype_means
      sub2_SAD$P_var <- phenotype_vars
      
      # add inner sub super data frame to outer sub super dataframe
      sub1_SAD <- rbind(sub1_SAD, sub2_SAD)
      
      # print(y)
      # 10 percent updates
      # write to progress text file, every 10% done
      if ((y/burn_in*100) %% 10 == 0) {
        TIME2 <- format(lubridate::now())
        
        if (evolve == TRUE) {
          update2 <- paste(TIME2, ' - evolution - ', trait, ' - ', rate, ' - ', 
                           trt, ' - beta ', beta, ' - ', nsims, ' sims - ', 
                           burn_in, ' burn-in years - sim ', i, ' out of ', 
                           nsims, ' - ', y/burn_in*100, '% done!', sep = '')
        } else {
          update2 <- paste(TIME2, ' - ', trt, ' - beta ', beta, ' - ', nsims, 
                           ' sims - ', burn_in, ' burn-in years - sim ', i, 
                           ' out of ', nsims, ' - ', y/burn_in*100, 
                           '% done!', sep = '')
        }
        
        write(update2, file = '../output/SAD_progress.txt', append = TRUE)
        
      }
      
      # break out of loop if there are zero males at any age
      if (sum(N_new[c(2, 4), ], na.rm = TRUE) < 1) { break }
      
      # print(paste(format(lubridate::now()), 
      #             ' - year ', y, ' out of ', burn_in, sep = ''))
      
      # new N
      
      if (sum(N_new) >= 2e6) {
        
        N <- round(N_new * 0.10) 
        
      } else { N <- N_new }
      
    }
    
    SAD <- rbind(SAD, sub1_SAD)
    
  }
  
  # add rolling 10 and 100 year median proportion values
  SADdf <- SAD %>%
    group_by(Sex, Age) %>%
    mutate(Prop_10yr_median = slide_dbl(Proportion, 
                                        median, 
                                        .before = 9, 
                                        .complete = TRUE)) %>%
    mutate(Prop_100yr_median = slide_dbl(Proportion, 
                                         median, 
                                         .before = 99, 
                                         .complete = TRUE)) 
  
  # update progress text file with total time it took to run the thing
  time3 <- lubridate::now()
  total_time <- format(round(time3 - time1, 3))
  
  if (evolve == TRUE) {
    update3 <- paste('evolution - ', trait, ' - ', rate, ' - ', 
                     trt, ' - beta ', beta, ' - ', nsims, ' sims - ', 
                     burn_in, ' burn-in years - ', nsims, 
                     ' sims - total time: ', total_time, '\n', sep = '')
  } else {
    update3 <- paste(trt, ' - beta ', beta, ' - ', nsims, ' sims - ', 
                     burn_in, ' burn-in years - ', nsims, 
                     ' sims - total time: ', total_time, '\n', sep = '')
  }
  
  write(update3, file = '../output/SAD_progress.txt', append = TRUE)
  
  if (evolve == TRUE) {
    
    f1 <- paste('_evolution', trait, '_', rate, sep = '')
    
  } else { f1 <- '_base' }
  
  # save object
  save(SADdf, 
       file = paste('../output/SAD/', date_to_use, f1, '_n', nsims, '_b', 
                    burn_in, '_', trt, '_beta', beta, '.Rdata', sep = ''))
  
}
