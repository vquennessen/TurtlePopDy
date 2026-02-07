# run_base_model <- function(scenario) {

# run_base_model <- function(scenarios, num_sims, betas) {
run_base_model <- function(arguments) {
  
  ###### model inputs ##########################################################
  
  # function arguments
  noise     <- arguments$Var1
  trt       <- as.character(arguments$Var2)
  scenario  <- arguments$Var3
  beta      <- arguments$Var4
  yrs       <- arguments$Var5
  nsims     <- arguments$Var6
  evolve    <- arguments$Var7
  trait     <- as.character(arguments$Var8) 
  rate      <- as.character(arguments$Var9)
  conserve  <- arguments$Var10
  intensity <- arguments$Var11
  frequency <- arguments$Var12
  folder    <- arguments$Var13
  
  # load in init age distribution (can't pass it as an argument anymore???)
  load('../output/init_age_distribution.Rdata')
  
  # write to progress text file
  TIME2 <- lubridate::now()
  time2 <- format(TIME2)
  
  if (evolve == TRUE) {
    
    EV <- paste('evolution - ', trait, ' - ', rate, ' - ', sep = '')
    
  } else { EV <- '' }
  
  if (conserve == TRUE) {
    
    CO <- paste('conservation - i', intensity, ' - F', frequency, ' - ', 
                sep = '')
    
  } else { CO <- '' }
  
  update1 <- paste(time2, ' - ', EV, CO, trt, ' - ', scenario, 'C - beta ', 
                   beta, ' - ', nsims, ' sims - ', yrs, ' years', sep = '')
  
  write(update1, file = '../output/progress.txt', append = TRUE)
  
  # model parameters to modulate
  temp_mu <- 31.8                         # base incubation temp mean
  climate_stochasticity <- TRUE           # whether or not to add in
  season_temp_sd <- 0.364                 # variance in temp at season level
  clutch_temp_sd <- 0.790                 # variance in temp at clutch level
  AC <- 0.5                               # autocorrelation coefficient
  
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
  M <- round(pnorm(q = 1:max_age,           # maturity ogive
                   mean = age_maturity_mu,
                   sd = age_maturity_sd),
             3)
  
  F_remigration_int <- 3.87                 # remigration interval - females
  M_remigration_int <- 1.82                 # remigration interval - males
  clutches_mu <- 4.95                       # mean # of clutches/F/season
  clutches_sd <- 2.09                       # sd of # of clutches/F/season
  eggs_mu <- 100.58                         # mean number of eggs/clutch - 100.58
  eggs_sd <- 22.61                          # sd of number of eggs/clutch - 22.61
  emergence_success_A <- 0.86               # logistic by temp - A
  emergence_success_k <- -1.7               # logistic by temp - beta
  emergence_success_t0 <- 32.7              # logistic by temp - t0
  T_threshold <- 35                         # lethal temperature threshold
  k_piv <- ifelse(trt == 'narrow',
                  -1.54,
                  -0.77)  
  T_piv <- 29.4                             # thermal reaction norm midpoint
  
  # stable age distributions to start with
  if (evolve == TRUE) {
    
    IAD <- init_age_distribution %>%
      filter(Beta == beta) %>%
      filter(TRT == trt) %>%
      filter(Trait == trait) %>%
      filter(Rate == rate)
    
  } else {
    
    IAD <- init_age_distribution %>%
      filter(Beta == beta) %>%
      filter(TRT == trt) 
    
  }
  
  IF_init <- IAD %>%
    ungroup() %>%
    filter(Sex == 'IF') %>%
    select(Abundance, G_mean, G_var, P_mean, P_var)
  
  IM_init <- IAD %>%
    ungroup() %>%
    filter(Sex == 'IM') %>%
    select(Abundance, G_mean, G_var, P_mean, P_var)
  
  MF_init <- IAD %>%
    ungroup() %>%
    filter(Sex == 'MF') %>%
    select(Abundance, G_mean, G_var, P_mean, P_var)
  
  MM_init <- IAD %>%
    ungroup() %>%
    filter(Sex == 'MM') %>%
    select(Abundance, G_mean, G_var, P_mean, P_var)
  
  if (evolve == TRUE) {

    # probabilities of females mating with 1-10 males
    male_probs <- c(0.188, 0.280, 0.236, 0.150, 0.080, 0.038, 0.017, 0.007, 0.003, 0.001)

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
    
    # genotype and phenotype summary stats, dimensions sex * age * year * # stats
    # stats: median, Q25, Q75, mean, var
    sims_G_stats <- array(rep(NA, times = 4 * max_age * yrs * 5 * nsims), 
                          dim = c(4, max_age, yrs, 5, nsims))
    
    sims_P_stats <- array(rep(NA, times = 4 * max_age * yrs * 5 * nsims), 
                          dim = c(4, max_age, yrs, 5, nsims))
    
  } else {
    
    value                 <- NULL
    h2                    <- NULL
    varGenetic            <- NULL
    varPhenotypic         <- NULL
    G                     <- NULL
    P                     <- NULL
    male_probs            <- NULL
    contributions         <- NULL
    
  }
  
  ##### conservation #############################################################
  
  # if conservation is TRUE    
  effect_size <- ifelse(conserve == TRUE, 1.3, NA)
  
  ##### initialize output ######################################################
  
  # initialize population size array by sex/maturity, age, years, sims
  sims_N <- array(rep(NA, times = 4 * max_age * yrs * nsims),
                  dim = c(4, max_age, yrs, nsims))
  
  sims_OSR <- array(rep(NA, times = yrs * nsims),
                    dim = c(yrs, nsims))
  
  sims_ES <- array(rep(NA, times = yrs * nsims),
                   dim = c(yrs, nsims))
  
  ##### run sims and save output #############################################
  
  # run the model for each simulation
  for (i in 1:nsims) {
    
    output <- base_model(scenario, beta, yrs, max_age,
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
                         conserve, frequency, intensity, effect_size)
    
    # save the output arrays
    sims_N[, , , i]             <- output[[1]]
    sims_OSR[, i]               <- output[[2]]
    sims_ES[, i]                <- output[[3]]
    
    if (evolve == TRUE) {
      sims_G_stats[, , , , i]       <- output[[4]]
      sims_P_stats[, , , , i]       <- output[[5]]
    }
    
    # write to progress text file
    if ((i/nsims*100) %% 10 == 0) {
      
      time2.5 <- format(lubridate::now())
      
      update2 <- paste(time2.5, ' - ', EV, CO, trt, ' - ', scenario, 
                       'C - beta ', beta, ' - ', nsims, ' sims - ', yrs, 
                       ' years - ', i/nsims*100, '% done!', sep = '')
      
      write(update2, file = '../output/progress.txt', append = TRUE)
      
    }
    
  }
  
  # get filepaths to save objects to
  filepath1 = paste('../output/', folder, '/',  trt, '/', scenario, 'C/beta', 
                    beta, '/', nsims, '_N.Rda', sep = '')
  filepath2 = paste('../output/', folder, '/',  trt, '/', scenario, 'C/beta', 
                    beta, '/', nsims, '_OSR.Rda', sep = '')
  filepath3 = paste('../output/', folder, '/',  trt, '/', scenario, 'C/beta', 
                    beta, '/', nsims, '_ES.Rda', sep = '')
  
  # save objects
  save(sims_N, file = filepath1)
  save(sims_OSR, file = filepath2)
  save(sims_ES, file = filepath3)
  
  if (evolve == TRUE) {
    
    filepath4 = paste('../output/', folder, '/',  trt, '/', scenario, 'C/beta', 
                      beta, '/', nsims, '_G_stats.Rda', sep = '')
    filepath5 = paste('../output/', folder, '/',  trt, '/', scenario, 'C/beta', 
                      beta, '/', nsims, '_P_stats.Rda', sep = '')
    
    save(sims_G_stats, file = filepath4)
    save(sims_P_stats, file = filepath5)
    
  }
  
  # update progress text file with total time it took to run the thing
  TIME3 <- lubridate::now()
  total_time <- format(round(TIME3 - TIME2, 3))
  
  update3 <- paste(EV, CO, trt, ' - ', scenario, 'C - beta ', beta, ' - ', 
                   nsims, ' sims - ', yrs, ' years - total time: ', 
                   total_time, '\n', sep = '')
  
  write(update3,
        file = '../output/progress.txt', append = TRUE)
  
}
