# function to get beta axis labels

OSRs_to_betas <- function(OSRs) {
  
  # load libraries
  library(ggplot2)
  library(dplyr)
  
  # values of x to plot
  # OSRs <- seq(from = 0.01, to = 0.5, by = 0.01)
  # OSRs <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)
  
  # beta values to cycle through
  betas <- seq(from = 1, to = 100, by = 0.01)
  # betas <- c(1, 1.35, 1.94, 3.1, 6.57, 8.31, 11.19, 16.94, 34.14)
  
  # initialise DF
  DF <- data.frame(Operational_Sex_Ratio = rep(OSRs, each = length(betas)), 
                   Beta = rep(betas, times = length(OSRs)), 
                   Reproductive_Success = NA)
  
  # CDF function
  for (a in 1:length(OSRs)) {
    
    start <- (a - 1)*length(betas) + 1
    stop <- start + length(betas) - 1
    DF$Reproductive_Success[start:stop] <- pbeta(2 * OSRs[a], 
                                                 shape1 = 1, 
                                                 shape2 = betas)
    
  }
  
  # make Beta factor
  #DF$Beta <- as.factor(DF$Beta)
  
  acceptable_loss_of_males <- DF %>% 
    group_by(Operational_Sex_Ratio) %>%
    filter(Reproductive_Success <= 0.99) %>%
    filter(Beta == max(as.numeric(Beta)))
  
  output <- as.character(acceptable_loss_of_males$Beta)
  
  return(output)
  
}
