# table and/or heatmap of ideal incubation temperatures

# load libraries
library(ggplot2)
library(tidyverse)
library(readr)
library(ggpattern)
library(patchwork)

##### calculate ideal temperatures (adjusted) ##################################

# source functions
source('~/Projects/iliketurtles3/code/mating function/OSRs_to_betas.R')

# operational proportion male
OSRs <- c(0.49, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05)
# xFs <- round(1/OSRs - 1, 2)
betas <- as.numeric(OSRs_to_betas(OSRs))

# hatchling sex ratios * emergence success
temps <- seq(from = 25, to = 35, by = 0.01)
emergence <- 0.86 / (1 + exp(1.7 * (temps - 32.7)))

# narrow TRT proportions male (patricio et al. 2017)
t_piv1 <- 29.4
k1 <- -1.54
narrow_TRT_pM <- round(1/(1 + exp(-k1*(temps - t_piv1))), 5)

# wide TRT proportions male (godfrey & mrosovsky 2006)
t_piv2 <- 29.4
k2 <- -0.77
wide_TRT_pM <- round(1/(1 + exp(-k2*(temps - t_piv2))), 5)

# TRTs
TRTs <- c('Narrow transitional range', 'Wide transitional range')
M_remigration <- 1.47
F_remigration <- 3.87

hatchlings <- data.frame(TRT = rep(TRTs, each = length(betas) * length(temps)), 
                         OSR = rep(OSRs, 
                                   each = length(temps), times = length(TRTs)),
                         Beta = rep(rep(betas, each = length(temps)), 
                                    times = length(TRTs)), 
                         Temp = rep(temps, times = length(TRTs) * length(betas)), 
                         Emergence = rep(emergence, 
                                         times = length(TRTs) * length(betas)), 
                         PSR = c(rep(narrow_TRT_pM, times = length(betas)), 
                                 rep(wide_TRT_pM, times = length(betas)))) %>%
  # mutate(xF = round(1/OSR - 1, 2)) %>%
  mutate(nM = 100 * PSR / M_remigration) %>%
  mutate(nF = 100 * (1 - PSR) / F_remigration) %>%
  mutate(breeding_success = pbeta(q = 2*(nM / (nM + nF)), 
                                    shape1 = 1, 
                                    shape2 = Beta)) %>%
  mutate(nEggs = nF * breeding_success * 4.95 * 100.58) %>%
  mutate(nHatchlings = nEggs * emergence)

ideals <- hatchlings %>%
  group_by(TRT, Beta) %>%
  filter(nHatchlings == max(nHatchlings))
  # mutate(xF = round(1/OSR - 1, 2)) 

# make values factors
ideals$OSR <- factor(ideals$OSR, levels = OSRs)
ideals$TRT <- factor(ideals$TRT, levels = TRTs)
# ideals$xF <- factor(ideals$xF, levels = unique(ideals$xF))

# save as table
save(ideals, 
     file = '~/Projects/iliketurtles3/output/ideals.Rdata')

ideals_without_emergence <- hatchlings %>%
  group_by(TRT, Beta) %>%
  filter(nEggs == max(nEggs))

# make values factors
ideals_without_emergence$OSR <- factor(ideals_without_emergence$OSR, levels = OSRs)
ideals_without_emergence$TRT <- factor(ideals_without_emergence$TRT, levels = TRTs)

# save as table
save(ideals_without_emergence, 
     file = '~/Projects/iliketurtles3/output/ideals_without_emergence.Rdata')

##### make the figure ##########################################################

# load ideals object
# ideals <- read.csv('output/ideals.csv')
load("~/Projects/iliketurtles3/output/ideals.Rdata")

# # xaxis labels
# xlabs <- paste(rep(unique(ideals_adjusted$xF), times = 2), 
#                '\n (', 
#                unique(ideals_adjusted$iPSR), 
#                ')', 
#                sep = '')

# apply to factor
# ideals_adjusted$xlabs <- xlabs

# adjust dataframe to get other useful columns
to_plot <- ideals %>%
  # mutate(Above_init_temp = as.character(Temp > 31.8)) %>%
  # mutate(temps_below = replace(Temp, Temp <= 31.8, '')) %>%
  # mutate(temps_above = replace(Temp, Temp > 31.8, '')) %>%
  mutate(labs = paste(round(Temp, 1), '\n (', round(PSR, 2), ')', sep = ''))

# actually do the heatmap thing
ideal_temps_heatmap <- ggplot(data = to_plot, 
                              aes(x = OSR, 
                                  y = TRT, 
                                  fill = Temp)) +
  geom_tile(color = 'white') +
  labs(fill = "Incubation \n temperature \n (\u00B0C)") +
  xlab("Minimum OSR required for 99% female reproductive success \n (associated hatchling sex ratio)") +
  scale_y_discrete(labels = c("Narrow \n transitional \n range", 
                              "Wide \n transitional \n range"
  )) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(0.5, 0.25, 1, 1), units = 'cm')) +
  theme(axis.title.x = element_text(size = 11, vjust = -3)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text = element_text(size = 9)) +
  theme(strip.text = element_text(size = 9)) +
  theme(legend.title = element_text(size = 9)) +
  geom_text(aes(label = labs), 
            size = 3, 
            fontface = 'bold')
  # geom_text(aes(label = round(as.numeric(temps_above), 2)), 
  #           size = 3)

# guides(pattern = 'none')

ideal_temps_heatmap

# ggsave(ideal_temps_heatmap,
#        file = '~/Projects/iliketurtles3/figures/ideal_temps_heatmap.png',
#        width = 8, height = 3)

##### ideal temps (not adjusted by emergence success) ##########################

# load in object
load("~/Projects/iliketurtles3/output/ideals_without_emergence.Rdata")

# adjust dataframe to get other useful columns
to_plot2 <- ideals_without_emergence %>%
  mutate(labs = paste(round(Temp, 1), '\n (', round(PSR, 2), ')', sep = ''))

# actually do the heatmap thing
ideal_temps_without_emergence_heatmap <- ggplot(data = to_plot2, 
                              aes(x = OSR, 
                                  y = TRT, 
                                  fill = Temp)) +
  geom_tile(color = 'white') +
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', 
                       midpoint = 29.4) +
  labs(fill = "Incubation \n temperature \n (\u00B0C)") +
  xlab("Minimum OSR required for 99% female breeding success \n (associated hatchling sex ratio)") +
  scale_y_discrete(labels = c("Narrow \n transitional \n range", 
                              "Wide \n transitional \n range"
  )) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(0.5, 0.25, 1, 1), units = 'cm')) +
  theme(axis.title.x = element_text(size = 11, vjust = -3)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text = element_text(size = 9)) +
  theme(strip.text = element_text(size = 9)) +
  theme(legend.title = element_text(size = 9)) +
  geom_text(aes(label = labs), 
            size = 3, 
            fontface = 'bold')

# guides(pattern = 'none')

ideal_temps_without_emergence_heatmap

# ggsave(ideal_temps_without_emergence_heatmap,
#        file = '~/Projects/iliketurtles3/figures/ideal_temps_without_emergence_heatmap.png',
#        width = 8, height = 3)

##### joint figure #############################################################

A <- ideal_temps_heatmap +
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', 
                       midpoint = 29.4, limits = c(29, 35)) + 
  guides(fill = 'none')  + 
  labs(tag = '(A)') +
  theme(plot.tag.position = c(0, 1), 
        plot.tag = element_text(hjust = 0, vjust = 1, size = 12, face = 'bold'))  +
  xlab('')

A

B <- ideal_temps_without_emergence_heatmap +
  theme() +
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', 
                       midpoint = 29.4, limits = c(29, 35)) +
  labs(tag = '(B)') +
  theme(axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.tag.position = c(0, 1), 
        plot.tag = element_text(hjust = 0, vjust = 1, size = 12, face = 'bold'))

joint_fig <- A / plot_spacer() / B +
  plot_layout(ncol = 1, heights = c(1, -0.25, 1), guides = 'collect')

joint_fig

# ggsave(joint_fig,
#        file = '~/Projects/iliketurtles3/figures/ideal_temps_joined_heatmap.png',
#        width = 7.5, height = 5)

##### joint fig + difference between the two ###################################

joint_to_plot <- to_plot %>%
  ungroup() %>%
  mutate(WE_Temp = to_plot2$Temp) %>%
  mutate(WE_PSR = to_plot2$PSR) %>%
  mutate(diff_Temp = round(Temp - WE_Temp, 1)) %>%
  mutate(diff_PSR = round(PSR - WE_PSR, 2)) %>%
  mutate(diff_labs = paste(diff_Temp, '\n (', diff_PSR, ')', sep = ''))


joint_fig_diff <- ggplot(data = joint_to_plot, 
         aes(x = OSR, 
             y = TRT, 
             fill = diff_Temp)) +
  geom_tile(color = 'white') +
  scale_fill_gradient2(high = 'red', mid = 'white', low = 'dodgerblue3', 
                       midpoint = 0) +
  labs(fill = "\n Difference in \n incubation \n temperature \n (\u00B0C)") +
  xlab("Minimum OSR required for 99% female breeding success \n (associated hatchling sex ratio)") +
  scale_y_discrete(labels = c("Narrow \n transitional \n range", 
                              "Wide \n transitional \n range"
  )) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(0.5, 0.25, 1, 1), units = 'cm')) +
  theme(axis.title.x = element_text(size = 11, vjust = -3)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text = element_text(size = 9)) +
  theme(strip.text = element_text(size = 9)) +
  theme(legend.title = element_text(size = 9)) +
  geom_text(aes(label = diff_labs), size = 3, fontface = 'bold')

joint_fig_diff

B2 <- B +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())

A2 <- A  +
  theme(axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

joint_fig2 <- A2 / plot_spacer() / B2 +
  plot_layout(ncol = 1, heights = c(1, -0.5, 1), guides = 'collect')

final_fig <- joint_fig2 / plot_spacer() / joint_fig_diff +
  plot_layout(ncol = 1, heights = c(1, -0.45, 1, -0.45, 1)) +
  labs(tag = '(C)')  +
  theme(plot.tag.position = c(0, 1), 
        plot.tag = element_text(hjust = 0, vjust = 1, size = 12, face = 'bold'))

final_fig

ggsave(final_fig,
       file = '~/Projects/iliketurtles3/figures/ideal_temps_all_heatmap.png',
       width = 8, height = 6)
