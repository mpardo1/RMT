rm(list = ls())
library("tidyverse")
library("deSolve")
library("ggplot2")

sol_k_patch <- read.csv(file = "~/Documents/1_vs_k_2022-03-08.csv")
head(sol_k_patch)

sol_k_patch <- sol_k_patch[-1,]
sol_k_patch$diff_outl <- sol_k_patch$outk - sol_k_patch$out1
sol_k_patch$ms_outl <- ((sol_k_patch$outk - sol_k_patch$out1)/sol_k_patch$out1)^2
ggplot(sol_k_patch) + 
  geom_point(aes(alpha,ms_outl, colour = N.patches))

ggplot(sol_filt1) +
  geom_point(aes(alpha,ms_outl, colour = N.patches))

sol_filt_beta$N.patches <- factor(sol_filt_beta$N.patches)
sol_filt_beta$alpha <- factor(sol_filt_beta$alpha)

ggplot(sol_filt_beta, aes(x = N.patches, y = alpha, fill = ms_outl)) +
  geom_tile(colour="white", size=0.25) +
  scale_fill_gradientn(colours = c("#D9ED92", "#76C893", "#34A0A4", "#1A759F", "#184E77"),
                       na.value = "gray90", limits = c(min(sol_filt_beta$ms_outl),max(sol_filt_beta$ms_outl))) +
  theme(panel.border=element_blank(), panel.background=element_blank()) +
  labs(x="", y="") +
  coord_fixed()

err_outl <- mean(sol_k_patch$ms_outl)

sol_k_patch$diff_inf <- sol_k_patch$max.inf.k - sol_k_patch$max.inf.1
sol_k_patch$ms_inf <- (sol_k_patch$max.inf.k - sol_k_patch$max.inf.1)^2/sol_k_patch$max.inf.1
err_inf <- sum(sol_k_patch$diff_inf)
