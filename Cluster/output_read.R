rm(list = ls())
library("tidyverse")
library("deSolve")
library("ggplot2")

sol_k_patch <- read.csv(file = "~/Documents/1_vs_k_2022-03-01.csv")
head(sol_k_patch)

sol_k_patch <- sol_k_patch[-1,]
sol_k_patch$diff_outl <- sol_k_patch$outk - sol_k_patch$out1
sol_k_patch$ls_outl <- (sol_k_patch$outk - sol_k_patch$out1)^2/sol_k_patch$out1
err_outl <- sum(sol_k_patch$diff_outl)

sol_k_patch$diff_inf <- sol_k_patch$max.inf.k - sol_k_patch$max.inf.1
sol_k_patch$ls_inf <- (sol_k_patch$max.inf.k - sol_k_patch$max.inf.1)^2/sol_k_patch$max.inf.1
err_inf <- sum(sol_k_patch$diff_inf)
