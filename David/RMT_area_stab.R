####### RANDOM MATRICES FOR METAPOPULATION MODELS #######
### 
### parent script
###
### generate, plot and integrate metapopulation
### epidemiological models
### 
rm(list = ls())
source("~/RMT/David/RMT_genrandom.R")
source("~/RMT/David/RMT_plotmobility.R")
source("~/RMT/David/d_functions_eigen_int.R")
####### GENERATE JACOBIAN ###############################

# number of patches
N <- 10

# epidemiological
#all rates must lie in (0,1) except for betas

Deltas <- rep(0.6, N) # birth rate
mub <- 0.2
sb <- 0.001
betas <- rep(mub, N) # transmission rates
# betas <- rgamma(N, shape = (mub/sb)^2, rate = mub/(sb^2))
thetas <- rep(0.3, N) # loss of immunity rates
mud <- 0.6
deaths <- rep(mud, N) # not disease-related death rates
mua <- 0.2
alphas <- rep(mua, N) # recovery rates
mudel <- 0
deltas <- rep(mudel, N) # disease-related death rates
gammas = deaths + alphas + deltas

# mobility
#commuting and migration networks
muw <- 0.1
sw <- 0.05
rhow <- 0 #original rho (Gamma of baron et al)
Gammaw <- 0 #gamma of baron et al
rw <- 0
cw <- 0

muc <- 0.01
sc <- 0.00001
rhoc <- 0
Gammac <- 0
rc <- 0
cc <- 0

step <- 0.01
beta_vec <- seq(0.01,0.9,step)
muw_vec <- seq(0.01,0.9,step)
df_sol <- data.frame(beta = 0, gamma = 0, N = 0, muw = 0, state = FALSE)
N = 100
for(i in c(1:length(beta_vec))){
  print(paste0(" i : ", i))
  for(j in c(1:length(beta_vec))){
    if(is.na(muw_vec[j]) | muw_vec[j] > 1 | muw_vec[j] <0  ){
      print(paste0("Problem muw: ", muw_vec[j]))
    }else{
      COMMUTING <- rand_mat(N, muw_vec[j], sw, distrib = "beta")
      diag(COMMUTING) <- 0
      MIGRATION <- rand_mat(N, muc, sc, distrib = "beta")
      diag(MIGRATION) <- 0
      
      # EPI param:
      betas <- rep(beta_vec[i], N)
      deltas <- rep(mudel, N)
      deaths <- rep(mud, N) 
      alphas <- rep(mua, N)
      
      jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
        diag(deaths + alphas + deltas + colSums(MIGRATION))
      eigen <- eigen_mat(jacobian)
      max_eig <- max(eigen$re)
      state <- ifelse(max_eig >0, FALSE, TRUE)
      df_sol[nrow(df_sol) + 1,1:4] <- list(betas[1], gammas[1], N, muw_vec[j])
      df_sol[nrow(df_sol) ,5] <- state
    }
  }
}

df_sol <- df_sol[-1,]

path <- paste0("~/RMT/David/OUTPUT/area_gen_",Sys.Date(),".csv")
write.csv(df_sol, path,row.names = TRUE)

# path <- "~/RMT/David/OUTPUT/area_gen_2022-03-22.csv"
# df_sol <- read.csv(file = path)
# df_sol$Stability <- ifelse(df_sol$state == TRUE, "Stable", "Unstable")
# 
# vec <- seq(0,1,0.01)
# df_sol <- df_sol[which(df_sol$beta == vec),]
# library(latex2exp)
# ggplot(df_sol) +
#   geom_point(aes(beta,muw, colour = Stability)) + theme_bw()  +
#   scale_color_manual(values=c("#909CC2", "#084887")) +
#   ylab(TeX("$\\mu_c$")) +
#   xlab(TeX("$\\beta$")) +
#   # ggtitle(""*gamma/beta~": 4")
#   ggtitle(paste0("N: ",N)) +
#   coord_fixed() +
#   theme(text = element_text(size = 15)) +
#   geom_point(aes(0.07,0.07), colour= "#F58A07", size = 0.1) +
#   geom_point(aes(0.07,0.1), colour= "#F58A07") +
#   geom_point(aes(0.1,0.07), colour= "#F58A07") 
