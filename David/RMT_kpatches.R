rm(list=ls())
source("~/RMT/David/RMT_genrandom.R")
source("~/RMT/David/RMT_plotmobility.R")
source("~/RMT/David/d_functions_eigen_int.R")

####### GENERATE JACOBIAN ###############################

# number of patches
N <- 100

# epidemiological
#all rates must lie in (0,1) except for betas

Deltas <- rep(0.3, N) # birth rate
mub <- 0.5
sb <- 0.001
betas <- rep(mub, N) # transmission rates
# betas <- rgamma(N, shape = (mub/sb)^2, rate = mub/(sb^2))
thetas <- rep(0.3, N) # loss of immunity rates
mud <- 0.3
deaths <- rep(mud, N) # not disease-related death rates
mua <- 0.2
alphas <- rep(mua, N) # recovery rates
mudel <- 0
deltas <- rep(mudel, N) # disease-related death rates
gammas = deaths + alphas + deltas

# mobility
#commuting and migration networks
muw <- 0.1 
sw <- 0.08
rhow <- 0 #original rho (Gamma of baron et al)
Gammaw <- 0 #gamma of baron et al
rw <- 0
cw <- 0

muc <- 0.001
sc <- 0.0001
rhoc <- 0
Gammac <- 0
rc <- 0
cc <- 0

COMMUTING <- rand_mat(N, muw, sw, distrib = "beta")
diag(COMMUTING) <- 0
# COMMUTING <- rand_mat_ell(N, muw, sw, rhow, distrib = "beta")
# COMMUTING[sample.int(N^2, round(p*N^2))] <- 0

MIGRATION <- rand_mat(N, muc, sc, distrib = "beta")
diag(MIGRATION) <- 0

step <- 0.001
k_vec <- seq(1,N,1)
df_sol <- data.frame(k = 0, alp = 0, max_eig_1 = 0, max_eig_k = 0)
N = 100
mug <- gammas[1]
alp <-  2
ind <- sample(1:N,1)
for(i in c(1:length(k_vec))){
  betas <- rep(mub, N)
  betas[ind] <- betas[ind] + alp
  jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
    diag(deaths + alphas + deltas + colSums(MIGRATION))
  eig_jac <- eigen_mat(jacobian)
  max_eig_1 <- max(eig_jac$re)
  
  betas <- rep(mub, N)
  betas[1:k_vec[i]] <- betas[1:k_vec[i]] + alp/k_vec[i]
  jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
    diag(deaths + alphas + deltas + colSums(MIGRATION))
  eig_jac <- eigen_mat(jacobian)
  max_eig_k <- max(eig_jac$re)
  
  df_sol[nrow(df_sol)+1,] <- c(k_vec[i], alp, max_eig_1, max_eig_k)
  
  print(paste0("i:",i))
}

df_sol <- df_sol[-1,]
df_sol <-  df_sol[,c(1,3,4)]
df_sol$err <- (df_sol$max_eig_1 - df_sol$max_eig_k)^2/df_sol$max_eig_1

gg_err <- ggplot(df_sol) + 
  geom_point(aes(k, err), color ="#2C2C54", size = 0.9) +
  geom_line(aes(k, err), color ="#A40E4C", size = 0.4) +
  ylab("Mean square error") + 
  theme_bw()

# Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/Gen/"
# path <- paste0(Path,"kpatch_err_g0,5_b0,5_muc_0,001_sc0,0001_muw0,1_sw0,05_",Sys.Date(), ".png")
# ggsave(path,
#        plot = gg_err, device = "png")

df_sol <-  df_sol[,c(1,2,3)]
df_plot <-  reshape2::melt(df_sol, id.vars = c("k"))

ggplot(df_plot) + 
  geom_line(aes(k, value, colour = variable))

Path <- "~/RMT/David/OUTPUT/"
path <- paste0(Path,"Areaepi_g0,5_muc_0,001_sc0,0001_muw0,1_sw0,05_",Sys.Date(), ".csv")
write.csv(df_sol, path,row.names = TRUE)

sol_k_patch <- read.csv(file = "~/Documents/1_vs_k_2022-03-08.csv")
head(sol_k_patch)

sol_k_patch <- sol_k_patch[-1,]

gammas = deaths + alphas + deltas
mub <- 0.1
mug <- 1.8
sb <- 1.5
alp_vec <- seq(0.5,2,0.1)
df_sol_k <- data.frame(k = 0, alp = 0, max_eig_1 = 0, max_eig_k = 0)
for(i in c(1:N)){
  ind <- sample(1:length(alp_vec),1)
  alp <- alp_vec[ind]
  ind <- sample(1:N,i)
  betas <- rep(mub, N) 
  betas[ind[1]] <- betas[ind[1]] + alp
  jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
    diag(gammas + colSums(MIGRATION))
  # plot_eigen(jacobian) +
  #   geom_vline(xintercept = 0, color = "blue", linetype = "longdash") 
    
  eig_jac <- eigen_mat(jacobian)
  max_eig_1 <- max(eig_jac$re)
  
  betas <- rep(mub, N)
  betas[ind] <- betas[ind] + alp/i
  jacobian_k <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
    diag(gammas + colSums(MIGRATION))
  # plot_eigen(jacobian) +
  #   geom_vline(xintercept = 0, color = "blue", linetype = "longdash") 
  #   
  eig_jac_k <- eigen_mat(jacobian_k)
  max_eig_k <- max(eig_jac_k$re)
  
  df_sol_k[nrow(df_sol_k)+1,] <- c(i, alp, max_eig_1, max_eig_k)
  
  print(paste0("i:",i))
}
df_sol_k <- df_sol_k[-1,]
df_plot <- df_sol_k[,c(1,3,4)]

colnames(df_plot) <- c("k","1 patch","k patches")
df_plot <- reshape2::melt(df_plot, id.vars = "k")
gg_1_vs_k <- ggplot(df_plot) + 
  geom_line(aes(k,value, colour = variable)) + 
  ylab("Right most eigenvalue") + 
  xlab("Number of patches") +
  labs(color='') +
  theme_bw()

Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/Gen/"
path <- paste0(Path,"1_vs_k_","N", N,
               "muw",format(muw,decimal.mark=","),
               "sw",format(sw,decimal.mark=","),
               "muc",format(muc,decimal.mark=","),
               "b",format(betas[1],decimal.mark=","),
               "d",format(deltas[1],decimal.mark=","), 
               "D",format(Deltas[1],decimal.mark=","),
               "a",format(alphas[1],decimal.mark=","),
               "t",format(thetas[1],decimal.mark=","),".pdf")
ggsave(path,
       plot = gg_1_vs_k, device = "png")

alp <- 1



outl1 <- (1/2)*(N*a + b + sqrt((N*a)^2 - (2*N-4)*a*b + (4*N-4)*a*c + b^2))
outl1 <- outl1 + (mub*(1-muw) - N*muc - mug)

gammas = deaths + alphas + deltas
mub <- 0.1
mug <- 1.8
sb <- 1.5

a <- mub*muw + muc
b <- alp
c <- alp*muw

alp_vec <- seq(0.5,2,0.1)
df_sol_k_1 <- data.frame(k = 0, alp = 0, 
                         max_eig_1 = 0, max_eig_k = 0 ,
                         pred_eig_1 = 0, pred_eig_k = 0)
for(i in c(1:N)){
  betas <- rep(mub, N) 
  betas[1] <- betas[1] + alp
  jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
    diag(gammas + colSums(MIGRATION))
  # plot_eigen(jacobian) +
  #   geom_vline(xintercept = 0, color = "blue", linetype = "longdash") 
  
  eig_jac <- eigen_mat(jacobian)
  max_eig_1 <- max(eig_jac$re)
  
  betas <- rep(mub, N)
  betas[1:i] <- betas[1:i] + alp/i
  jacobian_k <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
    diag(gammas + colSums(MIGRATION))
  # plot_eigen(jacobian) +
  #   geom_vline(xintercept = 0, color = "blue", linetype = "longdash") 
  #   
  eig_jac_k <- eigen_mat(jacobian_k)
  max_eig_k <- max(eig_jac_k$re)
  
  K <- i
  outlk <- (1/2)*(N*a + b + (K-1)*c + sqrt((N*a)^2 - (2*N-4*K)*a*b +
                                             (2*(K+1)*N-4*K)*a*c + b^2 + (2*K-2)*b*c + (K-1)^2*c^2))
  outlk <- outlk + (mub*(1-muw) - N*muc - mug)
  
  df_sol_k_1[nrow(df_sol_k_1)+1,] <- c(i, alp,
                                   max_eig_1, max_eig_k,
                                   outl1, outlk)
  
  print(paste0("i:",i))
}

df_sol_k_1 <- df_sol_k_1[-1,]
df_plot <- df_sol_k_1[,c(1,3,4,5,6)]

colnames(df_plot) <- c("k","real 1 patch","real k patches",
                       "pred 1 patch","pred k patches")
df_plot <- reshape2::melt(df_plot, id.vars = "k")
gg_1_vs_k <- ggplot(df_plot) + 
  geom_line(aes(k,value, colour = variable)) + 
  ylab("Right most eigenvalue") + 
  xlab("Number of patches") +
  labs(color='') +
  theme_bw()


