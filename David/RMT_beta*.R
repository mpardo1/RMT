####### RANDOM MATRICES FOR METAPOPULATION MODELS #######
### 
### article plots
###
### plots for the article
### works over "RMT_parent.R"
###
rm(list = ls())
source("~/RMT/David/RMT_genrandom.R")
source("~/RMT/David/RMT_plotmobility.R")
source("~/RMT/David/d_functions_eigen_int.R")
library("viridis")

####### COMMUTING: PANEL ################################
library("copula")
library("copula")
library("viridis")
library("ggsci")
library("grid")
library("reshape")
library("ggpubr")
library("cowplot")
library("ggforce")
library("latex2exp")
#----------------------------------------------------------------------------#
####### Parameter values #####

mub*muw*(N-1)+mub-mua-mud-mudel

col_unstab <- "#0D0C18"
colA <- "#F69A79"
colB <- "#F26430"
colC <- "#D2430F"
colD <- "#85DCFF"
colE <- "#0ABAFF"
colF <- "#0084B8"

# N <- 100
N <- 100
Deltas <- rep(0.4, N) # birth rate
mub <- 0.05
# sb <- 0.001
betas <- rep(mub, N) # transmission rates
thetas <- rep(0.1, N) # loss of immunity rates
mud <- 0.4
deaths <- rep(mud, N) # not disease-related death rates
mua <- 0.4
alphas <- rep(mua, N) # recovery rates
mudel <- 0
deltas <- rep(mudel, N) # disease-related death rates
gammas = deaths + alphas + deltas
mug <- gammas[1]
alp_bet <- 0.5
k <- 1
  
muw <- 0.05
sw <- 0.05
rhow <- 0             # original rho (Gamma of baron et al)

muc <- 0.001
sc <- 0.0005
rhoc <- .001


COMMUTING <- rand_mat(N, muw, sw, distrib = "beta")
MIGRATION <- rand_mat(N, muc, sc, distrib = "beta")
# MIGRATION <- matrix(0,N,N)
diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)

# para los plots
sus_init <- rep(100000, N) # initial susceptibles
inf_init <- runif(N, min = 50,100)  # initial infecteds
end_time <- 200
end_time_rand <- 20

################################################################
##### Change K:
stragBET <- function(k,N){
  a <- mub*muw +muc
  sj2 <- (N/2)*a +  (alp_bet/2)*(1+ (k-1)*muw) + mub*(1-muw) - mug - N*muc +
    (1/2)*sqrt(N^2*a^2 + alp_bet^2*(1+(k-1)*muw)^2 + 
                 2*alp_bet*a*(N*(1+(k-1)*muw) + 2*(N-k)*(muw-1)))
  sj1 <- mub*muw*(N-1)+mub-mua-mud-mudel
  if(sj2>sj1){
    sj <- sj2
  }else{
    sj <- sj1
  }
  return(sj)
} 

# COMMUTING <- rand_mat_cor_beta(N, muw*N, sw*N, rhow, Gammaw, rw, cw)
# MIGRATION <- rand_mat_cor_beta(N, muc*N, sc*N, rhoc, Gammac, rc, cc)
# diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)
betas[1:k] <- betas[1:k] + alp_bet
jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))
deaths[1] + alphas[1] + deltas[1]
plot_eigen(jacobian) + 
  geom_point(aes(stragBET(k,N),0), color = "red")

count = 0
max_out <- c()
while(count < 100){
  COMMUTING <- rand_mat(N, muw, sw, distrib = "beta")
  MIGRATION <- rand_mat(N, muc, sc, distrib = "beta")
  # MIGRATION <- matrix(0,N,N)
  diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)
  count = count + 1
  jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
    diag(deaths + alphas + deltas + colSums(MIGRATION))
  max_out[count] <- max(eigen_mat(jacobian)$re)
}

pred_old <- (N/2)*(mub*muw + muc) + (alp_bet/2)*(1 + (k-1)*muw) + mub*(1-muw) -
  mug - N*muc + (1/2)*sqrt(N^2*(mub*muw + muc)^2 + alp_bet^2*(1 + (k-1)*muw)^2 + 
                             2*alp_bet*(mub*muw + muc)*(N*(1+(k-1)*muw) + 2*(N-k)*(muw-1))) 

df <- data.frame(it=seq(1,100,1),out = max_out, 
                 out_pred = stragBET(k,N) , out_old = pred_old)
df_plot <- reshape2::melt(df, id.vars = c("it"))
ggplot(df_plot) + 
  geom_line(aes(it, value, color = variable))

####################################################################
######### CHANGING the PERTURBATION #############
# Eigenvalues:
stragBET <- function(muc){
  a <- mub*muw +muc
  sj2 <- (N/2)*a +  (alp_bet/2)*(1+ (k-1)*muw) + mub*(1-muw) - mug - N*muc +
    (1/2)*sqrt(N^2*a^2 + alp_bet^2*(1+(k-1)*muw)^2 + 
                 2*alp_bet*a*(N*(1+(k-1)*muw) + 2*(N-k)*(muw-1)))
  sj1 <- mub*muw*(N-1)+mub-mua-mud-mudel
  if(sj2>sj1){
    sj <- sj2
  }else{
    sj <- sj1
  }
  return(sj)
} 

vec <- seq(0,1,0.01)
outBETmuc <- sapply(vec, stragBET)

df_out <- data.frame(muc = vec, 
                     outBETmuc = outBETmuc)

mig_plot <- ggplot(df_out) + 
  geom_line(aes(muc,outBETmuc)) +
  ylab("s(J)") + xlab(TeX("$\\mu_m$")) + 
  theme_bw()

muc <- 0.001
# Eigenvalues:
stragBET <- function(muw){
  a <- mub*muw +muc
  sj2 <- (N/2)*a +  (alp_bet/2)*(1+ (k-1)*muw) + mub*(1-muw) - mug - N*muc +
    (1/2)*sqrt(N^2*a^2 + alp_bet^2*(1+(k-1)*muw)^2 + 
                 2*alp_bet*a*(N*(1+(k-1)*muw) + 2*(N-k)*(muw-1)))
  sj1 <- mub*muw*(N-1)+mub-mua-mud-mudel
  if(sj2>sj1){
    sj <- sj2
  }else{
    sj <- sj1
  }
  return(sj)
} 

vec <- seq(0,1,0.01)
outBETmuw <- sapply(vec, stragBET)

df_out <- data.frame(muw = vec, 
                     outBET = outBETmuw)

com_plot <- ggplot(df_out) + 
  geom_line(aes(muw,outBET)) +
  ylab("s(J)") + xlab(TeX("$\\mu_c$")) + 
  theme_bw()

jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))
plot_eigen(jacobian)

muw <- 0.05

######### CHANGING the PERTURBATION #############
# Eigenvalues:
stragBET <- function(mub){
  a <- mub*muw +muc
  sj2 <- (N/2)*a +  (alp_bet/2)*(1+ (k-1)*muw) + mub*(1-muw) - mug - N*muc +
    (1/2)*sqrt(N^2*a^2 + alp_bet^2*(1+(k-1)*muw)^2 + 
                 2*alp_bet*a*(N*(1+(k-1)*muw) + 2*(N-k)*(muw-1)))
  sj1 <- mub*muw*(N-1)+mub-mua-mud-mudel
  if(sj2>sj1){
    sj <- sj2
  }else{
    sj <- sj1
  }
  return(sj)
} 

vec <- seq(0,1,0.01)
outBETmub <- sapply(vec, stragBET)

df_out <- data.frame(mub = vec, 
                     outBET = outBETmub)

bet_plot <- ggplot(df_out) + 
  geom_line(aes(mub,outBET)) +
  ylab("s(J)") + xlab(TeX("$\\beta$")) + 
  theme_bw()

mub <- 0.05
##### Change K:
stragBET <- function(k){
  a <- mub*muw +muc
  sj2 <- (N/2)*a +  (alp_bet/2)*(1+ (k-1)*muw) + mub*(1-muw) - mug - N*muc +
    (1/2)*sqrt(N^2*a^2 + alp_bet^2*(1+(k-1)*muw)^2 + 
                 2*alp_bet*a*(N*(1+(k-1)*muw) + 2*(N-k)*(muw-1)))
  sj1 <- mub*muw*(N-1)+mub-mua-mud-mudel
  if(sj2>sj1){
    sj <- sj2
  }else{
    sj <- sj1
  }
  return(sj)
} 

vec <- seq(0,N/2,1)
out_BET <- sapply(vec, stragBET)

df_out <- data.frame(k = vec, 
                     outBET = out_BET)

k_plot <- ggplot(df_out) + 
  geom_line(aes(k,outBET), color = "#4C6085", size = 0.8) +
  ylab("s(J)") + xlab("k") + 
  theme_bw()

colA <- "#F69A79"
colB <- "#F26430"
colC <- "#D2430F"
colD <- "#85DCFF"
colE <- "#0ABAFF"
colF <- "#0084B8"

df_combmuc <- data.frame(x = seq(0,1,0.01), out_BET = outBETmuc, type = "muc")
df_combmub <- data.frame(x = seq(0,1,0.01), out_BET = outBETmub, type = "mub")
df_combmuw <- data.frame(x = seq(0,1,0.01), out_BET = outBETmuw, type = "muw")
df_comb <- rbind(df_combmuc, df_combmub, df_combmuw)
plot_3 <- ggplot(df_comb) + 
  geom_line(aes(x, out_BET, color = type), size = 0.8) +
  scale_color_manual(values = c(colA, colC, colD),
                     name = NULL,
                     labels = c(TeX("$\\beta$"), TeX("$\\mu_{m}$"),TeX("$\\mu_{c}$") )) + 
  theme_bw() + 
  theme(text = element_text(size = 15),legend.position = c(0.15, 0.8),
        legend.text.align = 0) + xlab("") +
  
  ylab("s(J)")


### panel
col <- "#5C5D8D"
plotf <- plot_grid(plot_3 + 
                     geom_hline(yintercept = 0, color = "red", linetype = "dashed") + 
                     theme(text = element_text(size = 16)),
          k_plot + 
            geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
            ylab("") + 
            scale_colour_manual(values = c(col)) + 
            theme(text = element_text(size = 16),legend.position = c(0.15, 0.8),
                  legend.text.align = 0))

plotf
Path <- "~/Documentos/PHD/2022/RMT_SIR/Plots/SM/"
muw
mub
N
mug
muc
ggsave(file=paste0(Path,"betaSJ50N0_05mum0_1muc0_1mub0_95mug.pdf"))


#############################Figura 8 SM############################################
mub <- 0.1
muw <- 0.05
COMMUTING <- rand_mat(N, muw, sw, distrib = "beta")
MIGRATION <- rand_mat(N, muc, sc, distrib = "beta")
# MIGRATION <- matrix(0,N,N)
diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)
stragBET <- function(k){
  a <- mub*muw + muc
  (N/2)*a + (alp_bet/2)*(1+ (k-1)*muw) + mub*(1-muw) - mug - N*muc + 
    (1/2)*sqrt((N^2)*(a^2) + (alp_bet^2)*(1+(k-1)*muw)^2 + 
                 2*alp_bet*a*(N*(1+(k-1)*muw) + 2*(N-k)*(muw-1)))
} 

outl_mean <- mub*muw*(N-1) + mub - mug 
outl_mean == stragBET(0)
k_max = 15
mat <- matrix(1,k_max+1,3)
for(i in c(0:k_max)){
  bet <- rep(mub,N)
  if(i>0){
    bet[1:i] <- bet[1:i] + alp_bet
  }
  mubk <- mean(bet)
  print(mubk)
  outl_mean <- mubk *(muw*(N-1)+1) - mug
  out_k <- stragBET(i)
  jacobian <- (COMMUTING + diag(N)) %*% diag(bet) + MIGRATION -
    diag(deaths + alphas + deltas + colSums(MIGRATION))
  reout <- max(eigen_mat(jacobian)$re)
  mat[i+1,] <- c(outl_mean, out_k, reout)
}

df_out <- data.frame( k = seq(0,k_max,1), outm = mat[,1], outk = mat[,2], outr = mat[,3] )
df_plot <- reshape2::melt(df_out, id.vars = "k")
df_plot$var2 <- "1"
df_plot[which(df_plot$variable == "outm")]

plot_pred_vs_real <- ggplot(df_plot) + 
  geom_line(aes(k, value, color = variable, linetype = ), size = 1) +
  scale_color_manual(values = c(colA, colC, colD),
                     name = NULL,
                     labels = c("Prediction CL", "Prediction LRP ", "Rightmost eigenvalue") ) +
  theme_bw() + 
  theme(text = element_text(size = 16),legend.position = c(0.3, 0.82),
        legend.text.align = 0) + 
  geom_hline(yintercept = 0, color = "red", linetype = "dashed")

plot_pred_vs_real

######## 1 VS K patches #######
alp <- alp_bet
alp_cte <- alp_bet
outl1 <- (N/2)*(mub*muw + muc) + (alp/2)*(1 + (k-1)*muw) + mub*(1-muw) -
  mug - N*muc + (1/2)*sqrt(N^2*(mub*muw + muc)^2 + alp^2*(1 + (k-1)*muw)^2 + 
                             2*alp*(mub*muw + muc)*(N*(1+(k-1)*muw) + 2*(N-k)*(muw-1)))

outlRMT <- mub - mug + mub*muw*(N-1)

alp_vec <- seq(0.5,2,0.1)
df_sol_k_1 <- data.frame(k = 0, alph = 0, 
                         max_eig_1 = 0, max_eig_k = 0 ,
                         pred_eig_1 = 0, pred_eig_k = 0, outl_RMT = outlRMT)

for(i in c(1:(N/2))){
  betas <- rep(mub, N) 
  alp <- i*alp_cte
  betas[1] <- betas[1] + alp
  jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
    diag(gammas + colSums(MIGRATION))
  plot_eigen(jacobian) +
    geom_vline(xintercept = 0, color = "blue", linetype = "longdash") +
    geom_point(aes(outl1,0), color = "purple")
  k <- 1
  outl1 <- (N/2)*(mub*muw + muc) + (alp/2)*(1 + (k-1)*muw) + mub*(1-muw) -
    mug - N*muc + (1/2)*sqrt(N^2*(mub*muw + muc)^2 + alp^2*(1 + (k-1)*muw)^2 + 
                               2*alp*(mub*muw + muc)*(N*(1+(k-1)*muw) + 2*(N-k)*(muw-1)))
  
  eig_jac <- eigen_mat(jacobian)
  max_eig_1 <- max(eig_jac$re)
  
  betas <- rep(mub, N)
  betas[1:i] <- betas[1:i] + alp_cte
  jacobian_k <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
    diag(gammas + colSums(MIGRATION))
  
  eig_jac_k <- eigen_mat(jacobian_k)
  max_eig_k <- max(eig_jac_k$re)
  
  k <- i
  alp <- alp_cte
  
  outlk <- (N/2)*(mub*muw + muc) + (alp/2)*(1 + (k-1)*muw) + mub*(1-muw) -
    mug - N*muc + (1/2)*sqrt(N^2*(mub*muw + muc)^2 + alp^2*(1 + (k-1)*muw)^2 + 
                               2*alp*(mub*muw + muc)*(N*(1+(k-1)*muw) + 2*(N-k)*(muw-1)))
  
  plot_eigen(jacobian) +
    geom_vline(xintercept = 0, color = "blue", linetype = "longdash") +
    geom_point(aes(outlk,0), color = "purple")
  
  df_sol_k_1[nrow(df_sol_k_1)+1,] <- c(i, alp,
                                       max_eig_1, max_eig_k,
                                       outl1, outlk, outlRMT)
  
  print(paste0("i:",i))
}

df_sol_k_1 <- df_sol_k_1[-1,]
df_plot <- df_sol_k_1[,c(1,3,4,5,6)]

colA <- "#F69A79"
colB <- "#F26430"
colC <- "#D2430F"
colD <- "#85DCFF"
colE <- "#0ABAFF"
colF <- "#0084B8"
# colnames(df_plot) <- c("k","real 1 patch","real k patches",
# "pred 1 patch","pred k patches")
# df_filt <- df_sol_k_1[,c(1,5,6)]
df_plot <- reshape2::melt(df_plot, id.vars = "k")
df_plot$type <- substr(df_plot$variable,1,3)
df_plot$numk <- substr(df_plot$variable,
                       nchar(as.character(df_plot$variable)),
                       nchar(as.character(df_plot$variable)))

df_plot$type[which(df_plot$type == "max")] <- "Empiric"
df_plot$type[which(df_plot$type == "pre")] <- "Prediction"
df_plot$kalp <- alp_cte*df_plot$k
gg_1_vs_k <- ggplot(df_plot) + 
  geom_line(aes(kalp,value, colour = numk, linetype = type), size = 1.4, alpha = 0.6) + 
  ylab("value") + 
  xlab("") +
  labs(color='') +
  scale_color_manual(values = c(colA, colD),
                     name = NULL,
                     labels = c("1", "k") ) +
  theme_bw() + 
  theme(text = element_text(size = 16),
        legend.position = c(0.2, 0.72),
        legend.text.align = 0)

gg_1_vs_k

plot_grid(plot_pred_vs_real  + ylab("s(J)") ,
          gg_1_vs_k  + ylab("") + xlab(""), nrow = 1)
