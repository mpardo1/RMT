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
Deltas <- rep(0.1, N) # birth rate
mub <- 0.35
# sb <- 0.001
betas <- rep(mub, N) # transmission rates
#betas <- rgamma(N, shape = (mub/sb)^2, rate = mub/(sb^2))
thetas <- rep(0.1, N) # loss of immunity rates
mud <- 0.2
deaths <- rep(mud, N) # not disease-related death rates
mua <- 0.4
alphas <- rep(mua, N) # recovery rates
mudel <- 0
deltas <- rep(mudel, N) # disease-related death rates
gammas = deaths + alphas + deltas
mug <- gammas[1]
alp_bet <- 1
k <- 5
  
muw <- 0.25
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
  (N/2)*a +  (alp_bet/2)*(1+ (k-1)*muw) + mub*(1-muw) - mug - N*muc +
    (1/2)*sqrt(N^2*a^2 + alp_bet^2*(1+(k-1)*muw)^2 + 2*alp_bet*a*(N*(1+(k-1)*muw) + 2*(N-k)*(muw-1)))
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
  (N/2)*a +  (alp_bet/2)*(1+ (k-1)*muw) + mub*(1-muw) - mug - N*muc +
    (1/2)*sqrt(N^2*a^2 + alp_bet^2*(1+(k-1)*muw)^2 + 2*alp_bet*a*(N*(1+(k-1)*muw) + 2*(N-k)*(muw-1)))
} 

vec <- seq(0,1,0.01)
outBETmuc <- sapply(vec, stragBET)

df_out <- data.frame(muc = vec, 
                     outBETmuc = outBETmuc)

mig_plot <- ggplot(df_out) + 
  geom_line(aes(muc,outBETmuc)) +
  ylab("s(J)") + xlab(TeX("$\\mu_m$")) + 
  theme_bw()


# Eigenvalues:
stragBET <- function(muw){
  a <- mub*muw +muc
  (N/2)*a +  (alp_bet/2)*(1+ (k-1)*muw) + mub*(1-muw) - mug - N*muc +
    (1/2)*sqrt(N^2*a^2 + alp_bet^2*(1+(k-1)*muw)^2 + 2*alp_bet*a*(N*(1+(k-1)*muw) + 2*(N-k)*(muw-1)))
} 

vec <- seq(0,1,0.01)
outBETmuw <- sapply(vec, stragBET)

df_out <- data.frame(muw = vec, 
                     outBET = outBETmuw)

com_plot <- ggplot(df_out) + 
  geom_line(aes(muw,outBET)) +
  ylab("s(J)") + xlab(TeX("$\\mu_c$")) + 
  theme_bw()

######### CHANGING the PERTURBATION #############
# Eigenvalues:
stragBET <- function(mub){
  a <- mub*muw +muc
  (N/2)*a +  (alp_bet/2)*(1+ (k-1)*muw) + mub*(1-muw) - mug - N*muc +
    (1/2)*sqrt(N^2*a^2 + alp_bet^2*(1+(k-1)*muw)^2 + 2*alp_bet*a*(N*(1+(k-1)*muw) + 2*(N-k)*(muw-1)))
} 

vec <- seq(0,1,0.01)
outBETmub <- sapply(vec, stragBET)

df_out <- data.frame(mub = vec, 
                     outBET = outBETmub)

bet_plot <- ggplot(df_out) + 
  geom_line(aes(mub,outBET)) +
  ylab("s(J)") + xlab(TeX("$\\beta$")) + 
  theme_bw()

##### Change K:
stragBET <- function(k){
  a <- mub*muw +muc
  (N/2)*a +  (alp_bet/2)*(1+ (k-1)*muw) + mub*(1-muw) - mug - N*muc +
    (1/2)*sqrt(N^2*a^2 + alp_bet^2*(1+(k-1)*muw)^2 + 2*alp_bet*a*(N*(1+(k-1)*muw) + 2*(N-k)*(muw-1)))
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
plotf <- plot_grid(plot_3 ,
          k_plot + ylab("")+ 
            scale_colour_manual(values = c(col)) + 
            theme(text = element_text(size = 15),legend.position = c(0.15, 0.8),
                  legend.text.align = 0)  )
Path <- "~/Documentos/PHD/2022/RMT_SIR/Plots/SM/"
muw
mub
N
mug
muc
ggsave(file=paste0(Path,"betaSJ50N0_05mum0_1muc0_1mub0_95mug.pdf"))


#########################################################################
outl_mean <- mub*(muw*(N-1)+1) - mug
k_max = 15
mat <- matrix(1,k_max,3)
for(i in c(1:k_max)){
  bet <- rep(mub,N)
  bet[1:i] <- bet[1:i] + alp_bet
  mubk <- mean(bet)
  outl_mean <- mubk *(muw*(N-1)+1) - mug
  out_k <- stragBET(i)
  jacobian <- (COMMUTING + diag(N)) %*% diag(bet) + MIGRATION -
    diag(deaths + alphas + deltas + colSums(MIGRATION))
  reout <- max(eigen_mat(jacobian)$re)
  mat[i,] <- c(outl_mean, out_k, reout)
}

df_out <- data.frame( k = seq(1,k_max,1), outm = mat[,1], outk = mat[,2], outr = mat[,3] )
df_plot <- reshape2::melt(df_out, id.vars = "k")
df_plot$var2 <- "1"
df_plot[which(df_plot$variable == "outm")]

plot_pred_vs_real <- ggplot(df_plot) + 
  geom_line(aes(k, value, color = variable, linetype = ), size = 1) +
  scale_color_manual(values = c(colA, colC, colD),
                     name = NULL,
                     labels = c("Prediction CL", "Prediction LRP ", "Right most eigenvalue") ) +
  theme_bw() + 
  theme(text = element_text(size = 15),legend.position = c(0.3, 0.8),
        legend.text.align = 0)

