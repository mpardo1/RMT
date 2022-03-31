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
N <- 100

# epidemiological
#all rates must lie in (0,1) except for betas

Deltas <- rep(0.3, N) # birth rate
mub <- 0.02
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
muw <- 0.2 
sw <- 0.05
rhow <- 0 #original rho (Gamma of baron et al)
Gammaw <- 0 #gamma of baron et al
rw <- 0
cw <- 0

muc <- 0.01
sc <- 0.001
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

# jacobian
betas <- rep(mub, N)
jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))

# plot the eigenvalues of the system
library("ggforce")
plot_stab <- plot_eigen_rmt(jacobian,
               N,mub,mug = mud + mua + mudel,
               muw,sw,rhow,Gammaw,
               muc,sc,rhoc,Gammac,
               tau = 0, alp = 0, K = 1) +
                scale_y_continuous( breaks=c(0)) 
# + xlim(c(-60,-50))
print(plot_eigen(jacobian))
eigen <-  eigen_mat(jacobian)
####### INTEGRATE SYSTEM ################################

# initial populations
# for constant populations, set deltas = 0, Deltas = deaths

sus_init <- rep(10000, N) # initial susceptibles
inf_init <- rep(100, N)    # initial infecteds

end_time <- 50
sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)

plot_stab <- plot_int1(N, sol, state = "INF") +
  theme_bw() +theme(legend.position="none") 
plot_stab

vec_col <-  vector(mode="character", length=N)
vec_col[1:N] <- "#A63446"

plot.inf.stab <- plot_int(N, sol, state = "INF") +
  scale_colour_manual(values = vec_col) +
  theme_bw() + theme(legend.position="none") 

plot.inf.stab
#### 1 PATCH ####
alp_bet <- 1.2

betas <- rep(mub, N) 
betas[1] <- alp_bet + betas[1]

jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))


plot_eigen_rmt(jacobian,
               N,mub,mug = mud + mua + mudel,
               muw,sw,rhow,Gammaw,
               muc,sc,rhoc,Gammac,
               tau = 0, alp = alp_bet, K = 1) +
  scale_y_continuous( breaks=c(0)) 


sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)

# plot SUS, INF, REC or TOT population
plot_unstab1 <- plot_int1(N, sol, state = "INF") +
  theme_bw() +theme(legend.position="none") 

plot_unstab1

vec_col <-  vector(mode="character", length=N)
vec_col[1:N] <- "#A63446"
vec_col[1] <- "#3066BE"

plot.inf <- plot_int1(N, sol, state = "INF") +
  scale_colour_manual(values = vec_col) +
  theme_bw() + theme(legend.position="none") 

plot.inf
muw
plot.inf.bet.cte <- plot.inf + ggtitle(""*mu[c]~": 0.1 , "*mu[w]~": 0.1 ") 

####### PERTURBATIONS ###################################
## 1 PATCH DIFFERENT MOBILITY ######
alp <- 7.4
K = 1
ind <- sample(1:N,K)
betas <- rep(mub, N)
betas[ind] <- alp + betas[ind]
vec_col <-  vector(mode="character", length=N)
vec_col[1:N] <- "royalblue3"
vec_col[ind] <- "red4"
# jacobian

jacobian1 <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))

# plot the eigenvalues of the system

plot_eigen_rmt(jacobian1,
               N,mub,mug = mud + mua + mudel,
               muw,sw,rhow,Gammaw,
               muc,sc,rhoc,Gammac,
               tau = 0, alp, K)

# integro el sistema con condiciones iniciales 
end_time <- 50
sol1 <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)

plot_int(N, sol1, state = "INF")
# plot SUS, INF, REC or TOT population
plot.inf1 <- plot_int(N, sol1, state = "INF") +
  scale_colour_manual(values = vec_col) +
  theme_bw() + theme(legend.position="none")

plot.inf.1 <- plot.inf1 + xlim(c(0,5)) + ylim(c(0,200000))

muc
muw
plot.inf.hm.lc.lim <- plot.inf.1 +
  ggtitle("a  "*mu[c]~": 0.1, "*mu[w]~": 0.01 ") 

##### DF SOL ####
sol_df1 <-  as.data.frame(sol1)
for(i in c(1:N)){
  colnames(sol_df)[i+1] <-  paste0("S",i)
  colnames(sol_df)[N+i+1] <-  paste0("I",i)
  colnames(sol_df)[2*N+i+1] <-  paste0("R",i)
}

####### SAVE FILE #####
gammas.cte <- format(round(gammas,2), decimal.mark = ',')
beta.cte <- format(round(mub,2), decimal.mark = ',')
# mu_w_w <- format(round(mu_w,2), decimal.mark = ',')
sw.cte <- format(round(sw,2), decimal.mark = ',')
# mu_c_w <- format(round(mu_c,2), decimal.mark = ',')
sc.cte <- format(round(sc,2), decimal.mark = ',')

Path <- "~/Documentos/PHD/2022/RMT_SIR/Plots/"
Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/Epi_param/1patch/"

path <- paste0(Path,"mob_epi","N",N,"b",beta.cte,"g",gammas.cte,
               "sw",sw.cte,"sc",sc.cte,".png")
ggfull <-   ggarrange(plot.inf.bet.cte, plot.inf.lm.hc.lim,
                      plot.inf.hm.lc.lim, plot.inf.hm.hc.lim)
ggsave(path,
       plot =ggfull, device = "png")
