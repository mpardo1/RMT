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
N <- 50

# epidemiological
#all rates must lie in (0,1) except for betas 

Deltas <- rep(0.6, N) # birth rate
mub <- 0.2
sb <- 0.001
betas <- rep(mub, N) # transmission rates
# betas <- rgamma(N, shape = (mub/sb)^2, rate = mub/(sb^2))
thetas <- rep(0.3, N) # loss of immunity rates
mud <- 0.2
deaths <- rep(mud, N) # not disease-related death rates
mua <- 0.2
alphas <- rep(mua, N) # recovery rates
mudel <- 0
deltas <- rep(mudel, N) # disease-related death rates
gammas = deaths + alphas + deltas

# mobility
#commuting and migration networks
muw <- 0.9
sw <- 0.05
rhow <- 0 #original rho (Gamma of baron et al)
Gammaw <- 0 #gamma of baron et al
rw <- 0
cw <- 0

muc <- 0.008
sc <- 0.005
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

plot_jac_inf <- plot_eigen_rmt(jacobian,
               N,mub,mug = mud + mua + mudel,
               muw,sw,rhow,Gammaw,
               muc,sc,rhoc,Gammac,
               tau = 0, alp = 0, K = 0) 
# + xlim(c(-60,-50))
eigen <-  eigen_mat(jacobian)
print(plot_eigen(jacobian))
####### INTEGRATE SYSTEM ################################

# initial populations
# for constant populations, set deltas = 0, Deltas = deaths

sus_init <- rep(200, N) # initial susceptibles
inf_init <- rep(50, N)    # initial infecteds

end_time <- 100

# integro el sistema con condiciones iniciales 
sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)

# plot SUS, INF, REC or TOT population
plot_int(N, sol, state = "INF")
plot_int(N, sol, state = "INF") +  xlim(c(0,0.5))
plot_inf <- plot_int(N, sol, state = "TOT") + theme_bw() 
# +
#   xlim(c(0,50)) + ylim(c(0,60))
# +ggtitle(""*mu[c]~": 0.5") 
muc
sc
plot_hm_ls <- plot_int(N, sol, state = "INF")+ theme_bw()  +
  ggtitle(TeX("$\\mu_{c}: 0.5, \\s_{c}: 0.0001$")) 
plot_eigen_hm_ls <- plot_jac_inf
#### FULL JACOBIAN MATRIX ####
full_jac <- full_mat(N,Deltas,betas,deaths,deltas,
                     thetas,alphas,COMMUTING, MIGRATION) 

eigen_full <-  eigen_mat(full_jac)
plot_full <- plot_eigen(full_jac)

ggarrange(plot_jac_inf,plot_full, plot_inf)

#### JACOBIAN MATRIX 1 PATCH####
mat_SIR_1 <- mat_SIR_1p(birth,betas,deaths,deltas,
                       thetas,alphas)
eigen_1 <-  eigen_mat(mat_SIR_1)
plot_1 <- plot_eigen(mat_SIR_1)
sol <- int1(1, Deltas[1],betas[1],deaths[1],thetas[1],alphas[1],deltas[1],
           COMMUTING[1,1],MIGRATION[1,1],
           sus_init[1],inf_init[1],end_time)
sol <-  as.data.frame(sol)
colnames(sol) <- c("time","S","I","R")
df_plot <- reshape2::melt(as.data.frame(sol), id.vars = c("time"))
plot  <- ggplot(df_plot,aes(time, value)) + 
  geom_line(aes(colour = variable))  +
  ylab("Number of individuals") 


#### SAVE FILE ##
alpha_ct_w <- format(round(alphas,2), decimal.mark = ',')
deaths_ct_w <- format(round(deaths,2), decimal.mark = ',')
delta_ct_w <- format(round(deltas,2), decimal.mark = ',')
birth_ct_w <- format(round(Deltas,2), decimal.mark = ',')
beta_ct_w <- format(round(betas,2), decimal.mark = ',')
theta_ct_w <- format(round(thetas,2), decimal.mark = ',')
mu_w_w <- format(round(muw,2), decimal.mark = ',')
s_w_w <- format(round(sw,2), decimal.mark = ',')
s_c_w <- format(round(sc,2), decimal.mark = ',')

Path <- "~/Documentos/PHD/2022/RMT_SIR/Plots/"
Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/Gen/"

path <- paste0(Path,"muc","N",N,"b",beta_ct_w,"a",alpha_ct_w,"d",deaths_ct_w,
               "del",delta_ct_w,"birt",birth_ct_w,"t",theta_ct_w,"mw",
               mu_w_w,"sw","0.01","sm",s_c_w,".png")
ggsave(path,
       plot =ggarrange(gg_inf, gg_eigen, ncol = 1, nrow=2), device = "png")
