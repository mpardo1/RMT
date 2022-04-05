####### RANDOM MATRICES FOR METAPOPULATION MODELS #######
### 
### parent script
###
### generate, plot and integrate metapopulation
### epidemiological models
### 
#########################################################
source("~/RMT/David/RMT_genrandom.R")
source("~/RMT/David/RMT_plotmobility.R")
source("~/RMT/David/d_functions_eigen_int.R")

####### GENERATE JACOBIAN ###############################

# number of patches
N <- 50

# epidemiological
#all rates must lie in (0,1) except for betas

Deltas <- rep(0.6, N) # birth rate
mub <- 0.6
sb <- 0.001
betas <- rep(mub, N) # transmission rates
# betas <- rgamma(N, shape = (mub/sb)^2, rate = mub/(sb^2))
thetas <- rep(0.3, N) # loss of immunity rates
mud <- 0.6
deaths <- rep(mud, N) # not disease-related death rates
mua <- 0.5
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

muc <- 0.1
sc <- 0.05
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

jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))

# plot the eigenvalues of the system

# plot_eigen_rmt(jacobian,
               # N,mub,mug = mud + mua + mudel,
               # muw,sw,rhow,Gammaw,
               # muc,sc,rhoc,Gammac,
               # tau = 0, alp = 0, K = 0)

####### INTEGRATE SYSTEM ################################

# initial populations
# for constant populations, set deltas = 0, Deltas = deaths
mu_sus <- 10000
s_sus <- 100
# alphag <- alphagamma(mu_sus,s_sus)
# betag <- betagamma(mu_sus,s_sus)
# sus_init <- rgamma(N,alphag,betag) 
sus_init <- rep(10000, N) # initial susceptibles

mu_inf <- 200
s_inf <- 150
# alphag <- alphagamma(mu_inf,s_inf)
# betag <- betagamma(mu_inf,s_inf)
# inf_init <- rgamma(N,alphag,betag) 
df_inf <- data.frame(inf_init=0,end_inf=0,sumrow_max=0,time_max=0)
mui <- 2000
si <- 1000
inf_init <- rgamma(N, shape = (mui/si)^2, rate = mui/(si^2)) 
list_sol <- list()
vec_init_inf <- seq(1,1000,20)
for(i in c(1:length(vec_init_inf))){
  inf_init <- rep(vec_init_inf[i], N)    # initial infecteds
  end_time <- 20
  Deltas <- rep(0.6, N) # birth rate
  print(paste0("i:", i))
  sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
              COMMUTING,MIGRATION,
              sus_init,inf_init,end_time)
  sol <- sol[,c(1,(N+2):(2*N+1))]
  end_inf <- sum(sol[nrow(sol),(2:ncol(sol))])
  rsum <- rowSums(sol[,(2:ncol(sol))])
  max_inf <- max(rsum)
  time_max_inf <- sol[which(rowSums(sol[,(2:ncol(sol))]) == max_inf),1]
  df_inf[nrow(df_inf) +1,] <- c(vec_init_inf[i], end_inf, max_inf, time_max_inf)
}

# df_inf <- df_inf[-1,]
# plot_int(N,sol, "INF") + ylim(c(0,10000))
path <- paste0("~/RMT/David/OUTPUT/ci_inf_",Sys.Date(),".csv")
write.csv(df_inf, path, row.names = TRUE)

path <- "~/RMT/David/OUTPUT/ci_inf_2022-04-05.csv"
df_inf_ci <- read.csv(file = path)
# df_inf_ci <- read.csv(file = "~/RMT/Integration/ci_inf2022-03-17.csv")
# head(df_inf_ci)

df_inf <- df_inf_ci[-1,]

size_let <- 13 

plot_inf_ci <- ggplot(df_inf) +
  geom_line(aes(inf_init, end_inf))  +
  xlab("Initial infected individuals")+
  ylab("Infected individuals at equilibrium")+
  theme_bw()  +
  theme(text = element_text(size = size_let))

Path <- "~/Documentos/PHD/2022/RMT_SIR/Plots/Gen/"
path <- paste0(Path,"Plot_inf_ci_b0,1_g0,95_muc_0,01_sc0,002_muw0,6_sw0,05.png")
ggsave(path,
       plot = plot_inf_ci, device = "png")

ggplot(df_inf) +
  geom_line(aes(inf_init, max_inf))  +
  xlab("Number of initial infected individuals")+
  ylab("Number of infected individuals at equilibrium")+
  theme_bw()

ggplot(df_inf) +
  geom_line(aes(inf_init, time_max))  +
  xlab("Initial number of infected individuals")+
  ylab("Time")+
  theme_bw()


####### Rand or cte IC ########
mui <- 4000
si <- 2000
inf_init <- rgamma(N, shape = (mui/si)^2, rate = mui/(si^2)) 
# inf_init <- rep(mui, N)

sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)

plot_int(N, sol, state = "INF") +
  theme_bw() 

library("latex2exp")
plot_cte_4 <- plot_int(N, sol, state = "INF") +
  theme_bw() + theme(legend.position="none") +
  ggtitle(TeX("$\\mu_{inf} = 4000, \\sigma_{inf} = 0$"))

plot_rand_1 <- plot_rand_1 +
  ggtitle(TeX("$\\mu_{inf} = 1000, \\sigma_{inf} = 2000$"))

library("ggpubr")
gg1 <- ggarrange(plot_rand_1  + ylim(c(0,9000)) + xlim(c(0,5)),
                 plot_rand_4 + ylim(c(0,9000)) + xlim(c(0,5)) + ylab(""),
                 ncol = 2, nrow = 1)

gg2 <- ggarrange(plot_inf_ci,
                     plot_cte_4 + ylim(c(0,9000)) + xlim(c(0,5)),
                      nrow = 1, ncol = 2)
gg_full <- ggarrange(plot_inf_ci,
                     plot_cte_4 + 
                       ylim(c(0,9000)) + xlim(c(0,5)) + xlab("") +
                       theme(text = element_text(size = size_let)),
                     plot_rand_1  +
                       ylim(c(0,9000)) + xlim(c(0,5)) +  
                       theme(text = element_text(size = size_let)),
                     plot_rand_4 + 
                       ylim(c(0,9000)) + xlim(c(0,5)) + ylab("") +
                       theme(text = element_text(size = size_let)),
                     ncol = 2, nrow = 2,
                     labels = c("a", "b", "c", "d"),
                     align = "h")

Path <- "~/Documentos/PHD/2022/RMT_SIR/Plots/Gen/"
path <- paste0(Path,"Rand_inf_ci_b0,1_g0,95_muc_0,01_sc0,002_muw0,6_sw0,05.png")
ggsave(path,
       plot = gg_full, device = "png")
