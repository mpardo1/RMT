####### RANDOM MATRICES FOR METAPOPULATION MODELS #######
### 
### parent script
###
### generate, plot and integrate metapopulation
### epidemiological models
### 
# rm(list = ls())
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
muw <- 0.1 
sw <- 0.05
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

mub <- 0.5
df_epi <- data.frame(mean = mub, mean_arit = 0, sigma = 0,
                     max_rand = 0, max_mean = 0,
                     max_time_rand = 0, max_time_mean = 0)
sig_vec <- seq(0,1,0.01)
len <- length(sig_vec)
count = 0
while( count < 1000){
  for(i in c(1:len)){
    ##### Random beta ####
    sb <- sig_vec[i]
    betas <- rgamma(N, shape = (mub/sb)^2, rate = mub/(sb^2)) 
    
    sus_init <- rep(10000, N) # initial susceptibles
    inf_init <- rep(100, N)    # initial infecteds
    
    end_time <- 50
    sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
                    COMMUTING,MIGRATION,
                    sus_init,inf_init,end_time)
    sol_inf <- sol[,c(1,(N+2):(2*N+1))]
    sum_inf <- rowSums(sol_inf[,2:(N+1)])
    max_rand <- max(sum_inf)
    ind <- which(sum_inf == max_rand)
    time_max_rand <- sol_inf[ind,1]
    
    # Mean(rand(betas)):
    betas <- rep(mean(betas),N)
    sus_init <- rep(10000, N) # initial susceptibles
    inf_init <- rep(100, N)    # initial infecteds
    
    end_time <- 50
    sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
                    COMMUTING,MIGRATION,
                    sus_init,inf_init,end_time)
    sol_inf <- sol[,c(1,(N+2):(2*N+1))]
    sum_inf <- rowSums(sol_inf[,2:(N+1)])
    max_mean <- max(sum_inf)
    ind <- which(sum_inf == max_mean)
    time_max_mean <- sol_inf[ind,1]
    
    df_epi[nrow(df_epi)+1,] <- c(mub, mean(betas), sb,
                                 max_rand, max_mean,
                                 time_max_rand,time_max_mean)
  }
  count = count + 1
}

#### Write in a CSV ######
Path <- "~/RMT/David/OUTPUT/"
path <- paste0(Path,"rand_bet_int_g0,5_muc_0,001_sc0,0001_muw0,2_sw0,05_",
               Sys.Date(), ".csv")
write.csv(df_epi, path,row.names = TRUE)

#### Read from a CSV ######
Path <- "~/RMT/David/OUTPUT/"
path <- paste0(Path,"rand_bet_g0,5_muc_0,001_sc0,0001_muw0,2_sw0,05_2022-04-01.csv")
df_rand <- read.csv(file = path)
df_rand <- df_rand[-1,]

df_rand$diff_max <- df_rand$max_rand - df_rand$max_mean
df_rand_group <- df_rand %>%  group_by(sigma) %>% 
  summarise(mean = mean(sq_err), n = n())

err_max_inf <- ggplot(df_rand_group) + 
  geom_point(aes(sigma,mean), size = 1) +
  xlab(TeX("$\\sigma_{\\beta}$")) +
  ylab("Mean squared error") +
  theme_bw()
# df_rand$diff_time <- df_rand$time_max_rand - df_rand$time_max_mean
