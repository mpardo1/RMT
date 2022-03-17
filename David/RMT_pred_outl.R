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
N <- 100

# epidemiological
#all rates must lie in (0,1) except for betas

Deltas <- rep(0.6, N) # birth rate
mub <- 0.1
sb <- 0.001
betas <- rep(mub, N) # transmission rates
# betas <- rgamma(N, shape = (mub/sb)^2, rate = mub/(sb^2))
thetas <- rep(0.3, N) # loss of immunity rates
mud <- 0.9
deaths <- rep(mud, N) # not disease-related death rates
mua <- 0.5
alphas <- rep(mua, N) # recovery rates
mudel <- 0.15
deltas <- rep(mudel, N) # disease-related death rates
mug = deaths + alphas + deltas

# mobility
#commuting and migration networks

muw <- 0.1 
sw <- 0.05
rhow <- 0 #original rho (Gamma of baron et al)
Gammaw <- 0 #gamma of baron et al
rw <- 0
cw <- 0

muc <- 0.001
sc <- 0.0005
rhoc <- 0
Gammac <- 0
rc <- 0
cc <- 0

df_err_outl <- data.frame(N = 0, bet = 0, muw = 0, outl = 0, pred = 0)
# jacobian
d <- 10000
alphag <- alphagamma(0.5,0.2)
betag <- betagamma(0.5,0.2)
muw_vec <- rgamma(d,alphag,betag)
muw_vec[which(muw_vec >= 1)] <- trunc(muw_vec[which(muw_vec >= 1)],2,3 )
muw_vec[which(muw_vec == 1)] <- muw_vec[which(muw_vec == 1)] <- 0.2
alphag <- alphagamma(1,2)
betag <- betagamma(1,2)
bet_vec <- rgamma(d,alphag,betag)
for(j in c(1:d)){
  muw = muw_vec[j]
  for(i in c(30:100,10)){
    # print(paste0("i: ",i))
    print(paste0("j: ",j))
    betas <- rep(bet_vec[j],i)
    deaths <- rep(mud,i) # not disease-related death rates
    alphas <- rep(mua,i) # recovery rates
    deltas <- rep(mudel,i) # disease-related death rates
    mug = deaths + alphas + deltas
    mug <- mug[1]
    mub <- betas[1]
    # Generate random mat:
    COMMUTING <- rand_mat(i, muw, sw, distrib = "beta")
    diag(COMMUTING) <- 0
    MIGRATION <- rand_mat(i, muc, sc, distrib = "beta")
    diag(MIGRATION) <- 0
    
    # Create jacobian
    jacobian <- (COMMUTING + diag(i)) %*% diag(betas) + MIGRATION -
      diag(deaths + alphas + deltas + colSums(MIGRATION))
    outlier <- mub - mug + mub*muw*(i-1)
                      
    eigen <- eigen_mat(jacobian)
    df_err_outl[nrow(df_err_outl)+1,] <- c(i, betas[1], muw,
                                           max(eigen$re),  outlier)
  }
}

# plot_eigen_rmt(jacobian,
               # N,mub,mug = mud + mua + mudel,
               # muw,sw,rhow,Gammaw,
               # muc,sc,rhoc,Gammac,
               # tau = 0, 0, 0)

# df_err_outl$err <- ((df_err_outl$pred - df_err_outl$outl)/df_err_outl$outl)^2
# df_err_outl <- df_err_outl[-1,]
# df_err_outl_group <- df_err_outl %>% group_by(N) %>% 
#   summarise(mean = mean(err))

path <- paste0("~/RMT/David/outl_pred",Sys.Date(),".csv")
write.csv(df_err_outl,path, row.names = TRUE)
# 
# ggplot(df_err_outl_group) +
#   geom_point(aes(N,mean))
