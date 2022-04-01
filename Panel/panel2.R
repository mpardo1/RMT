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

#--------------------------------------------------------#
##### All betas the same
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
                            tau = 0, alp = 0, K = 0) +
  scale_y_continuous( breaks=c(0)) 
# + xlim(c(-60,-50))
print(plot_eigen(jacobian))
eigen <-  eigen_mat(jacobian)

# Integrate the system:
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

plot.inf.stab <- plot_int1(N, sol, state = "INF") +
  scale_colour_manual(values = vec_col) +
  theme_bw() +
  theme(text = element_text(size = 20), legend.position="none")

plot.inf.stab

Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/panel2/"
path <- paste0(Path,"Plot_inf_g0,5_muc_0,01_sc0,001_muw0,08_sw0,05.png")
ggsave(path,
       plot = plot.inf.stab, device = "png")
#### 1 PATCH modified ####
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

#  Change color 
vec_col <-  vector(mode="character", length=N)
vec_col[1:N] <- "#A63446"
vec_col[1] <- "#3066BE"

plot.inf.1 <- plot_int1(N, sol, state = "INF") +
  scale_colour_manual(values = vec_col) +
  theme_bw() + theme(legend.position="none") +
  theme(text = element_text(size = 20))

plot.inf.1


Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/panel2/"
path <- paste0(Path,"Plot_inf1_g0,5_muc_0,01_sc0,001_muw0,08_sw0,05.png")
ggsave(path,
       plot = plot.inf.1, device = "png")

##### Sum of all the infected by alpha ####
alp_bet_vec <- seq(0,2.5,0.01)
end_time <- 50
sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)
df_sum <- data.frame(time = sol[,1])
for(i in c(1:length(alp_bet_vec))){
  print(paste0("i: ", i))
  betas <- rep(mub, N) 
  betas[1] <- betas[1] + alp_bet_vec[i]
  sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
             COMMUTING,MIGRATION,
             sus_init,inf_init,end_time)
  sol_inf <- sol[,c(1,(N+2):(2*N+1))]
  sum_inf <- rowSums(sol_inf[,2:(N+1)])
  df_sum[,ncol(df_sum)+1] <- sum_inf
}

## DF for maximum number of infected individuals and time of maximum.
max_inf <- df_sum %>% summarise_if(is.numeric, max)
time_max_vec <- c()
for(i in c(1:251)){
  ind <- which(df_sum[,i+1] == as.numeric(max_inf[i+1]) )
  time_max_vec[i] <- df_sum[ind,1]
}

df_sum_group <- data.frame(alp <- alp_bet_vec[1:251], 
                           max_inf <- t(max_inf)[1:251],
                           time_max <- time_max_vec[1:251])

df_sum_group <-  df_sum_group[-1,]
colnames(df_sum_group) <-  c("alpha", "Max_inf", "Time_max")
library("latex2exp")
plot_inf_max <- ggplot(df_sum_group) + 
  geom_line(aes(alpha, Max_inf)) + 
  xlab(TeX("$\\alpha$")) +
  ylab("Max of infected individuals") +
  theme_bw() +
  theme(text = element_text(size = 20)) 
  


plot_time_max <-  ggplot(df_sum_group) + 
  geom_line(aes(alpha, Time_max)) + 
  xlab(TeX("$\\alpha$")) +
  ylab("Time of max of infected individuals") +
  theme_bw() +
  theme(text = element_text(size = 20)) 

colnames(df_sum) <- c("time",as.character(alp_bet_vec))
df_plot <- reshape2::melt(df_sum, id.vars = c("time"))

sum_inf <- ggplot(data = df_plot, aes(x = time, y = value,
                      color = as.numeric(as.character(variable)),
                      group = variable)) +
  geom_line() +
  scale_colour_gradient(name = ""*alpha~" ", 
                        low = "blue", high = "red") +
  ylab("Sum of infected individuals") +
  theme_bw() +
  theme(text = element_text(size = 20)) 

sum_inf

Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/panel2/"
path <- paste0(Path,"Sum_inf_g0,5_muc_0,01_sc0,001_muw0,08_sw0,05.png")
ggsave(path,
       plot = sum_inf, device = "png")


##### Area stability #####
step <- 0.001
beta_vec <- seq(0,0.1,step)
alp_vec <- seq(0,2.5,step)
df_sol <- data.frame(beta = 0, gamma = 0, N = 0, alp = 0, state = FALSE)
N = 100
mug <- gammas[1]
for(i in c(1:length(beta_vec))){
  print(paste0(" i : ", i))
  for(j in c(1:length(alp_vec))){
      a <- beta_vec[i]*muw + muc
      b <- alp_vec[j]
      c <- alp_vec[j]*muw
      outl <- (1/2)*(N*a + b + sqrt((N*a)^2 - (2*N-4)*a*b + (4*N-4)*a*c + b^2))
      outl <- outl + (mub*(1-muw) - N*muc - mug)
      state <- ifelse(outl < 0, TRUE, FALSE)
      df_sol[nrow(df_sol) + 1,1:4] <- list(beta_vec[i], gammas[1], N, alp_vec[j])
      df_sol[nrow(df_sol) ,5] <- state
  }
}

df_sol <- df_sol[-1,]

library("latex2exp")
df_sol$Stability <- ifelse(df_sol$state == TRUE, "Stable", "Unstable")
plot_area <- ggplot(df_sol) +
  geom_point(aes(beta,alp, colour = Stability)) + theme_bw()  +
  scale_color_manual(values=c("#3066BE", "#A63446")) +
  ylab(TeX("$\\alpha$")) +
  xlab(TeX("$\\beta$")) +
  # ggtitle(""*gamma/beta~": 4")
  ggtitle(paste0("N: ",N)) +
  theme(text = element_text(size = 30), legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size=5)))
plot_area

Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/panel2/"
path <- paste0(Path,"Area_g0,5_muc_0,01_sc0,001_muw0,08_sw0,05.png")
ggsave(path,
       plot = plot_area, device = "png")

path <- paste0("~/RMT/David/OUTPUT/area_gen_beta_alp_",Sys.Date(),".csv")
write.csv(df_sol, path,row.names = TRUE)
##### Random beta ####
mub <- 1
sb <- 2
betas <- rgamma(N, shape = (mub/sb)^2, rate = mub/(sb^2)) 

jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))

mub <- mean(betas)
plot_eigen_rmt(jacobian,
               N,mub,mug = mud + mua + mudel,
               muw,sw,rhow,Gammaw,
               muc,sc,rhoc,Gammac,
               tau = 0, alp = 0, K = 0) +
  scale_y_continuous( breaks=c(0)) 

print(plot_eigen(jacobian))

sus_init <- rep(10000, N) # initial susceptibles
inf_init <- rep(100, N)
end_time <- 20
sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)

plot_stab <- plot_int1(N, sol, state = "INF") +
  theme_bw() +theme(legend.position="none") 
plot_stab