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
library("viridis")
library("ggpubr")

###### COLORS & SIZE#####
color_stab <- "#084C61"
color_unstab <- "#DB504A"
size_let <- 15
####### GENERATE JACOBIAN ###############################
# number of patches
N <- 50

# epidemiological
#all rates must lie in (0,1) except for betas

Deltas <- rep(0.3, N) # birth rate
mub <- 0.2
sb <- 0.001
betas <- rep(mub, N) # transmission rates
# betas <- rgamma(N, shape = (mub/sb)^2, rate = mub/(sb^2))
thetas <- rep(0.3, N) # loss of immunity rates
mud <- 0.3
deaths <- rep(mud, N) # not disease-related death rates
mua <- 0.5
alphas <- rep(mua, N) # recovery rates
mudel <- 0
deltas <- rep(mudel, N) # disease-related death rates
gammas = deaths + alphas + deltas

# mobility
#commuting and migration networks
muw <- 0.05
sw <- 0.01
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

sus_init <- rep(1000, N) # initial susceptibles
inf_init <- rep(10, N)    # initial infecteds

end_time <- 50
sol.stab <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)

plot_stab <- plot_int1(N, sol.stab, state = "INF") +
  theme_bw() +theme(legend.position="none")
plot_stab


vec_col <-  vector(mode="character", length=N)
vec_col[1:N] <- color_stab

plot.inf.stab <- plot_int1(N, sol.stab, state = "INF") +
  scale_colour_manual(values = vec_col) +
  theme_bw() +
  theme(text = element_text(size = size_let),
        legend.position="none")

plot.inf.stab

Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/panel2/"
path <- paste0(Path,"Plot_inf_b0,02_g0,5_muc_0,01_sc0,001_muw0,08_sw0,05.png")
ggsave(path,
       plot = plot.inf.stab, device = "png")

#### 1 PATCH modified ####
alp_bet <- 0.8

betas <- rep(mub, N)
betas[1] <- alp_bet + betas[1]

jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))

plot_eigen(jacobian)
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
vec_col[1:N] <- color_stab
vec_col[1] <- color_unstab

plot.inf.1 <- plot_int1(N, sol, state = "INF") +
  scale_colour_manual(values = vec_col) +
  theme_bw() + theme(legend.position="none") +
  theme(text = element_text(size = size_let))

plot.inf.1


Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/panel2/"
path <- paste0(Path,"Plot_inf1_alp_1,1_b0,02_g0,5_muc_0,01_sc0,001_muw0,08_sw0,05.png")
ggsave(path,
       plot = plot.inf.1, device = "png")

##### Sum of all the infected by alpha ####
sus_init <- rep(1000, N) # initial susceptibles
inf_init <- rep(10, N)    # initial infecteds

end_time <- 500
alp_bet_vec <- seq(1.7,1.8,0.01)
# end_time <- 500
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

vec_eq <- c()
rm_ind <- c()
rm_ind[1] <- 0
for(i in c(1: (nrow(df_sum)-10))){
  for(j in c(1:(ncol(df_sum)-1))){
    if((floor(df_sum[i,j+1]) == floor(df_sum[i +30,j+1])) &
       length(which(j==rm_ind)) == 0){
      print("Cumple condición")
      vec_eq[length(vec_eq) + 1] <- df_sum[i,j+1]
      rm_ind[length(rm_ind)+1] <- j
    }
  }
}

vec_eq
rm_ind
## DF for maximum number of infected individuals and time of maximum.
max_inf <- df_sum %>% summarise_if(is.numeric, max)
max_inf_1 <- apply(df_sum, 2, max)
#[2:length(apply(df_sum, 2, max))]
df <- data.frame(max_inf,max_inf_1)
time_max_vec <- c()
ind_vec <- c()
len_vec <- length(alp_bet_vec) - 2
for(i in c(1:len_vec)){
  vec <- df_sum[,i+1]
  ind <- which(vec == vec_eq[i])[1]
  # time_max_vec[i] <- df_sum[ind,1]
  ind_vec[i] <- ind
}

time_max_vec <- df_sum[ind_vec,1]
df_sum_group <- data.frame(alp <- alp_bet_vec[1:len_vec], 
                           max_inf <- t(max_inf)[1:len_vec],
                           time_max <- time_max_vec[1:len_vec])

df_sum_group <-  df_sum_group[-1,]
colnames(df_sum_group) <-  c("alpha", "Max_inf", "Time_max")
library("latex2exp")
plot_inf_max <- ggplot(df_sum_group) + 
  geom_line(aes(alpha, Max_inf)) + 
  xlab(TeX("$\\beta^*$")) +
  ylab("Max of infected individuals") +
  theme_bw() +
  theme(text = element_text(size = size_let)) 
plot_inf_max

# Check why it happens the fluctuations
df_plot_filt <- df_sum[,c(1,65)]
colnames(df_plot_filt) <- c("Time", "I")
ggplot(df_plot_filt) + 
  geom_line(aes(Time, log10(I))) + 
  ylim(c(58675.2,58675.3))

plot_time_max <-  ggplot(df_sum_group) + 
  geom_line(aes(alpha, Time_max), colour = "#69995D" , size = 0.5) + 
  geom_point(aes(alpha, Time_max), colour = "#050223" , size = 0.8) + 
  xlab(TeX("$\\beta^*$")) +
  ylab("Time of max of infected individuals") +
  theme_bw() +
  theme(text = element_text(size = size_let)) 
plot_time_max
# + xlim(c(0.7,1.7)) + ylim(c(80,110))

seq <- seq(1,len_vec,5 )
colnames(df_sum) <- c("time",as.character(alp_bet_vec))
filt_sum <- df_sum[, seq]
df_plot <- reshape2::melt(filt_sum, id.vars = c("time"))
# df_plot <- reshape2::melt(df_sum, id.vars = c("time"))

size_let <- 13
sum_inf <- ggplot(data = df_plot, aes(x = time, y = value,
                      color = as.numeric(as.character(variable)),
                      group = variable)) +
  geom_line() +
  scale_colour_gradient(name = TeX("$\\beta^*$"),
                        low = "#F8F053", high = "#4C0EF6") +
  ylab("Sum of infected individuals") +
  # xlim( c(0,10) ) +
  theme_bw() +
  theme(text = element_text(size = size_let),
        legend.position = "right") 

sum_inf + xlim(c(0,40)) 

### Solve equation for alpha ####
k <-  1
mug <- gammas[1]
outl1 <- function(x){
  (N/2)*(mub*muw + muc) + (x/2)*(1 + (k-1)*muw) + mub*(1-muw) -
    mug - N*muc + (1/2)*sqrt(N^2*(mub*muw + muc)^2 + x^2*(1 + (k-1)*muw)^2 + 
                               2*x*(mub*muw + muc)*(N*(1+(k-1)*muw) + 2*(N-k)*(muw-1)))
} 
curve(outl1(x), 0,3)
abline(h = 0, lty = 3)
uni <- uniroot(outl1, c(0, 3))$root
print(paste0("i: ", i))

betas <- rep(mub, N) 
betas[1] <- betas[1] + uni
sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)
sol_inf <- sol[,c(1,(N+2):(2*N+1))]
sum_inf <- rowSums(sol_inf[,2:(N+1)])
df_root <- data.frame(time = as.data.frame(sol)$time , sum_inf = sum_inf)

ggplot() + geom_line(data = df_root, aes(time,sum_inf))
sum_inf <- ggplot(data = df_plot) +
  geom_line(aes(x = time, y = value,
                  color = as.numeric(as.character(variable)),
                  group = variable)) +
  geom_line(data = df_root, aes(time,sum_inf), linetype = "dashed") +
  scale_colour_gradient(name = TeX("$\\beta^*$"),
                        low = "#F8F053", high = "#4C0EF6") +
  ylab("Sum of infected individuals") +
  xlim( c(0,20) ) +
  theme_bw() +
  theme(text = element_text(size = size_let), legend.position = "right") 

sum_inf

#### GGARRANGE ####
require(grid) 
gg1 <- ggarrange(plot.inf.stab + rremove("ylab") + rremove("xlab"),
                 plot.inf.1 + rremove("ylab") + rremove("xlab"),
                 nrow=1)
gg1<- annotate_figure(gg1, left = textGrob("Infected individuals", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                bottom = textGrob("Time", gp = gpar(cex = 1.3)))

gg2 <- ggarrange(sum_inf  + rremove("xlab") + ylab("Sum Infected"), 
          gg1, ncol = 1, labels = c("a", "b"))

gg3 <- ggarrange(plot_inf_max  +
                   rremove("xlab") + 
                   ylab("Maximum Infected"), 
                 plot_time_max + 
                   ylab("Time for maximum"),
                 ncol = 1)
ggall <- ggarrange(gg2,
                   gg3,
                   widths = c(1.8,1),
                   labels = c("","c"))
ggall
## Save image
Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/panel2/"
path <- paste0(Path,"Sum_inf_b0,1_g0,5_muc_0,01_sc0,001_muw0,05_sw0,001.png")
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
  theme(text = element_text(size = size_let), legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size=5)))
plot_area

Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/panel2/"
path <- paste0(Path,"Area_g0,5_muc_0,01_sc0,001_muw0,08_sw0,05.png")
ggsave(path,
       plot = plot_area, device = "png")

path <- paste0("~/RMT/David/OUTPUT/area_gen_beta_alp_",Sys.Date(),".csv")
write.csv(df_sol, path,row.names = TRUE)

#### PLOTS #####
# Stability
vec_col <-  vector(mode="character", length=N)
vec_col[1:N] <- "#A63446"

size_let <- 13
plot.inf.stab <- plot_int1(N, sol.stab, state = "INF") +
  scale_colour_manual(values = vec_col) +
  theme_bw() +
  theme(text = element_text(size = size_let), legend.position="none")

plot.inf.stab

# 1 patch different transmission rate:
vec_col <-  vector(mode="character", length=N)
vec_col[1:N] <- "#A63446"
vec_col[1] <- "#3066BE"

plot.inf.1 <- plot_int1(N, sol, state = "INF") +
  scale_colour_manual(values = vec_col) +
  theme_bw() + theme(legend.position="none") +
  theme(text = element_text(size = size_let))

plot.inf.1

# Max infected
library("tidyverse")
Path <- "~/RMT/David/OUTPUT/"
path <- paste0(Path,"Suminf_g0,5_muc_0,001_sc0,0001_muw0,2_sw0,05_t2000_2022-04-04.csv")
alp_bet_vec <- seq(0,2.5,0.01)
df_sum <- read.csv(file = path)
df_sum <- df_sum[,-1]
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
  theme(text = element_text(size = size_let)) 


# Time max infected
plot_time_max <-  ggplot(df_sum_group) + 
  geom_line(aes(alpha, Time_max)) + 
  xlab(TeX("$\\alpha$")) +
  ylab("Time of max of infected individuals") +
  theme_bw() +
  theme(text = element_text(size = size_let)) 


colnames(df_sum) <- c("time",as.character(alp_bet_vec))
df_sum_filt <- df_sum[,c(1,seq(1,ncol(df_sum),7))]
df_plot <- reshape2::melt(df_sum_filt, id.vars = c("time"))

# Sum of infected
sum_inf <- ggplot(data = df_plot, aes(x = time, y = value,
                                      color = as.numeric(as.character(variable)),
                                      group = variable)) +
  geom_line() +
  scale_colour_gradient(name = ""*alpha~" ", 
                        low = "blue", high = "red",
                        limits = c(0,5),
                        breaks = c(0,1,5),
                        labels = c(0,1,5)) +
  ylab("Sum of infected individuals") +
  xlim(c(0,20)) +
  ggtitle("") + 
  theme_bw() +
  theme(text = element_text(size = size_let), legend.position = "top") 

sum_inf

# Area:
path <- "~/RMT/David/OUTPUT/Areaepi_g0,5_muc_0,001_sc0,0001_muw0,1_sw0,05_2022-04-01.csv"
df_sol <- read.csv(file = path)
df_sol <- df_sol[-1,]

library("latex2exp")
df_sol$Stability <- ifelse(df_sol$state == TRUE, "Stable", "Unstable")

plot_area <- ggplot(df_sol) +
  geom_point(aes(beta,alp, colour = Stability)) + theme_bw()  +
  scale_color_manual(values=c("#3066BE", "#A63446")) +
  ylab(TeX("$\\alpha$")) +
  xlab(TeX("$\\beta$")) +
  # ggtitle(""*gamma/beta~": 4")
  theme(text = element_text(size = size_let), legend.position = "top") +
  guides(colour = guide_legend(override.aes = list(size=3)))
plot_area

## join all plots
library("ggpubr")
gg_1 <- ggarrange(plot.inf.stab,
                  plot.inf.1,
                  plot_inf_max,
                  ncol = 3,nrow = 1,
                  labels = c("a", "b", "c"))

gg_2 <- ggarrange(plot_area,
                  sum_inf,plot_time_max,
                  ncol = 3,nrow = 1,
                  labels = c("d", "e", "f"))

gg_full <- ggarrange(gg_1,gg_2,
          nrow = 2, ncol = 1)
gg_full

Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/panel2/"
path <- paste0(Path,"panel2.png")
ggsave(path,
       plot = gg_full, device = "png")


