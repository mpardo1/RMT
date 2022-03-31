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
mub <- 0.2
sb <- 0.001
betas <- rep(mub, N) # transmission rates
# betas <- rgamma(N, shape = (mub/sb)^2, rate = mub/(sb^2))
thetas <- rep(0.1, N) # loss of immunity rates
mud <- 0.3
deaths <- rep(mud, N) # not disease-related death rates
mua <- 0.65
alphas <- rep(mua, N) # recovery rates
mudel <- 0
deltas <- rep(mudel, N) # disease-related death rates
gammas = deaths + alphas + deltas

# mobility
muw <- 0.07
sw <- 0.05
rhow <- 0 #original rho (Gamma of baron et al)
Gammaw <- 0 #gamma of baron et al
rw <- 0
cw <- 0

muc <- 0.001
sc <- 0.0002
rhoc <- 0
Gammac <- 0
rc <- 0
cc <- 0


MIGRATION <- rand_mat(N, muc, sc, distrib = "beta")
diag(MIGRATION) <- 0
# ----------------------------------------------------------------------#
##### Plot Stability #####
step <- 0.001
beta_vec <- seq(0.01,0.9,step)
muw_vec <- seq(0.01,0.9,step)
df_sol <- data.frame(beta = 0, gamma = 0, N = 0, muw = 0, state = FALSE)
N = 100
for(i in c(1:length(beta_vec))){
  print(paste0(" i : ", i))
  for(j in c(1:length(beta_vec))){
    if(is.na(muw_vec[j]) | muw_vec[j] > 1 | muw_vec[j] <0  ){
      print(paste0("Problem muw: ", muw_vec[j]))
    }else{
      # Computed by the numerically computed eigevalues;
      # COMMUTING <- rand_mat(N, muw_vec[j], sw, distrib = "beta")
      # diag(COMMUTING) <- 0
      # MIGRATION <- rand_mat(N, muc, sc, distrib = "beta")
      # diag(MIGRATION) <- 0
      
      # EPI param:
      # betas <- rep(beta_vec[i], N)
      # deltas <- rep(mudel, N)
      # deaths <- rep(mud, N) 
      # alphas <- rep(mua, N)
      
      # jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
      #   diag(deaths + alphas + deltas + colSums(MIGRATION))
      # eigen <- eigen_mat(jacobian)
      # max_eig <- max(eigen$re)
      
      # Computed by the prediction from RMT:
      state <- ifelse(((N-1)*muw_vec[j]) < ((gammas[1]/betas[1]) - 1), TRUE, FALSE)
      df_sol[nrow(df_sol) + 1,1:4] <- list(betas[1], gammas[1], N, muw_vec[j])
      df_sol[nrow(df_sol) ,5] <- state
    }
  }
}

df_sol <- df_sol[-1,]

path <- paste0("~/RMT/David/OUTPUT/area_gen_",Sys.Date(),".csv")
write.csv(df_sol, path,row.names = TRUE)

path <- "~/RMT/David/OUTPUT/area_gen_2022-03-30.csv"
df_sol <- read.csv(file = path)
df_sol$Stability <- ifelse(df_sol$state == TRUE, "Stable", "Unstable")

vec <- seq(0,nrow(df_sol),2)
df_sol <- df_sol[vec,]
library(latex2exp)

# Values for the points in the area graph:
stab_par <- 0.07
unstab_par <- 0.2

# Create annotate for labels at each point in the area graph:
annotation <- data.frame(
  x = c(stab_par + 0.5,stab_par + 0.5, unstab_par + 0.04),
  y = c(stab_par + 0.3,unstab_par + 0.02 , stab_par + 0.3),
  label = c("c", "d", "e")
)

library(ggstar)
plot_area <- ggplot(df_sol) +
  geom_point(aes(beta,muw, colour = Stability)) + theme_bw()  +
  scale_color_manual(values=c("#3066BE", "#A63446")) +
  ylab(TeX("$\\mu_w$")) +
  xlab(TeX("$\\beta$")) +
  # ggtitle(""*gamma/beta~": 4")
  ggtitle(paste0("N: ",N)) +
  coord_fixed() +
  geom_point(aes(stab_par,stab_par), colour= "#ADA544", size = 3) +
  geom_point(aes(stab_par,unstab_par), colour= "#ADA544", size = 3) +
  geom_point(aes(unstab_par,stab_par), colour= "#ADA544", size = 3) +
  geom_text(data=annotation, aes( x=x, y=y, label=label),
            color="#ADA544", 
            size=9 , angle=0, fontface="bold" ) + 
  theme(text = element_text(size = 30), legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size=5)))

plot_area <- plot_area + labs(title = "b")

# Save plot
# Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/Gen/"
path <- paste0(Path,"Area_g0,95_muc_0,01_sc0,00001_sw0,05.png")
ggsave(path,
       plot = plot_area, device = "png")

#------------------------------------------------------------------------#
#### Plots RMT and Integration ####
sus_init <- rep(100000, N) # initial susceptibles
inf_init <- rep(100, N)    # initial infecteds

end_time <- 20

### Stable:
mub <- stab_par
betas <- rep(mub, N)
muw <- stab_par
COMMUTING <- rand_mat(N, muw, sw, distrib = "beta")
diag(COMMUTING) <- 0
# COMMUTING <- rand_mat_ell(N, muw, sw, rhow, distrib = "beta")
# COMMUTING[sample.int(N^2, round(p*N^2))] <- 0

jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))

# Plot the eigenvalues of the system
library("ggforce")
eigen_stab <- plot_eigen_rmt(jacobian,
                             N,mub,mug = mud + mua + mudel,
                             muw,sw,rhow,Gammaw,
                             muc,sc,rhoc,Gammac,
                             tau = 0, alp = 0, K = 0) +
  scale_y_continuous( breaks=c(0)) 

eigen_stab

# Tntegro el sistema con condiciones iniciales 
sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)

# plot SUS, INF, REC or TOT population
plot_inf_stab <- plot_int(N, sol, state = "INF")  + theme_bw()+ theme(legend.position = "none")
plot_inf_stab

### Unstable by commuting:
mub <- stab_par
betas <- rep(mub, N)
muw <- unstab_par
COMMUTING <- rand_mat(N, muw, sw, distrib = "beta")
diag(COMMUTING) <- 0
# COMMUTING <- rand_mat_ell(N, muw, sw, rhow, distrib = "beta")
# COMMUTING[sample.int(N^2, round(p*N^2))] <- 0

jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))

# Plot the eigenvalues of the system
library("ggforce")
eigen_unstab_com <- plot_eigen_rmt(jacobian,
                             N,mub,mug = mud + mua + mudel,
                             muw,sw,rhow,Gammaw,
                             muc,sc,rhoc,Gammac,
                             tau = 0, alp = 0, K = 0) +
  scale_y_continuous( breaks=c(0)) 

eigen_unstab_com

# Tntegro el sistema con condiciones iniciales 
sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)

# plot SUS, INF, REC or TOT population
plot_inf_unstab_com <- plot_int(N, sol, state = "INF")  + theme_bw()+ theme(legend.position = "none")
plot_inf_unstab_com
# plot_int(N, sol, state = "TOT")

### Unstable by commuting:
mub <- unstab_par
betas <- rep(mub, N)
muw <- stab_par
COMMUTING <- rand_mat(N, muw, sw, distrib = "beta")
diag(COMMUTING) <- 0
# COMMUTING <- rand_mat_ell(N, muw, sw, rhow, distrib = "beta")
# COMMUTING[sample.int(N^2, round(p*N^2))] <- 0

jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))

# Plot the eigenvalues of the system
library("ggforce")
eigen_unstab_bet <- plot_eigen_rmt(jacobian,
                                   N,mub,mug = mud + mua + mudel,
                                   muw,sw,rhow,Gammaw,
                                   muc,sc,rhoc,Gammac,
                                   tau = 0, alp = 0, K = 0) +
  scale_y_continuous( breaks=c(0)) 

eigen_unstab_bet

# Tntegro el sistema con condiciones iniciales 
sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)

# plot SUS, INF, REC or TOT population
plot_inf_unstab_bet <- plot_int(N, sol, state = "INF")  +
  theme_bw() + theme(legend.position = "none")
plot_inf_unstab_bet
# plot_int(N, sol, state = "TOT")

#-----------------------------------------------------------------------------#
##### Contructi panel1 ####
library("ggpubr")
size_text <- 16
ggeigen <- ggarrange(eigen_stab  + xlab("")  + xlab("") + labs(title = "c") +
                       theme(text = element_text(size = size_text)) ,
                     eigen_unstab_com  + xlab("") + labs(title = "d") + 
                       theme(text = element_text(size = size_text)),
                     eigen_unstab_bet + labs(title = "e")  + 
                       theme(text = element_text(size = size_text)),
                     nrow = 3, ncol = 1)

gginf <- ggarrange(plot_inf_stab + labs(title = "c") +
                     ylab("Infected individuals") + 
                     xlab("") +
                     theme(text = element_text(size = size_text)),
                   plot_inf_unstab_com + labs(title = "d") +
                     ylab("Infected individuals") + xlab("") +
                     theme(text = element_text(size = size_text)),
                   plot_inf_unstab_bet + labs(title = "e") +
                     ylab("Infected individuals") +
                     theme(text = element_text(size = size_text)),
                   ncol = 1, nrow = 3)

ggarr <- ggarrange(ggeigen,gginf, ncol = 2)

path <- paste0(Path,"gg_g0,95_muc_0,001_sc0,0001_sw0,05.png")
ggsave(path,
       plot = ggarr, device = "png")


Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/diagram_c.png"
library("jpeg")
diagram <- readPNG(Path)
im_A <- ggplot() + 
  background_image(diagram) +
  # This ensures that the image leaves so me space at the edges
  theme(plot.margin = margin(t=3, l=4.8, r=4.8, b=3, unit = "cm"))

gg_izq <- ggarrange(plot_area,
                    im_A, 
                    ncol = 2)
ggarrange(gg_izq,ggall)