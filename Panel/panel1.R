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
thetas <- rep(0.1, N) # loss of immunity rates
mud <- 0.3
deaths <- rep(mud, N) # not disease-related death rates
mua <- 0.65
alphas <- rep(mua, N) # recovery rates
mudel <- 0
deltas <- rep(mudel, N) # disease-related death rates
gammas = deaths + alphas + deltas
mug = gammas[1]
# mobility
muw <- 0.07
sw <- 0.06
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


MIGRATION <- rand_mat(N, muc, sc, distrib = "beta")
diag(MIGRATION) <- 0
# ----------------------------------------------------------------------#
##### Plot Stability #####
step <- 0.0025
beta_vec <- seq(0.01,0.9,step)
muw_vec <- seq(0.01,0.9,step)
df_sol <- data.frame(beta = 0, gamma = 0, N = 0, muw = 0, state = FALSE)
N = 50
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
      betas <- rep(beta_vec[i], N)
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

# path <- paste0("~/RMT/David/OUTPUT/area_gen_",Sys.Date(),".csv")
# write.csv(df_sol, path,row.names = TRUE)

# path <- "~/RMT/David/OUTPUT/area_gen_2022-03-30.csv"
# df_sol <- read.csv(file = path)
df_sol$Stability <- ifelse(df_sol$state == TRUE, "Stable", "Unstable")

# df_sol <- df_sol[vec,]
library(latex2exp)

# Values for the points in the area graph:
stab_par <- 0.1
unstab_par <- 0.25

dif <- 0.02


betnew <- 0.35
muwnew <- ((betnew*stab_par*(N-1)) + betnew - stab_par)/(stab_par*(N-1))

betnew - mug + betnew*stab_par*(N-1)
betnew - mug + betnew*stab_par*(N-1) == stab_par - mug + stab_par*muwnew*(N-1)
# Create annotate for labels at each point in the area graph:
annotation <- data.frame(
  x = c(stab_par + dif, stab_par + dif, betnew + dif),
  y = c(stab_par + dif, muwnew + dif, stab_par + dif),
  label = c("c", "d", "e")
)

library(ggstar)
library("latex2exp")

color_stab <- "#48639C"
color_unstab <- "#BB4430"
# 
# color_stab <- "#034732"
# color_unstab <- "#008148"

plot_area <- ggplot(df_sol) +
  geom_point(aes(beta,muw, colour = Stability)) + theme_bw()  +
  scale_color_manual(values=c(color_stab, color_unstab)) +
  ylab(TeX("$\\mu_c$")) +
  xlab(TeX("$\\beta$")) +
  scale_x_continuous(breaks=c(0,0.25,0.50), limits = c(0, 0.5),
                     labels = c("0", "0.25", "0.50")) +
  scale_y_continuous(breaks=c(0,0.25,0.50, 0.75), limits = c(0, 0.75),
                     labels = c("0", "0.25", "0.50", "0.75")) +
  coord_fixed() +
  theme(text = element_text(size = 15), legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size=3))) 

data_points <- data.frame(mub = c(stab_par,stab_par,betnew), 
                          muw = c(stab_par,muwnew,stab_par))
color_points <- "#050505"
plot_area <-  ggplot(df_sol) +
  geom_point(aes(beta,muw, colour = Stability)) + theme_bw()  +
  scale_color_manual(values=c(color_stab, color_unstab)) +
  ylab(TeX("$\\mu_c$")) +
  xlab(TeX("$\\beta$")) +
  scale_x_continuous(breaks=c(0,0.25,0.50), limits = c(0, 0.5),
                     labels = c("0", "0.25", "0.50")) +
  scale_y_continuous(breaks=c(0,0.25,0.50, 0.75), limits = c(0, 0.75),
                     labels = c("0", "0.25", "0.50", "0.75")) +
  geom_point(data = data_points, aes(mub, muw),
             colour= color_points, size = 1.5) +
  # geom_point(aes(c(mub,muw)),
  #            colour= color_points, size = 1.5) +
  # geom_point(aes(mub,muwnew),
  #            colour= color_points, size = 1.5) +
  # geom_point(aes(betnew,muw),
  #            colour= color_points, size = 1.5) +
  geom_text(data=annotation, aes( x=x, y=y, label=label),
            color=color_points, 
            size=5.5 , angle=0, fontface="bold" ) +
  coord_fixed() +
  theme(text = element_text(size = 15), legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size=3)))  

plot_area
# plot_area <- plot_area + labs(title = "b")

# Save plot
# Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/Gen/"
path <- paste0(Path,"Area_g0,95_muc_0,01_sc0,00001_sw0,05.png")
ggsave(path,
       plot = plot_area, device = "png")
#-----------------------------------------------------------------------------#
######### PLOTS ##############

# #### Stability Area ####
# path <- "~/RMT/David/OUTPUT/area_gen_2022-03-30.csv"
# df_sol <- read.csv(file = path)
# df_sol$Stability <- ifelse(df_sol$state == TRUE, "Stable", "Unstable")
# 
# vec <- seq(0,nrow(df_sol),2)
# df_sol <- df_sol[vec,]
# library(latex2exp)
# 
# # Values for the points in the area graph:
# 
# # Create annotate for labels at each point in the area graph:
# text_size <- 15
# 
# annotation <- data.frame(
#   x = c(stab_par + 0.05, stab_par + 0.05, unstab_par + 0.05),
#   y = c(stab_par + 0.05,unstab_par + 0.05 , stab_par + 0.05),
#   label = c("c", "d", "e")
# )
# 
# library(ggstar)
# color_points <- "#FFFFFF"
# color_stab <- "#3066BE"
# color_unstab <- "#A63446"
# plot_area <- ggplot(df_sol) +
#   geom_point(aes(beta,muw, colour = Stability)) + theme_bw()  +
#   scale_color_manual(values=c(color_stab, color_unstab)) +
#   ylab(TeX("$\\mu_w$")) +
#   xlab(TeX("$\\beta$")) +
#   # ggtitle(""*gamma/beta~": 4")
#   ggtitle(paste0("N: ",N)) +
#   coord_fixed() +
#   geom_point(aes(stab_par,stab_par), colour= color_points, size = 2) +
#   geom_point(aes(stab_par,unstab_par), colour= color_points, size = 2) +
#   geom_point(aes(unstab_par,stab_par), colour= color_points, size = 2) +
#   geom_text(data=annotation, aes( x=x, y=y, label=label),
#             color=color_points, 
#             size=9 , angle=0, fontface="bold" ) + 
#   theme(text = element_text(size = text_size), legend.position = "bottom") +
#   guides(colour = guide_legend(override.aes = list(size=5)))
# 
# plot_area <- plot_area + labs(title = "b")
# 
# # Save plot
# Path <- "~/Documentos/PHD/2022/RMT_SIR/Plots/Gen/"
# path <- paste0(Path,"Area_g0,95_muc_0,01_sc0,00001_sw0,05.png")
# ggsave(path,
#        plot = plot_area, device = "png")
# 

#------------------------------------------------------------------------#
col_stab <- "#F3A712"
col_unst_c <- "#129490"
col_unst_b <- "#7A306C"

# col_stab <- "#C6C013"
# col_unst_c <- "#EF8A17"
# col_unst_b <- "#EF2917"
#### Plots RMT and Integration ####
sus_init <- rep(100000, N) # initial susceptibles
inf_init <- rep(100, N)    # initial infecteds

end_time <- 100

### Stable:
mub <- stab_par
betas <- rep(mub, N)
muw <- stab_par
COMMUTING <- rand_mat(N, muw, sw, distrib = "beta")
diag(COMMUTING) <- 0
# COMMUTING <- rand_mat_ell(N, muw, sw, rhow, distrib = "beta")
# COMMUTING[sample.int(N^2, round(p*N^2))] <- 0

jacobian_stab <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))

eig_stab <- eigen_mat(jacobian_stab)

# Plot the eigenvalues of the system
library("ggforce")
eigen_stab <- plot_eigen_rmt(jacobian_stab,
                             N,mub,mug = mud + mua + mudel,
                             muw,sw,rhow,Gammaw,
                             muc,sc,rhoc,Gammac,
                             tau = 0, alp = 0, K = 0) +
  scale_y_continuous( breaks=c(0)) 

eigen_stab

# Tntegro el sistema con condiciones iniciales 
sol_stab_f <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)

sol_stab <- sol_stab_f
for(i in c(1:N)){
  colnames(sol_stab)[i+1] <-  paste0("S",i)
  colnames(sol_stab)[N+i+1] <-  paste0("I",i)
  colnames(sol_stab)[2*N+i+1] <-  paste0("R",i)
}

sol_stab <- as.data.frame(sol_stab)
# z  <- z   %>% filter( z$time < 1) 
sol_stab <- reshape2::melt(sol_stab, id.vars = c("time"))

# Filter Infected:
sol_stab$type <- substr(sol_stab$variable,1,1)
sol_stab <- sol_stab  %>% filter( substr(sol_stab$variable,1,1) == "I")

# plot SUS, INF, REC or TOT population
plot_inf_stab <- plot_int(N, sol_stab_f, state = "INF")  + 
  theme_bw() + theme(legend.position = "none")

vec_col <-  vector(mode="character", length=N)
vec_col[1:N] <- col_stab

plot_inf_stab <- plot_inf_stab +
  xlim(c(0,20))  + 
  scale_colour_manual(values = vec_col) 

plot_inf_stab

### Unstable by commuting:
mub <- stab_par
betas <- rep(mub, N)
muw <- muwnew
COMMUTING <- rand_mat(N, muw, sw, distrib = "beta")
diag(COMMUTING) <- 0
# COMMUTING <- rand_mat_ell(N, muw, sw, rhow, distrib = "beta")
# COMMUTING[sample.int(N^2, round(p*N^2))] <- 0

jacobian_uns_com <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))

eig_uns_com <- eigen_mat(jacobian_uns_com)
# Plot the eigenvalues of the system
library("ggforce")
eigen_unstab_com <- plot_eigen_rmt(jacobian_uns_com,
                                   N,mub,mug = mud + mua + mudel,
                                   muw,sw,rhow,Gammaw,
                                   muc,sc,rhoc,Gammac,
                                   tau = 0, alp = 0, K = 0) +
  scale_y_continuous( breaks=c(0)) 

eigen_unstab_com

# Tntegro el sistema con condiciones iniciales 
sol_uns_com_f <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)

sol_uns_com <- sol_uns_com_f 
for(i in c(1:N)){
  colnames(sol_uns_com)[i+1] <-  paste0("S",i)
  colnames(sol_uns_com)[N+i+1] <-  paste0("I",i)
  colnames(sol_uns_com)[2*N+i+1] <-  paste0("R",i)
}

sol_uns_com <- as.data.frame(sol_uns_com)
# z  <- z   %>% filter( z$time < 1) 
sol_uns_com <- reshape2::melt(sol_uns_com, id.vars = c("time"))

# Filter Infected:
sol_uns_com$type <- substr(sol_uns_com$variable,1,1)
sol_uns_com <- sol_uns_com  %>% filter( substr(sol_uns_com$variable,1,1) == "I")

# plot SUS, INF, REC or TOT population
plot_inf_unstab_com <- plot_int(N, sol_uns_com_f, state = "INF")  +
  theme_bw() + 
  theme(legend.position = "none")

plot_inf_unstab_com
# plot_int(N, sol, state = "TOT")
vec_col <-  vector(mode="character", length=N)
vec_col[1:N] <- col_unst_c

plot_inf_unstab_com <- plot_inf_unstab_com +
  scale_colour_manual(values = vec_col) 
plot_inf_unstab_com <- plot_inf_unstab_com + xlim(c(0,20))


plot  <- ggplot(sol_stab,
                aes(time, value), show.legend = NA) + 
  geom_line(aes( group =variable, colour = type),size=0.5)  +
  # geom_line(aes( colour =variable),size=0.5)  +
  ylab("Number of infected individuals") 

### Unstable by commuting:
mub <- betnew
betas <- rep(mub, N)
muw <- stab_par
COMMUTING <- rand_mat(N, muw, sw, distrib = "beta")
diag(COMMUTING) <- 0
# COMMUTING <- rand_mat_ell(N, muw, sw, rhow, distrib = "beta")
# COMMUTING[sample.int(N^2, round(p*N^2))] <- 0

jacobian_uns_bet <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))

eig_uns_bet <- eigen_mat(jacobian_uns_bet)
# Plot the eigenvalues of the system
library("ggforce")
eigen_unstab_bet <- plot_eigen_rmt(jacobian_uns_bet,
                                   N,mub,mug = mud + mua + mudel,
                                   muw,sw,rhow,Gammaw,
                                   muc,sc,rhoc,Gammac,
                                   tau = 0, alp = 0, K = 0) +
  scale_y_continuous( breaks=c(0)) 

eigen_unstab_bet

# Tntegro el sistema con condiciones iniciales 
sol_uns_bet_f <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)


sol_uns_bet <- sol_uns_bet_f 
for(i in c(1:N)){
  colnames(sol_uns_bet)[i+1] <-  paste0("S",i)
  colnames(sol_uns_bet)[N+i+1] <-  paste0("I",i)
  colnames(sol_uns_bet)[2*N+i+1] <-  paste0("R",i)
}

sol_uns_bet <- as.data.frame(sol_uns_bet)
# z  <- z   %>% filter( z$time < 1) 
sol_uns_bet <- reshape2::melt(sol_uns_bet, id.vars = c("time"))

# Filter Infected:
sol_uns_bet$type <- substr(sol_uns_bet$variable,1,1)
sol_uns_bet <- sol_uns_bet  %>% filter( substr(sol_uns_bet$variable,1,1) == "I")

# plot SUS, INF, REC or TOT population
plot_inf_unstab_bet <- plot_int(N, sol_uns_bet_f, state = "INF")  +
  theme_bw() + theme(legend.position = "none")
plot_inf_unstab_bet <- plot_inf_unstab_bet + xlim(c(0,20))
# plot_int(N, sol, state = "TOT")
vec_col <-  vector(mode="character", length=N)
vec_col[1:N] <- col_unst_b

plot_inf_unstab_bet <- plot_inf_unstab_bet +
  scale_colour_manual(values = vec_col) 
plot_inf_unstab_bet

plot  <- ggplot(sol_uns_bet,
                aes(time, value), show.legend = NA) + 
  geom_line(aes( group =variable, colour = type),
            color = col_unst_c , size=0.5, linetype = "dashed")  +
  ylab("Number of infected individuals") +
  geom_line(data = sol_uns_com,
            aes( group =variable, colour = type),
            color = col_unst_b , size=0.5) + theme_bw()
  
###### ALL EIGENVALUE DIST ####
size_text <- 16
eig_uns_com$lab <- "unstab_com"
eig_uns_bet$lab <- "unstab_bet"
eig_stab$lab <- "stab"
eig_full <- rbind(rbind(eig_uns_com,eig_uns_bet),eig_stab)
ggplot(eig_full) + 
  geom_point(aes(re,im), size = 0.2) +
  theme_bw() + coord_fixed() +
  theme(text = element_text(size = size_text),
        legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size=3)))

mug <- gammas[1]
tau <- 0
c_stab <- stab_par*(1-stab_par) - mug - N*muc
r_stab <- sqrt(stab_par^2*sw^2 + 2*stab_par*tau + sc^2)*sqrt(N)
o_stab <- stab_par - mug + stab_par*stab_par*(N-1)

c_unstab_c <- stab_par*(1-muwnew) - mug - N*muc
r_unstab_c <- sqrt(stab_par^2*sw^2 + 2*stab_par*tau + sc^2)*sqrt(N)
o_unstab_c <- stab_par - mug + stab_par*muwnew*(N-1)

c_unstab_b <- betnew*(1-stab_par) - mug - N*muc
r_unstab_b <- sqrt(betnew^2*sw^2 + 2*betnew*tau + sc^2)*sqrt(N)
o_unstab_b <- betnew - mug + betnew*stab_par*(N-1)

colors <- c("c" = col_stab, 
            "d" = col_unst_c,
            "e" = col_unst_b)

eigen_full <- ggplot(eig_full) + 
  geom_point(aes(re,im), size = 0.2) +
  theme_bw() + coord_fixed() +
  theme(text = element_text(size = size_text),
        legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  geom_ellipse(aes(x0 = c_unstab_b, y0 = 0, 
                   a = r_unstab_b,
                   b = r_unstab_b, 
                   angle = 0, color = "e")) +
  geom_ellipse(aes(x0 = c_unstab_c, y0 = 0, 
                   a = r_unstab_c,
                   b = r_unstab_c, 
                   angle = 0, color = "d"))  +
  geom_ellipse(aes(x0 = c_stab, y0 = 0, 
                   a = r_stab,
                   b = r_stab, 
                   angle = 0, color = "c")) +
  geom_point(aes(o_stab,0), 
             color = col_stab, shape = 21) +
  geom_point(aes(o_unstab_b,0),
             color = col_unst_b, shape = 19, size = 2.4) +
  geom_point(aes(o_unstab_c,0), 
             color = col_unst_c, shape = 21, size = 0.6) +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed")  +
  scale_color_manual(values = colors,
                     name = " Scenario")
  
eigen_full  
  
##### Construct panel1 ####
library("ggpubr")

ggeigen <- ggarrange(eigen_stab  + rremove("xlab")  + rremove("ylab") + labs(title = "c") +
                       theme(text = element_text(size = size_text)) ,
                     eigen_unstab_com   + rremove("xlab")  + rremove("ylab") + labs(title = "d") + 
                       theme(text = element_text(size = size_text)),
                     eigen_unstab_bet  + rremove("xlab")  + rremove("ylab") + labs(title = "e")  + 
                       theme(text = element_text(size = size_text)),
                     nrow = 3, ncol = 1)
ggeigen  <- annotate_figure(ggeigen,
                         bottom = text_grob("Real part", color = "black",
                                            size = 15),
                         left = text_grob("Imaginary part",
                                          color = "black", rot = 90,  size = 15))
ggeigen

gginf <- ggarrange(plot_inf_stab + labs(title = "c") +
                     rremove("xlab")  + rremove("ylab") +
                     scale_y_continuous( breaks=c(0, 50, 100)) +
                     theme(text = element_text(size = size_text)),
                   plot_inf_unstab_com + labs(title = "d") +
                     rremove("xlab")  + rremove("ylab") +
                     scale_y_continuous( breaks=c(0,7500, 15000)) +
                     theme(text = element_text(size = size_text)),
                   plot_inf_unstab_bet + labs(title = "e") +
                      rremove("xlab")  + rremove("ylab") +
                     scale_y_continuous( breaks=c(0,7500, 15000)) +
                     theme(text = element_text(size = size_text)),
                   ncol = 1, nrow = 3)

### plot grid
gg1 <- plot_inf_stab + labs(title = "c") +
  scale_y_continuous(breaks=c(0, 50, 100),
                     labels = function(x) format(x, scientific = TRUE)) +
  theme(text = element_text(size = size_text),
        axis.title.y=element_blank(),
        axis.title.x=element_blank())

gg2 <-  plot_inf_unstab_com + labs(title = "d") +
  scale_y_continuous(breaks=c(0, 15000, 30000),
                     labels = function(x) format(x, scientific = TRUE)) +
  theme(text = element_text(size = size_text),
        axis.title.y=element_blank(),
        axis.title.x=element_blank())

gg3 <- plot_inf_unstab_bet + labs(title = "e") +
  scale_y_continuous(breaks=c(0,15000, 30000),
                     labels = function(x) format(x, scientific = TRUE)) +
  theme(text = element_text(size = size_text),
        axis.title.y=element_blank(),
        axis.title.x=element_blank())

library("cowplot")
library("ggpubr")
gg_grid <- plot_grid(gg1,
                     gg2,
                     gg3,
                     align = "v", ncol = 1)

gginf  <- annotate_figure(gg_grid,
                          bottom = text_grob("Time", color = "black",
                                             size = 15),
                          left = text_grob("Number of infected individuals",
                                           color = "black", rot = 90,  size = 15))

ggarr1 <- plot_grid(plot_area + ggtitle("a"), 
                    NULL,
                    gginf,
                    rel_widths = c(1,0.1,1),
                    ncol = 3) 
plot_grid(ggarr1,
          eigen_full  + ggtitle("b"),
          nrow = 2,
          ncol = 1,
          rel_heights = c(1.6,1))

########
gginf  <- annotate_figure(gginf,
                            bottom = text_grob("Time", color = "black",
                                               size = 15),
                            left = text_grob("Number of infected individuals",
                                             color = "black", rot = 90,  size = 15))

gginf

ggall <- ggarrange(eigen_full,gginf, ncol = 2)

ggall
path <- paste0(Path,"gg_g0,95_muc_0,001_sc0,0001_sw0,05.png")
ggsave(path,
       plot = ggarr, device = "png")


Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/panel1/diagram_c.png"
library("png")
diagram <- readPNG(Path)
im_A <- ggplot() + 
  background_image(diagram) +
  # This ensures that the image leaves so me space at the edges
  theme(plot.margin = margin(t=1, l=1, r=1, b=1, unit = "cm"))

gg_izq <- ggarrange(plot_area,
                    im_A, 
                    widths = c(1.2,2),
                    ncol = 2, common.legend = TRUE)
gg_izq

gg_tot <- ggarrange(gg_izq,ggall, 
                    heights = c(1.2,2),
                    nrow = 2)
gg_tot
