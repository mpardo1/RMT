####### RANDOM MATRICES FOR METAPOPULATION MODELS #######
### 
### parent script
###
### generate, plot and integrate metapopulation
### epidemiological models
### 
#########################################################
rm(list=ls())
source("~/RMT/David/RMT_genrandom.R")
source("~/RMT/David/RMT_plotmobility.R")
source("~/RMT/David/d_functions_eigen_int.R")
plotmobility <- function(mob, color1 = "#0000FF", cmin = "N", cmax = "N"){
  
  # cmin <- ifelse(cmin == "N", min(mob), cmin)
  # cmax <- ifelse(cmax == "N", max(mob), cmax)
  cmin = min(mob)
  cmax = max(mob)
  collow <- "white"
  
  N <- nrow(mob)
  if (is.matrix(mob) & (nrow(mob) == ncol(mob))) {
    
    mob <- t(mob)
    
    diag(mob) <- rep(0,N)
    
    #mob <- mob[seq(N,1,-1),]
    
    mob <- as.data.frame(mob)
    names(mob) <- c(1:N)
    mob$x <- c(N:1)
    mob <- mob %>% pivot_longer(as.character(c(1:N)), names_to = "y", values_to = "mob")
    
    mob$x <- factor(mob$x, c(N:1))
    mob$y <- factor(mob$y, c(N:1))
    
    ggplot(mob, aes(x,y)) +
      geom_tile(aes(fill = mob)) +
      scale_fill_gradient(low = collow, high = color1, na.value = "yellow", limits = c(cmin,cmax)) +
      theme_void() +
      theme(legend.position="none",
            panel.background=element_blank(),
            plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm")) +
      coord_fixed()
    
  } else {
    print("mob needs to be a square matrix")
  }
}
####### GENERATE JACOBIAN ###############################

# number of patches
N <- 50

# epidemiological
#all rates must lie in (0,1) except for betas

Deltas <- rep(0.3, N) # birth rate
mub <- 0.1
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
mug <- gammas[1]
# mobility
#commuting and migration networks

muw <- 0.6
sw <- 0.05
rhow <- 0 #original rho (Gamma of baron et al)
Gammaw <- 0 #gamma of baron et al
rw <- 0
cw <- 0

muc <- 0.01
sc <- 0.002
rhoc <- 0
Gammac <- 0
rc <- 0
cc <- 0

#Colors:
col_n <- "#3186D0"
col_s <- "#D7612A"

COMMUTING <- rand_mat(N, muw, sw, distrib = "beta")
diag(COMMUTING) <- 0
# COMMUTING <- rand_mat_ell(N, muw, sw, rhow, distrib = "beta")
# COMMUTING[sample.int(N^2, round(p*N^2))] <- 0

MIGRATION <- rand_mat(N, muc, sc, distrib = "beta")
diag(MIGRATION) <- 0
# MIGRATION <- matrix(0, N, N)
# jacobian

jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))

### SPARSE MATRIX ####
p <- 0.4
vec_p <- which(rbinom(N^2, 1, p)==0)
COM_SPA <- COMMUTING
COM_SPA[vec_p] <- 0
MIGRATION[vec_p] <- 0

jacobian_sparse <- (COM_SPA + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))

# Plot the eigenvalues of the system
sigma <- sqrt(mub^2*sw^2  + sc^2)
radius <- sigma*sqrt(N)

library("ggpubr")
library("ggforce")
outl <- mub*(muw*(N-1)+1) - mug
center <- mub*(1-muw) - mug - N*muc

plot_eig <- plot_eigen(jacobian) + 
  geom_point(aes(outl,0), color = "#5448C8") + 
  geom_circle(aes(x0 = center, y0 = 0, r = radius), color = "#5448C8")
plot_eig

outl_spa <- mub*(p*muw*(N-1)+1) - mug
center_spa <- mub*(1-p*muw) - mug - N*p*muc
# Variance of the product of Bernouilli and the initial jacobian:
sigma_spa <- sqrt((p*(1-p)*sigma^2)+(p*(1-p)*(mub*muw + muc)^2) + sigma^2*p^2)
radius_spa <- sqrt(N)*sigma_spa

plot_eig_spa <- plot_eigen(jacobian_sparse) + 
  geom_point(aes(outl_spa,0), color = "#5448C8") + 
  geom_circle(aes(x0 = center_spa, y0 = 0, r = radius_spa), color = "#5448C8") + 
  theme_bw()
plot_eig_spa

ggarrange(plot_eig, plot_eig_spa)

### MOBILITY MATRIX ##

# plot the mobility network
#legend for plotmobility2 can be found in RMT_plotmobility
plotmobility(COM_SPA)
plot_mob <- plotmobility(COM_SPA, color1 = col_s)
plotmobility(MIGRATION, color1 = col_s)
#### MEAN = p*muw ###  
muw <- p*muw
muc <- p*muc
center_mean <- mub*(1-muw) - mug - N*muc
outl_mean <- mub*(muw*(N-1)+1) - mug

COMMUTING <- rand_mat(N, muw, sw, distrib = "beta")
diag(COMMUTING) <- 0
MIGRATION <- rand_mat(N, muc, sc, distrib = "beta")
diag(MIGRATION) <- 0

jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))

plot_eig_mean <- plot_eigen(jacobian) + 
  geom_point(aes(outl_mean,0), color = "#5448C8") + 
  geom_circle(aes(x0 = center_mean, y0 = 0, r = radius), color = "#5448C8") +
  theme_bw()

#### FULL EIGEN ###
df_eigen_n <- eigen_mat(jacobian)
df_eigen_n$mat <- "Fully connected"
df_eigen_s <- eigen_mat(jacobian_sparse)
df_eigen_s$mat <- "Sparse"
df_eigen <- rbind(df_eigen_n, df_eigen_s)


plot <- ggplot(df_eigen) + 
  geom_point(aes(re, im, color = mat), size = 0.2) +
  scale_color_manual(name = NULL,values = c(col_n,col_s )) + 
  geom_point(aes(outl_mean,0), color = "#000000") + 
  geom_circle(aes(x0 = center_mean, y0 = 0, r = radius), color = "#000000") + 
  geom_point(aes(outl_spa,0), color = "#000000") + 
  geom_circle(aes(x0 = center_spa, y0 = 0, r = radius_spa), color = "#000000") + 
  theme_bw() + guides(colour = guide_legend(override.aes = list(size=2))) + 
  theme(legend.position = "left")

leg <- get_legend(plot)
plot <- ggplot(df_eigen) + 
  geom_point(aes(re, im), size = 0.1)  + 
  geom_point(aes(outl_spa,0), color = col_s, size = 3, pch = 8) + 
  geom_point(aes(outl_mean,0), color = col_n, size = 2) + 
  geom_circle(aes(x0 = center_mean, y0 = 0, r = radius), color = col_n) + 
  geom_circle(aes(x0 = center_spa, y0 = 0, r = radius_spa), color = col_s) + 
  theme_bw() + guides(colour = guide_legend(override.aes = list(size=2)))### MOBILITY MATRIX ##

plot <- plot_grid(plot + rremove("xlab") + rremove("ylab"),
          leg, nrow = 1, rel_widths  = c(1,0.3))
# plot the mobility network
#legend for plotmobility2 can be found in RMT_plotmobility
plotmobility(COMMUTING)
plot_mob_mean <- plotmobility(COMMUTING, color1 = col_n)

library("latex2exp")
library("cowplot")

gg_arr <- plot_grid(plot_mob + theme(aspect.ratio = 1),
          plot_mob_mean+ theme(aspect.ratio = 1),
          nrow = 1, ncol = 2,
          align = "h")
# gg_arr <- plot_grid(plot,gg_arr , 
#           ncol = 1, rel_heights = c(0.7,1))
# gg_arr

Path <- "~/Documentos/PHD/2022/RMT_SIR/Plots/SM/"
path <- paste0(Path,"Plot1_spa_mean1_b0,1_g0,95_muc_0,004_sc0,002_muw0,24_sw0,05.pdf")
ggsave(path,
       plot = plot_f, device = "pdf")

#######INTEGRATION#######
end_time <- 100
sus_init <- rep(1000, N) # initial susceptibles
inf_init <- rep(10, N)
sol_n <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
                  COMMUTING,MIGRATION,
                  sus_init,inf_init,end_time)

sol_s <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
                  COM_SPA,MIGRATION,
                  sus_init,inf_init,end_time)

vec_col <-  vector(mode="character", length=N)
vec_col[1:N] <- col_n
plotint_n <- plot_int(N,sol_n, state ="INF")   + 
  rremove("ylab") +
  scale_colour_manual(values = vec_col) + theme_bw()  +
  theme(text = element_text(size = 15),legend.position = "none")
  
vec_col <-  vector(mode="character", length=N)
vec_col[1:N] <- col_s
plotint_s <- plot_int(N,sol_s, state ="INF")+ 
  ylab("Infected Individuals") +
  scale_colour_manual(values = vec_col) + theme_bw() +
  theme(text = element_text(size = 15),legend.position = "none")

plotint <- plot_grid(plotint_s + xlim(c(0,20))+
                       theme(plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm")) ,
                     plotint_n + ylab("")  + xlim(c(0,20)) +
                       theme(plot.margin = unit(c(0.1,0.3,0.1,0.1),"cm")),  ncol = 2)

ggarr <- plot_grid(plotint_s + theme(aspect.ratio = 1), 
                    plotint_n + ylab("")+ theme(aspect.ratio = 1),
                    plot_mob + theme(aspect.ratio = 1),
                   plot_mob_mean+ theme(aspect.ratio = 1),
                   nrow = 2, ncol = 2, align = "v")
plot_f <- plot_grid(gg_arr , plot ,plotint,  ncol = 1, rel_heights = c(0.9,0.6,1))
plot_f
# Compute the difference between the right most eigenvalue with sparse
# # and with the matix with mean p*muc, p*muw
# p_vec <- seq(0.1,1,0.01)
# dim <- length(p_vec)
# df_spa <- data_frame(p = 0, max_eig_m = 0, max_eig_spa = 0)
# for(i in c(1:dim)){
#   COMMUTING <- rand_mat(N, p_vec[i]*muw, sw, distrib = "beta")
#   diag(COMMUTING) <- 0
#   MIGRATION <- rand_mat(N, p_vec[i]*muc, sc, distrib = "beta")
#   diag(MIGRATION) <- 0
#   
#   jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
#     diag(deaths + alphas + deltas + colSums(MIGRATION))
#   
#   eig_mean <- eigen_mat(jacobian)
#   max_eig_mean <- max(eig_mean$re)
#   
#   #### Sparse matrix ###
#   COMMUTING <- rand_mat(N, muw, sw, distrib = "beta")
#   diag(COMMUTING) <- 0
#   MIGRATION <- rand_mat(N, muc, sc, distrib = "beta")
#   diag(MIGRATION) <- 0
#   vec_p <- which(rbinom(N^2, 1, p_vec[i])==0)
#   COMMUTING[vec_p] <- 0
#   MIGRATION[vec_p] <- 0
#   
#   ### MOBILITY MATRIX #
#   
#   # plot the mobility network
#   #legend for plotmobility2 can be found in RMT_plotmobility
#   plotmobility(COMMUTING)
#   plotmobility2(MIGRATION, COMMUTING)
#   jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
#     diag(deaths + alphas + deltas + colSums(MIGRATION))
#   
#   eig_spa <- eigen_mat(jacobian)
#   max_eig_spa <- max(eig_spa$re)
#   
#   df_spa[nrow(df_spa)+1,] <- list(p_vec[i], max_eig_mean, max_eig_spa)
# }
# 
# df_spa <- df_spa[-1,]
# df_plot <- reshape2::melt(df_spa, id.vars="p")
# gg_max_eig <- ggplot(df_plot) +
#   geom_line(aes(p, value, colour=variable)) +
#   xlab("Right most eigenvalue real part") +
#   theme_bw()
# 
Path <- "~/Documentos/PHD/2022/RMT_SIR/Plots/SM/"
path <- paste0(Path,"Plot1_spa_mean1_b0,1_g0,95_muc_0,004_sc0,002_muw0,24_sw0,05.pdf")
ggsave(path,
       plot = plot, device = "pdf")
