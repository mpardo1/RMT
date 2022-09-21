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
library("latex2exp")

###### COLORS & SIZE#####
color_stab <- "#FDC016"
color_unstab <- "#7A16FD"
size_let <- 15
####### GENERATE JACOBIAN ###############################
# number of patches
N <- 50

# epidemiological
#all rates must lie in (0,1) except for betas

Deltas <- rep(0.3, N) # birth rate
mub <- 0.15
sb <- 0.001
betas <- rep(mub, N) # transmission rates
# betas <- rgamma(N, shape = (mub/sb)^2, rate = mub/(sb^2))
thetas <- rep(0.3, N) # loss of immunity rates
mud <- 0.3
deaths <- rep(mud, N) # not disease-related death rates
mua <- 0.4
alphas <- rep(mua, N) # recovery rates
mudel <- 0
deltas <- rep(mudel, N) # disease-related death rates
gammas = deaths + alphas + deltas

# mobility
#commuting and migration networks
muw <- 0.05
sw <- 0.03
rhow <- 0 #original rho (Gamma of baron et al)
Gammaw <- 0 #gamma of baron et al
rw <- 0
cw <- 0

muc <- 0.002
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

#--------------------------------------------------------#
##### All betas the same
# jacobian
betas <- rep(mub, N)
jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))

# plot the eigenvalues of the system
library("ggforce")
plot_eig_stab <- plot_eigen_rmt(jacobian,
                                N,mub,mug = mud + mua + mudel,
                                muw,sw,rhow,Gammaw,
                                muc,sc,rhoc,Gammac,
                                tau = 0, alp = 0, K = 0) +
  scale_y_continuous( breaks=c(0))
# + xlim(c(-60,-50))
print(plot_eigen(jacobian))
eigen_stab <-  eigen_mat(jacobian)

# Integrate the system:
# initial populations
# for constant populations, set deltas = 0, Deltas = deaths

sus_init <- rep(10000, N) # initial susceptibles
inf_init <- rep(100, N)    # initial infecteds

end_time <- 500
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
alp_bet <- 0.7

betas <- rep(mub, N)
betas[1] <- alp_bet + betas[1]

jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))

plot_eigen(jacobian)
plot_eigen_uns <- plot_eigen_rmt(jacobian,
                                 N,mub,mug = mud + mua + mudel,
                                 muw,sw,rhow,Gammaw,
                                 muc,sc,rhoc,Gammac,
                                 tau = 0, alp = alp_bet, K = 1) +
  scale_y_continuous( breaks=c(0))

eigen_unst <- eigen_mat(jacobian)

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

ggarrange(plot_eig_stab, plot_eigen_uns)

Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/panel2/"
path <- paste0(Path,"Plot_inf1_alp_1,1_b0,02_g0,5_muc_0,01_sc0,001_muw0,08_sw0,05.png")
ggsave(path,
       plot = plot.inf.1, device = "png")

### Eigenvalue distribution combined ####
col_stab <- color_stab
col_unst_c <- color_unstab

mug <- gammas[1]
tau <- 0
c_stab <- mub*(1-muw) - mug - N*muc
r_stab <- sqrt(mub^2*sw^2 + 2*mub*tau + sc^2)*sqrt(N)
o_stab <- mub - mug + muw*mub*(N-1)

alp <- alp_bet
a <- mub*muw +muc
b <- alp
c <- alp*muw

outl <- (1/2)*(N*a + b + sqrt((N*a)^2 - (2*N-4)*a*b + (4*N-4)*a*c + b^2))
outl2 <- (1/2)*(N*a + b - sqrt((N*a)^2 - (2*N-4)*a*b + (4*N-4)*a*c + b^2))
outl1 <- outl + (mub*(1-muw) - N*muc - mug)
outl2 <- outl2 + (mub*(1-muw) - N*muc - mug)

colors <- c("c" = col_stab, 
            "d" = col_unst_c)

eig_full <- rbind(eigen_unst, eigen_stab)
size_text <- 16

eigen_full <- ggplot(eig_full) + 
  geom_point(aes(re,im), size = 0.2) +
  theme(text = element_text(size = size_text),
        legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  geom_ellipse(aes(x0 = c_stab, y0 = 0, 
                   a = r_stab,
                   b = r_stab, 
                   angle = 0, color = "c")) +
  geom_ellipse(aes(x0 = c_stab, y0 = 0, 
                   a = r_stab,
                   b = r_stab, 
                   angle = 0, color = "d"), linetype = "dashed")  +
  geom_point(aes(o_stab,0), 
             color = col_stab, shape = 21) +
  geom_point(aes(outl1,0), 
             color = col_unst_c, shape = 21) +
  geom_point(aes(outl2,0),
             color = col_unst_c, shape = 21, size = 2.4) +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed")  +
  scale_color_manual(values = colors,
                     name = " ", 
                     labels = c("Equal transmission","Perturbed case")) +
  theme_bw() + coord_fixed() 

eigen_full

plot.inf.stab
plot.inf.1
library("cowplot")
text_tit = 20
grid1 <- plot_grid(plot.inf.stab +
                     xlim(c(0,30)) + ggtitle("b")  + 
                     theme(plot.title = element_text(size = text_tit, face = "bold")),
                   plot.inf.1 + xlim(c(0,30)) + ggtitle("") + 
                     theme(plot.title = element_text(size = text_tit, face = "bold")) + 
                     rremove("ylab") + 
                     scale_y_continuous( breaks=c(0,500,1000,1400)),
                   nrow = 1)
unst_stab <- plot_grid(grid1 ,
                       eigen_full, 
                       ncol = 1)
unst_stab

##### Sum of all the infected by alpha ####
sus_init <- rep(1000, N) # initial susceptibles
inf_init <- rep(10, N)    # initial infecteds

end_time <- 500
alp_bet_vec <- seq(0,2,0.03)
# end_time <- 500
# sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           # COMMUTING,MIGRATION,
           # sus_init,inf_init,end_time)
# df_sum <- data.frame(time = sol[,1])
# for(i in c(1:length(alp_bet_vec))){
#   print(paste0("i: ", i))
#   betas <- rep(mub, N)
#   betas[1] <- betas[1] + alp_bet_vec[i]
#   sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
#              COMMUTING,MIGRATION,
#              sus_init,inf_init,end_time)
#   sol_inf <- sol[,c(1,(N+2):(2*N+1))]
#   sum_inf <- rowSums(sol_inf[,2:(N+1)])
#   df_sum[,ncol(df_sum)+1] <- sum_inf
# }

# save(df_sum,file="/home/marta/Documentos/PHD/2022/RMT_SIR/df_sum.Rda")
load(file="/home/marta/Documentos/PHD/2022/RMT_SIR/df_sum.Rda")

 alp_bet_vec1 <- seq(0,2,0.01)
 # end_time <- 500
 df_bet_out <- data.frame(alp = 0, re = 0)
 for(i in c(1:length(alp_bet_vec1))){
   alp <- alp_bet_vec1[i]
   a <- mub*muw +muc
   b <- alp
   c <- alp*muw
#
   outl <- (1/2)*(N*a + b + sqrt((N*a)^2 - (2*N-4)*a*b + (4*N-4)*a*c + b^2))
   outl1 <- outl + (mub*(1-muw) - N*muc - mug)
   outl2 <- (1/2)*(N*a + b - sqrt((N*a)^2 - (2*N-4)*a*b + (4*N-4)*a*c + b^2))
   outl2 <- outl2 + (mub*(1-muw) - N*muc - mug)

   df_bet_out[nrow(df_bet_out) + 1,] <- c(alp_bet_vec1[i], outl1)
   df_bet_out[nrow(df_bet_out) + 2,] <- c(alp_bet_vec1[i], outl2)
 }

 df_bet_out$im <- 0
 eigen_unst$alp = 0

 plot_alp <- ggplot(eigen_unst, aes(re,im)) +
   geom_point(size = 0.3) +
  geom_point(data = df_bet_out, aes(re,im, colour = alp),
              size = 0.8) +
  scale_colour_gradient(name = TeX("$\\beta^*$"),
                         low = "#F8F053", high = "#4C0EF6",
                        limits=c(0,2),
                        breaks = c(0, 0.5, 1, 1.5,2)) +
   scale_y_continuous( breaks=c(-0.02,0,0.02)) +
  geom_ellipse(aes(x0 = c_stab, y0 = 0,    
                   a = r_stab,
                   b = r_stab,
                  angle = 0))  +
  geom_point(aes(o_stab,0),
             color = color_stab, shape = 8, size = 2.8) +
  geom_point(aes(outl1,0),
             color = color_unstab, shape = 8, size = 2.8) +
 geom_point(aes(outl2,0),
             color = color_unstab, shape = 8, size = 2.8) +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
   theme_bw() +
   theme(text = element_text(size = size_text),
         legend.position = "none", legend.box = "horizontal",
         plot.title = element_text(size = text_tit, face = "bold"))  + 
   ggtitle("c") + rremove("xlab") + rremove("ylab")
 plot_alp

 leg_plot <- get_legend(ggplot(df_bet_out) +
                          geom_point(aes(re,im, colour = alp)) +
                          scale_colour_gradient(name = TeX("$\\beta^*$"),
                                                low = "#F8F053", high = "#4C0EF6",
                                                limits=c(0,2),
                                                breaks = c(0, 0.5, 1, 1.5,2)) + 
                          theme(text = element_text(size = size_text),
                                legend.key.height = unit(0.2, 'cm'),
                                legend.key.width = unit(0.8, 'cm'),
                                legend.key.size = unit(0.1, 'cm'),
                                legend.position = "bottom", legend.box = "horizontal"))
 
 plot_alp1 <- ggdraw() +
   draw_plot(plot_alp) +
   draw_plot(leg_plot, x = 0.65, y = .52, width = .25, height = .25)
 
 plot_alp1

 ##### MAx infected individuals ######
dec <- 0
max_vec <- c()
vec_eq <- c()
for(i in c(1:(ncol(df_sum)-1))){
  max_vec[i] <- max(df_sum[,1+i])
}
for(i in c(1: (ncol(df_sum)-1))){
  # ind <- which(floor(max_inf[i+1]) == floor(df_sum[,i+1]))[1]
  ind <- which(max_vec[i] == df_sum[,i+1])[1]
  if(is.na(which(max_vec[i] == df_sum[,i+1])[1])){
    print("NA found")
    break
  }
  vec_eq[length(vec_eq) + 1] <- df_sum[ind,1]
}

len_vec <- length(alp_bet_vec)
time_max_vec <- vec_eq
df_sum_group <- data.frame(alp <- alp_bet_vec[1:len_vec], 
                           max_inf <- t(max_vec)[1:len_vec],
                           time_max <- time_max_vec[1:len_vec])

colnames(df_sum_group) <-  c("alpha", "Max_inf", "Time_max")
library("latex2exp")
plot_inf_max <- ggplot(df_sum_group) + 
  geom_line(aes(alpha, Max_inf, color= alpha), size = 1.4) + 
  xlab(TeX("$\\beta^*$")) +
  ylab("Max of infected ind.") +
  scale_colour_gradient(name = TeX("$\\beta^*$"),
                        low = "#F8F053", high = "#4C0EF6") +
  # scale_y_continuous( breaks=c(0,1000,2000,2900)) +
  theme_bw() +
  theme(text = element_text(size = size_let), legend.position = "none") 
plot_inf_max

plot_time_max <-  ggplot(df_sum_group) + 
  geom_line(aes(alpha, Time_max), size = 0.5) + 
  geom_point(aes(alpha, Time_max, color= alpha) , size = 1) +
  # scale_y_continuous( breaks=c(0,50,100,150, 190)) +
  xlab(TeX("$\\beta^*$")) +
  ylab("Time to max of infected ind.") +
  scale_colour_gradient(name = TeX("$\\beta^*$"),
                        low = "#F8F053", high = "#4C0EF6") +
  theme_bw() +
  theme(text = element_text(size = size_let), legend.position = "none") 
plot_time_max 
# + xlim(c(0,1.75))
# + xlim(c(0.7,1.7)) + ylim(c(80,110))

seq <- seq(1,len_vec,5 )
colnames(df_sum) <- c("time",as.character(alp_bet_vec))
# filt_sum <- df_sum[, seq]
# df_plot <- reshape2::melt(filt_sum, id.vars = c("time"))
df_plot <- reshape2::melt(df_sum, id.vars = c("time"))

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

sum_inf + xlim(c(0,20)) 

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
                        low = "#F8F053", high = "#4C0EF6", breaks = c(0,0.5,1,1.5,2)) +
  ylab("Sum of infected individuals") +
  xlim( c(0,30) ) +
  theme_bw() +
  theme(text = element_text(size = size_let), legend.position = "right") 

sum_inf

Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/panel2/"
path <- paste0(Path,"panel2.png")
ggsave(path,
       plot = gg_full, device = "png")

unst_stab
plot_inf_max + xlim(c(uni,1))
plot_time_max + xlim(c(uni,1))
sum_inf 
plot_alp

# Create Panel 2:
text_tit = 17
grid1 <- plot_grid(plot.inf.stab +
                     xlim(c(0,30)) + ggtitle("b")  + 
                     theme(plot.title = element_text(size = text_tit, face = "bold")),
                   plot.inf.1 + xlim(c(0,30)) + ggtitle("") + 
                     theme(plot.title = element_text(size = text_tit, face = "bold")) + 
                     rremove("ylab") ,
                   nrow = 1)


plotsum <- plot_grid(sum_inf + theme(legend.position = "none") + ggtitle("a")+ 
                       theme(plot.title = element_text(size = text_tit, face = "bold")),
                     grid1 ,
           ncol = 2, nrow = 1,
          rel_widths = c(0.8,1,1))

plot_sum1 <- plot_grid(plotsum ,
                       plot_alp1,
          nrow = 2,
          rel_heights = c(1.1,0.7))

grid2 <- plot_grid(plot_inf_max + ggtitle("d") + 
                     theme(plot.title = element_text(size = text_tit, face = "bold")),
                    plot_time_max  + ggtitle(" ")+ 
                     theme(plot.title = element_text(size = text_tit, face = "bold")), nrow = 1)
plot_grid(plot_sum1,
          grid2 , nrow = 2, ncol = 1, 
          rel_heights  = c(1.9,1) )
