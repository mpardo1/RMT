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
library("latex2exp")
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

size_text <- 16

MIGRATION <- rand_mat(N, muc, sc, distrib = "beta")
diag(MIGRATION) <- 0
# ----------------------------------------------------------------------#
##### Plot Stability #####

beta_vec <- seq(0, 0.5, 0.001)
stab_func <- function(bet){
 (mug - bet)/(bet*(N-1)) 
}
line_thresh <- sapply(beta_vec, stab_func)
df_thres <- data.frame(bet = beta_vec, muc = line_thresh)
df_thres[which(df_thres$muc>0.75), 2] = 0.75


# Data points for the values beta and muc for integration and eigenvalues plots:
stab_par = 0.1
betnew <- 0.35
muwnew <- ((betnew*stab_par*(N-1)) + betnew - stab_par)/(stab_par*(N-1))

betnew - mug + betnew*stab_par*(N-1) == stab_par - mug + stab_par*muwnew*(N-1)

data_points <- data.frame(mub = c(stab_par,stab_par,betnew), muc = c(stab_par,muwnew,stab_par))
df_stab <- data.frame(mub = c(stab_par), muc = c(stab_par))
df_unst_c <- data.frame(mub = c(stab_par), muc = c(muwnew))
df_unst_b <- data.frame(mub = c(betnew), muc = c(stab_par))

# col_stab <- "#F3A712"
# col_unst_c <- "#129490"
# col_unst_b <- "#7A306C"
col_stab <- "#F2C911"
col_unst_c <- "#22A884FF"
col_unst_b <- "#440154FF"

# col_stab <- "#414487FF"
# col_unst_c <- "#7AD151FF"
# col_unst_b <- "#440154FF"


col_stab_r = "#484C4E"
col_unstab_r = "#BCBBC9"

plot_area <- ggplot(df_thres, aes(bet, muc)) + 
  geom_ribbon(aes(x = bet, ymin = muc, ymax = 0.75, fill = col_unstab_r), alpha= 0.7) + 
  geom_ribbon(aes(x = bet, ymin = 0, ymax = muc, fill = col_stab_r), alpha= 0.7) + 
  scale_fill_manual(values = c(col_stab_r,col_unstab_r), name = NULL,
                    labels = c("Stable", "Unstable"), position = "bottom") + 
  xlab(TeX("$\\beta$")) + ylab(TeX("$\\mu_c$")) +
  geom_point(data = df_stab, aes(mub, muc),
             colour= col_stab, size = 2.7) +
  geom_point(data = df_unst_c, aes(mub, muc),
             colour= col_unst_c, size = 2.7) +
  geom_point(data = df_unst_b, aes(mub, muc),
             colour= col_unst_b, size = 2.7)  +
  scale_x_continuous(breaks=c(0,0.25,0.50), limits = c(0, 0.5),
                     labels = c("0", "0.25", "0.50")) +
  scale_y_continuous(breaks=c(0,0.25,0.50, 0.75), limits = c(0, 0.75),
                     labels = c("0", "0.25", "0.50", "0.75")) +
  theme_bw() +
  theme(text = element_text(size = 15),legend.position = c(0.7, 0.85)) 
  

# Save plot
# Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/Gen/"
path <- paste0(Path,"Area_g0,95_muc_0,01_sc0,00001_sw0,05.png")
ggsave(path,
       plot = plot_area, device = "png")
#-----------------------------------------------------------------------------#
######### PLOTS ##############
#------------------------------------------------------------------------
# col_stab <- "#C6C013"
# col_unst_c <- "#EF8A17"
# col_unst_b <- "#EF2917"
#### Plots RMT and Integration ####
sus_init <- rep(50, N) # initial susceptibles
inf_init <- rep(10, N)    # initial infecteds

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
sol_stab_f <- reshape2::melt(as.data.frame(sol_stab_f[,c(1,52:101)]), id.vars = c("time"))
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

# Plot integrations together
sol_uns_bet$type <- "uns_bet"
sol_uns_com$type <- "uns_com"
sol_stab_f$type <- "stab"

sol_uns_bet$variable1 <- paste0("UB",sol_uns_bet$variable)
sol_uns_com$variable1 <- paste0("UC",sol_uns_bet$variable)
sol_stab_f$variable1 <- paste0("S",sol_uns_bet$variable)

sol_tot <- rbind(sol_uns_bet,sol_uns_com,sol_stab_f)

plot_integration  <- ggplot(sol_tot,
                            aes(time, value), show.legend = NA) + 
  geom_line(aes( group =variable1, colour = type), size=0.5)  +
  ylab("Number of infected individuals") +
  scale_color_manual(values = c(col_stab,col_unst_c,col_unst_b),
                     name = NULL, labels = c("Stable        ", 
                                           TeX("Unstable $\\mu_c$"), 
                                           TeX("Unstable $\\beta$ "))) +
  xlim(c(0,15)) + 
  theme_bw() +
  theme(text = element_text(size = size_text),
        legend.position = c(0.75, 0.85),
        legend.text.align = 0)

###### ALL EIGENVALUE DIST ####

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
  theme_bw() + 
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
let_size <- 15
##### Construct panel1 ####
plot1 <- plot_grid(plot_area  + 
                     # ggtitle("a") +
            theme(plot.title = element_text(size = let_size), 
                  plot.margin = unit(c(0.2, 0.4, 0.2, 0.2), "cm")) +
              scale_x_continuous( expand = c(0, 0)) +
              scale_y_continuous( expand = c(0, 0)),
          plot_integration +
            # ggtitle("c") +
            scale_x_continuous(breaks=c(0,2,4,6,8,10,12),
                               limits = c(0,12), expand = c(0, 0)) +
            theme(plot.title = element_text(size = let_size), 
                  plot.margin = unit(c(0.2, 0.4, 0.2, 0.2), "cm")) +
            scale_y_continuous( expand = c(0, 0)) +
            ylab("Infected Individuals"),
          rel_widths = c(1,1), nrow = 1, labels = c("a","c"),
          label_fontfamily = "Helvetica",
          label_size = size_text)

plot_full <- plot_grid(plot1,
                       eigen_full +
                         theme(plot.title = element_text(size = let_size),
                               legend.position = "none", 
                               plot.margin = unit(c(1.3, 0.2, 1.2, 0.2), "cm")) +
                         rremove("xlab") + rremove("ylab") +
                         scale_y_continuous(breaks=c(-0.1,0,0.1), expand = c(0, 0.01))  +
                         scale_x_continuous( expand = c(0, 0.07)) ,
          ncol = 1, rel_heights = c(2.2,1), axis = "v", scale = c(1,1),
          labels = c("","b"),
          label_fontfamily = "Helvetica",
          label_size = size_text)
plot_full


