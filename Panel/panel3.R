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
##### Random beta ####
mub <- 0.5
sb <- 0.8
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

sus_init <- rep(10000, N) # Initial susceptibles
inf_init <- rep(100, N)
end_time <- 20
sol.rand <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)

plot_stab.rand <- plot_int(N, sol.rand, state = "INF") +
  theme_bw() + theme(legend.position="none") 
plot_stab.rand

for(i in c(1:N)){
  colnames(sol.rand)[i+1] <-  paste0("S",i)
  colnames(sol.rand)[N+i+1] <-  paste0("I",i)
  colnames(sol.rand)[2*N+i+1] <-  paste0("R",i)
}

sol.rand <- as.data.frame(sol.rand)
df_plot <- reshape2::melt(sol.rand, id.vars = c("time"))
df_plot$type <- substr(df_plot$variable,1,1)
df_inf <- df_plot  %>% 
  filter( substr(df_plot$variable,1,1) == "I")
df_inf$bet <- 0
for(i in c(1:nrow(df_inf))){
  df_inf$bet[i] <- betas[as.numeric(substr(df_plot$variable[i],2,3))]
}

plot_stab.rand  <- ggplot(df_inf) + 
  geom_line(aes(time, value, colour = bet, group = variable))  +
  # geom_line(aes( colour =variable),size=0.5)  +
  ylab("Number of infected individuals")+
  scale_colour_gradient(name = ""*beta~" ", 
                        low = "blue", high = "red",
                        breaks=c(0,1,4),
                        labels=c(0,1,4),
                        limits=c(0,5))

plot_stab.rand

###### Mean(rand(betas)) ######
betas_cte <- rep(mean(betas),N)
jacobian <- (COMMUTING + diag(N)) %*% diag(betas_cte) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))

mub <- mean(betas_cte)
plot_eigen_rmt(jacobian,
               N,mub,mug = mud + mua + mudel,
               muw,sw,rhow,Gammaw,
               muc,sc,rhoc,Gammac,
               tau = 0, alp = 0, K = 0) +
  scale_y_continuous( breaks=c(0)) 

eig <- eigen_mat(jacobian)
max(eig$re)

# print(plot_eigen(jacobian))

sol.mean <- int(N, Deltas,betas_cte,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)

plot_stab.mean <- plot_int(N, sol.mean, state = "INF") +
  theme_bw() +theme(legend.position="none") 
plot_stab.mean


library("ggpubr")
plot1 <- plot_stab.rand +
  xlim(c(0,10))  +
         scale_y_continuous(breaks=c(0, 2000, 4000,6000,8000))  + 
  ylab("Infected infividuals") +
  theme_bw()+
  theme(text = element_text(size = 15)) 
 
vec_col <-  vector(mode="character", length=N)
vec_col[1:N] <- "#7018D5"
 
plot2 <- plot_stab.mean +
  xlim(c(0,10))  + 
  scale_colour_manual(values = vec_col) +
  ylab("Infected infividuals") +
  scale_y_continuous(breaks=c(0, 2000, 4000,6000,8000)) +
  theme(text = element_text(size = 15)) 

gg_comp <- ggarrange(plot1,plot2,
                     common.legend = TRUE, legend = "bottom",
                     labels = c("a", "b"))
gg_comp
##### Numerical comparison rand and mean(bet) right most eigen #####
# Read from a CSV 
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
  theme_bw() +
  theme(text = element_text(size = 15)) 
err_max_inf

Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/panel3/diagrams.svg"
svg_image <- SVGMapping::loadSVG()
