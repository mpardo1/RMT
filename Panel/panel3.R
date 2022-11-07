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
library("latex2exp")
library("cowplot")
library("ggpubr")
####### GENERATE JACOBIAN ###############################

# number of patches
N <- 50
# epidemiological
#all rates must lie in (0,1) except for betas

Deltas <- rep(0.6, N) # birth rate
# betas <- rgamma(N, shape = (mub/sb)^2, rate = mub/(sb^2))
thetas <- rep(0.4, N) # loss of immunity rates
mud <- 0.6
deaths <- rep(mud, N) # not disease-related death rates
mua <- 0.7
alphas <- rep(mua, N) # recovery rates
mudel <- 0
deltas <- rep(mudel, N) # disease-related death rates
gammas = deaths + alphas + deltas

# mobility
#commuting and migration networks
muw <- 0.6
sw <- 0.3
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

COMMUTING <- rand_mat(N, muw, sw, distrib = "beta")
diag(COMMUTING) <- 0
# COMMUTING <- rand_mat_ell(N, muw, sw, rhow, distrib = "beta")
# COMMUTING[sample.int(N^2, round(p*N^2))] <- 0

MIGRATION <- rand_mat(N, muc, sc, distrib = "beta")
diag(MIGRATION) <- 0
##### Random beta ####
mub <- 0.6
sb <- 0.6
betas <- rgamma(N, shape = (mub/sb)^2, rate = mub/(sb^2)) 

jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))

mub <- mean(betas)

print(plot_eigen(jacobian)) + coord_fixed()
library("ggforce")
plot_eigen_rmt(jacobian,
               N,mub = mean(betas),mug = mud + mua + mudel,
               muw,sw,rhow,Gammaw,
               muc,sc,rhoc,Gammac,
               tau = 0, alp = 0, K = 0)

sus_init <- rep(50, N) # Initial susceptibles
inf_init <- sample(1:100,N)
end_time <- 40
sol.rand <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)

plot_stab.rand <- plot_int(N, sol.rand, state = "INF") +
  theme_bw() + theme(legend.position="none") 
plot_stab.rand + xlim(c(0,10))

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
  df_inf$bet[i] <- betas[as.numeric(substr(df_inf$variable[i],2,3))]
}

library(viridis)

plot_stab.rand  <- ggplot(df_inf) + 
  geom_line(aes(time, value, colour = bet, group = variable))  +
  # geom_line(aes( colour =variable),size=0.5)  +
  ylab("Number of infected individuals")+
  scale_colour_gradient(low = "#F8F053", high = "#4C0EF6",
                        name = ""*beta~" ", 
                        # low = "blue",
                        # # mid = "yellow",
                        # high = "green",
                        # midpoint = 1,
                        breaks=c(0,0.5,1),
                        labels=c(0,0.5,1),
                        limits=c(0,1) )  +  theme_bw()

plot_stab.rand + xlim(c(0,3))

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

print(plot_eigen(jacobian))

sol.mean <- int(N, Deltas,betas_cte,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)

plot_stab.mean <- plot_int(N, sol.mean, state = "INF") +
  theme_bw() +theme(legend.position="none") 
plot_stab.mean


library("ggpubr")
max_time <- 1
text_size <- 15
plot1 <- plot_stab.rand +
  xlim(c(0,max_time))  +
         # scale_y_continuous(breaks=c(0, 200, 400,600,800))  + 
  ylab("Infected infividuals") +
  xlab("Time") +
  theme_bw()+
  theme(text = element_text(size = text_size), legend.position = "none",
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))+
  scale_x_continuous( expand = c(0, 0), limits = c(0,max_time)) +
  scale_y_continuous( expand = c(0, 0))
 
leg_plot <- get_legend(plot_stab.rand  + 
                         theme(text = element_text(size = text_size),
                               legend.key.height = unit(0.2, 'cm'),
                               legend.key.width = unit(0.8, 'cm'),
                               legend.key.size = unit(0.1, 'cm'),
                               legend.position = "bottom", legend.box = "horizontal"))
plot_stab.rand1 <- ggdraw() +
  draw_plot(plot1) +
  draw_plot(leg_plot, x = 0.5, y = .15, width = .25, height = .25)

plot_stab.rand1

vec_col <-  vector(mode="character", length=N)
vec_col[1:N] <- "#A25BC1"
 
plot2 <- plot_stab.mean +
  scale_colour_manual(values = vec_col) +
  ylab("Infected infividuals") +
  scale_x_continuous( expand = c(0, 0), limits = c(0,max_time)) +
  scale_y_continuous( expand = c(0, 0))+
  # scale_y_continuous(breaks=c(0, 200, 400,600,800)) +
  theme(text = element_text(size = text_size),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) 

gg_comp <- plot_grid(plot_stab.rand1 + xlab("Time") ,
                     plot2 + xlab("Time"), nrow = 1,
                     rel_heights = c(0.9,1))
gg_comp

# #######Plot lines eigenvalues#########
# # df_lines <- data.frame(CV = numeric(0),
# #                        eig_cte = numeric(0),
# #                        eig_rand = numeric(0))
# df_lines <- data.frame(mub = numeric(0),
#                        muab = numeric(0),
#                        sb = numeric(0),
#                        CV = numeric(0),
#                        eig_cte = numeric(0),
#                        eig_rand = numeric(0))
# # count = 0
# mub = 0.5
# sb = 0.5
# sb_vec <- seq(0.001,2,0.01)
# mu_vec <- seq(0.001,2,0.01)
# count = 0
# # while(count < 10000){
# for(i in c(1:length(sb_vec))){
#   # sb <- sb_vec[i]
#   # mub <- sb/CV
#   sb <- sb_vec[i]
#   count = 1
#   while(count <= 100){
#     # Rand beta
#     betas <- rgamma(N, shape = (mub/sb)^2, rate = mub/(sb^2)) 
#     jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
#       diag(deaths + alphas + deltas + colSums(MIGRATION))
#     outrand <- max(eigen_mat(jacobian)$re)
#     # Cte beta
#     mean_b <- mean(betas)
#     betas <- rep(mean_b,N)
#     jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
#       diag(deaths + alphas + deltas + colSums(MIGRATION))
#     outcte <- max(eigen_mat(jacobian)$re)
#     
#     df_lines[nrow(df_lines)+1,] <- c(mub,mean_b,sb,sb/mub,outcte,outrand)
#     print(paste0("count:",count))
#     count = count + 1
#   }
# }
# 
# df_linesg <- df_lines %>%  group_by(mub,sb) %>%
#   summarise(muab = mean(muab), CV = min(CV), 
#             eig_cte = mean(eig_cte),
#             eig_rand = mean(eig_rand), n = n())
# 
# df_linesg$diff <- abs(df_linesg$eig_rand - df_linesg$eig_cte)
# 
# gg_CVmu <- ggplot(df_linesg) + 
#   geom_line(aes(CV,diff)) +
#   ggtitle("mu = 0.5") +
#   xlim(c(0,5)) + 
#   theme_bw()
# 
# plot_grid(gg_sig,gg_CV,gg_CVmu)
# 
# gg_sig <- gg_sig + 
#   xlab(TeX("$\\sigma_{\\beta}$")) + 
#   ylab("Error") + ggtitle("")
# 
# df_lines$diff <- abs(df_lines$eig_rand - df_lines$eig_cte)
# 
# ggplot(df_lines) + 
#   geom_line(aes(sb,diff)) +
#   ggtitle("CV = 0.5") +
#   theme_bw()
# 
# df_linesp <- reshape2::melt(df_lines, id.vars = c("CV","sb","mub", "muab"))
# 
# betas <- rep(mub,N)
# jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
#   diag(deaths + alphas + deltas + colSums(MIGRATION))
# outRe <- max(eigen_mat(jacobian)$re)
# 
# plot_lines <- ggplot(df_linesp) +
#   geom_line(aes(CV,value, color = variable)) + 
#   geom_point(aes(CV,value, color = variable)) + 
#   ylab("Right-most eigenvalue") + 
#   xlab(TeX("$C_{V}$")) +
#   geom_hline(yintercept = outRe, linetype = "dashed") +
#   scale_color_manual(name = NULL, labels = c("Constant", "Random"),
#                        values = c("#331491", "#6CB7BA")) +
#   theme_bw()
# plot_lines
# 
# ### READ SVG ######
# library(grImport2)
# library(rsvg)
# library(ggimage)
# library(pdftools)
# Path <- "/home/marta/Documentos/PHD/2022/RMT_SIR/Plots/panel3/metap_model_new.svg"
# ggdraw() +
#   draw_image(Path)
# 
# ggdraw() + draw_image(magick::image_read_pdf(Path, density = 600))
# # rsvg_svg(Path)
# diag <- readPicture(Path)
# # diag <- grid.picture(diag)
# ggplot(diag) + 
#   geom_image()
# 
# ##### Read PNG #####
# library("png")
# 
# Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/panel3/diagram_cte_bet.png"
# diagram_cte <- readPNG(Path)
# margin_ud <- 0.1
# margin_lr <- 1.2
# im_cte <- ggplot() + 
#   background_image(diagram_cte) +
#   # This ensures that the image leaves so me space at the edges
#   theme(plot.margin = margin(t=margin_ud+0.3 , l=margin_lr+0.3, r=margin_lr+0.3,
#                              b=margin_ud+0.3, unit = "cm"))
# 
# Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/panel3/diagram_rand_bet.png"
# Path <- "~/Documentos/PHD/2022/RMT_SIR/Plots/panel3/circ-1.png"
# diagram_rand <- readPNG(Path)
# im_rand <- ggplot() + 
#   background_image(diagram_rand) +
#   # This ensures that the image leaves so me space at the edges
#   theme(plot.margin = margin(t=margin_ud+0.3, l=margin_lr, r=margin_lr, b=margin_ud, unit = "cm"))
# 
# 
# library("ggpubr")
# gg1 <- ggarrange(plot1, plot2, nrow = 1, ncol = 2, common.legend = TRUE)
# gg2 <- ggarrange( im_rand,im_cte, nrow = 1, ncol = 2)
# gg3 <- ggarrange(gg1,gg2,
#                  ncol = 1, nrow = 2,
#                  labels= c("a", "b"),
#                  heights = c(1.3, 1) )
# ggfull <- ggarrange(gg3,err_max_inf, 
#                     ncol = 2, nrow =1,
#                     widths = c(1.5, 1),
#                     labels = c("", "c"))
# 
# 
# path <- paste0(Path,"~/Documents/PHD/2022/RMT_SIR/Plots/panel3/panel3.png")
# ggsave(path,
#        plot = ggfull, device = "png")
# 
# 
# 
#           