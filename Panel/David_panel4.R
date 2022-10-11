##### RANDOM MATRICES FOR METAPOPULATION MODELS #########
### 
### panel 4. works over "RMT_parent.R"
###
rm(list = ls())
source("~/RMT/David/RMT_genrandom.R")
source("~/RMT/David/RMT_plotmobility.R")
source("~/RMT/David/d_functions_eigen_int.R")
library("viridis")
library("latex2exp")
library(cowplot)

##### FUNCTIONS #########################################


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

plot_panelc <- function(N, sol, COMMUTING, IMIN = "N", IMAX = "N",CMMIN = "N", CMMAX = "N",
                        colormean = "green", colorcomm = "#0000FF") {
  
  comm_flows <- data.frame(variable = paste0(rep("S",N),as.character(c(1:N))),
                           incoming = rowSums(COMMUTING)/(N-1),
                           outgoing = colSums(COMMUTING)/(N-1)) %>%
    mutate(both = (incoming+outgoing)/2)
  
  sol_df <- as.data.frame(sol) %>% select(c(1,(N+2):(2*N+1)))
  names(sol_df) <- c("time", paste0(rep("S",N),as.character(c(1:N))))
  sol_df <- mutate(sol_df, Smean = rowSums(sol_df[,-1])/N) %>%
    reshape2::melt(id.vars = c("time")) %>%
    left_join(comm_flows, by = "variable")
  
  IMIN <- ifelse(IMIN == "N", min(sol_df$value), IMIN)
  IMAX <- ifelse(IMAX == "N", max(sol_df$value), IMAX)
  CMMIN <- ifelse(CMMIN == "N", min(sol_df$both), CMMIN)
  CMMAX <- ifelse(CMMAX == "N", max(sol_df$both), CMMAX)
  
  ggplot(filter(sol_df, variable != "Smean"), aes(time, value, group = variable)) + 
    geom_line(aes(colour = both), size = .5) +
    scale_color_gradient(low = "white", high = colorcomm, na.value = "gray90", limits = c(CMMIN,CMMAX)) +
    #geom_line(data = filter(sol_df, variable == "Smean"), aes(y = value), color = colormean, size = 1.4) +
    ylab("Infected individuals") +
    ylim(c(IMIN,IMAX)) +
    theme(legend.position = "none") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(color = "none") +
    #theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.ticks = element_blank())
}

stragAB <- function(k, nu){
  (N/2)*(mub*muw+muc) + ((k-1)/2)*(mub*nu) + mub*(1-muw) - mug - N*muc + 
    (1/2)*sqrt(N^2*(mub*muw + muc)^2 + ((k-1)*mub*nu)^2 + 2*(N+N*k-4*k)*(mub*muw + muc)*mub*nu)
} 

stragC <- function(k, nu){
  N*(mub*muw+muc)/2 + (k/2-1)*mub*nu/2 +
    sqrt((N*(mub*muw+muc))^2+mub*nu*(2*mub*muw+muc)*(N+3*N*k/2-2*(k/2+1))+
           ((mub*nu)^2)*(2*k*(N-k/2)+(k/2-1)^2))/2 +
    mub*(1-muw) - mug - N*muc
} 

stragDEF <- function(k, muwstar){
  mub*(muw*(N-k)+muwstar*k/N)*(N-1)/N + mub - mug
} 


stragstrucC <- function(N,km) {
  
  N*(mub*muw+muc)/2 + (km-1)*mub*nu/2 +
    sqrt((N*(mub*muw+muc))^2+2*mub*nu*(mub*muw+muc)*(N+km*(3*N-2*km-2))+
           ((mub*nu)^2)*(4*km*(N-km)+(km-1)^2))/2
  
}

##### PARAMETER VALUES ##################################

N <- 50

sus_init <- rep(50000, N) # initial susceptibles
inf_init <- runif(N, min = 50,100)  # initial infecteds
tolinf <- 5
end_time <- 25
end_time_rand <- 20

col_unstab <- "#0D0C18"
colA <- "#F69A79"
colB <- "#F26430"
colC <- "#D2430F"
colD <- "#85DCFF"
colE <- "#0ABAFF"
colF <- "#0084B8"

Deltas <- rep(0.1, N) # birth rate
mub <- 0.1
sb <- 0.001
betas <- rep(mub, N) # transmission rates
#betas <- rgamma(N, shape = (mub/sb)^2, rate = mub/(sb^2))
thetas <- rep(0.1, N) # loss of immunity rates
mud <- 0.1
deaths <- rep(mud, N) # not disease-related death rates
mua <- 0.45
alphas <- rep(mua, N) # recovery rates
mudel <- 0
deltas <- rep(mudel, N) # disease-related death rates
gammas = deaths + alphas + deltas

muw <- 0.12
sw <- 0.07/3
rhow <- 0 #original rho (Gamma of baron et al)

muc <- 0.0001
sc <- 0.00001
rhoc <- .001

mub*muw*(N-1)+mub-mua-mud-mudel

COMMUTING <- rand_mat_ell(N, muw, sw, rhow, distrib = "beta")
MIGRATION <- rand_mat_ell(N, muc, sc, rhoc, distrib = "beta")
diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)

sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)

##### SCENARIOS #########################################

# UNSTABLE

plot_stab_mob <- plotmobility(COMMUTING, color1 = col_unstab)
plot_stab_int <- plot_panelc(N, sol, COMMUTING, colormean = col_unstab, colorcomm = col_unstab)
plot_grid(plot_stab_mob + ggtitle("Unstable Scenario"),
          plot_stab_int , 
          nrow  =2 ,
          vjust = 1.1)

nodes <- 4
muwstar <- 0
COMMUTINGA <- COMMUTINGB <- COMMUTINGC <- COMMUTINGD <- COMMUTINGE <- COMMUTINGF <- COMMUTING

# STRATEGY A
fs <- sample(c(1:N), nodes)
COMMUTINGA[fs,] <- muwstar
solA <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTINGA,MIGRATION,
           sus_init,inf_init,end_time)

# STRATEGY B
cs <- sample(c(1:N), nodes)
COMMUTINGB[,fs] <- muwstar
solB <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
            COMMUTINGB,MIGRATION,
            sus_init,inf_init,end_time)

# STRATEGY C
fcs <- sample(c(1:N), nodes/2)
fcs <- fs[c(1:(nodes/2))]
COMMUTINGC[fcs,] <- COMMUTINGC[,fcs] <- muwstar
solC <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
            COMMUTINGC,MIGRATION,
            sus_init,inf_init,end_time)

# STRATEGY D
inds <- sample(c(1:N^2), nodes*N)
COMMUTINGD[inds] <- muwstar
solD <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
            COMMUTINGD,MIGRATION,
            sus_init,inf_init,end_time)

# STRATEGY E
rinds <- inds[c(1:(nodes*N/2))]
COMMUTINGE[rinds] <- muwstar
COMMUTINGE <- t(COMMUTINGE)
COMMUTINGE[rinds] <- muwstar
COMMUTINGE <- t(COMMUTINGE)
solE <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
            COMMUTINGE,MIGRATION,
            sus_init,inf_init,end_time)

# STRATEGY F
COMMUTINGF <- COMMUTINGF*(1-nodes/N)
solF <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
            COMMUTINGF,MIGRATION,
            sus_init,inf_init,end_time)

##### MOB NETWORKS ######################################

CMAXDET <- max(c(COMMUTINGA,COMMUTINGB,COMMUTINGC))
CMAXRAND <- max(c(COMMUTINGD,COMMUTINGE, COMMUTINGF))
CMAXRAND <- CMAXDET

IMAX <- max(unlist(c(select(as.data.frame(solA), c((N+2):(2*N+1))),
            select(as.data.frame(solB), c((N+2):(2*N+1))),
            select(as.data.frame(solC), c((N+2):(2*N+1))),
            select(as.data.frame(solD), c((N+2):(2*N+1))),
            select(as.data.frame(solE), c((N+2):(2*N+1))),
            select(as.data.frame(solF), c((N+2):(2*N+1)))))) + tolinf
IMIN <- min(unlist(c(select(as.data.frame(solA), c((N+2):(2*N+1))),
              select(as.data.frame(solB), c((N+2):(2*N+1))),
              select(as.data.frame(solC), c((N+2):(2*N+1))),
              select(as.data.frame(solD), c((N+2):(2*N+1))),
              select(as.data.frame(solE), c((N+2):(2*N+1))),
              select(as.data.frame(solF), c((N+2):(2*N+1)))))) - tolinf
CMMAX <- max(c(rowSums(COMMUTINGA)/(N-1),colSums(COMMUTINGA)/(N-1),
               rowSums(COMMUTINGB)/(N-1),colSums(COMMUTINGB)/(N-1),
               rowSums(COMMUTINGC)/(N-1),colSums(COMMUTINGC)/(N-1),
               rowSums(COMMUTINGD)/(N-1),colSums(COMMUTINGD)/(N-1),
               rowSums(COMMUTINGE)/(N-1),colSums(COMMUTINGE)/(N-1),
               rowSums(COMMUTINGF)/(N-1),colSums(COMMUTINGF)/(N-1)))
CMMIN <- min(c(rowSums(COMMUTINGA)/(N-1),colSums(COMMUTINGA)/(N-1),
               rowSums(COMMUTINGB)/(N-1),colSums(COMMUTINGB)/(N-1),
               rowSums(COMMUTINGC)/(N-1),colSums(COMMUTINGC)/(N-1),
               rowSums(COMMUTINGD)/(N-1),colSums(COMMUTINGD)/(N-1),
               rowSums(COMMUTINGE)/(N-1),colSums(COMMUTINGE)/(N-1),
               rowSums(COMMUTINGF)/(N-1),colSums(COMMUTINGF)/(N-1)))              

mobstab <- plotmobility(COMMUTING, cmax = CMAXDET, color1 = col_unstab)
mobA <- plotmobility(COMMUTINGA , cmax = CMAXDET, color1 = colA)
mobB <- plotmobility(COMMUTINGB, cmax = CMAXDET, color1 = colB)
mobC <- plotmobility(COMMUTINGC, cmax = CMAXDET, color1 = colC)
mobD <- plotmobility(COMMUTINGD, cmax = CMAXRAND, color1 = colD)
mobE <- plotmobility(COMMUTINGE, cmax = CMAXRAND, color1 = colE)
mobF <- plotmobility(COMMUTINGF, cmax = CMAXRAND, color1 =colF)

##### INFECTEDS #########################################

intu <- plot_panelc(N, sol, COMMUTING, colormean = col_unstab, colorcomm = col_unstab)
intA <- plot_panelc(N, solA, COMMUTINGA, colormean = col_stab, colorcomm = colA,
                    IMIN = IMIN, IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX) + ylim(c(0,100))
intB <- plot_panelc(N, solB, COMMUTINGB, colormean = col_stab, colorcomm = colB,
                    IMIN = IMIN, IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX) + ylim(c(0,100))
intC <- plot_panelc(N, solC, COMMUTINGC, colormean = col_stab, colorcomm = colC,
                    IMIN = IMIN, IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX)  + ylim(c(0,100))
intD <- plot_panelc(N, solD, COMMUTINGD, colormean = col_stab, colorcomm = colD,
                    IMIN = IMIN, IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX) + ylim(c(0,100))
intE <- plot_panelc(N, solE, COMMUTINGE, colormean = col_stab, colorcomm = colE,
                    IMIN = IMIN, IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX) + ylim(c(0,100))
intF <- plot_panelc(N, solF, COMMUTINGF, colormean = col_stab, colorcomm = colF,
                    IMIN = IMIN, IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX)  + ylim(c(0,100))

##### COMPARISON ########################################

eig_comp <- data.frame(k = integer(0),
                       rmt = numeric(0),
                       lrpAB = numeric(0),
                       lrpC = numeric(0),
                       AA = numeric(0),
                       BB = numeric(0),
                       CC = numeric(0),
                       DD = numeric(0),
                       EE = numeric(0),
                       FF = numeric(0))

nu <- muwstar - muw
mug <- mud + mua + mudel
max_nodes <- floor(N/4)
for (nodes in seq(0,max_nodes,by = 1)) {
  
  COMMUTINGA <- COMMUTINGB <- COMMUTINGC <- COMMUTINGD <- COMMUTINGE <- COMMUTINGF <- COMMUTING
  fs <- sample(c(1:N), nodes)
  COMMUTINGA[fs,] <- muwstar
  COMMUTINGB[,fs] <- muwstar
  inds <- sample(c(1:N^2), nodes*N)
  COMMUTINGD[inds] <- muwstar
  rinds <- inds[c(1:(nodes*N/2))]
  COMMUTINGE[rinds] <- muwstar
  COMMUTINGE <- t(COMMUTINGE)
  COMMUTINGE[rinds] <- muwstar
  COMMUTINGE <- t(COMMUTINGE)
  COMMUTINGF <- COMMUTINGF*(1-nodes/N)
  
  jacA <- (COMMUTINGA + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  jacB <- (COMMUTINGB + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  jacD <- (COMMUTINGD + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  jacE <- (COMMUTINGE + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  jacF <- (COMMUTINGF + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  
  k <- nodes
  
  if ((k %% 2) == 0) {
    # k even
    fcs <- sample(c(1:N), nodes/2)
    fcs <- fs[c(1:(nodes/2))]
    COMMUTINGC[fcs,] <- COMMUTINGC[,fcs] <- muwstar
    jacC <- (COMMUTINGC + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
    km <- k/2
    eig_comp[k+1,] <- c(k,
                        mub*(muw*(N-k)+muwstar*k/N)*(N-1)/N + mub - mug,
                        N*(mub*muw+muc)/2 + (k-1)*mub*nu/2 +
                          sqrt((N*(mub*muw+muc))^2+((k-1)*mub*nu)^2+2*(mub*nu)*(mub*muw+muc)*(N+N*k-4*k))/2  +
                          mub*(1-muw) - mug - N*muc,
                        #N*(mub*muw+muc)/2 + (km-1)*mub*nu/2 +
                        #  sqrt((N*(mub*muw+muc))^2+2*mub*nu*(mub*muw+muc)*(N+km*(3*N-2*km-2))+
                        #         ((mub*nu)^2)*(4*km*(N-km)+(km-1)^2))/2 +
                        #  mub*(1-muw) - mug - N*muc,
                        stragstrucC(N,km)+ mub*(1-muw) - mug - N*muc,
                        max(eigen_mat(jacA)$re),
                        max(eigen_mat(jacB)$re),
                        max(eigen_mat(jacC)$re),
                        max(eigen_mat(jacD)$re),
                        max(eigen_mat(jacE)$re),
                        max(eigen_mat(jacF)$re))
  } else {
    # k odd
    eig_comp[k+1,] <- c(k,
                        mub*(muw*(N-k)+muwstar*k/N)*(N-1)/N + mub - mug,
                        N*(mub*muw+muc)/2 + (k-1)*mub*nu/2 +
                          sqrt((N*(mub*muw+muc))^2+((k-1)*mub*nu)^2+2*(mub*nu)*(mub*muw+muc)*(N+N*k-4*k))/2  +
                          mub*(1-muw) - mug - N*muc,
                        stragstrucC(N,k/2)+ mub*(1-muw) - mug - N*muc,
                        max(eigen_mat(jacA)$re),
                        max(eigen_mat(jacB)$re),
                        NA,
                        max(eigen_mat(jacD)$re),
                        max(eigen_mat(jacE)$re),
                        max(eigen_mat(jacF)$re))
  }
}

eig_comp_real <- select(eig_comp,k,AA,BB,CC,DD,EE,FF) %>%
  pivot_longer(cols = (2:ncol(.)), names_to = "prediction", values_to = "eigenvalue")
eig_comp_pred <- select(eig_comp, k, rmt, lrpAB, lrpC) %>%
  pivot_longer(cols = (2:ncol(.)), names_to = "prediction", values_to = "eigenvalue")

plot_lines <-  eig_comp_pred %>% ggplot(aes(x = k, y = eigenvalue, color = prediction)) +
  #geom_line(size = 0.8) + 
  geom_line() +
  geom_point(data = eig_comp_real, aes(x = k, y = eigenvalue, color = prediction)) +
  xlab("Number of affected nodes") +
  ylab("Epidemiological threshold") +
  scale_color_manual(values=c(colA, colB, colC,
                              colD, colE, colF,
                              colA, colC, colD)) + 
  scale_x_continuous(breaks = c(0,3,6,9), minor_breaks = seq(0,10,1)) +
  theme_bw() +
  theme(legend.position = "none") +
  #theme(text = element_text(size = 12)) +
  annotate("rect", xmin = 3.5, xmax = 4.5,
           ymin = -0.035, ymax = 0,
           fill = "yellow", alpha = .2) +
  annotate("text", x = 6, y = 0.005, label = "In panel") +
  theme(aspect.ratio = 2)
plot_lines

##### PANEL 4 ###########################################

relh <- 1.2

gu <- plot_grid(mobstab,
                intu + theme(aspect.ratio = 1),
                nrow = 2,
                align = "v",
                labels = NULL,
                #labels = c("Unstable scenario", NULL),
                rel_heights = c(1,relh))
###PPT###
plot_grid(mobstab,
          intu + 
            rremove("ylab") + rremove("xlab") + 
            theme(text = element_text(size = size_text),aspect.ratio = 1),
          nrow = 2,
          align = "v",
          labels = NULL,
          #labels = c("Scenario A", NULL),
          rel_heights = c(1,relh))

########

ga <- plot_grid(mobA,
                intA + theme(aspect.ratio = 1),
                nrow = 2,
                align = "v",
                labels = NULL,
                #labels = c("Scenario A", NULL),
                rel_heights = c(1,relh))

###PPT###
plot_grid(mobA,
          intA + 
            rremove("ylab") + rremove("xlab") + 
            theme(text = element_text(size = size_text),aspect.ratio = 1),
          nrow = 2,
          align = "v",
          labels = NULL,
          #labels = c("Scenario A", NULL),
          rel_heights = c(1,relh))

########

gb <- plot_grid(mobB,
                intB + theme(aspect.ratio = 1) + ylab(""),
                nrow = 2,
                align = "v",
                labels = NULL,
                #labels = c("Scenario B", NULL),
                rel_heights = c(1,relh))

###PPT###
plot_grid(mobB,
          intB + 
            rremove("ylab") + rremove("xlab") + 
            theme(text = element_text(size = size_text),aspect.ratio = 1),
          nrow = 2,
          align = "v",
          labels = NULL,
          #labels = c("Scenario A", NULL),
          rel_heights = c(1,relh))

########

gc <- plot_grid(mobC,
                intC + theme(aspect.ratio = 1) + ylab(""),
                nrow = 2,
                align = "v",
                labels = NULL,
                #labels = c("Scenario C", NULL),
                rel_heights = c(1,relh))

###PPT###
plot_grid(mobC,
          intC + 
            rremove("ylab") + rremove("xlab") + 
            theme(text = element_text(size = size_text),aspect.ratio = 1),
          nrow = 2,
          align = "v",
          labels = NULL,
          #labels = c("Scenario A", NULL),
          rel_heights = c(1,relh))

########

gd <- plot_grid(mobD,
                intD + theme(aspect.ratio = 1),
                nrow = 2,
                align = "v",
                labels = NULL,
                #labels = c("Scenario D", NULL),
                rel_heights = c(1,relh))

ge <- plot_grid(mobE,
                intE + theme(aspect.ratio = 1) + ylab(""),
                nrow = 2,
                align = "v",
                labels = NULL,
                #labels = c("Scenario E", NULL),
                rel_heights = c(1,relh))

gf <- plot_grid(mobF,
                intF + theme(aspect.ratio = 1) + ylab(""),
                nrow = 2,
                align = "v",
                labels = NULL,
                #labels = c("Scenario F", NULL),
                rel_heights = c(1,relh))

panel4 <- plot_grid(gu, ga, gb, gc,
          plot_lines, gd, ge, gf,
          nrow = 2,
          align = "h",
          #axis = "r",
          #labels = NULL,
          labels = c("Unstable scenario",
                     "Strategy A", "Strategy B", "Strategy C",
                     "Performance", "Strategy D", "Strategy E", "Strategy F"),
          label_y = 1.05,
          hjust = 0,
          label_x = 0.15,
          rel_widths = c(1,1,1,1)) +
  theme(aspect.ratio = 1) +
  theme(plot.margin = margin(20,20,20,20))
panel4

ggsave(filename = "panel4.jpg", panel4, width = 12, height = 12, dpi = 450)

# si no en el save, meter el aspect ratio deseado (~1)

######### CORRELATED ########
####### COMMUTING: DIRECTED CORRELATED PERTURBATIONS ####

nu <- -muw
mug <- mud + mua + mudel

eig_comp <- data.frame(k = integer(0),
                       random = numeric(0),
                       correlated = numeric(0),
                       real_r = numeric(0),
                       real_c = numeric(0))


for (k in seq(2,(N^2)/10,by = 2)) {
  
  COMM_R <- COMM_C <- COMMUTING
  COMM_R[sample(c(1:N^2),k)] <- muw + nu
  
  ind <- vector()
  for (K in c(2:N)) {
    ind <- c(ind, seq(((K-2)*N+K),(K-1)*N))
  }
  cor_ind <- sample(ind, k/2)
  COMM_C[cor_ind] <- muw + nu
  COMM_C <- t(COMM_C)
  COMM_C[cor_ind] <- muw + nu
  COMM_C <- t(COMM_C)
  
  jacobian_r <- (COMM_R + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  jacobian_c <- (COMM_C + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  
  eig_comp[k,] <- c(k,
                    mub*(N-1)*(muw+nu*k/(N^2-N))+mub-mug + 1/(N*mug),
                    mub*(N-1)*(muw+nu*k/(N^2-N))+mub-mug,
                    max(eigen_mat(jacobian_r)$re),
                    max(eigen_mat(jacobian_c)$re))
}

eig_comp <- filter(eig_comp, !is.na(k))
eig_comp %>%
  pivot_longer(cols = (2:ncol(.)), names_to = "prediction", values_to = "eigenvalue") %>% 
  ggplot(aes(x = k, y = eigenvalue, color = prediction)) +
  geom_line(size = 1)
#scale_color_manual(values=c("#0000FF","#FF0000","#0000FF","#FF0000","#63B159","#63B159"))

####### CORRELATIONS: MONOTONICITY IN GAMMA #############

N <- 500

muw <- .05
sw <- .01
rhow <- .1
rw <- cw <- .7
COMMUTING <- rand_mat_ell(N, muw, sw, rhow, distrib = "beta")

muc <- 0.001
sc <- 0.0005
rhoc <- 0
Gammac <- 0
MIGRATION <- rand_mat_ell(N, muc, sc, rhoc, distrib = "beta")

outlier_df <- data.frame("Gammaw" = vector(), "realGammaw" = vector(), "bgpred" = vector(), "rmtpred" = vector(), "real" = vector())
for (Gammaw in seq(0.01,0.99,by = .005)) {
  
  if ((Gammaw/sqrt(rw*cw) < 1) & ((N*rhow-2*Gammaw)/(N-(rw+cw)) < 1)) {
    
    gCOMMUTING <- rand_mat_cor_beta(N, muw*N, sw*N, rhow, Gammaw, rw, cw)  
    
    COMMUTING <- gCOMMUTING[[1]]
    diag(COMMUTING) <- rep(0,N)
    
    diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)
    sum(colSums(MIGRATION) > 1)
    
    jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
    
    sigma <- sqrt(mub^2*sw^2 + sc^2)
    rho <- (mub^2*rhow*sw^2 + rhoc*sc^2)/(sigma^2)
    Gamma <- (mub^2*Gammaw*(sw^2)/N + Gammac*(sc^2)/N)/(sigma^2)
    
    outlier_df <- rbind(outlier_df, c(Gammaw,
                                      gCOMMUTING[[2]][4,2],
                                      mub - mud - mua - mudel + mub*muw*(N-1) +
                                        (1/2)*(N*mub*muw+N*muc)*(1+rho/Gamma)*(sqrt(1+(4*Gamma*sigma^2)/((mub*muw + muc)^2))-1),
                                      mub - mud - mua - mudel + mub*muw*(N-1),
                                      max(eigen_mat(jacobian)$re)))
    print(Gammaw)
  }
}

names(outlier_df) <- c("Gammaw", "realGammaw", "bgpred", "rmtpred", "real")

outlier_df <- pivot_longer(outlier_df, cols = (2:ncol(outlier_df)) , names_to = "prediction", values_to = "outlier")

filter(outlier_df, prediction %in% c("realGammaw")) %>%
  ggplot(aes(x = Gammaw, y = outlier, color = prediction)) +
  geom_line(size = 1)

filter(outlier_df, prediction %in% c("bgpred", "rmtpred", "real")) %>%
  ggplot(aes(x = Gammaw, y = outlier, color = prediction)) +
  geom_line(size = 1)

outlier_df <- pivot_wider(outlier_df, names_from = prediction, values_from = outlier) %>%
  mutate(bgpredreal = mub - mud - mua - mudel + mub*muw*(N-1) +
           (1/2)*(N*mub*muw+N*muc)*(1+rho/realGammaw)*(sqrt(1+(4*realGammaw*sigma^2)/((mub*muw + muc)^2))-1))

outlier_df %>%
  ggplot(aes(x = Gammaw, y = N*realGammaw)) +
  geom_point() +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  geom_line(aes(x = Gammaw, y = rmtpred)) +
  geom_line(aes(x = Gammaw, y = bgpred))

filter(outlier_df, prediction %in% c("bgpred", "rmtpred", "real", "bgpredreal")) %>%
  ggplot(aes(x = Gammaw, y = outlier, color = prediction)) +
  geom_line(size = 1)

####### CORRELATIONS: INSTABILITY #######################

N <- 200

mub <- 0.02
betas <- rep(mub, N) # transmission rates
mud <- 0.1282
deaths <- rep(mud, N)
mua <- 0.25
alphas <- rep(mua, N)
mudel <- 0
deltas <- rep(mudel, N)

# re-scaling only mu and sigma: RMT and BG are indistinguishable
# muw <- 0.09
# sw <- 0.06
# rhow <- 0.45
# Gammaw <- 0.3
# rw <- 0.4
# cw <- 0.4

# re-scaling all parameters: ellipse is lost, outlier is still accurate
muw <- 0.09
sw <- 0.06
rhow <- 0.45
Gammaw <- 0.3*N
rw <- 0.4*N
cw <- 0.4*N

(Gammaw/sqrt(rw*cw) < 1) & ((N*rhow-2*Gammaw)/(N-(rw+cw)) < 1)

muc <- muw/100
sc <- sw/120
rhoc <- rhow/100
Gammac <- Gammaw/100
rc <- rw/100
cc <- cw/100

(Gammac/sqrt(rc*cc) < 1) & ((N*rhoc-2*Gammac)/(N-(rc+cc)) < 1)

COMMUTING <- rand_mat_cor_beta(N, N*muw, sqrt(N)*sqrt(N)*sw, rhow, Gammaw, rw, cw)
MIGRATION <- rand_mat_cor_norm(N, N*muc, sqrt(N)*sqrt(N)*sc, rhoc, Gammac, rc, cc)
# COM[[2]]
# MIG[[2]]
# COMMUTING <- COM[[1]]
# MIGRATION <- MIG[[1]]
#MIGRATION <- matrix(rep(0,N^2), nrow = N)
#muc <- sc <- rhoc <- Gammac <- 0
diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)
jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
eigmayor <- eigen_mat(jacobian)
redmayor <- plotmobility(COMMUTING[1:80,1:80])

COMMUTING <- rand_mat_cor_norm(N, N*muw, sqrt(N)*sqrt(N)*sw, -rhow, -Gammaw, rw, cw)
MIGRATION <- rand_mat_cor_norm(N, N*muc, sqrt(N)*sqrt(N)*sc, -rhoc, -Gammac, rc, cc)
diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)
jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
eigmenor <- eigen_mat(jacobian)
redmenor <- plotmobility(COMMUTING[1:80,1:80])

COMMUTING <- rand_mat_cor_norm(N, N*muw, sqrt(N)*sqrt(N)*sw, rhow, 0.001, 0.001, 0.001)
MIGRATION <- rand_mat_ell(N, muc, sc, rhoc, distrib = "beta")
diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)
jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
eigRMT <- eigen_mat(jacobian)
reduncor <- plotmobility(COMMUTING[1:80,1:80])

sig <- sqrt((mub*sw)^2 + sc^2)
rho <- (rhow*((mub*sw)^2) + rhoc*sc^2)/sig^2
Gamma <- (Gammaw*((mub*sw)^2) + Gammac*sc^2)/sig^2

outlierREmayor <- max(eigmayor$re)
outlierREmenor <- max(eigmenor$re)
outlierRMT <- (N-1)*mub*muw + mub - mud - mua - mudel
outlierBGmayor <- (N-1)*mub*muw + mub - mud - mua - mudel +
  (N/2)*(mub*muw+muc)*(1+rho/Gamma)*(sqrt(1+(4*Gamma*sig^2/(N*(mub*muw+muc)^2)))-1)
outlierBGmenor <- (N-1)*mub*muw + mub - mud - mua - mudel +
  (N/2)*(mub*muw+muc)*(1+rho/Gamma)*(sqrt(1-(4*Gamma*sig^2/(N*(mub*muw+muc)^2)))-1)

center <- mub*(1-muw)-mua-mud-mudel-N*muc
radius <- sqrt(N)*sig

color0 <- "black"
color1 <- "red"
color2 <- "blue"
color3 <- "green"
color4 <- "purple"

BGvsRMT <- ggplot(data = eigmenor) + 
  geom_ellipse(aes(x0 = center, y0 = 0, a = (1+rho)*radius,
                   b = (1-rho)*radius, angle = 0), color = color1) +
  geom_vline(xintercept = 0, color = color2, linetype = "dashed") +
  geom_point(aes(x = outlierRMT, y = 0), color = color1) +
  geom_point(aes(x = outlierBGmenor, y = 0), color = color3) +
  geom_point(aes(x = outlierBGmayor, y = 0), color = color4) +
  geom_point(aes(re,im), size = 0.05, color = color3, alpha = .5) +
  geom_point(data = eigmayor, aes(re,im), size = 0.05, color = color4, alpha = .5) +
  geom_point(data = eigRMT, aes(re,im), size = 0.05, color = color0, alpha = .5) +#+ xlim(c(-0.7,-0.6))
  labs(x="Real part", y="Imaginary part") +
  theme_bw()

ggsave(file="BGvsRMT.svg", plot=BGvsRMT, width=10, height=8)

# prueba <- data.frame(param = numeric(0), eigen_dist = numeric(0))
# ind <- 0
# for (mub in seq(0,1,by = .05)) {
#   if ((Gammaw/sqrt(rw*cw) < 1) & ((N*rhow-2*Gammaw)/(N-(rw+cw)) < 1)) {
#   
#   ind <- ind+1
#   sig <- sqrt((mub*sw)^2 + sc^2)
#   rho <- (rhow*((mub*sw)^2) + rhoc*sc^2)/sig^2
#   Gamma <- (Gammaw*((mub*sw)^2) + Gammac*sc^2)/sig^2
#   
#   outlierRMT <- (N-1)*mub*muw + mub - mud - mua - mudel
#   outlierBG <- (N-1)*mub*muw + mub - mud - mua - mudel +
#     (N/2)*(mub*muw+muc)*(1+rho/Gamma)*(sqrt(1+(4*Gamma*sig^2/(N*(mub*muw+muc)^2)))-1)
#   
#   prueba[ind,] <- c(mub,abs(outlierRMT-outlierBG))
#   }
# }
# ggplot(data = prueba, aes(x = param, y = eigen_dist)) + geom_point() #+ coord_fixed()

####### CORRELATIONS: MOBILITY NETWORKS #################
col_neg <- "#2A9D8F"
col_pos <- "#7B287D"
col_non <- "#09071A"
N <- 50

# positive rhow vs negative rhow
rhow <- .85
rhow <- -.85
sus_init <- rep(100, N) # initial susceptibles
inf_init <- rep(10, N)    # initial infecteds

end_time <- 100
COMMUTING <- rand_mat_ell(N, muw, sw, rhow, distrib = "beta")
diag(COMMUTING) <- rep(0,N)
plotmobility(COMMUTING, color1 = col_non)

sol.stab <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
                COMMUTING,MIGRATION,
                sus_init,inf_init,end_time)

plot_stab <- plot_int1(N, sol.stab, state = "INF") +
  theme_bw() +theme(legend.position="none")
plot_stab

vec_col <-  vector(mode="character", length=N)
vec_col[1:N] <- col_non

plot.inf.stab <- plot_int1(N, sol.stab, state = "INF") +
  scale_colour_manual(values = vec_col) +
  theme_bw() +
  theme(text = element_text(size = 20),
        legend.position="none") + rremove("xlab") + rremove("ylab")

plot.inf.stab + xlim(c(0,20))

# high rw/cw (exchange values)
muw <- 0.2
sw <- 0.06
rhow <- 0 #original rho (Gamma of baron et al)
Gammaw <- .15 #gamma of baron et al
rw <- .7*N
cw <- .1*N
(Gammaw/sqrt(rw*cw) < 1) & ((N*rhow-2*Gammaw)/(N-(rw+cw)) < 1)

COMMUTING <- rand_mat(N, muw, sw, distrib = "beta")
COMMUTING <- rand_mat_cor_beta(N, muw*N, sw*N, rhow, Gammaw, rw, cw)
diag(COMMUTING) <- rep(0,N)
plotmobility(COMMUTING,color1 = col_pos)

sus_init <- rep(100, N) # initial susceptibles
inf_init <- rep(10, N)    # initial infecteds

end_time <- 100
sol.stab <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
                COMMUTING,MIGRATION,
                sus_init,inf_init,end_time)

plot_stab <- plot_int1(N, sol.stab, state = "INF") +
  theme_bw() +theme(legend.position="none")
plot_stab

vec_col <-  vector(mode="character", length=N)
vec_col[1:N] <- col_pos

plot.inf.stab <- plot_int1(N, sol.stab, state = "INF") +
  scale_colour_manual(values = vec_col) +
  theme_bw() +
  theme(text = element_text(size = 20),
        legend.position="none") + rremove("xlab") + rremove("ylab")

plot.inf.stab + xlim(c(0,20))

# positive vs negative Gammaw
muw <- 0.25
sw <- 0.1
rhow <- 0 #original rho (Gamma of baron et al)
Gammaw <- .22*N #gamma of baron et al
Gammaw <- -.22*N #gamma of baron et al
rw <- .25*N
cw <- .25*N
(Gammaw/sqrt(rw*cw) < 1) & ((N*rhow-2*Gammaw)/(N-(rw+cw)) < 1)

COMMUTING <- rand_mat(N, muw, sw, distrib = "beta")
COMMUTING <- rand_mat_cor_beta(N, muw*N, sw*N, rhow, Gammaw, rw, cw)
diag(COMMUTING) <- rep(0,N)
plotmobility(COMMUTING, color1 = col_neg)

sol.stab <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
                COMMUTING,MIGRATION,
                sus_init,inf_init,end_time)

plot_stab <- plot_int1(N, sol.stab, state = "INF") +
  theme_bw() +theme(legend.position="none")
plot_stab

vec_col <-  vector(mode="character", length=N)
vec_col[1:N] <- col_neg

plot.inf.stab <- plot_int1(N, sol.stab, state = "INF") +
  scale_colour_manual(values = vec_col) +
  theme_bw() +
  theme(text = element_text(size = size_let),
        legend.position="none") + rremove("xlab") + rremove("ylab")

plot.inf.stab + xlim(c(0,20))

# "high" Gammaw and rhow

muw <- 0.25
sw <- 0.1
rhow <- 0.85 #original rho (Gamma of baron et al)
Gammaw <- .4*N #gamma of baron et al
rw <- .45*N #positive
cw <- .45*N #positive
(abs(Gammaw)/sqrt(rw*cw) < 1) & (rw+cw<N) & (abs(N*rhow-2*Gammaw)/(N-(rw+cw)) < 1)

# negative values generate a network with the same properties

muw <- 0.25
sw <- 0.1
rhow <- -0.85 #original rho (Gamma of baron et al)
Gammaw <- -.4*N #gamma of baron et al
rw <- .45*N #positive
cw <- .45*N #positive
(abs(Gammaw)/sqrt(rw*cw) < 1) & (rw+cw<N) & (abs(N*rhow-2*Gammaw)/(N-(rw+cw)) < 1)

COMMUTING <- rand_mat(N, muw, sw, distrib = "beta")
COMMUTING <- rand_mat_cor_beta(N, muw*N, sw*N, rhow, Gammaw, rw, cw)
diag(COMMUTING) <- rep(0,N)
plotmobility(COMMUTING,color1 = col_neg)
end_time <- 200
sol.stab <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
                COMMUTING,MIGRATION,
                sus_init,inf_init,end_time)

plot_stab <- plot_int1(N, sol.stab, state = "INF") +
  theme_bw() +theme(legend.position="none")
plot_stab

vec_col <-  vector(mode="character", length=N)
vec_col[1:N] <- col_neg

plot.inf.stab <- plot_int1(N, sol.stab, state = "INF") +
  scale_colour_manual(values = vec_col) +
  theme_bw() +
  theme(text = element_text(size = size_let),
        legend.position="none") + rremove("xlab") + rremove("ylab")

plot.inf.stab + xlim(c(0,100))

# preventing disease

muw <- 0.25
sw <- 0.1
rhow <- 0.35 #original rho (Gamma of baron et al)
Gammaw <- -.1*N #gamma of baron et al
rw <- .15*N #positive
cw <- .15*N #positive
(abs(Gammaw)/sqrt(rw*cw) < 1) & (rw+cw<N) & (abs(N*rhow-2*Gammaw)/(N-(rw+cw)) < 1)

COMMUTING <- rand_mat(N, muw, sw, distrib = "beta")
COMMUTING <- rand_mat_cor_beta(N, muw*N, sw*N, rhow, Gammaw, rw, cw)
diag(COMMUTING) <- rep(0,N)
plotmobility(COMMUTING, color1 =  col_neg)
sol.stab <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
                COMMUTING,MIGRATION,
                sus_init,inf_init,end_time)

plot_stab <- plot_int1(N, sol.stab, state = "INF") +
  theme_bw() +theme(legend.position="none")
plot_stab

vec_col <-  vector(mode="character", length=N)
vec_col[1:N] <- col_neg

plot.inf.stab <- plot_int1(N, sol.stab, state = "INF") +
  scale_colour_manual(values = vec_col) +
  theme_bw() +
  theme(text = element_text(size = size_let),
        legend.position="none") + rremove("xlab") + rremove("ylab")

plot.inf.stab + xlim(c(0,100))
####### CORRELATIONS: PANEL #################

N <- 200

mub <- 0.02
betas <- rep(mub, N) # transmission rates
mud <- 0.1282
deaths <- rep(mud, N)
mua <- 0.25
alphas <- rep(mua, N)
mudel <- 0
deltas <- rep(mudel, N)

# re-scaling all parameters: ellipse is lost, outlier is still accurate
muw <- 0.09
sw <- 0.06
rhow <- 0.45
Gammaw <- 0.3*N
rw <- 0.4*N
cw <- 0.4*N

muc <- muw/100
sc <- sw/120
rhoc <- rhow/100
Gammac <- Gammaw/100
rc <- rw/100
cc <- cw/100

COMMUTING <- rand_mat_cor_beta(N, N*muw, sqrt(N)*sqrt(N)*sw, rhow, Gammaw, rw, cw)[[1]]
MIGRATION <- rand_mat_cor_beta(N, N*muc, sqrt(N)*sqrt(N)*sc, rhoc, Gammac, rc, cc)[[1]]
diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)
jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
eigmayor <- eigen_mat(jacobian)
redmayor <- plotmobility(COMMUTING[1:60,1:60])

COMMUTING <- rand_mat_cor_beta(N, N*muw, sqrt(N)*sqrt(N)*sw, -rhow, -Gammaw, rw, cw)[[1]]
MIGRATION <- rand_mat_cor_beta(N, N*muc, sqrt(N)*sqrt(N)*sc, -rhoc, -Gammac, rc, cc)[[1]]
diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)
jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
eigmenor <- eigen_mat(jacobian)
redmenor <- plotmobility(COMMUTING[1:60,1:60])

COMMUTING <- rand_mat_cor_beta(N, N*muw, sqrt(N)*sqrt(N)*sw, rhow, 0.001, 0.001, 0.001)[[1]]
MIGRATION <- rand_mat_cor_beta(N, N*muc, sqrt(N)*sqrt(N)*sc, rhoc, 0.00001, 0.00001, 0.00001)[[1]]
diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)
jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
eigRMT <- eigen_mat(jacobian)
reduncor <- plotmobility(COMMUTING[1:60,1:60])

sig <- sqrt((mub*sw)^2 + sc^2)
rho <- (rhow*((mub*sw)^2) + rhoc*sc^2)/sig^2
Gamma <- (Gammaw*((mub*sw)^2) + Gammac*sc^2)/sig^2

outlierREmayor <- max(eigmayor$re)
outlierREmenor <- max(eigmenor$re)
outlierRMT <- (N-1)*mub*muw + mub - mud - mua - mudel
outlierBGmayor <- (N-1)*mub*muw + mub - mud - mua - mudel +
  (N/2)*(mub*muw+muc)*(1+rho/Gamma)*(sqrt(1+(4*Gamma*sig^2/(N*(mub*muw+muc)^2)))-1)
outlierBGmenor <- (N-1)*mub*muw + mub - mud - mua - mudel +
  (N/2)*(mub*muw+muc)*(1+rho/Gamma)*(sqrt(1-(4*Gamma*sig^2/(N*(mub*muw+muc)^2)))-1)

center <- mub*(1-muw)-mua-mud-mudel-N*muc
radius <- sqrt(N)*sig

color0 <- "black"
color1 <- "red"
color2 <- "blue"
color3 <- "green"
color4 <- "purple"

BGvsRMT <- ggplot(data = eigmenor) + 
  geom_ellipse(aes(x0 = center, y0 = 0, a = (1+rho)*radius,
                   b = (1-rho)*radius, angle = 0), color = color2) +
  geom_vline(xintercept = 0, color = color2, linetype = "dashed") +
  geom_point(aes(x = outlierRMT, y = 0), size = 3, color = color1, shape = 15) +
  geom_point(aes(x = outlierBGmenor, y = 0), size = 3, color = color3, shape = 16) +
  geom_point(aes(x = outlierBGmayor, y = 0), size = 3, color = color4, shape = 17) +
  geom_point(aes(re,im), size = 1.7, color = color3, shape = 16) +  #, alpha = .4
  geom_point(data = eigmayor, aes(re,im), size = 1.7, color = color4, shape = 17) +
  geom_point(data = eigRMT, aes(re,im), size = 1.7, color = color1, shape = 15) +#+ xlim(c(-0.7,-0.6))
  labs(x="Real part", y="Imaginary part") +
  theme_bw()
BGvsRMT
#ggsave(file="BGvsRMT.svg", plot=BGvsRMT)

#ggsave(file="redmenor.svg", plot=redmenor)
#ggsave(file="redmayor.svg", plot=redmayor)
#ggsave(file="reduncor.svg", plot=reduncor)

corpanel <-ggarrange(ggarrange(redmenor, reduncor, redmayor, labels = c("A", "B", "C"),
                               ncol = 3, nrow = 1),
                     BGvsRMT,
                     ncol = 1, nrow = 2)
corpanel

ggsave(file="corpanel.svg", plot=corpanel)
ggsave(corpanel, file="corpanel.eps", device = cairo_ps, fallback_resolution = 600)

####### CORRELATIONS: PANEL ROUND 2 #########

library(cowplot)

N <- 200

mub <- 0.02
betas <- rep(mub, N) # transmission rates
mud <- 0.1282
deaths <- rep(mud, N)
mua <- 0.25
alphas <- rep(mua, N)
mudel <- 0
deltas <- rep(mudel, N)

Deltas <- rep(0.1, N) # birth rate
thetas <- rep(0.1, N) # loss of immunity rates

# re-scaling all parameters: ellipse is lost, outlier is still accurate
muw <- 0.09
sw <- 0.06
rhow <- 0.45
Gammaw <- 0.3*N
rw <- 0.4*N
cw <- 0.4*N

muc <- muw/1000
sc <- sw/1200
rhoc <- rhow/1000
Gammac <- Gammaw/1000
rc <- rw/1000
cc <- cw/1000

(Gammaw/sqrt(rw*cw) < 1) & ((N*rhow-2*Gammaw)/(N-(rw+cw)) < 1)

MIGRATION <- rand_mat_cor_beta(N, N*muc, sqrt(N)*sqrt(N)*sc, rhoc, 0.001, 0.001, 0.001)[[1]]
MIGRATION <- matrix(rep(0, N^2),nrow = N)
muc <- sc <- rhoc <- Gammac <- 0

COMMUTINGUNC <- rand_mat_cor_beta(N, N*muw, sqrt(N)*sqrt(N)*sw, rhow, 0.001, 0.001, 0.001)[[1]]
jacobianunc <- (COMMUTINGUNC + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
eigRMT <- eigen_mat(jacobianunc)

COMMUTINGMEN <- rand_mat_cor_beta(N, N*muw, sqrt(N)*sqrt(N)*sw, -rhow, -Gammaw, rw, cw)[[1]]
jacobianmen <- (COMMUTINGMEN + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
eigmenor <- eigen_mat(jacobianmen)

COMMUTINGMAY <- rand_mat_cor_beta(N, N*muw, sqrt(N)*sqrt(N)*sw, rhow, Gammaw, rw, cw)[[1]]
#MIGRATION <- rand_mat_cor_beta(N, N*muc, sqrt(N)*sqrt(N)*sc, rhoc, Gammac, rc, cc)[[1]]
jacobianmay <- (COMMUTINGMAY + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
eigmayor <- eigen_mat(jacobianmay)

titlesize <- 22
margin <- 20
netsize <- 50
CMAX <- max(c(max(COMMUTINGMAY),max(COMMUTINGMEN),max(COMMUTINGUNC)))
color0 <- "black"
color1 <- "red"
color2 <- "blue"
color3 <- "#009528" #menor
color4 <- "#FF0000" #mayor
color5 <- "#0000FF" #uncor
rescale <- 2

redmayor <- plotmobility(rescale*COMMUTINGMAY[1:netsize,1:netsize], cmax = CMAX, color1 = color4) +
  ggtitle(expression("C: "~Gamma>0)) +
  #theme(plot.title = element_text(size=titlesize)) +
  theme(plot.title = element_text(margin=margin(0,0,5,0))) +
  theme(plot.margin = margin(margin,margin,margin,margin, "pt"))

redmenor <- plotmobility(rescale*COMMUTINGMEN[1:netsize,1:netsize], cmax = CMAX, color1 = color3) +
  ggtitle(expression("A: "~Gamma<0)) +
  #theme(plot.title = element_text(size=titlesize)) +
  theme(plot.title = element_text(margin=margin(0,0,5,0))) +
  theme(plot.margin = margin(margin,margin,margin,margin, "pt"))

reduncor <- plotmobility(rescale*COMMUTINGUNC[1:netsize,1:netsize], cmax = CMAX, color1 = color5) +
  ggtitle(expression("B: "~Gamma==0)) +
  #theme(plot.title = element_text(size=titlesize)) +
  theme(plot.title = element_text(margin=margin(0,0,5,0))) +
  theme(plot.margin = margin(margin,margin,margin,margin, "pt"))

sig <- sqrt((mub*sw)^2 + sc^2)
rho <- (rhow*((mub*sw)^2) + rhoc*sc^2)/sig^2
Gamma <- (Gammaw*((mub*sw)^2) + Gammac*sc^2)/sig^2

outlierREmayor <- max(eigmayor$re)
outlierREmenor <- max(eigmenor$re)
outlierRMT <- (N-1)*mub*muw + mub - mud - mua - mudel
outlierBGmayor <- (N-1)*mub*muw + mub - mud - mua - mudel +
  (N/2)*(mub*muw+muc)*(1+rho/Gamma)*(sqrt(1+(4*Gamma*sig^2/(N*(mub*muw+muc)^2)))-1)
outlierBGmenor <- (N-1)*mub*muw + mub - mud - mua - mudel +
  (N/2)*(mub*muw+muc)*(1+rho/Gamma)*(sqrt(1-(4*Gamma*sig^2/(N*(mub*muw+muc)^2)))-1)

center <- mub*(1-muw)-mua-mud-mudel-N*muc
radius <- sqrt(N)*sig

just <- 0.7

BGvsRMT <- ggplot(data = eigmenor) + 
  #geom_ellipse(aes(x0 = center, y0 = 0, a = (1+rho)*radius,
  #                 b = (1-rho)*radius, angle = 0), color = color1) +
  geom_vline(xintercept = 0, color = color2, linetype = "dashed") +
  geom_point(aes(x = outlierRMT, y = 0), size = 3, color = color5, shape = 15) +
  geom_text(aes(x = outlierRMT, y = 0, label = "B", vjust = just, hjust = -just)) +
  geom_text(aes(x = outlierBGmenor, y = 0, label = "A", vjust = just, hjust = -just)) +
  geom_text(aes(x = outlierBGmayor, y = 0, label = "C", vjust = just, hjust = -just)) +
  geom_point(aes(x = outlierBGmenor, y = 0), size = 3, color = color3, shape = 16) +
  geom_point(aes(x = outlierBGmayor, y = 0), size = 3, color = color4, shape = 17) +
  geom_point(aes(re,im), size = 1.7, color = color3, shape = 16, alpha = .6) +
  geom_point(data = eigmayor, aes(re,im), size = 1.7, color = color4, shape = 17, alpha = .6) +
  geom_point(data = eigRMT, aes(re,im), size = 1.7, color = color5, shape = 15, alpha = .6) +#+ xlim(c(-0.7,-0.6))
  labs(x="Real part", y="Imaginary part") +
  theme_bw()
BGvsRMT
#ggsave(file="BGvsRMT.svg", plot=BGvsRMT)

plot_grid(redmenor, reduncor, redmayor, nrow = 1, labels = NULL)

#ggsave(file="redmenor.svg", plot=redmenor)
#ggsave(file="redmayor.svg", plot=redmayor)
#ggsave(file="reduncor.svg", plot=reduncor)

corpanel <-ggarrange(ggarrange(redmenor, reduncor, redmayor,
                               ncol = 3, nrow = 1),
                     BGvsRMT,
                     ncol = 1, nrow = 2)
corpanel

ggsave(file="corpanel.svg", plot=corpanel)
ggsave(file="corpanel.png", plot=corpanel)
#ggsave(corpanel, file="corpanel.eps", device = cairo_ps, fallback_resolution = 600)

sus_init <- rep(100000, N) # initial susceptibles
inf_init <- rep(100, N)    # initial infecteds

end_time <- 10

# integro el sistema con condiciones iniciales 
sol_unc <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
               COMMUTINGUNC,MIGRATION,
               sus_init,inf_init,end_time) 
sol_may <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
               COMMUTINGMAY,MIGRATION,
               sus_init,inf_init,end_time)
sol_men <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
               COMMUTINGMEN,MIGRATION,
               sus_init,inf_init,end_time)

inf_unc <- sol_unc[,c(1,(N+2):(2*N+1))]%>%
  as.data.frame() %>%
  reshape2::melt(id.vars = c("time"))
inf_may <- sol_may[,c(1,(N+2):(2*N+1))]%>%
  as.data.frame() %>%
  reshape2::melt(id.vars = c("time"))
inf_men <- sol_men[,c(1,(N+2):(2*N+1))]%>%
  as.data.frame() %>%
  reshape2::melt(id.vars = c("time"))

mininf <- min(c(inf_unc$value,inf_may$value,inf_men$value))
maxinf <- max(c(inf_unc$value,inf_may$value,inf_men$value))

tlinea <- .15

unc_sol_plot <- inf_unc %>% ggplot(aes(time, value)) + 
  geom_line(colour = color5, size = tlinea)  +
  ylab("") + 
  ylim(c(mininf,maxinf)) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(aspect.ratio=1)

may_sol_plot <- inf_may %>% ggplot(aes(time, value)) + 
  geom_line(colour = color4, size = tlinea)  +
  ylab("") + 
  ylim(c(mininf,maxinf)) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(aspect.ratio=1)

men_sol_plot <- inf_men %>% ggplot(aes(time, value)) + 
  geom_line(colour = color3, size = tlinea)  +
  ylab("Infected individuals") + 
  ylim(c(mininf,maxinf)) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(aspect.ratio=1)

#plot_grid(men_sol_plot, unc_sol_plot, may_sol_plot, nrow = 1, labels = NULL)

plot_grid(
  plot_grid(
    plot_grid(redmenor, reduncor, redmayor, nrow = 1, labels = NULL),
    plot_grid(men_sol_plot, unc_sol_plot, may_sol_plot, nrow = 1, labels = NULL),
    nrow = 2,
    labels = NULL,
    rel_heights = c(1.3,1),
    align = "v"),
  plot_grid(NULL, BGvsRMT,
            nrow = 1,
            labels = NULL,
            align = "h",
            rel_widths = c(1,2)),
  nrow = 2,
  labels = NULL,
  align = "h",
  rel_widths = c(1,1),
  rel_heights = c(2,1))
ggsave(file = "panel.svg")

plot_grid(redmenor, reduncor, redmayor, nrow = 1, labels = NULL)
ggsave(file = "networks.svg")
plot_grid(men_sol_plot, unc_sol_plot, may_sol_plot, nrow = 1, labels = NULL)
ggsave(file = "infecteds.svg")
BGvsRMT + theme(aspect.ratio=1/2)
ggsave(file = "BGvsRMT.svg")