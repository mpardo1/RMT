####### RANDOM MATRICES FOR METAPOPULATION MODELS #######
### 
### article plots
###
### plots for the article
### works over "RMT_parent.R"
###
rm(list = ls())
source("~/RMT/David/RMT_genrandom.R")
source("~/RMT/David/RMT_plotmobility.R")
source("~/RMT/David/d_functions_eigen_int.R")
library("viridis")

####### COMMUTING: PANEL ################################
library("copula")
library("copula")
library("viridis")
library("ggsci")
library("grid")
library("reshape")
library("ggpubr")
library("cowplot")
library("ggforce")

######### FUNCTION FOR PLOTS ######################

if( getwd() == "/home/marta"){
  Path <- "~/Documentos/PHD/2022/RMT_SIR/Plots/panel4/"
}else{
  Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/panel4/"
}

plot_panelc <- function(N, sol, COMMUTING, colormean, color_low = "#FFFFFF", color_high, IMIN = "N", IMAX = "N",CMMIN = "N", CMMAX = "N",
                        size_text = 15) {
  
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
    scale_color_gradient(low = color_low, high = color_high,
                         na.value = "white", limits = c(CMMIN,CMMAX)) +
    # scale_fill_gradient(low="grey90", high="purple") +
    geom_line(data = filter(sol_df, variable == "Smean"), aes(y = value), 
              color = colormean, size = 0.8, linetype = "dashed") +
    ylab("Infected individuals") +
    ylim(c(IMIN,IMAX)) +
    guides(color = "none") +
    theme_bw() + 
    theme(text = element_text(size = 12), legend.position = "none",
          plot.margin = unit(c(0,0,0,0),"cm")) 
  
}


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

### NETWORK
plotmobilityMPA <- function(mat, high_col, low_col = "#FFFFFF"){
  diag(mat) <- NA
  longData <- melt(mat)
  colnames(longData) <- c("X1", "X2", "value")
  ggplot(longData, aes(x = X2, y = X1)) + 
    geom_raster(aes(fill=value)) + 
    # scale_color_brewer(palette = "Dark2")+
    scale_fill_gradient(low=low_col, high=high_col, na.value = "white") +
    theme_bw() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position = "none") + 
    theme_void()  +
    coord_fixed() +
    ylim(c(N,0)) 
}

# Eigenvalues:

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

#----------------------------------------------------------------------------#

# N <- 100
N <- 40
# para los plots
sus_init <- rep(100000, N) # initial susceptibles
inf_init <- runif(N, min = 50,100)  # initial infecteds
end_time <- 200
end_time_rand <- 20

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

muw <- 0.1
sw <- 0.07/3
rhow <- 0 #original rho (Gamma of baron et al)

muc <- 0.0001
sc <- 0.00001
rhoc <- .001

mub*muw*(N-1)+mub-mua-mud-mudel

col_unstab = "#0D0C18"   # negro
col_stabA = "#06D622"  # verde
col_stabB = "#E96F1D"  # naranja
col_stabC = "#F8DC22"  # amarillo
col_stabD = "#B62F2F"  # rojo
col_stabE = "#2C5CE1"  # azul
col_stabF = "#AA2FB5"  # violeta
col_circ = "#010C0C" # negro

color_stab <- "#4464AD"
color_unstab <- "#A4B0F5"

COMMUTING <- rand_mat_ell(N, muw, sw, rhow, distrib = "beta")
MIGRATION <- rand_mat_ell(N, muc, sc, rhoc, distrib = "beta")
# MIGRATION <- matrix(0,N,N)
diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)
sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)

col_low <- "#90F0CE"
col_high <- "#3F83E9"
col_mean <- "#000205"
plot_stab_mob <- plotmobility(COMMUTING, col_low, col_high) + 
  theme(legend.position = "left") 
ggsave(file=paste0(Path,"unstnet.svg"))
plot_stab_int <- plot_panelc(N, sol, COMMUTING,col_mean, col_low, col_high )
ggsave(file=paste0(Path,"unstsol.svg"))

plot_grid(plot_stab_mob + ggtitle("Unstable Scenario"),
          plot_stab_int, 
          nrow  =2 ,
          vjust = 1.1)

ggsave(file=paste0(Path,"unstable.svg"))


#### Eigenvalues
mug <- gammas[1]
tau <- 0
c_unstab <- mub*(1-muw) - mug - N*muc
r_unstab <- sqrt(mub^2*sw^2 + 2*mub*tau + sc^2)*sqrt(N)
o_unstab <- mub - mug + muw*mub*(N-1)

jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))
plot_eigen(jacobian) + 
  geom_ellipse(aes(x0 = c_unstab, y0 = 0, 
                   a = r_unstab,
                   b = r_unstab, 
                   angle = 0), color = col_circ) + geom_point(aes(o_unstab,0), 
                                                              color = col_unstab, shape = 21) +
  geom_vline(xintercept = 0) +
  theme_bw()

### CONTROL
nodes <- 4
muwstar <- 1
MIGRATIONA <- MIGRATIONB <- MIGRATIONC <- MIGRATIOND <- MIGRATIONE <- MIGRATIONF <- MIGRATION

color_stab <- "#4464AD"
color_unstab <- "#8AEA92"

# STRATEGY A
fs <- sample(c(1:N), nodes)
MIGRATIONA[fs,] <- muwstar
solA <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
            COMMUTING,MIGRATIONA,
            sus_init,inf_init,end_time)

jacobianA <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATIONA -
  diag(deaths + alphas + deltas + colSums(MIGRATIONA))

vec_colA <-  vector(mode="character", length=N)
vec_colA[1:N] <- color_unstab
vec_colA[fs] <- color_stab

plot.infA <- plot_int1(N, solA, state = "INF") +
  scale_colour_viridis_d() +
  # scale_colour_manual(values = vec_colA) +
  theme_bw() +
  theme(text = element_text(size = 15),
        legend.position="none") + xlim(c(0,20))
plot.infA

# STRATEGY B
# cs <- sample(c(1:N), nodes)
MIGRATIONB[,fs] <- muwstar
solB <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
            COMMUTING,MIGRATIONB,
            sus_init,inf_init,end_time)

jacobianB <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATIONB -
  diag(deaths + alphas + deltas + colSums(MIGRATIONB))

vec_colB <-  vector(mode="character", length=N)
vec_colB[1:N] <- color_unstab
vec_colB[fs] <- color_stab

plot.infB <- plot_int1(N, solB, state = "INF") +
  scale_colour_manual(values = vec_colB) +
  theme_bw() +
  theme(text = element_text(size = 15),
        legend.position="none") + xlim(c(0,20))
plot.infB

# STRATEGY C
fcs <- sample(c(1:N), nodes/2)
fcs <- fs[c(1:(nodes/2))]
MIGRATIONC[fcs,] <- MIGRATIONC[,fcs] <- muwstar
solC <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
            COMMUTING,MIGRATIONC,
            sus_init,inf_init,end_time)

jacobianC <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATIONC -
  diag(deaths + alphas + deltas + colSums(MIGRATIONC))

vec_colC <-  vector(mode="character", length=N)
vec_colC[1:N] <- col_stabA
vec_colC[fcs] <- color_stab

plot.infC <- plot_int1(N, solC, state = "INF") +
  scale_colour_manual(values = vec_colC) +
  theme_bw() + theme(text = element_text(size = 15),
        legend.position="none") + xlim(c(0,20))
plot.infC

plot_grid(plot.infA,plot.infB,plot.infC)
# RAND STRATEGIES #
end_time <- 200

# STRATEGY D
inds <- sample(c(1:N^2), nodes*N)
MIGRATIOND[inds] <- muwstar

solD <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
            COMMUTING,MIGRATIOND,
            sus_init,inf_init,end_time)

jacobianD <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATIOND -
  diag(deaths + alphas + deltas + colSums(MIGRATIOND))

# STRATEGY E
rinds <- inds[c(1:(nodes*N/2))]
MIGRATIONE[rinds] <- muwstar
MIGRATIONE <- t(MIGRATIONE)
MIGRATIONE[rinds] <- muwstar
MIGRATIONE <- t(MIGRATIONE)

solE <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
            COMMUTING,MIGRATIONE,
            sus_init,inf_init,end_time)

jacobianE <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATIONE -
  diag(deaths + alphas + deltas + colSums(MIGRATIONE))

# STRATEGY F
MIGRATIONF <- MIGRATIONF*(1-nodes/N)
solF <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
            COMMUTING,MIGRATIONF,
            sus_init,inf_init,end_time)

jacobianF <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATIONF -
  diag(deaths + alphas + deltas + colSums(MIGRATIONF))


#### MOB HEAT MAP########

###########################################################


CMAXDET <- max(c(MIGRATIONA,MIGRATIONB,MIGRATIONC))
CMAXRAND <- max(c(MIGRATIOND,MIGRATIONE, MIGRATIONF))

IMAX <- max(unlist(c(select(as.data.frame(solA), c((N+2):(2*N+1))),
                     select(as.data.frame(solB), c((N+2):(2*N+1))),
                     select(as.data.frame(solC), c((N+2):(2*N+1))),
                     select(as.data.frame(solD), c((N+2):(2*N+1))),
                     select(as.data.frame(solE), c((N+2):(2*N+1))),
                     select(as.data.frame(solF), c((N+2):(2*N+1)))))) + 50
IMIN <- min(unlist(c(select(as.data.frame(solA), c((N+2):(2*N+1))),
                     select(as.data.frame(solB), c((N+2):(2*N+1))),
                     select(as.data.frame(solC), c((N+2):(2*N+1))),
                     select(as.data.frame(solD), c((N+2):(2*N+1))),
                     select(as.data.frame(solE), c((N+2):(2*N+1))),
                     select(as.data.frame(solF), c((N+2):(2*N+1)))))) - 50
CMMAX <- max(c(rowSums(MIGRATIONA)/(N-1),colSums(MIGRATIONA)/(N-1),
               rowSums(MIGRATIONB)/(N-1),colSums(MIGRATIONB)/(N-1),
               rowSums(MIGRATIONC)/(N-1),colSums(MIGRATIONC)/(N-1),
               rowSums(MIGRATIOND)/(N-1),colSums(MIGRATIOND)/(N-1),
               rowSums(MIGRATIONE)/(N-1),colSums(MIGRATIONE)/(N-1),
               rowSums(MIGRATIONF)/(N-1),colSums(MIGRATIONF)/(N-1)))
CMMIN <- min(c(rowSums(MIGRATIONA)/(N-1),colSums(MIGRATIONA)/(N-1),
               rowSums(MIGRATIONB)/(N-1),colSums(MIGRATIONB)/(N-1),
               rowSums(MIGRATIONC)/(N-1),colSums(MIGRATIONC)/(N-1),
               rowSums(MIGRATIOND)/(N-1),colSums(MIGRATIOND)/(N-1),
               rowSums(MIGRATIONE)/(N-1),colSums(MIGRATIONE)/(N-1),
               rowSums(MIGRATIONF)/(N-1),colSums(MIGRATIONF)/(N-1)))              
# 
# commpanel <-ggarrange(plotmobility(COMMUTINGA, cmax = CMAX), plot_panelc(N, solA, COMMUTINGA, IMIN = IMIN, IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX),
#                      plotmobility(COMMUTINGB, cmax = CMAX), plot_panelc(N, solB, COMMUTINGB, IMIN = IMIN, IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX),
#                      plotmobility(COMMUTINGC, cmax = CMAX), plot_panelc(N, solC, COMMUTINGC, IMIN = IMIN, IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX),
#                      plotmobility(COMMUTINGD, cmax = CMAX), plot_panelc(N, solD, COMMUTINGD, IMIN = IMIN, IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX),
#                      plotmobility(COMMUTINGE, cmax = CMAX), plot_panelc(N, solE, COMMUTINGE, IMIN = IMIN, IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX),
#                      plotmobility(COMMUTINGF, cmax = CMAX), plot_panelc(N, solF, COMMUTINGF, IMIN = IMIN, IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX),
#                      ncol = 2, nrow = 6)
col_lstab <- "#ACA8C4"
col_stab <- "#615E71"
col_ldet <- "#B1C2B6"
col_det <- "#40B565"
col_lrand <- "#7F9BD5"
col_rand <- "#567DD1"
mobstab <- plotmobility(MIGRATION, cmax = CMAXDET, color1 =col_stab )
mobstab
mobA <- plotmobility(MIGRATIONA , cmax = CMAXDET, color1 =col_det )
# mobA <- plotmobilityMPA(reescale*COMMUTINGA - 1, col_det )
mobA

mobB <- plotmobility(MIGRATIONB, cmax = CMAXDET, color1 =col_det )
mobB
mobC <- plotmobility(MIGRATIONC, cmax = CMAXDET, color1 =col_det )
mobC
mobD <- plotmobility(MIGRATIOND, cmax = CMAXRAND, color1 =col_rand )
mobD
mobE <- plotmobility(MIGRATIONE, cmax = CMAXRAND, color1 =col_rand )
mobE
mobF <- plotmobility(MIGRATIONF, cmax = CMAXRAND, color1 =col_rand )
mobF

####INTEGRATION####
int <- plot_panelc(N, sol, COMMUTING,col_mean, col_lstab, col_stab )
int
intA <- plot_panelc(N, solA, MIGRATIONA, col_mean, col_lstab, col_stab )
intA <- plot_int(N, solA, state = "INF") 
intA
intB <- plot_int(N, solB, state = "INF")
intB
intC <- plot_int(N, solC, state = "INF")
intC
intD <- plot_int(N, solD, state = "INF")
intD
intE <- plot_int(N, solE, state = "INF")
intE
intF <- plot_int(N, solF, state = "INF")
intF
text_tit = 15

### COMPARISON

eig_comp <- data.frame(k = integer(0),
                       RMT = numeric(0),
                       LRP_in = numeric(0),
                       LRP_out = numeric(0),
                       LRP_inout = numeric(0),
                       RMT_allr = numeric(0),
                       Real_in = numeric(0),
                       Real_out = numeric(0),
                       Real_inout = numeric(0),
                       Real_sparse = numeric(0),
                       Real_corr = numeric(0),
                       Real_dec = numeric(0)
)

nu <- muwstar - muw
mug <- mud + mua + mudel
max_nodes <- floor(N/4)
max_nodes <- 20
for (nodes in seq(0,max_nodes,by = 2)) {
  
  COMMUTINGA <- COMMUTINGB <- COMMUTINGC <- COMMUTINGD <- COMMUTINGE <- COMMUTINGF <- COMMUTING
  fs <- sample(c(1:N), nodes)
  COMMUTINGA[fs,] <- muwstar
  COMMUTINGB[,fs] <- muwstar
  fcs <- sample(c(1:N), nodes/2)
  fcs <- fs[c(1:(nodes/2))]
  COMMUTINGC[fcs,] <- COMMUTINGC[,fcs] <- muwstar
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
  jacC <- (COMMUTINGC + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  jacD <- (COMMUTINGD + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  jacE <- (COMMUTINGE + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  jacF <- (COMMUTINGF + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  
  k <- nodes
  eig_comp[k+1,] <- c(k,
                      mub - mug + muw*mub*(N-1),
                      N*(mub*muw+muc)/2 + (k-1)*mub*nu/2 +
                        sqrt((N*(mub*muw+muc))^2+((k-1)*mub*nu)^2+2*(mub*nu)*(mub*muw+muc)*(N+N*k-4*k))/2  +
                        mub*(1-muw) - mug - N*muc,
                      N*(mub*muw+muc)/2 + (k-1)*mub*nu/2 +
                        sqrt((N*(mub*muw+muc))^2+((k-1)*mub*nu)^2+2*(mub*nu)*(mub*muw+muc)*(N+N*k-4*k))/2  +
                        mub*(1-muw) - mug - N*muc,
                      N*(mub*muw+muc)/2 + (k/2-1)*mub*nu/2 +
                        sqrt((N*(mub*muw+muc))^2+mub*nu*(2*mub*muw+muc)*(N+3*N*k/2-2*(k/2+1))+
                               ((mub*nu)^2)*(2*k*(N-k/2)+(k/2-1)^2))/2 +
                        mub*(1-muw) - mug - N*muc,
                      mub*(muw*(N-k)+muwstar*k/N)*(N-1)/N + mub - mug,
                      max(eigen_mat(jacA)$re),
                      max(eigen_mat(jacB)$re),
                      max(eigen_mat(jacC)$re),
                      max(eigen_mat(jacD)$re),
                      max(eigen_mat(jacE)$re),
                      max(eigen_mat(jacF)$re)
  )
}

eig_comp <- filter(eig_comp, !is.na(k))

eig_comp_real <- select(eig_comp, k, Real_in, Real_out, Real_inout, Real_sparse, Real_corr, Real_dec) %>%
  pivot_longer(cols = (2:ncol(.)), names_to = "prediction", values_to = "eigenvalue")
eig_comp_pred <- select(eig_comp, k, LRP_in, LRP_inout, RMT_allr) %>%
  pivot_longer(cols = (2:ncol(.)), names_to = "prediction", values_to = "eigenvalue")

# plot_eig_lines <- 
plot_lines <-  eig_comp_pred %>% ggplot(aes(x = k, y = eigenvalue, color = prediction)) +
  geom_line(size = 0.8) + 
  geom_point(data = eig_comp_real, 
             aes(x = k, y = eigenvalue, color = prediction)) +
  scale_color_manual(values=c("#0000FF","#FF0000","#63B159",
                              "#0000FF","#0000FF","#FF0000",
                              "#63B159","#63B159","#63B159")) + 
  theme_bw() + theme(text = element_text(size = 12))
ggsave(file=paste0(Path,"commparison.svg"))


######PANEL 4#####
gstab <- plot_grid(NULL,
                   mobstab + ggtitle("Unstable Scenario\n"),
                   int + theme(aspect.ratio = 1) ,
                   NULL,
                   nrow = 4, align = "v")


gdet <- plot_grid(mobA + ylab("") + ggtitle("Scenario A \n") ,
                  mobB + ggtitle("Scenario B \n"),
                  mobC + ggtitle("Scenario C \n"),
                  intA + theme(aspect.ratio = 1) ,
                  intB  + theme(aspect.ratio = 1) + rremove("ylab") ,
                  intC  + theme(aspect.ratio = 1) + rremove("ylab"), 
                  mobD + ylab("") + ggtitle("Scenario D \n"),
                  mobE + ggtitle("Scenario E \n"),
                  mobF + ggtitle("Scenario F \n"),
                  intD + theme(aspect.ratio = 1)  ,
                  intE + theme(aspect.ratio = 1) + rremove("ylab") ,
                  intF + theme(aspect.ratio = 1) + rremove("ylab") ,
                  ncol = 3, nrow = 4, align = "v")

plotb <- plot_grid(gstab,
                   NULL,
                   gdet,
                   ncol = 3,
                   rel_widths = c(1,0.2,3))

# plotu <- plot_grid( plot_eig_lines + ggtitle("a") + 
#                       theme(legend.position = "none", aspect.ratio = 1),
#                   gg_eigen_full + ggtitle("b"),
#                   ncol = 2, nrow = 1,
#                   rel_widths = c(1,3))

plotm <- plot_grid( plotb,
                    NULL,
                    gg_eigen_full,
                    nrow = 3, rel_heights = c(4,0.3,1))

# This is to get the legend for the eigenvalues vs k plot.
legend <- get_legend(
  # create some space to the left of the legend
  plot_eig_lines +
    theme(legend.box.margin = margin(0, 0, 0, 12), legend.position = "bottom")
)

plot_grid(uplot, legend, nrow = 2, rel_heights = c(3, .4))


####### MIGRATION: DIRECTED VS RANDOM EXODUS ############

# base parameters + reduced muw (to induce stability with control)
muc <- 0.001
sc <- 0.0005

mucstar <- 0.005
MIGRATION <- rand_mat_cor_beta(N, muc*N, sc*N, rhoc, Gammac, rc, cc)[[1]]

# control: two rows
f1 <- 30
f2 <- 55
MIGRATION[f1,c(1:(f1-1),(f1+1):N)] <- MIGRATION[f2,c(1:(f2-1),(f2+1):N)] <- mucstar

# control: two columns
c1 <- 45
c2 <- 80
MIGRATION[c(1:(c1-1),(c1+1):N),c1] <- MIGRATION[c(1:(c2-1),(c2+1):N),c2] <- mucstar

# control: one row and one column
fc <- 65
MIGRATION[fc,c(1:(fc-1),(fc+1):N)] <- MIGRATION[c(1:(fc-1),(fc+1):N),fc] <- mucstar

# control: random
MIGRATION[sample(c(1:N^2),2*(N-1))] <- mucstar


diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)
sum(colSums(MIGRATION) > 1)
mobplot <- plotmobility(MIGRATION)
mobplot
#print(mobplot, vp = viewport(angle = -90))
jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)

plot_eigen_rmt(jacobian,
               N,mub,mug = mud + mua + mudel,
               muw,sw,rhow,Gammaw,
               muc,sc,rhoc,Gammac,
               tau = 0) #+ xlim(c(-215,-210))

sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)
plot_int(N, sol, state = "INF")

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

COMMUTING <- rand_mat_cor_beta(N, N*muw, sqrt(N)*sqrt(N)*sw, rhow, Gammaw, rw, cw)[[1]]
MIGRATION <- rand_mat_cor_norm(N, N*muc, sqrt(N)*sqrt(N)*sc, rhoc, Gammac, rc, cc)[[1]]
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

COMMUTING <- rand_mat_cor_norm(N, N*muw, sqrt(N)*sqrt(N)*sw, -rhow, -Gammaw, rw, cw)[[1]]
MIGRATION <- rand_mat_cor_norm(N, N*muc, sqrt(N)*sqrt(N)*sc, -rhoc, -Gammac, rc, cc)[[1]]
diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)
jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
eigmenor <- eigen_mat(jacobian)
redmenor <- plotmobility(COMMUTING[1:80,1:80])

COMMUTING <- rand_mat_cor_norm(N, N*muw, sqrt(N)*sqrt(N)*sw, rhow, 0.001, 0.001, 0.001)[[1]]
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

N <- 30

# positive rhow vs negative rhow
rhow <- .85
rhow <- -.85
COMMUTING <- rand_mat_ell(N, muw, sw, rhow, distrib = "beta")
diag(COMMUTING) <- rep(0,N)
plotmobility(COMMUTING)

ggsave(file="mobnet_negrho.svg")

# high rw/cw (exchange values)
muw <- 0.2
sw <- 0.06
rhow <- 0 #original rho (Gamma of baron et al)
Gammaw <- .15 #gamma of baron et al
rw <- .7*N
cw <- .1*N
(Gammaw/sqrt(rw*cw) < 1) & ((N*rhow-2*Gammaw)/(N-(rw+cw)) < 1)

COMMUTING <- rand_mat(N, muw, sw, distrib = "beta")
COMMUTING <- rand_mat_cor_beta(N, muw*N, sw*N, rhow, Gammaw, rw, cw)[[1]]
diag(COMMUTING) <- rep(0,N)
plotmobility(COMMUTING)

ggsave(file="mobnet_posr.svg")

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
COMMUTING <- rand_mat_cor_beta(N, muw*N, sw*N, rhow, Gammaw, rw, cw)[[1]]
diag(COMMUTING) <- rep(0,N)
plotmobility(COMMUTING)

ggsave(file="mobnet_negG.svg")

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
COMMUTING <- rand_mat_cor_beta(N, muw*N, sw*N, rhow, Gammaw, rw, cw)[[1]]
diag(COMMUTING) <- rep(0,N)
plotmobility(COMMUTING)

# preventing disease

muw <- 0.25
sw <- 0.1
rhow <- 0.35 #original rho (Gamma of baron et al)
Gammaw <- -.1*N #gamma of baron et al
rw <- .15*N #positive
cw <- .15*N #positive
(abs(Gammaw)/sqrt(rw*cw) < 1) & (rw+cw<N) & (abs(N*rhow-2*Gammaw)/(N-(rw+cw)) < 1)

COMMUTING <- rand_mat(N, muw, sw, distrib = "beta")
COMMUTING <- rand_mat_cor_beta(N, muw*N, sw*N, rhow, Gammaw, rw, cw)[[1]]
diag(COMMUTING) <- rep(0,N)
plotmobility(COMMUTING)

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
ggsave(file="BGvsRMT.svg", plot=BGvsRMT)

ggsave(file="redmenor.svg", plot=redmenor)
ggsave(file="redmayor.svg", plot=redmayor)
ggsave(file="reduncor.svg", plot=reduncor)

corpanel <-ggarrange(ggarrange(redmenor, reduncor, redmayor, labels = c("A", "B", "C"),
                               ncol = 3, nrow = 1),
                     BGvsRMT,
                     ncol = 1, nrow = 2)
corpanel

ggsave(file="corpanel.svg", plot=corpanel)
ggsave(corpanel, file="corpanel.eps", device = cairo_ps, fallback_resolution = 600)
# 
# ##########################################################################################
# out_mig <- function(x){
#   (1/2)*(N*mub*muw+muc) + sqrt((mub*muw+muc)*((N^2*(mub*muw+muc) + (2*(N-1)^2)*x))) + mub*(1-muw) - mug - (N-1)*muc - x
# }
# 
# N = 50
# Deltas <- rep(0.3, N) # birth rate
# mub <- 0.2
# sb <- 0.001
# betas <- rep(mub, N) # transmission rates
# # betas <- rgamma(N, shape = (mub/sb)^2, rate = mub/(sb^2))
# thetas <- rep(0.3, N) # loss of immunity rates
# mud <- 0.3
# deaths <- rep(mud, N) # not disease-related death rates
# mua <- 0.4
# alphas <- rep(mua, N) # recovery rates
# mudel <- 0
# deltas <- rep(mudel, N) # disease-related death rates
# gammas = deaths + alphas + deltas
# 
# seq <- seq(0,2,0.1)
# out <- sapply(seq,out_mig)
# df_mig <- data.frame(seq, out)
# ggplot(df_mig) + 
#   geom_line(aes(seq,out))
# 
# COMMUTING <- rand_mat(N, muw, sw, distrib = "beta")
# diag(COMMUTING) <- 0
# # COMMUTING <- rand_mat_ell(N, muw, sw, rhow, distrib = "beta")
# # COMMUTING[sample.int(N^2, round(p*N^2))] <- 0
# 
# MIGRATION <- rand_mat(N, muc, sc, distrib = "beta")
# diag(MIGRATION) <- 0
# 
# jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
#   diag(deaths + alphas + deltas + colSums(MIGRATION))
# 
# print(plot_eigen(jacobian))
# 
# sus_init <- rep(10000, N) # initial susceptibles
# inf_init <- rep(100, N)    # initial infecteds
# 
# end_time <- 500
# sol.stab <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
#                 COMMUTING,MIGRATION,
#                 sus_init,inf_init,end_time)
# 
# plot_stab <- plot_int1(N, sol.stab, state = "INF") +
#   theme_bw() +theme(legend.position="none")
# plot_stab
# 
# 
# MIGRATION[,c(2:N)] <- 2
# jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
#   diag(deaths + alphas + deltas + colSums(MIGRATION))
# 
# print(plot_eigen(jacobian))
# 
# sus_init <- rep(10000, N) # initial susceptibles
# inf_init <- rep(100, N)    # initial infecteds
# 
# end_time <- 500
# sol.unstab <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
#                 COMMUTING,MIGRATION,
#                 sus_init,inf_init,end_time)
# 
# plot_unstab <- plot_int1(N, sol.unstab, state = "INF") +
#   theme_bw() +theme(legend.position="none")
# plot_unstab
#  
# N = 3
# 
# ac = 1
# bc = 2
# a <- matrix(ac,N,N)
# a[,c(2:N)] <- ac + bc
# max(eigen_mat(a)$re) == (1/2)*(N*ac + sqrt(ac*(N^2*ac + 2*(N-1)*bc)))
