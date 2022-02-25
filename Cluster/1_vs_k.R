rm(list = ls())
library("tidyverse")
library("deSolve")
library("ggplot2")

#----------------------------------------------------------------------------#
source("~/RMT/Integration/functions_eigen_int.R")
#----------------PARAMETERS-----------------
#-------------------LOGIC----------------------
# Mobility parameter, 0: Just commuting, 1: just migration 2: migration & commuting.
MOB <- 2
# Integration parameter, 0: No integration, 1: integration.
INT <- 1
# Parameter for initial population. 0: No cte, 1: cte.
CTE_POP <- 1
# Parameter for constant population at each patch. 0: No cte, 1: cte.
CTE_POP_Patch <- 0
# Parameter for transmission rate. 0: No cte, 1: cte.
BETA_CTE <- 0
# Parameter for initial infected ind. 0: No cte , 1: cte.
CTE_INF <- 1
# Parameter for distribution, "normal", "beta", "gamma":
DIST <-  "beta"

#-------------------EPIDEMIOLOGICAL----------------------
N = 100 # Number of patches
# CTE parameters:
del_N <- rep(0.7, N) # Birth rate

#-------------------- MOBILITY ---------------------------
### Migration:
mu_c <- 0.01
s_c <- 0.00001
alp_c <- beta_a_b(mu_c, s_c)[1]
bet_c <- beta_a_b(mu_c, s_c)[2]
# Compute mean and sd:
print(paste0("mu :", mu_c))
print(paste0("sigma :", s_c))

migrate_mat <- mat_conect(N,alp_c,bet_c,DIST)
### Commuting
mu_w <- 0.4
s_w <- 0.2
alp_w <- beta_a_b(mu_w, s_w)[1]
bet_w <- beta_a_b(mu_w, s_w)[2]

# Compute mean and sd:
print(paste0("mu :", mu_w))
print(paste0("sigma :",s_w))

commut_mat <- mat_conect(N,alp_w,bet_w,DIST)

tau_ct <- 0
# Initial populations:
init_pop <- matrix(100, nrow = N)

#-----------------POPULATION INIT----------------------#
# Number of initial individuals by compartments:
SUS_INIT <- 10000 #Susceptible
INF_INIT <- 100    #Infected

# End time integration:
end_time <- 100

#-------------------------------------------------------------------------#
# ------------------COMPARISON RIGHT MOST----------------------------
d_vec <- rep(0.7, N) # Natural mortality rate
alp.vec <- rep(0.72, N) # Rate of disease overcome
delt <- rep(0, N) # Diseases related mortality rate
# Check if the right most eigenvalue is the same for many iterations
thet <- rep(0.6, N) # Rate of loss of immunity
gamma_ct <-  alp.vec[1] + delt[1] + d_vec[1]

d <- 100
alphag <- alphagamma(2,3)
betag <- betagamma(2,3)
bet_vec <- rgamma(d,alphag,betag) 
alp_vec <- rgamma(d,alphag,betag) 
df.comp <- data.frame(N.patches=0,out1=0,outk=0,diff=0,
                      out.pred.1=0, out.pred.k=0, max.inf.1=0,max.inf.k=0)
for(j in c(1:d)){
  ind <-  sample(1:N,1)
  bet_cte <- bet_vec[j]
  bet_new <- alp_vec[j]
  alp <- bet_new
  dim_l <- N
  mat.comp <-  matrix(0, ncol = 8, nrow = dim_l)
  for(i in c(1:dim_l)){
    # Compute  right most eigenvalue for 1 patch:
    bet <-  rep(bet_cte,N)
    bet[ind] <- bet[ind] + alp
    # Compute the jacobian matrix:
    jac <- jacobian(N, bet, gamma_ct, commut_mat, migrate_mat, mu_c, MOB)
    eig <- eigen_mat(jac)
    out1 <-  max(eig$re)
    out.pred.1 <- max(outl_1patch(alp, bet_cte, N, mu_w, mu_c, gamma_ct))
    
    sol <- int(N, del_N,bet,d_vec,thet,alp.vec,delt,
               commut_mat,migrate_mat,5,
               MOB, CTE_POP, CTE_INF,SUS_INIT, INF_INIT,init_pop)
    
    # Compute  right most eigenvalue for K patches:
    bet <-  rep(bet_cte,N)
    bet[1:i] <- bet[ind] + alp/i
    # Compute the jacobian matrix:
    jac <- jacobian(N,bet,gamma_ct, commut_mat, migrate_mat,mu_c, MOB)
    eig <- eigen_mat(jac)
    outk <-  max(eig$re)
    out.pred.k <- outl_Kpatch(i, alp/i, bet_cte, N, mu_w, mu_c, gamma_ct)
    
    sol.k <- int(N, del_N,bet,d_vec,thet,alp.vec,delt,
                 commut_mat,migrate_mat,5,
                 MOB, CTE_POP, CTE_INF,SUS_INIT, INF_INIT,init_pop)
    
    print(paste0("i:",i))
    mat.comp[i,] = c(i, out1, outk,
                     abs((outk -out1)/outk),out.pred.1,
                     out.pred.k, max(sol[,c((N+2):(2*N+1))]),
                     max(sol.k[,c((N+2):(2*N+1))]))
  }
  print(paste0("mat.comp",mat.comp[1:2,1:5]))
  print("Uno los dos df")
  df.comp1 <-  as.data.frame(mat.comp)
  head(df.comp1)
  colnames(df.comp1) <-  c("N.patches","out1","outk","diff",
                           "out.pred.1", "out.pred.k", "max.inf.1","max.inf.k")
  df.comp <- rbind(df.comp, df.comp1)
}
path <- paste0("~/RMT/Integration/1_vs_k_",Sys.Date(),".csv")
write.csv(df.comp,path, row.names = TRUE)
