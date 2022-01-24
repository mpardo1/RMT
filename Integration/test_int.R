rm(list = ls())

library("deSolve")
library("ggplot2")
library("gifski")
# Carga el fichero de las funciones:
source("~/RMT/Integration/functions_eigen_int.R")

# NOTE: In the system the variables are ordered as susceptible all patches, infected all patches
# and recovered all patches, so N*3 + 1 are the output variables of the integration, the +1 is
# the time, which is the first column after the rest.

#----------------PARAMETERS-----------------------#
#-------------------LOGIC-----------------#
# Mobility parameter, 0: Just commuting, 1: just migration 2: migration & commuting.
MOB <- 2
# Integration parameter, 0: No integration, 1: integration.
INT <- 1
# Parameter for initial population. 0: No cte, 1: cte.
CTE_POP <- 0
# Parameter for transmission rate. 0: No cte, 1: cte.
BETA_CTE <- 0
# Parameter for initial infected ind. 0: No cte , 1: cte.
CTE_INF <- 1


#-------------------EPIDEMIOLOGICAL------------------------#
N = 25 # Number of patches
del_N <- rep(0.6, N) # Birth rate
# bet_cte <-  0.001
# bet <- rep(bet_cte, N)  # Transmission rate cte
bet <- abs(rnorm(N,1,1))  # Transmission rate random
d_vec <- rep(0.8, N) # Natural mortality rate
thet <- rep(0.1, N) # Rate of loss of immunity
alp <- rep(0.02, N) # Rate of disease overcome
delt <- rep(0, N) # Diseases related mortality rate

#-------------------- MOBILITY ---------------------#
### Migration:
alp_m <- 0.01
bet_m <- 0.1

# Compute mean and sd:
mu_m <- alp_m/(alp_m + bet_m)
s_m <-  sqrt((alp_m*bet_m)/(((alp_m + bet_m)^2)*(1+alp_m+bet_m)))
print(paste0("mu :", mu_m))
print(paste0("sigma :", s_m))

migrate_mat <- mat_conect(N,alp_m,bet_m,MOB)
### Commuting
alp_c <- 0.1
bet_c <- 0.1

# Compute mean and sd:
mu_w <- alp_c/(alp_c + bet_c)
s_w <- sqrt((alp_c*bet_c)/((alp_c + bet_c)^2*(1+alp_c+bet_c)))
print(paste0("mu :", mu_w))
print(paste0("sigma :",s_w))

commut_mat <- mat_conect(N,alp_c,bet_c,MOB)

tau_ct <- 0

print(paste0("beta - gamma:", bet[1] - (alp[1] + delt[1] + d_vec[1])))
#-----------------POPULATION INIT----------------------#
# Number of initial individuals by compartments:
SUS_INIT <- 100000 #Susceptible
INF_INIT <- 100    #Infected

#-------------TEMPORAL---------#
# End time integration:
end_time <- 30


#--------------------------INTEGRATION---------------------------------#
# Integro el sistema con condiciones iniciales 
sol <- int(N, del_N,bet,d_vec,thet,alp,delt,
           commut_mat,migrate_mat,end_time,
           MOB, CTE_POP, CTE_INF,SUS_INIT, INF_INIT)

sol_df <-  as.data.frame(sol)
# Extract infected individuals:
inf_df <- sol_df[, c(1,(N+2):(2*N+1))]

# Plot the susceptible (SUS), infected (INF) and recovered (REC) o todos (TOT):
state <- "INF"
plot_inf <- plot_int(N, sol, state)

# Rename columns:
for(i in c(1:N)){
  colnames(sol_df)[i+1] <-  paste0("S",i)
  colnames(sol_df)[N+i+1] <-  paste0("I",i)
  colnames(sol_df)[2*N+i+1] <-  paste0("R",i)
}


