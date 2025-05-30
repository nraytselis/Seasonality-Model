# Loading packages
library(Matrix)
library(deSolve)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(deSolve)


Pond_ODE =function(t, y, parameters) {
  
  with(as.list(parameters),{
    N=y[1]; J=y[2]; A=y[3]; Es = y[4:(4+latent_stages - 1)]; I = y[4+latent_stages]; Preds = y[5+latent_stages]; L3F = y[6+latent_stages]
    VOL = 1
    
    Pred_A = f*Preds/(1 + f*h*(A + f_J*h_J*J + f_N*h_N*N)/VOL  + i_P*max(Preds-1, 0)/VOL)
    Pred_J = f*f_J*Preds/(1 + f*h*(A + f_J*h_J*J + f_N*h_N*N)/VOL  + i_P*max(Preds-1, 0)/VOL)
    Pred_N = f*f_N*Preds/(1 + f*h*(A + f_J*h_J*J + f_N*h_N*N)/VOL  + i_P*max(Preds-1, 0)/VOL)
    
    d_A_c = d_A*exp(comp_d / VOL * (c_N * N + c_J * J + A + sum(Es))) #density dependence in deaths
    
    d_J_c = d_J*exp(comp_d / VOL * (c_N * N + c_J * J + A + sum(Es))) #density dependence in deaths
    
    d_N_c = d_N*exp(comp_d / VOL * (c_N * N + c_J * J + A + sum(Es))) #density dependence in deaths
    
    m_N_c = m_N*exp(-comp_M/VOL*(c_N*N + c_J*J + A + sum(Es))) #density dependence in maturation
    
    m_J_c = m_J*exp(-comp_M/VOL*(c_N*N + c_J*J + A + sum(Es))) #density dependence in maturation
    
    dNdt = b_M*(A + sum(Es))/2*exp(-comp_b/VOL*(c_N*N + c_J*J + A + sum(Es))) - (m_N_c+d_N_c)*N - can*(A + I + sum(Es))*N - Pred_N*N
    
    dJdt = m_N_c*N - (m_J_c+d_J_c)*J -Pred_J*J
    
    dAdt = m_J_c*J - d_A_c*A - lambda*A - Pred_A*A
    
    # development of all stages
    latent_progression = latent_rate*Es
    # lost to next stage   #death      #gained from last stage
    dEsdt = -latent_progression - d_A_c*Es + c(lambda*A, latent_progression[1:(latent_stages - 1)]) - Pred_A*Es
    
    dIdt = as.numeric(latent_progression[latent_stages]) - d_A_c*I - Pred_A*I
    
    dPredsdt = Pred_N*N + Pred_J*J + Pred_A*A + Pred_A*sum(Es) - d_F*Preds   
    
    dL3Fdt <- Pred_A*I - d_W*L3F - d_F*L3F 
    
    result = c(dNdt,dJdt,dAdt, dEsdt, dIdt,dPredsdt, dL3Fdt)
    
    
    return(list(result))
  }
  )
}

ReboundParams <- readRDS(file = "Rebound_parameters3.RDA") #bring in best parameters
parameters = as.list(signif(ReboundParams$samples[which.max(ReboundParams$log.p),], 3))
parameters["latent_stages"] = 60
parameters["latent_rate"] = 4.3
parameters["lambda"] = 1.3
parameters["comp_M"] = 0.01
parameters["comp_b"] = 0.001
parameters["comp_d"] = 0.001
parameters["can"] = 0
parameters["f"] = 0.4
parameters["f_N"] = 4.704984e-02
parameters["f_J"] = 3.471900e-02
parameters["h"] = 0.3
parameters["h_N"] = 4.675860e-01
parameters["h_J"] = 3.951459e-01
parameters["i_P"] = 1.911824e-01
parameters["d_W"] = 0.05
parameters["d_F"] = 0.1

parameters = unlist(parameters)

Exposed_names = paste0("E", 1:parameters["latent_stages"])
Exposed_values = rep(0, times=parameters["latent_stages"])
names(Exposed_values) = Exposed_names
Exposed_values

Initial_conditions = c(N = 7500, J = 6000, A = 700, Exposed_values, I = 0, Preds = 15,L3F = 0)/15
timespan = 365*5

#fish added during rainy season
intro_fish = expand.grid(day = c(1:180, 270:365), year = 1:5) 
intro_fish$time <- intro_fish$day + 365 * intro_fish$year
n <- nrow(intro_fish)
event_data_fish = data.frame(var = "Preds", time = intro_fish$time, value = 5, method = "add")

#fishing during dry season (fish and L3Fs)
remove_fishing <- expand.grid(day = 181:269, year = 1:5)
remove_fishing$time <- remove_fishing$day + 365 * remove_fishing$year
nfishing <- nrow(remove_fishing)
var_choices <- sample(c("Preds", "L3F"), size = nfishing, replace = TRUE, prob = c(0.8, 0.2)) 
event_data_fishing = data.frame(var = var_choices, time = remove_fishing$time, value = -10, method = "add")

#combined list of events
all_events = rbind(event_data_fish, event_data_fishing)


#run simulation
PondSim = data.frame(ode(y = Initial_conditions, times=1:timespan, parms=parameters, hmax=1,
                            method="lsoda", func=Pond_ODE, events = list(data = all_events))) 



#reformat sim data frame
PondSim[,"Es"] = rowSums(PondSim) - PondSim[,"N"] - PondSim[,"J"]- PondSim[,"A"] - PondSim[,"I"] - PondSim[,"time"]
PondSim[,"Exposed"] = apply(X=PondSim[,which(str_detect(colnames(PondSim),"E"))],MARGIN=1,FUN = sum) 

#check to make sure that events worked, NEED TO CHANGE EVENT SO THAT CANNOT GO NEGATIVE 
ggplot(PondSim, aes(x = time)) +
  geom_line(aes(y = Preds), color = "blue") +
  geom_line(aes(y = L3F), color = "red") +
  labs(y = "Fish (blue = total, red = infected)", title = "Predator Fish Dynamics Over Time") +
  theme_minimal()

#reformat for plotting
PondSim = PondSim %>% select(time,N,J,A,Exposed,I,Preds, L3F)
PondSim = PondSim %>% pivot_longer(cols = c(N,J,A,Exposed,I,Preds, L3F)) 

#plot sim 
p1 = ggplot(PondSim, aes(x=time, y = value + 0.01 , group = name, color = name)) + 
  geom_line() + ylab("density per L") + theme_minimal() + scale_y_log10()
p1

#need to make sure L3Fs cannot go negative, as they currently can