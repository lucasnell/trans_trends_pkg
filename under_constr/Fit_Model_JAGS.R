#==================================================
#Preliminaries
#==================================================

# Load packages
library(tidyverse)
library(coda)
library(jagsUI)

#Generate_Data_JP.R





#==================================================
#Prepare Data
#==================================================

data_list = list(
    # Rearrange array to have dimensions site x species x time
    n_obs = cube %>% aperm(c(3,2,1)),
    # Initial population sizes
    n_init = cube[1,,] %>% aperm(c(2,1)),
    # Environmental variables
    x = X_ %>% aperm(c(2,1)),
    # Count dimensions
    n_loc = dim(cube)[3],n_spp = dim(cube)[2],n_time = dim(cube)[1])

# Define dimensions in environment
n_loc = data_list$n_loc
n_spp = data_list$n_spp
n_time = data_list$n_spp



#==================================================
#Fit Model
#==================================================
# Fix parameters
data_list$sigma_obs = 1e-10

# Parameters to monitor
params = c("mean_b0","mean_b1","mean_rho","sigma_b0","sigma_b1","sigma_rho","sigma_eps","b0","b1","rho","eps","n")

# Initial values
inits = function(){
    list(mean_b0 = runif(1,-10,10),
         mean_b1 = runif(1,-1,1),
         mean_rho = runif(1,-10,10),
         sigma_b0 = runif(1,0,1),
         sigma_b1 = runif(1,0,1),
         sigma_rho = runif(1,0,1),
         sigma_eps = runif(1,0,1))
}

# MCMC specifiations
n.ad = 100 # number of iterations for adaptation
n.it = 4000 # number of iterations for Markov Chain
n.burn = n.it/2 # number steps to omit for burnin
n.thin = 20 # thinning interval
n.chain = 4 # number of independent chains

# Fit model using MCMC via JAGS
output = jags.basic(data=data_list, model="under_constr/rep_time_model.txt", 
                    parameters.to.save=params, n.adapt=n.ad, n.burnin=n.burn, 
                    n.iter=n.it, n.thin=n.thin, n.chains=n.chain,parallel=T, DIC=F)



# Parameters to check convergence
conv_pars = c("mean_b0","mean_b1","sigma_b0","sigma_b1","rho[4]")
output_conv = output[,conv_pars]
gelman.diag(output_conv)






#==================================================
#Clean output
#==================================================
samples_fun = function(y) {sapply(1:length(y), function(x){
    y[[x]] %>% tbl_df %>%
        mutate(chain = x, step = 1:nrow(y[[1]]))
}, simplify=F, USE.NAMES = T) %>% bind_rows %>%
        mutate(
            chain = factor(chain)
        ) 
}

samples = samples_fun(output)

samples_long = samples %>%
    gather(name, value, -matches('chain'), -matches('step')) %>%
    mutate(
        var = strsplit(name, "\\[|\\]|,") %>% map_chr(~.x[1]),
        v1 = strsplit(name, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
        v2 = strsplit(name, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
        v3 = strsplit(name, "\\[|\\]|,") %>% map_int(~as.integer(.x[4]))
    ) %>% select(-name)




#==================================================
#Plot output
#==================================================

# by species

# b1
b1s = data.frame(species=1:n_spp,value=b1_mat_[,1])
samples_long %>% 
    filter(var=="b1") %>%
    rename(species=v1) %>%
    ggplot(aes(value, color=factor(species)))+
    geom_density()+
    geom_vline(xintercept = b1s$value)+
    theme_classic()

# b0
b0s = data.frame(species=1:n_spp,value=b0_mat_[,1])
samples_long %>% 
    filter(var=="b0") %>%
    rename(species=v1) %>%
    ggplot(aes(value, color=factor(species)))+
    geom_vline(xintercept = b0s$value)+
    geom_density()+
    theme_classic()

# rho0
samples_long %>% 
    filter(var=="rho") %>%
    rename(species=v1) %>%
    ggplot(aes(value, color=factor(species)))+
    geom_density()+
    theme_classic()
