#==================================================
#Preliminaries
#==================================================

# Load packages
library(tidyverse)
library(coda)
library(jagsUI)

#simulations.R





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
    nsites = dim(cube)[3],nspecies = dim(cube)[2],ntimes = dim(cube)[1])

# Define dimensions in environment
nsites = data_list$nsites
nspecies = data_list$nspecies
ntimes = data_list$nspecies



#==================================================
#Fit Model
#==================================================
# Fix observation error to 0.1
data_list$sigma_obs = 0.1

# Parameters to monitor
params = c("mean_b0","mean_b1","mean_rho","sigma_b0","sigma_b1","sigma_rho","sigma_e","b0","b1","rho","e","n")

# Initial values
inits = function(){
    list(mean_b0 = runif(1,-10,10),
         mean_b1 = runif(1,-1,1),
         mean_rho = runif(1,-10,10),
         sigma_b0 = runif(1,0,10),
         sigma_b1 = runif(1,0,1),
         sigma_rho = runif(1,0,1))
}

# MCMC specifiations
n.ad = 100 # number of iterations for adaptation
n.it = 5000 # number of iterations for Markov Chain
n.burn = n.it/2 # number steps to omit for burnin
n.thin = 20 # thinning interval
n.chain = 4 # number of independent chains

# Fit model using MCMC via JAGS
output = jags.basic(data=data_list, model="under_constr/rep_time_model.txt", 
                    parameters.to.save=params, n.adapt=n.ad, n.burnin=n.burn, 
                    n.iter=n.it, n.thin=n.thin, n.chains=n.chain,parallel=T, DIC=F)



# Parameters to check convergence
conv_pars = c("mean_b0","mean_b1","mean_rho","sigma_b0","sigma_b1","sigma_rho","sigma_e")
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
samples_long %>% 
    filter(var=="rho") %>%
    rename(species=v1) %>%
    ggplot(aes(value, color=factor(species)))+
    geom_density()+
    theme_classic()



