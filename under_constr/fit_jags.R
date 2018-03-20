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

# Define dimensions
n_loc = dim(cube)[3]
n_spp = dim(cube)[2]
n_time = dim(cube)[1]

# Function to center and scale
cent_scale = function(x){z = (x - mean(x,na.rm=T))/sd(x)}
cent = function(x){z = (x - mean(x,na.rm=T))}

# Define function to bundle data
bundle_fn = function(data){
    
    # Create list of day x hour data frames
    y = lapply(data$loc %>% unique, 
               function(x){data %>%
                       filter(loc==x) %>%
                       select(-n,-loc) %>%
                       arrange(time) %>%
                       spread(species, n_z) %>%
                       select(-time)})
    
    # Unlist and store in day x hour x year array
    y %>% 
        unlist %>% 
        array(dim=c(y[[1]] %>% dim, length(y))) %>%
        # Rearrange array to have dimensions year x day x hour
        aperm(c(3,2,1))
} 

# Center and scale data
data = cube_df %>% gather(species, n, -loc, -time) %>%
    group_by(species) %>%
    mutate(n_z = cent(n)) %>%
    ungroup

# Scale and bundle n_obs
n_obs = bundle_fn(data)

# Scale X
x = array(cent_scale(gen_data$X), c(n_time,n_loc)) %>% aperm(c(2,1))

data_list = list(
    # Rearrange array to have dimensions site x species x time
    n_obs = n_obs,
    # Initial population sizes
    n_init = n_obs[,,1],
    # Environmental variables
    x = x,
    # Count dimensions
    n_loc = n_loc, 
    n_spp = n_spp,
    n_time = n_time)







#==================================================
#Fit Model
#==================================================
# Fix parameters
data_list$sigma_obs = 0.5

# Parameters to monitor
params = c("mean_b1","sigma_b1","sigma_eps","b1","rho","eps","n")

# Initial values
inits = function(){
    list(mean_b1 = rnorm(1,0,1),
         sigma_b1 = runif(1,0,5),
         sigma_eps = runif(1,0,5),
         rho = runif(n_spp,0.2,0.8))
}

# MCMC specifiations
n.ad = 100 # number of iterations for adaptation
n.it = 6000 # number of iterations for Markov Chain
n.burn = n.it/2 # number steps to omit for burnin
n.thin = 20 # thinning interval
n.chain = 4 # number of independent chains

# Fit model using MCMC via JAGS
output = jags.basic(data=data_list, model="under_constr/rep_ts_jags.txt", 
                    parameters.to.save=params, n.adapt=n.ad, n.burnin=n.burn, 
                    n.iter=n.it, n.thin=n.thin, n.chains=n.chain,parallel=T, DIC=F)



# Parameters to check convergence
# Define function to define parameters by species
thread_f = function(var,n_spp){sapply(1:n_spp, function(x){paste0(var,"[",x,"]")})}
conv_pars = c("mean_b1","sigma_b1","sigma_eps",sapply(c("b1","rho"), thread_f, n_spp))
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
#Plot output: Parameters by Species
#==================================================

spec_par = data_frame(species = 1:ncol(cube),
                      b1 = gen_data$b1_mat[,1],
                      rho = gen_data$rho_mat[,1]) %>%
    gather(var,value,-species)

# Select variable v
v = "rho"
samples_long %>% 
    filter(var==v) %>%
    rename(species=v1) %>%
    ggplot(aes(value, color=factor(species)))+
    geom_density()+
    geom_vline(data = spec_par %>% filter(var==v), aes(xintercept = value, 
               color=factor(species)),
               linetype=2)+
    theme_classic()



# sigma_eps
samples_long %>% 
    filter(var=="sigma_eps") %>%
    rename(species=v1) %>%
    ggplot(aes(value))+
    geom_density()+
    theme_classic()

# sigma_b1
samples_long %>% 
    filter(var=="sigma_eps") %>%
    rename(species=v1) %>%
    ggplot(aes(value))+
    geom_density()+
    theme_classic()

