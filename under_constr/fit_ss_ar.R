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

# Function to center and scale
cent_scale = function(x){z = (x - mean(x,na.rm=T))/sd(x)}

# Scale observations
n_obs = cent_scale(cube[,,1])
x = cent_scale(X_[,1])

# Define time dimension
n_time = length(n_obs)

# Package data
data_list = list(
    # Rearrange array to have dimensions site x species x time
    n_obs = n_obs,
    # Initial population sizes
    n_init = n_obs[1],
    # Environmental variables
    x = x,
    # Count dimensions
    n_time = n_time)







#==================================================
#Fit Model
#==================================================

# Fixed values
data_list$sigma_obs = 0.5

# Parameters to monitor
params = c("b1","rho","sigma_eps","n")

# Initial values
inits = function(){
    list(b1 = rnorm(1,0,1),
         rho = runif(1,0.2,0.8),
         sigma_eps = runif(1,0,3),
         sigma_obs = runif(1,0,3))}

# MCMC specifiations
n.ad = 100 # number of iterations for adaptation
n.it = 5000 # number of iterations for Markov Chain
n.burn = n.it/2 # number steps to omit for burnin
n.thin = 20 # thinning interval
n.chain = 4 # number of independent chains

# Fit model using MCMC via JAGS
output = jags.basic(data=data_list, model="under_constr/ss_ar_jags.txt", 
                    parameters.to.save=params, n.adapt=n.ad, n.burnin=n.burn, 
                    n.iter=n.it, n.thin=n.thin, n.chains=n.chain,parallel=T, DIC=F)

# Parameters to check convergence
gelman.diag(output)






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
        time = strsplit(name, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))
    ) %>% select(-name)




#==================================================
#Plot output
#==================================================

# b1
samples_long %>% 
    filter(var=="b1") %>%
    ggplot(aes(value, color=factor(chain)))+
    geom_density()+
    theme_classic()


# rho
samples_long %>% 
    filter(var=="rho") %>%
    ggplot(aes(value, color=factor(chain)))+
    geom_density()+
    theme_classic()

# sigma_eps
samples_long %>% 
    filter(var=="sigma_eps") %>%
    ggplot(aes(value, color=factor(chain)))+
    geom_density()+
    theme_classic()

# population size
pop_sum = samples_long %>%
    filter(var=="n") %>%
    group_by(time) %>%
    summarize(upper = quantile(value, 0.95),
              lower = quantile(value, 0.05),
              value = median(value))

pop_sum %>%
    ggplot(aes(time, value))+
    geom_line(size=0.8, alpha=0.7)+
    geom_ribbon(aes(ymin=lower, ymax=upper),size=0.75, alpha=0.3)+
    geom_point(aes(y=n_obs), size=2)+
    theme_classic()


