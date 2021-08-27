########################################################################
##### RStan Bayesian GLM estimation and diagnostics 
##### Nearshore fish abundance and covariates
##### Written by Genoa Sullaway and Tessa Francis
##### Last edited August 23, 2021

# load libraries
library(tidyverse)
library(rstan)
library(rstanarm)
library(bayesplot)
library(tidybayes)
library(modelr)
library(LaplacesDemon)
library(here)
library(loo)
library(cowplot)
library(RColorBrewer)
library(testthat)
library(emmeans)
library(modelbased)
library(ggpubr)
library(bridgesampling)
library(stats)

# install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
# remotes::install_github("stan-dev/rstanarm")

set.seed(50)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

net<- here::here("data","net.csv")
net<- read.csv(net)

#snorkel<- here::here("data","snorkel.csv")
#snorkel<- read.csv(snorkel)


#set up dataframe for count data using net data
unique_df<- unique(net[c("year", "month","site","ipa")])
species <- c("Chinook", "Chum", "Surf Smelt", "Herring" )
all<- expand_grid(unique_df, species) 


df<- net %>%
  select(-1) %>%
  filter(species %in% c("Chinook", "Chum", "Surf Smelt", "Herring"))%>%
  group_by(year,month,site,ipa,species) %>%
  summarise(total= sum(species_count)) %>% 
  ungroup( ) %>%
  mutate(total = as.numeric(total)) %>%
  mutate_if(is.character, as.factor) %>%
  merge(all, all=TRUE) %>%
  mutate(total=replace_na(total, 0)) %>%
  mutate(veg = recode(site, 'COR' = "Present",
                      'TUR'="Present",
                      'FAM'="Absent",
                      'DOK'="Absent",
                      'EDG'="Absent",
                      'SHR'="Present")) %>%
  mutate(age = recode(site, 'COR' = 7,
                      'TUR'= 4,
                      'FAM'= 4,
                      'DOK'= 6,
                      'EDG'= 3,
                      'SHR'= 5)) %>% 
# to include the effect of period
# mutate(period = case_when(
# month > 6 & year > 2018 ~ "period4",
# month < 7  & year > 2018 ~ "period3",
# month > 6 & year < 2019 ~ "period2",
# TRUE ~ "period1")) %>%
filter(!total%in% c(188,357,592))  



# prepare to estimate models
col.filters <- unique(df$species) 

lapply(seq_along(col.filters), function(x) {
  filter(df, species == col.filters[x])
}
) -> df_list

names(df_list) <- col.filters


# List of models to estimate 
# also versions of each model (except a-c/1-3) excluding season
# model0  count~period
# modela  count~site
# modela_2  count~site + dock
# modelb  count~veg
# modelb_2  count~veg + dock
# modelc  count~ipa
# model4  count~ipa + veg
# model5  count~site + veg
# model6  count~site + ipa
# model6  count~site + ipa + ipa*dock
# model7  count~site + ipa + veg
# model7_2  count~site + ipa + veg + dock
# model8  count~ipa + veg + veg*ipa
# model9  count~site + ipa + site*ipa
# model10  count~site + ipa + veg + veg*ipa
# model11  count~site + ipa + veg + site*ipa
# model11_2  count~site + ipa + veg + site*ipa + dock
# model12  count~site + ipa + veg + site*ipa + veg*ipa
# model13 count~dock

# Best model for Chinook: modela
# Best model for chum: modelb
# Best model for herring: model7
# Best model for smelt: model11


model0 = list()
for(i in 1: length(df_list)) {
  model0[[i]] <- stan_glm(total ~ period,
                          data = df_list[[i]],
                          family = neg_binomial_2(link="log"), init_r=1,
                          adapt_delta = 0.999, iter=5000, warmup=3000)
                          # prior_intercept = normal(0,10), prior = normal(0,10)
}


modela = list()
for(i in 1: length(df_list)) {
  modela[[i]] <- stan_glm(total ~ site,
                          data = df_list[[i]],
                          family = neg_binomial_2(link="log"), init_r=1,
                          adapt_delta = 0.999, iter=5000, warmup=3000)
                          # prior_intercept = normal(0,10), prior = normal(0,10)
  
}

modela_2 = list()
for(i in 1: length(df_list)) {
  modela_2[[i]] <- stan_glm(total ~ site + dock,
                          data = df_list[[i]],
                          family = neg_binomial_2(link="log"), init_r=1,
                          adapt_delta = 0.999, iter=5000, warmup=3000)
  # prior_intercept = normal(0,10), prior = normal(0,10)
  
}

modelb = list()
for(i in 1: length(df_list)) {
  modelb[[i]] <- stan_glm(total ~ veg,
                          data = df_list[[i]],
                          family = neg_binomial_2(link="log"), init_r=1,
                          adapt_delta = 0.999, iter=5000, warmup=3000)
                          # prior_intercept = normal(0,10), prior = normal(0,10)

}

modelb_2 = list()
for(i in 1: length(df_list)) {
  modelb_2[[i]] <- stan_glm(total ~ veg + dock,
                          data = df_list[[i]],
                          family = neg_binomial_2(link="log"), init_r=1,
                          adapt_delta = 0.999, iter=5000, warmup=3000)
  # prior_intercept = normal(0,10), prior = normal(0,10)
  
}

modelc = list()
for(i in 1: length(df_list)) {
  modelc[[i]] <- stan_glm(total ~ ipa,
                          data = df_list[[i]],
                          family = neg_binomial_2(link="log"), init_r=1,
                          adapt_delta = 0.999, iter=5000, warmup=3000)
                          # prior_intercept = normal(0,10), prior = normal(0,10)

}

modelc_2 = list()
for(i in 1: length(df_list)) {
  modelc_2[[i]] <- stan_glm(total ~ ipa + ipa*dock,
                          data = df_list[[i]],
                          family = neg_binomial_2(link="log"), init_r=1,
                          adapt_delta = 0.999, iter=5000, warmup=3000)
  # prior_intercept = normal(0,10), prior = normal(0,10)
  
}

#model 1 not needed if we exclude season from our models
model1 = list()
for(i in 1: length(df_list)) {
  model1[[i]] <- stan_glm(total ~ period + veg,
                          data = df_list[[i]],
                          family = neg_binomial_2(link="log"), init_r=1,
                          adapt_delta = 0.999,
                          diagnostic_file = file.path(tempdir(), "df.csv")
                          # prior_intercept = normal(0,10), prior = normal(0,10)
  )
}
print(model1 )


#model 2 not needed if we exclude season from our models
model2 = list()
for(i in 1: length(df_list)) {
  model2[[i]] <- stan_glm(total ~ period + ipa,
                          data = df_list[[i]],
                          family = neg_binomial_2(link="log"), init_r=1,
                          adapt_delta = 0.999
                          # prior_intercept = normal(0,10), prior = normal(0,10)
  )
}
print(model2)

#model 3 not needed if we exclude season from our models
model3 = list()
for(i in 1: length(df_list)) {
  model3[[i]] <- stan_glm(total ~ period + site,
                          data = df_list[[i]],
                          family = neg_binomial_2(link="log"), init_r=1,
                          adapt_delta = 0.999
                          # prior_intercept = normal(0,10), prior = normal(0,10)
  )
}
print(model3)

# version with season
model4 = list()
for(i in 1: length(df_list)) {
  model4[[i]] <- stan_glm(total ~ period + ipa + veg,
                          data = df_list[[i]],
                          family = neg_binomial_2(link="log"), init_r=1,
                          adapt_delta = 0.999
                          # prior_intercept = normal(0,10), prior = normal(0,10)
  )
}
print(model4)

# version without season
model4 = list()
for(i in 1: length(df_list)) {
  model4[[i]] <- stan_glm(total ~ ipa + veg,
                          data = df_list[[i]],
                          family = neg_binomial_2(link="log"), init_r=1,
                          adapt_delta = 0.999
                          # prior_intercept = normal(0,10), prior = normal(0,10)
                          , 
                          iter=5000, warmup=3000)
}
print(model4)

#version with season
model5 = list()
for(i in 1: length(df_list)) {
  model5[[i]] <- stan_glm(total ~ period + site + veg,
                          data = df_list[[i]],
                          family = neg_binomial_2(link="log"), init_r=1,
                          adapt_delta = 0.999
                          # prior_intercept = normal(0,10), prior = normal(0,10)
  )
}
print(model5)

#version without season
model5 = list()
for(i in 1: length(df_list)) {
  model5[[i]] <- stan_glm(total ~ site + veg,
                          data = df_list[[i]],
                          family = neg_binomial_2(link="log"), init_r=1,
                          adapt_delta = 0.999, 
                          iter=5000, warmup=3000
                          # prior_intercept = normal(0,10), prior = normal(0,10)
  )
}
print(model5)

#version with season
model6 = list()
for(i in 1: length(df_list)) {
  model6[[i]] <- stan_glm(total ~ period + site + ipa,
                          data = df_list[[i]],
                          family = neg_binomial_2(link="log"), init_r=1,
                          adapt_delta = 0.999
                          # prior_intercept = normal(0,10), prior = normal(0,10)
  )
}
print(model6)

#version without season 
model6 = list()
for(i in 1: length(df_list)) {
  model6[[i]] <- stan_glm(total ~ site + ipa,
                          data = df_list[[i]],
                          family = neg_binomial_2(link="log"), init_r=1,
                          adapt_delta = 0.999, 
                          iter=5000, warmup=3000
                          # prior_intercept = normal(0,10), prior = normal(0,10)
  )
}
print(model6)

model6_b = list()
for(i in 1: length(df_list)) {
  model6_b[[i]] <- stan_glm(total ~ site + ipa + ipa*dock,
                          data = df_list[[i]],
                          family = neg_binomial_2(link="log"), init_r=1,
                          adapt_delta = 0.999, 
                          iter=5000, warmup=3000
                          # prior_intercept = normal(0,10), prior = normal(0,10)
  )
}

#version with season
model7 = list()
for(i in 1: length(df_list)) {
  model7[[i]] <- stan_glm(total ~ period + site + ipa + veg,
                          data = df_list[[i]],
                          family = neg_binomial_2(link="log"), init_r=1,
                          adapt_delta = 0.999
                          # prior_intercept = normal(0,10), prior = normal(0,10)
  )
}
print(model7)

#version without season
model7 = list()
for(i in 1: length(df_list)) {
  model7[[i]] <- stan_glm(total ~ site + ipa + veg,
                          data = df_list[[i]],
                          family = neg_binomial_2(link="log"), init_r=1,
                          adapt_delta = 0.999, 
                          iter=5000, warmup=3000
                          # prior_intercept = normal(0,10), prior = normal(0,10)
  )
}
print(model7)

model7_b = list()
for(i in 1: length(df_list)) {
  model7_b[[i]] <- stan_glm(total ~ site + ipa + veg + ipa*dock,
                          data = df_list[[i]],
                          family = neg_binomial_2(link="log"), init_r=1,
                          adapt_delta = 0.999, 
                          iter=5000, warmup=3000
                          # prior_intercept = normal(0,10), prior = normal(0,10)
  )
}

#version wtih season
model8 = list()
for(i in 1: length(df_list)) {
  model8[[i]] <- stan_glm(total ~ period + ipa + veg + ipa*veg,
                          data = df_list[[i]],
                          family = neg_binomial_2(link="log"), init_r=1,
                          adapt_delta = 0.999
                          # prior_intercept = normal(0,10), prior = normal(0,10)
  )
}
print(model8)

#version wtihout season
model8 = list()
for(i in 1: length(df_list)) {
  model8[[i]] <- stan_glm(total ~ ipa + veg + ipa*veg,
                          data = df_list[[i]],
                          family = neg_binomial_2(link="log"), init_r=1,
                          adapt_delta = 0.999, 
                          iter=5000, warmup=3000
                          # prior_intercept = normal(0,10), prior = normal(0,10)
  )
}
print(model8)

#version with season
model9 = list()
for(i in 1: length(df_list)) {
  model9[[i]] <- stan_glm(total ~ period + site + ipa + site*ipa,
                          data = df_list[[i]],
                          family = neg_binomial_2(link="log"), init_r=1,
                          adapt_delta = 0.999
                          # prior_intercept = normal(0,10), prior = normal(0,10)
  )
}
print(model9)

#version without season
model9 = list()
for(i in 1: length(df_list)) {
  model9[[i]] <- stan_glm(total ~ site + ipa + site*ipa,
                          data = df_list[[i]],
                          family = neg_binomial_2(link="log"), init_r=1,
                          adapt_delta = 0.999, 
                          iter=5000, warmup=3000
                          # prior_intercept = normal(0,10), prior = normal(0,10)
  )
}
print(model9)

#version with season
model10 = list()
for(i in 1: length(df_list)) {
  model10[[i]] <- stan_glm(total ~ period + site + ipa + veg + veg*ipa,
                          data = df_list[[i]],
                          family = neg_binomial_2(link="log"), init_r=1,
                          adapt_delta = 0.999
                          # prior_intercept = normal(0,10), prior = normal(0,10)
  )
}
print(model10)

#version without season
model10 = list()
for(i in 1: length(df_list)) {
  model10[[i]] <- stan_glm(total ~ site + ipa + veg + veg*ipa,
                           data = df_list[[i]],
                           family = neg_binomial_2(link="log"), init_r=1,
                           adapt_delta = 0.999, 
                           iter=5000, warmup=3000
                           # prior_intercept = normal(0,10), prior = normal(0,10)
  )
}
print(model10)

# version with season
model11 = list()
for(i in 1: length(df_list)) {
  model11[[i]] <- stan_glm(total ~ period + site + ipa + veg + site*ipa,
                           data = df_list[[i]],
                           family = neg_binomial_2(link="log"), init_r=1,
                           adapt_delta = 0.999
                           # prior_intercept = normal(0,10), prior = normal(0,10)
  )
}
print(model11)

# version without season
model11 = list()
for(i in 1: length(df_list)) {
  model11[[i]] <- stan_glm(total ~ site + ipa + veg + site*ipa,
                           data = df_list[[i]],
                           family = neg_binomial_2(link="log"), init_r=1,
                           adapt_delta = 0.999, 
                           iter=5000, warmup=3000
                           # prior_intercept = normal(0,10), prior = normal(0,10)
  )
}
print(model11)

model11_b = list()
for(i in 1: length(df_list)) {
  model11_b[[i]] <- stan_glm(total ~ site + ipa + veg + site*ipa + dock,
                           data = df_list[[i]],
                           family = neg_binomial_2(link="log"), init_r=1,
                           adapt_delta = 0.999, 
                           iter=5000, warmup=3000
                           # prior_intercept = normal(0,10), prior = normal(0,10)
  )
}

#version with season
model12 = list()
for(i in 1: length(df_list)) {
  model12[[i]] <- stan_glm(total ~ period + site + ipa + veg + site*ipa + veg*ipa,
                           data = df_list[[i]],
                           family = neg_binomial_2(link="log"), init_r=1,
                           adapt_delta = 0.999
                           # prior_intercept = normal(0,10), prior = normal(0,10)
  )
}
print(model12)

#version without season
model12 = list()
for(i in 1: length(df_list)) {
  model12[[i]] <- stan_glm(total ~ site + ipa + veg + site*ipa + veg*ipa,
                           data = df_list[[i]],
                           family = neg_binomial_2(link="log"), init_r=1,
                           adapt_delta = 0.999, 
                           iter=5000, warmup=3000
                           # prior_intercept = normal(0,10), prior = normal(0,10)
  )
}
print(model12)


# Model comparison using LOO
# and Bayesian model weights with stacking
loo_model0 <- list()
loo_model1 <- list()
loo_model2 <- list()
loo_model3 <- list()
loo_model4 <- list()
loo_model5 <- list()
loo_model6 <- list()
loo_model6_b <- list()
loo_model7 <- list()
loo_model7_b <- list()
loo_model8 <- list()
loo_model9 <- list()
loo_model10 <- list()
loo_model11 <- list()
loo_model11_b <- list()
loo_model12 <- list()
loo_modela <- list()
loo_modela_2 <- list()
loo_modelb <- list()
loo_modelb_2 <- list()
loo_modelc <- list()
loo_modelc_2 <- list()
loo_model13 <- list()

for(i in 1: length(df_list)) {
#for(i in 4) {
loo_model0[[i]] <- loo(model0[[i]], save_psis = TRUE, k_threshold = 0.7) #}#loo = leave-one-out
loo_model1[[i]] <- loo(model1[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
loo_model2[[i]] <- loo(model2[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
loo_model3[[i]] <- loo(model3[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
loo_model4[[i]] <- loo(model4[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
loo_model5[[i]] <- loo(model5[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
loo_model6[[i]] <- loo(model6[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
loo_model6_b[[i]] <- loo(model6_b[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
loo_model7[[i]] <- loo(model7[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
loo_model7_b[[i]] <- loo(model7_b[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
loo_model8[[i]] <- loo(model8[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
loo_model9[[i]] <- loo(model9[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
loo_model10[[i]] <- loo(model10[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
loo_model11[[i]] <- loo(model11[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
loo_model11_b[[i]] <- loo(model11_b[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
loo_model12[[i]] <- loo(model12[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
loo_modela[[i]] <- loo(modela[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
loo_modela_2[[i]] <- loo(modela_2[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
loo_modelb[[i]] <- loo(modelb[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
loo_modelb_2[[i]] <- loo(modelb_2[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
loo_modelc[[i]] <- loo(modelc[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
loo_modelc_2[[i]] <- loo(modelc_2[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
loo_model13[[i]] <- loo(model13[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
}

# # MODELS w SEASON
# loo_list <- list()
# for(i in 4) { #length(df_list)) {
# loo_list[[i]] <- list(loo_model0[[i]], loo_model1[[i]], loo_model2[[i]], loo_model3[[i]], loo_model4[[i]], loo_model5[[i]],
#                   loo_model6[[i]], loo_model7[[i]], loo_model8[[i]], loo_model9[[i]], loo_model10[[i]],
#                   loo_model11[[i]], loo_model12[[i]])
# }


# MODELS EXCLUDING SEASON
loo_list <- list()
for(i in 1: length(df_list)) {
#for(i in 4) {
  loo_list[[i]] <- list(loo_modela[[i]], loo_modelb[[i]], loo_modelc[[i]], loo_model4[[i]], loo_model5[[i]],
                        loo_model6[[i]], loo_model7[[i]], loo_model7_b[[i]],loo_model8[[i]], loo_model9[[i]], loo_model10[[i]],
                        loo_model11[[i]], loo_model12[[i]],
                        loo_model13[[i]], loo_modela_2[[i]], loo_modelb_2[[i]],
                        loo_modelc_2[[i]], loo_model6_b[[i]], 
                        loo_model11_b[[i]])
}



# Bayesian stacked weights
stacking_wts <- list()
for(i in 1: length(df_list)) {
#  for(i in 4) {
  stacking_wts[[i]]<- loo_model_weights(loo_list[[i]])
}

capture.output(stacking_wts, file = "Bayesian model weights.txt")


# Compare loos with ELPD   
w <- loo_compare(loo_model0, loo_model1, loo_model2, loo_model3, loo_model4, loo_model5, loo_model6, 
                 loo_model7, loo_model8, loo_model9, loo_model10, loo_model11, loo_model12)
w

# Plot PSIS
loo_modelf <- loo_model12
plot(loo_modelf[[4]], label_points=TRUE)


# DIAGNOSTICS START HERE
# Diagnostic order for code:
#   Traceplots
#   MCMC pairs
#   PSIS
#   Over dispersion
#   Credible intervals
#   MCMC areas
#   Zeros
#   Kernel density estimates


# Assign model for diagnostics
modelf <- model12

# assign model-specific parameters
#model1
base_pars = c("(Intercept)",
              "vegAbsent")

#model3
base_pars = c("(Intercept)",
              "siteDOK", "siteEDG","siteFAM",
              "siteSHR","siteTUR")

#model6
base_pars = c("(Intercept)",
              "siteDOK", "siteEDG","siteFAM",
              "siteSHR","siteTUR", "ipaNatural", "ipaRestored")

#model9
base_pars = c("(Intercept)",
              "siteDOK", "siteEDG","siteFAM",
              "siteSHR","siteTUR", "ipaNatural", "ipaRestored", "siteDOK:ipaNatural",
              "siteEDG:ipaNatural", "siteFAM:ipaNatural", "siteSHR:ipaNatural",
              "siteTUR:ipaNatural", "siteDOK:ipaRestored", "siteEDG:ipaRestored",
              "siteFAM:ipaRestored","siteSHR:ipaRestored", "siteTUR:ipaRestored")

traceplot1 <- list()
for(i in 1: length(df_list))
{
traceplot1[[i]] <- traceplot(modelf[[i]]$stanfit, pars = base_pars) #test model divergence
plot(traceplot1[[i]])
}

# plot credible intervals
ci.plot <- list()
posterior1 <- modelf
for(i in 1:length(df_list)) {
ci.plot[[i]] <- plot(posterior1[[i]])
}

# Density overlay of predicted vs observed
y1 = list()
y_rep1 = list()
plot_ppc1 = list()

for(i in 1: length(df_list))
{
  df = as.data.frame(df_list[[i]])
  NAME <- unique(df$species)

    # define actual observations:
  y1[[i]] <- df$total #species count
  #define draws from posterior predictive
  y_rep1[[i]] <- posterior_predict(modelf[[i]], draws = 1000)

  color_scheme_set("brightblue")
  plot_ppc1[[i]] <- ppc_dens_overlay(y1[[i]], y_rep1[[i]][1:100, ]) +
    coord_cartesian(xlim = c(0, 100)) +
    ggtitle(paste(as.character("Smelt")))

  plot(plot_ppc1[[i]]) 
}

# Histogram of obs vs predicted
ppc_hist(y1[[i]], y_rep1[[i]][1:5, ])


# Posterior probability distributions
mcmc_areas <- list()
for(i in 1: length(df_list))
{
mcmc_areas[[i]] <- mcmc_areas(as.matrix(modelf[[i]], prob_outer=.95))
plot(mcmc_areas[[i]])
}

# Plot histogram posterior probability distributions
mcmc_histo <- list()
for(i in 1: length(df_list))
{
mcmc_histo <- mcmc_hist(as.matrix(modelf[[4]]), pars=c("(Intercept)", "siteDOK", "siteEDG",
                                                        "siteFAM", "siteSHR", "siteTUR"), 
                                    binwidth=1)
  plot(mcmc_histo)
}

  
mcmc_histo <- mcmc_hist(as.matrix(model6[[3]]), pars=c("(Intercept)", "ipaNatural", 
                          "ipaRestored"), binwidth=1/2) + ggtitle("Herring") +
                          geom_vline(aes(xintercept = 0), linetype = "dashed")
  
mcmc_histo2 <- mcmc_hist(as.matrix(model9[[4]]), pars=c("(Intercept)", 
                          "ipaNatural", "ipaRestored"), binwidth=1/2) + ggtitle("Surf Smelt")+
                          geom_vline(aes(xintercept = 0), linetype = "dashed")
  
    ggarrange(mcmc_histo, mcmc_histo2, ncol=2, nrow=1,common.legend = TRUE 
              )

# Cross Validation Checking
# Pareto-smoothed importance sampling  leave-one-out cross-validation as a model checking tool.
# ID highly influential observations, indicating model misspecification (if too many)

# use above-estimated loo
loo_modelf <- list()
for(i in 1: length(df_list)) 
{
  loo_modelf[[i]] <- loo(modelf[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
  plot(loo_modelf[[i]], label_points=TRUE)
}


# posterior distribution of over-dispersion parameter
mcmc_areas2 <- list()
for(i in 1: length(df_list))
{
mcmc_areas2[[i]] <- mcmc_areas(as.matrix(modelf[[i]]), prob_outer = .999,
           pars = c("reciprocal_dispersion"))
}
print(mcmc_areas2)


# test posterior dependencies
#predictors should not be correlated
mcmc_pairs <- list()
for(i in 1: length(df_list))
{
  mcmc_pairs[[i]] <- mcmc_pairs(as.matrix(modelf[[i]]), pars=base_pars)
  plot(mcmc_pairs[[i]])
}


# Evaluate ability of model to predict the number of zeros in the data
prop_zero <- function(x) mean(x==0)

prop_zero_test <- list()
for(i in 1: length(df_list))
{
prop_zero_test[[i]] <- pp_check(modelf[[i]], plotfun= "stat", stat = "prop_zero")
}
prop_zero_test


# Evaluate ability of model to predict the maximum value(s)
max_test_nb <- pp_check(modelf[[i]], plotfun = "stat", stat = "max", binwidth=200) +
  coord_cartesian(xlim = c(-1, 1000))
max_test_nb


# Evaluate ability of model to predict the minimum value(s)
min_test_nb <- pp_check(modelf[[i]], plotfun = "stat", stat = "min", binwidth=200) +
  coord_cartesian(xlim = c(-1, 1000))
min_test_nb


# Check proportion of posteriors that are positive
prop_pos <- function(x) mean(x>0)
prop_pos(h.df$ipaRestored)


#LOO_PIT
loo_pit_me <- list()
loo_modelf <- loo_model12

for(i in 1: length(df_list))
{
  # define actual observations:
  df = as.data.frame(df_list[[i]])
  y1 <- df$total #species count
 
  # define draws from posterior predictive
  y_rep <- posterior_predict(modelf[[i]])
  
  # define loo weights
  lw <- weights(loo_modelf[[i]]$psis_object)

  loo_pit_me[[i]] <- ppc_loo_pit_overlay(y1, y_rep, lw)
}
loo_pit_me
