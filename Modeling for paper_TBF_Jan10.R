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
library(googlesheets4)
library(googledrive)
library(testthat)
library(emmeans)
library(modelbased)
library(ggpubr)
library(bridgesampling)
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
 # mutate(season = recode(month, '4' = "Early",
  #                       '5'="Early",
   #                      '6'="Early",
    #                     '7'="Late",
     #                    '8'="Late",
      #                   '9'="Late")) %>%
  mutate(period = case_when(
    month > 6 & year > 2018 ~ "period4",
    month < 7  & year > 2018 ~ "period3",
    month > 6 & year < 2019 ~ "period2",
    TRUE ~ "period1")) %>%
  filter(!total%in% c(188,357,592))  




# prepare to estimate models
col.filters <- unique(df$species) 

lapply(seq_along(col.filters), function(x) {
  filter(df, species == col.filters[x])
}
) -> df_list

names(df_list) <- col.filters

# #if want to limit data to only include sites where species was observed
# # smelt: COR, DOK, EDG, FAM
# df_sm <- df_list[[4]]
# df_sm <- df_sm[!(df_sm$site %in% c("SHR","TUR")),]
# df_list[[4]] <- df_sm
# 
# #  herring: COR, DOK, TUR, FAM
# df_he <- df_list[[3]]
# df_he <- df_he[!(df_he$site %in% c("SHR","EDG")),]
# df_list[[3]] <- df_he
# 
# #  herring: COR, FAM
# df_he2 <- df_list[[3]]
# df_he2 <- df_he2[!(df_he2$site %in% c("SHR","EDG","DOK", "TUR")),]
# df_list[[3]] <- df_he2


# List of models to estimate 
# also versions of each model (except a-c/1-3) excluding season
# model0  count~period
# modela  count~site
# modelb  count~veg
# modelc  count~ipa
# model4  count~ipa + veg
# model5  count~site + veg
# model6  count~site + ipa
# model7  count~site + ipa + veg
# model8  count~ipa + veg + veg*ipa
# model9  count~site + ipa + site*ipa
# model10  count~site + ipa + veg + veg*ipa
# model11  count~site + ipa + veg + site*ipa
# model12  count~site + ipa + veg + site*ipa + veg*ipa


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

modelb = list()
for(i in 1: length(df_list)) {
  modelb[[i]] <- stan_glm(total ~ veg,
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

#model 1 not needed if we exclude season from our models
# model1 = list()
# for(i in 1: length(df_list)) {
#   model1[[i]] <- stan_glm(total ~ period + veg,
#                           data = df_list[[i]],
#                           family = neg_binomial_2(link="log"), init_r=1,
#                           adapt_delta = 0.999,
#                           diagnostic_file = file.path(tempdir(), "df.csv")
#                           # prior_intercept = normal(0,10), prior = normal(0,10)
#   )
# }
# print(model1 )
# 
# 
# bayes_factor(bridge_sampler(model0[[1]]), bridge_sampler(model1[[1]]))


# #model 2 not needed if we exclude season from our models
# model2 = list()
# for(i in 1: length(df_list)) {
#   model2[[i]] <- stan_glm(total ~ period + ipa,
#                           data = df_list[[i]],
#                           family = neg_binomial_2(link="log"), init_r=1,
#                           adapt_delta = 0.999
#                           # prior_intercept = normal(0,10), prior = normal(0,10)
#   )
# }
# print(model2)
# 
# #model 3 not needed if we exclude season from our models
# model3 = list()
# for(i in 1: length(df_list)) {
#   model3[[i]] <- stan_glm(total ~ period + site,
#                           data = df_list[[i]],
#                           family = neg_binomial_2(link="log"), init_r=1,
#                           adapt_delta = 0.999
#                           # prior_intercept = normal(0,10), prior = normal(0,10)
#   )
# }
# print(model3)

# # version with season
# model4 = list()
# for(i in 1: length(df_list)) {
#   model4[[i]] <- stan_glm(total ~ period + ipa + veg,
#                           data = df_list[[i]],
#                           family = neg_binomial_2(link="log"), init_r=1,
#                           adapt_delta = 0.999
#                           # prior_intercept = normal(0,10), prior = normal(0,10)
#   )
# }
# print(model4)

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

# #version with season
# model5 = list()
# for(i in 1: length(df_list)) {
#   model5[[i]] <- stan_glm(total ~ period + site + veg,
#                           data = df_list[[i]],
#                           family = neg_binomial_2(link="log"), init_r=1,
#                           adapt_delta = 0.999
#                           # prior_intercept = normal(0,10), prior = normal(0,10)
#   )
# }
# print(model5)

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

# #version with season
# model6 = list()
# for(i in 1: length(df_list)) {
#   model6[[i]] <- stan_glm(total ~ period + site + ipa,
#                           data = df_list[[i]],
#                           family = neg_binomial_2(link="log"), init_r=1,
#                           adapt_delta = 0.999
#                           # prior_intercept = normal(0,10), prior = normal(0,10)
#   )
# }
# print(model6)

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

# #version with season
# model7 = list()
# for(i in 1: length(df_list)) {
#   model7[[i]] <- stan_glm(total ~ period + site + ipa + veg,
#                           data = df_list[[i]],
#                           family = neg_binomial_2(link="log"), init_r=1,
#                           adapt_delta = 0.999
#                           # prior_intercept = normal(0,10), prior = normal(0,10)
#   )
# }
# print(model7)

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

# #version wtih season
# model8 = list()
# for(i in 1: length(df_list)) {
#   model8[[i]] <- stan_glm(total ~ period + ipa + veg + ipa*veg,
#                           data = df_list[[i]],
#                           family = neg_binomial_2(link="log"), init_r=1,
#                           adapt_delta = 0.999
#                           # prior_intercept = normal(0,10), prior = normal(0,10)
#   )
# }
# print(model8)

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

# #version with season
# model9 = list()
# for(i in 1: length(df_list)) {
#   model9[[i]] <- stan_glm(total ~ period + site + ipa + site*ipa,
#                           data = df_list[[i]],
#                           family = neg_binomial_2(link="log"), init_r=1,
#                           adapt_delta = 0.999
#                           # prior_intercept = normal(0,10), prior = normal(0,10)
#   )
# }
# print(model9)

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

# #version with season
# model10 = list()
# for(i in 1: length(df_list)) {
#   model10[[i]] <- stan_glm(total ~ period + site + ipa + veg + veg*ipa,
#                           data = df_list[[i]],
#                           family = neg_binomial_2(link="log"), init_r=1,
#                           adapt_delta = 0.999
#                           # prior_intercept = normal(0,10), prior = normal(0,10)
#   )
# }
# print(model10)

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

# # version with season
# model11 = list()
# for(i in 1: length(df_list)) {
#   model11[[i]] <- stan_glm(total ~ period + site + ipa + veg + site*ipa,
#                            data = df_list[[i]],
#                            family = neg_binomial_2(link="log"), init_r=1,
#                            adapt_delta = 0.999
#                            # prior_intercept = normal(0,10), prior = normal(0,10)
#   )
# }
# print(model11)

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

# #version with season
# model12 = list()
# for(i in 1: length(df_list)) {
#   model12[[i]] <- stan_glm(total ~ period + site + ipa + veg + site*ipa + veg*ipa,
#                            data = df_list[[i]],
#                            family = neg_binomial_2(link="log"), init_r=1,
#                            adapt_delta = 0.999
#                            # prior_intercept = normal(0,10), prior = normal(0,10)
#   )
# }
# print(model12)

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


# 
# capture model output of models with most Bayesian weight
capture.output(summary(model0[[4]], digits=2, prob=c(.025, .05, .055, .1, .9, .945, .95, .975)), file = "summary_BestModelCoeffs_Jan11_mean_allCIs_All2.txt")
capture.output(summary(modela[[4]], digits=2, prob=c(.025, .05, .055, .1, .9, .945, .95, .975)), file = "summary_BestModelCoeffs_Jan11_mean_allCIs_All2.txt", append = TRUE)
capture.output(summary(modelb[[4]], digits=2, prob=c(.025, .05, .055, .1, .9, .945, .95, .975)), file = "summary_BestModelCoeffs_Jan11_mean_allCIs_All2.txt", append = TRUE)
capture.output(summary(modelc[[4]], digits=2, prob=c(.025, .05, .055, .1, .9, .945, .95, .975)), file = "summary_BestModelCoeffs_Jan11_mean_allCIs_All2.txt", append = TRUE)
capture.output(summary(model4[[4]], digits=2, prob=c(.025, .05, .055, .1, .9, .945, .95, .975)), file = "summary_BestModelCoeffs_Jan11_mean_allCIs_All2.txt", append = TRUE)
capture.output(summary(model5[[4]], digits=2, prob=c(.025, .05, .055, .1, .9, .945, .95, .975)), file = "summary_BestModelCoeffs_Jan11_mean_allCIs_All2.txt", append = TRUE)
capture.output(summary(model6[[4]], digits=2, prob=c(.025, .05, .055, .1, .9, .945, .95, .975)), file = "summary_BestModelCoeffs_Jan11_mean_allCIs_All2.txt", append = TRUE)
capture.output(summary(model7[[4]], digits=2, prob=c(.025, .05, .055, .1, .9, .945, .95, .975)), file = "summary_BestModelCoeffs_Jan11_mean_allCIs_All2.txt", append = TRUE)
capture.output(summary(model8[[4]], digits=2, prob=c(.025, .05, .055, .1, .9, .945, .95, .975)), file = "summary_BestModelCoeffs_Jan11_mean_allCIs_All2.txt", append = TRUE)
capture.output(summary(model9[[4]], digits=2, prob=c(.025, .05, .055, .1, .9, .945, .95, .975)), file = "summary_BestModelCoeffs_Jan11_mean_allCIs_All2.txt", append = TRUE)
capture.output(summary(model10[[4]], digits=2, prob=c(.025, .05, .055, .1, .9, .945, .95, .975)), file = "summary_BestModelCoeffs_Jan11_mean_allCIs_All2.txt", append = TRUE)
capture.output(summary(model11[[4]], digits=2, prob=c(.025, .05, .055, .1, .9, .945, .95, .975)), file = "summary_BestModelCoeffs_Jan11_mean_allCIs_All2.txt", append = TRUE)
capture.output(summary(model12[[4]], digits=2, prob=c(.025, .05, .055, .1, .9, .945, .95, .975)), file = "summary_BestModelCoeffs_Jan11_mean_allCIs_All2.txt", append = TRUE)


#save 89% credible intervals for all models, all species
capture.output(posterior_interval(model12[[1]], prob=0.89), file = "posterior_interval_ModelCoeffs_89_model0.txt", append = TRUE)
capture.output(posterior_interval(model12[[2]], prob=0.89), file = "posterior_interval_ModelCoeffs_89_model0.txt", append = TRUE)
capture.output(posterior_interval(model12[[3]], prob=0.89), file = "posterior_interval_ModelCoeffs_89_model0.txt", append = TRUE)
capture.output(posterior_interval(model12[[4]], prob=0.89), file = "posterior_interval_ModelCoeffs_89_model0.txt", append = TRUE)


# capture.output(emmip(model3[[1]], site ~ period, CIs=TRUE, plotit=FALSE),file = "Emmips_Nov20.txt")
# capture.output(emmip(model1[[2]], veg ~ period , CIs=TRUE, plotit=FALSE),file = "Emmips_Nov20.txt", append=TRUE)
# capture.output(emmip(model6[[3]], ipa ~ site, CIs=TRUE, plotit=FALSE),file = "Emmips_Nov20.txt", append=TRUE)
# capture.output(emmip(model9[[4]], ipa ~ site, CIs=TRUE, plotit=FALSE),file = "Emmips_Nov20.txt", append=TRUE)


# DIAGNOSTICS START HERE
# ORDER, FROM SUPPLEMENT
# [Note: diagnostic order for code:
#   Traceplots
#   MCMC pairs
#   PSIS
#   Over dispersion
#   Credible intervals
#   MCMC areas
#   Zeros
#   Kernel density estimates
  


# Model comparison using LOO
# and Bayesian model weights with stacking

loo_model0 <- list()
#loo_model1 <- list()
#loo_model2 <- list()
#loo_model3 <- list()
loo_model4 <- list()
loo_model5 <- list()
loo_model6 <- list()
loo_model7 <- list()
loo_model8 <- list()
loo_model9 <- list()
loo_model10 <- list()
loo_model11 <- list()
loo_model12 <- list()
loo_modela <- list()
loo_modelb <- list()
loo_modelc <- list()

for(i in 1: length(df_list)) {
#for(i in 4) {
loo_model0[[i]] <- loo(model0[[i]], save_psis = TRUE, k_threshold = 0.7) #}#loo = leave-one-out
#loo_model1[[i]] <- loo(model1[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
#loo_model2[[i]] <- loo(model2[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
#loo_model3[[i]] <- loo(model3[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
loo_model4[[i]] <- loo(model4[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
loo_model5[[i]] <- loo(model5[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
loo_model6[[i]] <- loo(model6[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
loo_model7[[i]] <- loo(model7[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
loo_model8[[i]] <- loo(model8[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
loo_model9[[i]] <- loo(model9[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
loo_model10[[i]] <- loo(model10[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
loo_model11[[i]] <- loo(model11[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
loo_model12[[i]] <- loo(model12[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
loo_modela[[i]] <- loo(modela[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
loo_modelb[[i]] <- loo(modelb[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
loo_modelc[[i]] <- loo(modelc[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
}

# ALL MODELS w SEASON
# loo_list <- list()
# for(i in 4) { #length(df_list)) {
# loo_list[[i]] <- list(loo_model0[[i]], loo_model1[[i]], loo_model2[[i]], loo_model3[[i]], loo_model4[[i]], loo_model5[[i]],
#                   loo_model6[[i]], loo_model7[[i]], loo_model8[[i]], loo_model9[[i]], loo_model10[[i]],
#                   loo_model11[[i]], loo_model12[[i]])
# }


# LOO LIST EXCLUDING SEASON
loo_list <- list()
for(i in 1: length(df_list)) {
#for(i in 4) {
  loo_list[[i]] <- list(loo_model0[[i]], loo_modela[[i]], loo_modelb[[i]], loo_modelc[[i]], loo_model4[[i]], loo_model5[[i]],
                        loo_model6[[i]], loo_model7[[i]], loo_model8[[i]], loo_model9[[i]], loo_model10[[i]],
                        loo_model11[[i]], loo_model12[[i]])
}



# Bayesian stacked weights
stacking_wts <- list()
for(i in 1: length(df_list)) {
#  for(i in 4) {
  stacking_wts[[i]]<- loo_model_weights(loo_list[[i]])
}


capture.output(stacking_wts, file = "Bayesian model weights_AllSp_noseason_Jan11_setseed2.txt")


# Compare with ELPD   
w <- loo_compare(loo_model0, loo_model1, loo_model2, loo_model3, loo_model4, loo_model5, loo_model6, 
                 loo_model7, loo_model8, loo_model9, loo_model10, loo_model11, loo_model12)
w

# Plot PSIS
loo_modelf <- loo_model12
plot(loo_modelf[[4]], label_points=TRUE)


# DIAGNOSTICS START HERE
# ORDER, FROM SUPPLEMENT
# [Note: diagnostic order for code:
#   Traceplots
#   MCMC pairs
#   PSIS
#   Over dispersion
#   Credible intervals
#   MCMC areas
#   Zeros
#   Kernel density estimates


# CHECK MODEL FITS
modelf <- model12


#TRACEPLOT
# adjust for model of interest
# Best models 9/4/2020
# Chinook, herring: model9
# Chum: model1
# Surf smelt: model11

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
#
#model9
base_pars = c("(Intercept)",
              "siteDOK", "siteEDG","siteFAM",
              "siteSHR","siteTUR", "ipaNatural", "ipaRestored", "siteDOK:ipaNatural",
              "siteEDG:ipaNatural", "siteFAM:ipaNatural", "siteSHR:ipaNatural",
              "siteTUR:ipaNatural", "siteDOK:ipaRestored", "siteEDG:ipaRestored",
              "siteFAM:ipaRestored","siteSHR:ipaRestored", "siteTUR:ipaRestored")
#model11
base_pars = c("(Intercept)",
              "siteDOK", "siteEDG","siteFAM", "siteSHR",
             "siteTUR", "ipaNatural", "ipaRestored", "vegAbsent","siteDOK:ipaNatural",
              "siteEDG:ipaNatural", "siteFAM:ipaNatural", "siteSHR:ipaNatural",
              "siteTUR:ipaNatural", "siteDOK:ipaRestored", "siteEDG:ipaRestored",
              "siteFAM:ipaRestored","siteSHR:ipaRestored", "siteTUR:ipaRestored")

#model12
base_pars = c("(Intercept)",
              "siteDOK", "siteEDG","siteFAM", "siteSHR",
              "siteTUR", "ipaNatural", "ipaRestored", "vegAbsent","siteDOK:ipaNatural",
              "siteEDG:ipaNatural", "siteFAM:ipaNatural", "siteSHR:ipaNatural",
              "siteTUR:ipaNatural", "siteDOK:ipaRestored", "siteEDG:ipaRestored",
              "siteFAM:ipaRestored","siteSHR:ipaRestored", "siteTUR:ipaRestored",
              "ipaNatural:vegAbsent", "ipaRestored:vegAbsent")

modelf<-model12
traceplot1 <- list()
for(i in 1: length(df_list))
{
traceplot1[[i]] <- traceplot(modelf[[i]]$stanfit, pars = base_pars) #test model divergence
plot(traceplot1[[i]])
}#[[i]]) # view plots


# plot credible intervals
posterior1 <- list()
posterior1 <- modelf
for(i in 1:length(df_list)) {
plot(posterior1[[i]])
}

#predicted vs observed
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
  plot_ppc1[[i]] <- ppc_dens_overlay(y1[[i]], y_rep1[[i]][1:1000, ]) +
    coord_cartesian(ylim = c(0, 1), xlim = c(0, 100000)) +
    ggtitle(paste(as.character("Smelt")))

  plot(plot_ppc1[[i]]) # view plots
}

# Analyze posterior probability distributions
mcmc_areas <- list()
for(i in 1: length(df_list))
{
mcmc_areas[[i]] <- mcmc_areas(as.matrix(modelf[[i]], prob_outer=.95))
plot(mcmc_areas[[i]])
}

# Plot histogram posterior probability distributions
modelf <- model6
#mcmc_histo <- list()
#for(i in 1: length(df_list))
#{
#  df = as.data.frame(df_list[[i]])
#  NAME <- unique(df$species)
mcmc_histo <- mcmc_hist(as.matrix(modelf[[4]]), pars=c("(Intercept)", "siteDOK", "siteEDG",
                                                        "siteFAM", "siteSHR", "siteTUR"), 
                                    binwidth=1)
  plot(mcmc_histo)
#}

  
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
# ID highly influencial observations, indicating model misspecification (if too many)

# use above-estimated loo
loo_modelf <- list()
for(i in 1: length(df_list)) 
{
  loo_modelf[[i]] <- loo(modelf[[i]], save_psis = TRUE, k_threshold = 0.7) #loo = leave-one-out
  plot(loo_modelf[[i]], label_points=TRUE)
}

# # other options
# pp_check(modelf, plotfun = "boxplot", nreps = 10, notch = FALSE)
# pp_check(modelf, plotfun = "stat_grouped", stat = "median", group = "site")
# pp_check(modelf, plotfun = "stat_grouped", stat = "median", group = "veg")
# pp_check(modelf, plotfun = "scatter_avg")
# pp_check(modelf, plotfun = "scatter", nreps = 3)
# 


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

 
# See if model can predict the number of zeros in the data
prop_zero_test <- list()
for(i in 1: length(df_list))
{
prop_zero_test[[i]] <- pp_check(modelf[[i]], plotfun= "stat", stat = "prop_zero")
}
prop_zero_test

 # Evaluate ability of model to predict the maximum value(s)
max_test_nb <- pp_check(modelf[[i]], plotfun = "stat", stat = "max", binwidth=10)+
  coord_cartesian(xlim = c(0, 100))
max_test_nb


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



# CONSIDER MOVING TO NEW R SCRIPT - PLOTS, MAYBE EVEN ONE PER PLOT 
# #######
 # PLOTS
 ```{r NEG BINOMIAL PLOTS}
# First, create lists/dataframes for use in building plots

site_perd = list() # site + period 
veg_perd = list() # veg + period 
site_ipa_perd = list() #site + ipa + period -
smelt_model = list() #site + ipa + veg + period + site*ipa 
site = list()
site2 = list()
site_plus_ipa = list()
site_plus_ipa2 = list()
ipa = list()
veg_plus_ipa=list()
veg_ipa=list()
site_ipa = list()
site_ipa2 = list()
veg = list()
veg2 = list()

#site + period
i=1             # if single species
#fits <- c(1,2) # if multiple species
#for(i in fits) # if multiple species
#{
  df = as.data.frame(df_list[[i]]) # to id each species
  NAME <- unique(df$species)
  
  site_perd[[i]] <-  as.data.frame(model3[[i]]) %>%
    transmute(DOK.Per1 = `(Intercept)` + siteDOK,
              DOK.Per2 = `(Intercept)` + siteDOK + periodperiod2,
              DOK.Per3 = `(Intercept)` + siteDOK + periodperiod3,
              DOK.Per4 = `(Intercept)` + siteDOK + periodperiod4,
              EDG.Per1 = `(Intercept)` + siteEDG,
              EDG.Per2 = `(Intercept)` + siteEDG + periodperiod2,
              EDG.Per3 = `(Intercept)` + siteEDG + periodperiod3,
              EDG.Per4 = `(Intercept)` + siteEDG + periodperiod4,
              FAM.Per1 = `(Intercept)` + siteFAM,
              FAM.Per2 = `(Intercept)` + siteFAM + periodperiod2,
              FAM.Per3 = `(Intercept)` + siteFAM + periodperiod3,
              FAM.Per4 = `(Intercept)` + siteFAM + periodperiod4,
              SHR.Per1 = `(Intercept)` + siteSHR,
              SHR.Per2 = `(Intercept)` + siteSHR + periodperiod2,
              SHR.Per3 = `(Intercept)` + siteSHR + periodperiod3,
              SHR.Per4 = `(Intercept)` + siteSHR + periodperiod4,
              TUR.Per1 = `(Intercept)` + siteTUR,
              TUR.Per2 = `(Intercept)` + siteTUR + periodperiod2,
              TUR.Per3 = `(Intercept)` + siteTUR + periodperiod3,
              TUR.Per4 = `(Intercept)` + siteTUR + periodperiod4,
              COR.Per1 = `(Intercept)`,
              COR.Per2 = `(Intercept)` + periodperiod2,
              COR.Per3 = `(Intercept)` + periodperiod3,
              COR.Per4 = `(Intercept)` + periodperiod4)%>%
    gather()%>%
    group_by(key) %>%
    mean_hdi(.width=0.89) %>% # hdi = 0.89 standard Bayesian CIs
    mutate_if(is.double, round, digits = 2) %>%
    mutate(median = exp(value)) %>%
    mutate(.lower = exp(.lower)) %>%
    mutate(.upper = exp(.upper)) %>%
    separate(key, c("site", "period")) %>%
    select(-c("value")) %>%
    mutate(species = paste(as.character("Chinook")))
  
write.csv(site_perd[[1]], "Chinook_bestmodel_transmuteoutput.csv")
  
# veg + period
# best model for Chum
#   fits <- c(1,2) #to select certain species
#  for(i in fits)
#  {
i=2
    df = as.data.frame(df_list[[i]]) # to id each species
    NAME <- unique(df$species)
    veg_perd[[i]] <-  as.data.frame(model1[[i]]) %>%
      transmute(Present.Period1 = `(Intercept)`,
                Present.Period2 = `(Intercept)` + periodperiod2,
                Present.Period3 = `(Intercept)` + periodperiod3,
                Present.Period4 = `(Intercept)` + periodperiod4,
                Absent.Period1 = `(Intercept)` + vegAbsent,
                Absent.Period2 = `(Intercept)` + vegAbsent + periodperiod2,
                Absent.Period3 = `(Intercept)` + vegAbsent + periodperiod3,
                Absent.Period4 = `(Intercept)` + vegAbsent + periodperiod4)%>%
      gather() %>%
      group_by(key) %>%
      mean_hdi(.width=0.89) %>%  # hdi = 0.89 standard Bayesian CIs
      mutate_if(is.double, round, digits = 2) %>%
      #mutate(key)
      mutate(median = exp(value)) %>%
      mutate(.lower = exp(.lower)) %>%
      mutate(.upper = exp(.upper)) %>%
      separate(key, c("veg", "period")) %>%
      select(-c("value")) %>%
      mutate(species = paste(as.character(NAME)))

write.csv(veg_perd[[2]], file = "Chum_bestmodel_transmuteoutput.csv")
    

# site + ipa + period
# best model for Herring
# VERSION EXCLUDING EDG & SHR WHERE HERRING ARE NEVER FOUND
NAME = c("Herring")
i=3
site_ipa_perd[[i]] <-  as.data.frame(model6[[i]]) %>%
  transmute(Per1.DOK.Natural = `(Intercept)` + siteDOK + ipaNatural,
            Per1.DOK.Armored = `(Intercept)` + siteDOK ,
            Per1.DOK.Restored = `(Intercept)` + siteDOK +ipaRestored,
            Per1.FAM.Natural = `(Intercept)` + siteFAM +ipaNatural,
            Per1.FAM.Armored = `(Intercept)` + siteFAM ,
            Per1.FAM.Restored = `(Intercept)` + siteFAM +ipaRestored,
            Per1.TUR.Natural = `(Intercept)` + siteTUR +ipaNatural,
            Per1.TUR.Armored = `(Intercept)` + siteTUR ,
            Per1.TUR.Restored = `(Intercept)` + siteTUR +ipaRestored,
            Per1.COR.Natural = `(Intercept)` + ipaNatural ,
            Per1.COR.Armored = `(Intercept)` ,
            Per1.COR.Restored = `(Intercept)` +ipaRestored,
            Per2.DOK.Natural = `(Intercept)` + siteDOK + ipaNatural + periodperiod2,
            Per2.DOK.Armored = `(Intercept)` + siteDOK + periodperiod2,
            Per2.DOK.Restored = `(Intercept)` + siteDOK +ipaRestored + periodperiod2,
            Per2.FAM.Natural = `(Intercept)` + siteFAM +ipaNatural + periodperiod2,
            Per2.FAM.Armored = `(Intercept)` + siteFAM + periodperiod2,
            Per2.FAM.Restored = `(Intercept)` + siteFAM +ipaRestored,
            Per2.TUR.Natural = `(Intercept)` + siteTUR +ipaNatural + periodperiod2,
            Per2.TUR.Armored = `(Intercept)` + siteTUR + periodperiod2,
            Per2.TUR.Restored = `(Intercept)` + siteTUR + ipaRestored + periodperiod2,
            Per2.COR.Natural = `(Intercept)` + ipaNatural + periodperiod2,
            Per2.COR.Armored = `(Intercept)` + periodperiod2,
            Per2.COR.Restored = `(Intercept)` +ipaRestored + periodperiod2,
            Per3.DOK.Natural = `(Intercept)` + siteDOK + ipaNatural + periodperiod3,
            Per3.DOK.Armored = `(Intercept)` + siteDOK + periodperiod3,
            Per3.DOK.Restored = `(Intercept)` + siteDOK +ipaRestored + periodperiod3,
            Per3.FAM.Natural = `(Intercept)` + siteFAM + ipaNatural + periodperiod3,
            Per3.FAM.Armored = `(Intercept)` + siteFAM + periodperiod3,
            Per3.FAM.Restored = `(Intercept)` + siteFAM + ipaRestored + periodperiod3,
            Per3.TUR.Natural = `(Intercept)` + siteTUR + ipaNatural + periodperiod3,
            Per3.TUR.Armored = `(Intercept)` + siteTUR + periodperiod3,
            Per3.TUR.Restored = `(Intercept)` + siteTUR + ipaRestored + periodperiod3,
            Per3.COR.Natural = `(Intercept)` + ipaNatural + periodperiod3,
            Per3.COR.Armored = `(Intercept)` + periodperiod3,
            Per3.COR.Restored = `(Intercept)` + ipaRestored + periodperiod3,
            Per4.DOK.Natural = `(Intercept)` + siteDOK + ipaNatural + periodperiod4,
            Per4.DOK.Armored = `(Intercept)` + siteDOK + periodperiod4,
            Per4.DOK.Restored = `(Intercept)` + siteDOK + ipaRestored + periodperiod4,
            Per4.FAM.Natural = `(Intercept)` + siteFAM + ipaNatural + periodperiod4,
            Per4.FAM.Armored = `(Intercept)` + siteFAM + periodperiod4,
            Per4.FAM.Restored = `(Intercept)` + siteFAM + ipaRestored + periodperiod4,
            Per4.TUR.Natural = `(Intercept)` + siteTUR + ipaNatural + periodperiod4,
            Per4.TUR.Armored = `(Intercept)` + siteTUR + periodperiod4,
            Per4.TUR.Restored = `(Intercept)` + siteTUR + ipaRestored + periodperiod4,
            Per4.COR.Natural = `(Intercept)` + ipaNatural + periodperiod4,
            Per4.COR.Armored = `(Intercept)` + periodperiod4,
            Per4.COR.Restored = `(Intercept)` + ipaRestored + periodperiod4) %>%
  gather()%>%
  group_by(key) %>%
  mean_hdi(.width=0.89) %>% # hdi = 0.89 is Bayesian standard
  mutate_if(is.double, round, digits = 2) %>%
  mutate(median = exp(value)) %>%
  mutate(.lower = exp(.lower)) %>%
  mutate(.upper = exp(.upper)) %>%
  separate(key, c("period", "site", "shoreline")) %>%
  select(-c("value")) %>%
  mutate(species = paste(as.character(NAME)))

write.csv(site_ipa_perd[[3]], file = "Herring_bestmodel_transmuteoutput.csv")


# period + site + ipa + site*ipa
# best model for smelt
#for(i in 1: length(df_list))
  #or to exclude species where model doesn't converge
  #fits <- c()
  #for(i in fits)
#{

  # WITH TUR & SHR REMOVED
  i=4 
  df = as.data.frame(df_list[[i]]) # to id each species
  #  NAME <- unique(df$species)
  NAME <- "Surf Smelt"   
  smelt_model[[i]] <-  as.data.frame(model9[[i]]) %>%
    transmute(Per1.COR.Natural = `(Intercept)` + ipaNatural ,
              Per1.COR.Armored = `(Intercept)` ,
              Per1.COR.Restored = `(Intercept)` +ipaRestored,
              Per2.COR.Natural = `(Intercept)` + ipaNatural + periodperiod2,
              Per2.COR.Armored = `(Intercept)`+ periodperiod2,
              Per2.COR.Restored = `(Intercept)` +ipaRestored + periodperiod2,
              Per3.COR.Natural = `(Intercept)` + ipaNatural + periodperiod3,
              Per3.COR.Armored = `(Intercept)` + periodperiod3,
              Per3.COR.Restored = `(Intercept)` +ipaRestored + periodperiod3,
              Per4.COR.Natural = `(Intercept)` + ipaNatural + periodperiod4,
              Per4.COR.Armored = `(Intercept)` + periodperiod4,
              Per4.COR.Restored = `(Intercept)` +ipaRestored + periodperiod4,
              Per1.DOK.Natural = `(Intercept)` + siteDOK + ipaNatural + `siteDOK:ipaNatural`,
              Per1.DOK.Armored = `(Intercept)` + siteDOK,
              Per1.DOK.Restored = `(Intercept)` + siteDOK + ipaRestored + `siteDOK:ipaRestored`,
              Per1.EDG.Natural = `(Intercept)` + siteEDG + ipaNatural + `siteEDG:ipaNatural`,
              Per1.EDG.Armored = `(Intercept)` + siteEDG,
              Per1.EDG.Restored = `(Intercept)` + siteEDG +ipaRestored + `siteEDG:ipaRestored`,
              Per1.FAM.Natural = `(Intercept)` + siteFAM +ipaNatural + `siteFAM:ipaNatural`,
              Per1.FAM.Armored = `(Intercept)` + siteFAM,
              Per1.FAM.Restored = `(Intercept)` + siteFAM +ipaRestored + `siteFAM:ipaRestored`,
              Per2.DOK.Natural = `(Intercept)` + siteDOK + ipaNatural + `siteDOK:ipaNatural` + periodperiod2,
              Per2.DOK.Armored = `(Intercept)` + siteDOK + periodperiod2,
              Per2.DOK.Restored = `(Intercept)` + siteDOK + ipaRestored + `siteDOK:ipaRestored` + periodperiod2,
              Per2.EDG.Natural = `(Intercept)` + siteEDG + ipaNatural + `siteEDG:ipaNatural` + periodperiod2,
              Per2.EDG.Armored = `(Intercept)` + siteEDG + periodperiod2,
              Per2.EDG.Restored = `(Intercept)` + siteEDG +ipaRestored + `siteEDG:ipaRestored` + periodperiod2,
              Per2.FAM.Natural = `(Intercept)` + siteFAM +ipaNatural + `siteFAM:ipaNatural` + periodperiod2,
              Per2.FAM.Armored = `(Intercept)` + siteFAM + periodperiod2,
              Per2.FAM.Restored = `(Intercept)` + siteFAM +ipaRestored + `siteFAM:ipaRestored` + periodperiod2,
              Per3.DOK.Natural = `(Intercept)` + siteDOK + ipaNatural + `siteDOK:ipaNatural` + periodperiod3,
              Per3.DOK.Armored = `(Intercept)` + siteDOK + periodperiod3,
              Per3.DOK.Restored = `(Intercept)` + siteDOK + ipaRestored + `siteDOK:ipaRestored` + periodperiod3,
              Per3.EDG.Natural = `(Intercept)` + siteEDG + ipaNatural + `siteEDG:ipaNatural` + periodperiod3,
              Per3.EDG.Armored = `(Intercept)` + siteEDG + periodperiod3,
              Per3.EDG.Restored = `(Intercept)` + siteEDG +ipaRestored + `siteEDG:ipaRestored` + periodperiod3,
              Per3.FAM.Natural = `(Intercept)` + siteFAM +ipaNatural + `siteFAM:ipaNatural` + periodperiod3,
              Per3.FAM.Armored = `(Intercept)` + siteFAM + periodperiod3,
              Per3.FAM.Restored = `(Intercept)` + siteFAM +ipaRestored + `siteFAM:ipaRestored` + periodperiod3,
              Per4.DOK.Natural = `(Intercept)` + siteDOK + ipaNatural + `siteDOK:ipaNatural` + periodperiod4,
              Per4.DOK.Armored = `(Intercept)` + siteDOK + periodperiod4,
              Per4.DOK.Restored = `(Intercept)` + siteDOK + ipaRestored + `siteDOK:ipaRestored` + periodperiod4,
              Per4.EDG.Natural = `(Intercept)` + siteEDG + ipaNatural + `siteEDG:ipaNatural` + periodperiod4,
              Per4.EDG.Armored = `(Intercept)` + siteEDG + periodperiod4,
              Per4.EDG.Restored = `(Intercept)` + siteEDG +ipaRestored + `siteEDG:ipaRestored` + periodperiod4,
              Per4.FAM.Natural = `(Intercept)` + siteFAM +ipaNatural + `siteFAM:ipaNatural` + periodperiod4,
              Per4.FAM.Armored = `(Intercept)` + siteFAM + periodperiod4,
              Per4.FAM.Restored = `(Intercept)` + siteFAM +ipaRestored + `siteFAM:ipaRestored` + periodperiod4) %>%
    gather()%>%
    group_by(key) %>%
    mean_hdi(.width=0.89) %>%
    mutate_if(is.double, round, digits = 2) %>%
    mutate(median = exp(value)) %>%
    mutate(.lower = exp(.lower)) %>%
    mutate(.upper = exp(.upper)) %>%
    separate(key, c("period", "site", "shoreline")) %>%
    select(-c("value")) %>%
    mutate(species = paste(as.character(NAME))) 

write.csv(smelt_model[[4]], file = "Smelt_bestmodel_transmuteoutput.csv")
  
  
# 
# # site + ipa + period
# # best model for Herring
# NAME = c("Herring")
# i=3
# site_ipa_perd[[i]] <-  as.data.frame(model6[[i]]) %>%
#   transmute(Per1.DOK.Natural = `(Intercept)` + siteDOK + ipaNatural,
#             Per1.DOK.Armored = `(Intercept)` + siteDOK ,
#             Per1.DOK.Restored = `(Intercept)` + siteDOK +ipaRestored,
#             Per1.EDG.Natural = `(Intercept)` + siteEDG +ipaNatural,
#             Per1.EDG.Armored = `(Intercept)` + siteEDG ,
#             Per1.EDG.Restored = `(Intercept)` + siteEDG +ipaRestored,
#             Per1.FAM.Natural = `(Intercept)` + siteFAM +ipaNatural,
#             Per1.FAM.Armored = `(Intercept)` + siteFAM ,
#             Per1.FAM.Restored = `(Intercept)` + siteFAM +ipaRestored,
#             Per1.SHR.Natural = `(Intercept)` + siteSHR +ipaNatural,
#             Per1.SHR.Armored = `(Intercept)` + siteSHR ,
#             Per1.SHR.Restored = `(Intercept)` + siteSHR +ipaRestored,
#             Per1.TUR.Natural = `(Intercept)` + siteTUR +ipaNatural,
#             Per1.TUR.Armored = `(Intercept)` + siteTUR ,
#             Per1.TUR.Restored = `(Intercept)` + siteTUR +ipaRestored,
#             Per1.COR.Natural = `(Intercept)` + ipaNatural ,
#             Per1.COR.Armored = `(Intercept)` ,
#             Per1.COR.Restored = `(Intercept)` +ipaRestored,
#             Per2.DOK.Natural = `(Intercept)` + siteDOK + ipaNatural + periodperiod2,
#             Per2.DOK.Armored = `(Intercept)` + siteDOK + periodperiod2,
#             Per2.DOK.Restored = `(Intercept)` + siteDOK +ipaRestored + periodperiod2,
#             Per2.EDG.Natural = `(Intercept)` + siteEDG +ipaNatural + periodperiod2,
#             Per2.EDG.Armored = `(Intercept)` + siteEDG + periodperiod2,
#             Per2.EDG.Restored = `(Intercept)` + siteEDG +ipaRestored + periodperiod2,
#             Per2.FAM.Natural = `(Intercept)` + siteFAM +ipaNatural + periodperiod2,
#             Per2.FAM.Armored = `(Intercept)` + siteFAM + periodperiod2,
#             Per2.FAM.Restored = `(Intercept)` + siteFAM +ipaRestored,
#             Per2.SHR.Natural = `(Intercept)` + siteSHR +ipaNatural + periodperiod2,
#             Per2.SHR.Armored = `(Intercept)` + siteSHR + periodperiod2,
#             Per2.SHR.Restored = `(Intercept)` + siteSHR +ipaRestored + periodperiod2,
#             Per2.TUR.Natural = `(Intercept)` + siteTUR +ipaNatural + periodperiod2,
#             Per2.TUR.Armored = `(Intercept)` + siteTUR + periodperiod2,
#             Per2.TUR.Restored = `(Intercept)` + siteTUR + ipaRestored + periodperiod2,
#             Per2.COR.Natural = `(Intercept)` + ipaNatural + periodperiod2,
#             Per2.COR.Armored = `(Intercept)` + periodperiod2,
#             Per2.COR.Restored = `(Intercept)` +ipaRestored + periodperiod2,
#             Per3.DOK.Natural = `(Intercept)` + siteDOK + ipaNatural + periodperiod3,
#             Per3.DOK.Armored = `(Intercept)` + siteDOK + periodperiod3,
#             Per3.DOK.Restored = `(Intercept)` + siteDOK +ipaRestored + periodperiod3,
#             Per3.EDG.Natural = `(Intercept)` + siteEDG + ipaNatural + periodperiod3,
#             Per3.EDG.Armored = `(Intercept)` + siteEDG + periodperiod3,
#             Per3.EDG.Restored = `(Intercept)` + siteEDG +ipaRestored + periodperiod3,
#             Per3.FAM.Natural = `(Intercept)` + siteFAM + ipaNatural + periodperiod3,
#             Per3.FAM.Armored = `(Intercept)` + siteFAM + periodperiod3,
#             Per3.FAM.Restored = `(Intercept)` + siteFAM + ipaRestored + periodperiod3,
#             Per3.SHR.Natural = `(Intercept)` + siteSHR + ipaNatural + periodperiod3,
#             Per3.SHR.Armored = `(Intercept)` + siteSHR + periodperiod3,
#             Per3.SHR.Restored = `(Intercept)` + siteSHR + ipaRestored + periodperiod3,
#             Per3.TUR.Natural = `(Intercept)` + siteTUR + ipaNatural + periodperiod3,
#             Per3.TUR.Armored = `(Intercept)` + siteTUR + periodperiod3,
#             Per3.TUR.Restored = `(Intercept)` + siteTUR + ipaRestored + periodperiod3,
#             Per3.COR.Natural = `(Intercept)` + ipaNatural + periodperiod3,
#             Per3.COR.Armored = `(Intercept)` + periodperiod3,
#             Per3.COR.Restored = `(Intercept)` + ipaRestored + periodperiod3,
#             Per4.DOK.Natural = `(Intercept)` + siteDOK + ipaNatural + periodperiod4,
#             Per4.DOK.Armored = `(Intercept)` + siteDOK + periodperiod4,
#             Per4.DOK.Restored = `(Intercept)` + siteDOK + ipaRestored + periodperiod4,
#             Per4.EDG.Natural = `(Intercept)` + siteEDG + ipaNatural + periodperiod4,
#             Per4.EDG.Armored = `(Intercept)` + siteEDG + periodperiod4,
#             Per4.EDG.Restored = `(Intercept)` + siteEDG + ipaRestored + periodperiod4,
#             Per4.FAM.Natural = `(Intercept)` + siteFAM + ipaNatural + periodperiod4,
#             Per4.FAM.Armored = `(Intercept)` + siteFAM + periodperiod4,
#             Per4.FAM.Restored = `(Intercept)` + siteFAM + ipaRestored + periodperiod4,
#             Per4.SHR.Natural = `(Intercept)` + siteSHR + ipaNatural + periodperiod4,
#             Per4.SHR.Armored = `(Intercept)` + siteSHR + periodperiod4,
#             Per4.SHR.Restored = `(Intercept)` + siteSHR +ipaRestored + periodperiod4,
#             Per4.TUR.Natural = `(Intercept)` + siteTUR + ipaNatural + periodperiod4,
#             Per4.TUR.Armored = `(Intercept)` + siteTUR + periodperiod4,
#             Per4.TUR.Restored = `(Intercept)` + siteTUR + ipaRestored + periodperiod4,
#             Per4.COR.Natural = `(Intercept)` + ipaNatural + periodperiod4,
#             Per4.COR.Armored = `(Intercept)` + periodperiod4,
#             Per4.COR.Restored = `(Intercept)` + ipaRestored + periodperiod4) %>%
#   gather()%>%
#   group_by(key) %>%
#   mean_hdi(.width=0.89) %>% # hdi = 0.89 is Bayesian standard
#   mutate_if(is.double, round, digits = 2) %>%
#   mutate(median = exp(value)) %>%
#   mutate(.lower = exp(.lower)) %>%
#   mutate(.upper = exp(.upper)) %>%
#   separate(key, c("period", "site", "shoreline")) %>%
#   select(-c("value")) %>%
#   mutate(species = paste(as.character(NAME)))

# 
# 
 ```{r NEG BINOMIAL FISH PLOT wsn}
#site + period
# best model for Chinook
#TO DO ON THIS FIGURE: (1) reverse the period order; (2) Remove "Chinook" title; (3) remove legend;
#(4) Reorder the sites (my mutate() code is not working).

  cols =c(brewer.pal(6,"Dark2"))
siteper <- bind_rows(site_perd) %>%
  mutate(site = as.factor(site)) %>%
  mutate(period = as.factor(period)) %>%
  mutate(site = recode(site, level= rev(c("FAM", "TUR", "COR", "SHR", "DOK", "EDG"))))
siteper_plot <- ggplot(siteper, aes(median,period)) +
  geom_point(aes(color = period)) +
  geom_errorbarh(aes(xmin = .lower, xmax=.upper,color=period),size =0.3) +
  scale_color_manual(values=rev(cols), name=" ")+
  guides(color = guide_legend(reverse = TRUE))+
  facet_grid(site~species, scales = "free_x")+
  # scale_x_continuous(limits = c(0,1))+
  xlab("Estimated Median Fish Count") +
  ylab("Period")+
  theme_classic() +
  theme(panel.spacing.y = unit(0,"line"), legend.position = "left")
siteper_plot



# veg + period
# best model for Chum
# TO DO FOR THIS PLOT: I don't think this wants to be a bar plot, b/c the 0.89 HDI calculates
# credible intervals, which aren't really the same as error bars. We want to be able to see the full
# 89% credible interval and compare those between eelgrass/no eelgrass sites.
my_color = c("snow3","black","snow3","black","snow3","black","snow3","black")

veg_only <- bind_rows(veg_perd)
veg_only_plot <- veg_only %>%
  mutate(Eelgrass = factor(veg, levels = c("Absent", "Present"))) %>%
  ggplot(aes(x=period, y=median, fill = Eelgrass)) +
  geom_bar(colour ="black",position="dodge", stat="identity")+
  scale_fill_manual(values = my_color)+
  geom_errorbar(aes(ymin = .lower, ymax=.upper, group = Eelgrass), width = 0.01, position=position_dodge(.9)) +
  facet_wrap(~species, scales = "free_y", ncol = 5) +
  ylab("Mean Fish Count") +
  xlab(" ")+
  theme_classic() +
  theme(strip.background =element_rect(fill="white"), axis.text.x = element_text(angle = 90)) 

veg_only_plot

#site + ipa + period
# best model for Herring 
cols =c(brewer.pal(2,"Paired"))
siteipa <- bind_rows(smelt_model[[4]])
#siteipa <- bind_rows(site_ipa_perd)
#siteipa_per1 <- siteipa[siteipa$period=="Per1",]
siteipa_plot <- ggplot(siteipa, aes(median, shoreline)) +
  geom_point(aes(color = shoreline)) +
  geom_errorbarh(aes(xmin = .lower, xmax=.upper,color=shoreline),width=0.1,size =0.3) +
  scale_color_manual(values=rev(cols), name=" "
                     #wes_palette(n=3, name="GrandBudapest1")
  )+
  guides(color = guide_legend(reverse = TRUE))+
  facet_grid(site~period, scales = "free_x")+
  # scale_x_continuous(limits = c(0,1))+
  xlab("Estimated Median Fish Count") +
  ylab("Shoreline Type")+
  theme_classic() +
  theme(panel.spacing.y = unit(0,"line"), legend.position = "left")
siteipa_plot

#period + site + ipa + site*ipa
# best model for Smelt
cols =c(brewer.pal(2,"Paired"))
siteipa_sm <- bind_rows(smelt_model[[4]])
siteipa_sm_plot <- ggplot(siteipa_sm, aes(median, shoreline)) +
  geom_point(aes(color = shoreline)) +
  geom_errorbarh(aes(xmin = .lower, xmax=.upper,color=shoreline),width=0.1,size =0.3) +
  scale_color_manual(values=rev(cols), name=" "
                     #wes_palette(n=3, name="GrandBudapest1")
  )+
  guides(color = guide_legend(reverse = TRUE))+
  facet_grid(site~period, scales = "free_x")+
  # scale_x_continuous(limits = c(0,1))+
  xlab("Estimated Median Fish Count") +
  ylab("Shoreline Type")+
  theme_classic() +
  theme(panel.spacing.y = unit(0,"line"), legend.position = "left")
siteipa_sm_plot
# 

# 
#
#
#
#
#
#
##########################
# CODE BELOW HERE NOT USED
# 
# Transmute code not used
# # site*ipa 
# for(i in 1: length(df_list))
# #or to exclude species where model doesn't converge
# #fits <- c()
# #for(i in fits)
# {
# # i=4 
#   df = as.data.frame(df_list[[i]]) # to id each species
# #  NAME <- unique(df$species)
# NAME <- "Surf Smelt"   
#   site_ipa2[[i]] <-  as.data.frame(model9[[i]]) %>%
#     transmute(Site5.Natural = `(Intercept)` + siteDOK + ipaNatural + `siteDOK:ipaNatural`,
#               Site5.Armored = `(Intercept)` + siteDOK,
#               Site5.Restored = `(Intercept)` + siteDOK + ipaRestored + `siteDOK:ipaRestored`,
#               Site6.Natural = `(Intercept)` + siteEDG + ipaNatural + `siteEDG:ipaNatural`,
#               Site6.Armored = `(Intercept)` + siteEDG,
#               Site6.Restored = `(Intercept)` + siteEDG +ipaRestored + `siteEDG:ipaRestored`,
#               Site1.Natural = `(Intercept)` + siteFAM +ipaNatural + `siteFAM:ipaNatural`,
#               Site1.Armored = `(Intercept)` + siteFAM ,
#               Site1.Restored = `(Intercept)` + siteFAM +ipaRestored + `siteFAM:ipaRestored`,
#               Site4.Natural = `(Intercept)` + siteSHR +ipaNatural + `siteSHR:ipaNatural`,
#               Site4.Armored = `(Intercept)` + siteSHR ,
#               Site4.Restored = `(Intercept)` + siteSHR +ipaRestored + `siteSHR:ipaRestored`,
#               Site2.Natural = `(Intercept)` + siteTUR +ipaNatural + `siteTUR:ipaNatural`,
#               Site2.Armored = `(Intercept)` + siteTUR ,
#               Site2.Restored = `(Intercept)` + siteTUR +ipaRestored + `siteTUR:ipaRestored`,
#               Site3.Natural = `(Intercept)` + ipaNatural ,
#               Site3.Armored = `(Intercept)` ,
#               Site3.Restored = `(Intercept)` +ipaRestored) %>%
#     gather()%>%
#     group_by(key) %>%
#     mean_hdi(.width=0.9) %>%
#     mutate_if(is.double, round, digits = 2) %>%
#     mutate(median = exp(value)) %>%
#     mutate(.lower = exp(.lower)) %>%
#     mutate(.upper = exp(.upper)) %>%
#     separate(key, c("site", "shoreline")) %>%
#     select(-c("value")) %>%
#       mutate(species = paste(as.character(NAME))) 
#   }

# # WITH TUR & SHR REMOVED
# i=4 
# df = as.data.frame(df_list[[i]]) # to id each species
# #  NAME <- unique(df$species)
# NAME <- "Surf Smelt"   
# smelt_model[[i]] <-  as.data.frame(model11[[i]]) %>%
#   transmute(Per1.COR.Natural.Pres = `(Intercept)` + ipaNatural ,
#             Per1.COR.Armored.Pres = `(Intercept)` ,
#             Per1.COR.Restored.Pres = `(Intercept)` +ipaRestored,
#             Per2.COR.Natural.Pres = `(Intercept)` + ipaNatural + periodperiod2,
#             Per2.COR.Armored.Pres = `(Intercept)`+ periodperiod2,
#             Per2.COR.Restored.Pres = `(Intercept)` +ipaRestored + periodperiod2,
#             Per3.COR.Natural.Pres = `(Intercept)` + ipaNatural + periodperiod3,
#             Per3.COR.Armored.Pres = `(Intercept)` + periodperiod3,
#             Per3.COR.Restored.Pres = `(Intercept)` +ipaRestored + periodperiod3,
#             Per4.COR.Natural.Pres = `(Intercept)` + ipaNatural + periodperiod4,
#             Per4.COR.Armored.Pres = `(Intercept)` + periodperiod4,
#             Per4.COR.Restored.Pres = `(Intercept)` +ipaRestored + periodperiod4,
#             Per1.DOK.Natural.Abs = `(Intercept)` + siteDOK + ipaNatural + `siteDOK:ipaNatural` + vegAbsent,
#             Per1.DOK.Armored.Abs = `(Intercept)` + siteDOK + vegAbsent,
#             Per1.DOK.Restored.Abs = `(Intercept)` + siteDOK + ipaRestored + `siteDOK:ipaRestored` + vegAbsent,
#             Per1.EDG.Natural.Abs = `(Intercept)` + siteEDG + ipaNatural + `siteEDG:ipaNatural` + vegAbsent,
#             Per1.EDG.Armored.Abs = `(Intercept)` + siteEDG + vegAbsent,
#             Per1.EDG.Restored.Abs = `(Intercept)` + siteEDG +ipaRestored + `siteEDG:ipaRestored` + vegAbsent,
#             Per1.FAM.Natural.Abs = `(Intercept)` + siteFAM +ipaNatural + `siteFAM:ipaNatural` + vegAbsent,
#             Per1.FAM.Armored.Abs = `(Intercept)` + siteFAM + vegAbsent,
#             Per1.FAM.Restored.Abs = `(Intercept)` + siteFAM +ipaRestored + `siteFAM:ipaRestored` + vegAbsent,
#             Per1.DOK.Natural.Abs = `(Intercept)` + siteDOK + ipaNatural + `siteDOK:ipaNatural` + periodperiod2 + vegAbsent,
#             Per2.DOK.Armored.Abs = `(Intercept)` + siteDOK + periodperiod2 + vegAbsent,
#             Per2.DOK.Restored.Abs = `(Intercept)` + siteDOK + ipaRestored + `siteDOK:ipaRestored` + periodperiod2 + vegAbsent,
#             Per2.EDG.Natural.Abs = `(Intercept)` + siteEDG + ipaNatural + `siteEDG:ipaNatural` + periodperiod2 + vegAbsent,
#             Per2.EDG.Armored.Abs = `(Intercept)` + siteEDG + periodperiod2 + vegAbsent,
#             Per2.EDG.Restored.Abs = `(Intercept)` + siteEDG +ipaRestored + `siteEDG:ipaRestored` + periodperiod2 + vegAbsent,
#             Per2.FAM.Natural.Abs = `(Intercept)` + siteFAM +ipaNatural + `siteFAM:ipaNatural` + periodperiod2 + vegAbsent,
#             Per2.FAM.Armored.Abs = `(Intercept)` + siteFAM + periodperiod2 + vegAbsent,
#             Per2.FAM.Restored.Abs = `(Intercept)` + siteFAM +ipaRestored + `siteFAM:ipaRestored` + periodperiod2 + vegAbsent,
#             Per3.DOK.Natural.Abs = `(Intercept)` + siteDOK + ipaNatural + `siteDOK:ipaNatural` + periodperiod3 + vegAbsent,
#             Per3.DOK.Armored.Abs = `(Intercept)` + siteDOK + periodperiod3 + vegAbsent,
#             Per3.DOK.Restored.Abs = `(Intercept)` + siteDOK + ipaRestored + `siteDOK:ipaRestored` + periodperiod3 + vegAbsent,
#             Per3.EDG.Natural.Abs = `(Intercept)` + siteEDG + ipaNatural + `siteEDG:ipaNatural` + periodperiod3 + vegAbsent,
#             Per3.EDG.Armored.Abs = `(Intercept)` + siteEDG + periodperiod3 + vegAbsent,
#             Per3.EDG.Restored.Abs = `(Intercept)` + siteEDG +ipaRestored + `siteEDG:ipaRestored` + periodperiod3 + vegAbsent,
#             Per3.FAM.Natural.Abs = `(Intercept)` + siteFAM +ipaNatural + `siteFAM:ipaNatural` + periodperiod3 + vegAbsent,
#             Per3.FAM.Armored.Abs = `(Intercept)` + siteFAM + periodperiod3 + vegAbsent,
#             Per3.FAM.Restored.Abs = `(Intercept)` + siteFAM +ipaRestored + `siteFAM:ipaRestored` + periodperiod3 + vegAbsent,
#             Per4.DOK.Natural.Abs = `(Intercept)` + siteDOK + ipaNatural + `siteDOK:ipaNatural` + periodperiod4 + vegAbsent,
#             Per4.DOK.Armored.Abs = `(Intercept)` + siteDOK + periodperiod4 + vegAbsent,
#             Per4.DOK.Restored.Abs = `(Intercept)` + siteDOK + ipaRestored + `siteDOK:ipaRestored` + periodperiod4 + vegAbsent,
#             Per4.EDG.Natural.Abs = `(Intercept)` + siteEDG + ipaNatural + `siteEDG:ipaNatural` + periodperiod4 + vegAbsent,
#             Per4.EDG.Armored.Abs = `(Intercept)` + siteEDG + periodperiod4 + vegAbsent,
#             Per4.EDG.Restored.Abs = `(Intercept)` + siteEDG +ipaRestored + `siteEDG:ipaRestored` + periodperiod4 + vegAbsent,
#             Per4.FAM.Natural.Abs = `(Intercept)` + siteFAM +ipaNatural + `siteFAM:ipaNatural` + periodperiod4 + vegAbsent,
#             Per4.FAM.Armored.Abs = `(Intercept)` + siteFAM + periodperiod4 + vegAbsent,
#             Per4.FAM.Restored.Abs = `(Intercept)` + siteFAM +ipaRestored + `siteFAM:ipaRestored` + periodperiod4 + vegAbsent) %>%
#   gather()%>%
#   group_by(key) %>%
#   mean_hdi(.width=0.89) %>%
#   mutate_if(is.double, round, digits = 2) %>%
#   mutate(median = exp(value)) %>%
#   mutate(.lower = exp(.lower)) %>%
#   mutate(.upper = exp(.upper)) %>%
#   separate(key, c("period", "site", "shoreline", "veg")) %>%
#   select(-c("value")) %>%
#   mutate(species = paste(as.character(NAME))) 
# #}
# 

#   #ipa
# #   ipa[[i]] <-  as.data.frame(model2[[i]]) %>%
# #     transmute(Armored = `(Intercept)`,
# #               Natural = `(Intercept)` + ipaNatural,
# #               Restored = `(Intercept)` + ipaRestored) %>%
# #     gather() %>%
# #     group_by(key) %>%
# #     median_hdi(.width=0.5) %>%
# #     mutate_if(is.double, round, digits = 2) %>%
# #     mutate(median = exp(value)) %>%
# #     mutate(.lower = exp(.lower)) %>%
# #     mutate(.upper = exp(.upper)) %>% 
# #     select(-c("value")) %>%
# #     mutate(species = paste(as.character(NAME)))  
# # }



# # veg
# fits <- c(1,2) #to select certain species
# for(i in fits)
# {
#   df = as.data.frame(df_list[[i]]) # to id each species
#   NAME <- unique(df$species)
#   
# veg[[i]] <-  as.data.frame(model1[[i]]) %>%
#   transmute(vegPresent = `(Intercept)`,
#             vegAbsent  = `(Intercept)` + vegAbsent)%>%
#   gather() %>%
#   group_by(key) %>%
#   mean_hdi(.width=0.95) %>%
#   mutate_if(is.double, round, digits = 2) %>%
#   #mutate(key)
#   mutate(median = exp(value)) %>%
#   mutate(.lower = exp(.lower)) %>%
#   mutate(.upper = exp(.upper)) %>%
#   select(-c("value")) %>%
#   mutate(species = paste(as.character(NAME)))
# }



# # site + ipa
# NAME = c("Herring")
# i=3
# site_plus_ipa <-  as.data.frame(model6[[3]]) %>%
#   transmute(Site5.Natural = `(Intercept)` + siteDOK + ipaNatural,
#             Site5.Armored = `(Intercept)` + siteDOK ,
#             Site5.Restored = `(Intercept)` + siteDOK +ipaRestored,
#             Site6.Natural = `(Intercept)` + siteEDG +ipaNatural,
#             Site6.Armored = `(Intercept)` + siteEDG ,
#             Site6.Restored = `(Intercept)` + siteEDG +ipaRestored,
#             Site1.Natural = `(Intercept)` + siteFAM +ipaNatural,
#             Site1.Armored = `(Intercept)` + siteFAM ,
#             Site1.Restored = `(Intercept)` + siteFAM +ipaRestored,
#             Site4.Natural = `(Intercept)` + siteSHR +ipaNatural,
#             Site4.Armored = `(Intercept)` + siteSHR ,
#             Site4.Restored = `(Intercept)` + siteSHR +ipaRestored,
#             Site2.Natural = `(Intercept)` + siteTUR +ipaNatural,
#             Site2.Armored = `(Intercept)` + siteTUR ,
#             Site2.Restored = `(Intercept)` + siteTUR +ipaRestored,
#             Site3.Natural = `(Intercept)` + ipaNatural ,
#             Site3.Armored = `(Intercept)` ,
#             Site3.Restored = `(Intercept)` +ipaRestored) %>%
#   gather()%>%
#   group_by(key) %>%
#   mean_hdi(.width=0.9) %>%
#   mutate_if(is.double, round, digits = 2) %>%
#   mutate(mean = exp(value)) %>%
#   mutate(.lower = exp(.lower)) %>%
#   mutate(.upper = exp(.upper)) %>%
#   separate(key, c("site", "shoreline")) %>%
#   select(-c("value")) %>%
#     mutate(species = paste(as.character(NAME)))



# #ipa + vegetation
# veg_plus_ipa = list()
# #for(i in 1: length(df_list))
# #{
# #  df = as.data.frame(df_list[[i]]) # to id each species
#       #NAME <- unique(df$species)
#       NAME <- c("Chum")
#         veg_plus_ipa <-  as.data.frame(model5[[1]]) %>%
#   transmute(Present.Armored = `(Intercept)` ,
#             Present.Restored = `(Intercept)` + ipaRestored,
#             Present.Natural = `(Intercept)`+ ipaNatural,
#             Absent.Armored = `(Intercept)` + vegAbsent,
#             Absent.Natural = `(Intercept)` + vegAbsent + ipaRestored,
#             Absent.Restored = `(Intercept)` + vegAbsent + ipaNatural)%>%
#       gather() %>%
#   group_by(key) %>%
#   median_hdi(.width=0.5) %>%
#   mutate_if(is.double, round, digits = 2) %>%
#   mutate(median = exp(value)) %>%
#   mutate(.lower = exp(.lower)) %>%
#   mutate(.upper = exp(.upper)) %>%
#   separate(key, c("Vegetation", "Shoreline")) %>%
#   select(-c("value")) %>%
#     mutate(species = paste(as.character(NAME)))
# #}



# # ipa * vegetation
# veg_times_ipa = list()
# for(i in 1: length(df_list))
# {
#   df = as.data.frame(df_list[[i]]) # to id each species
#       NAME <- unique(df$species)
# veg_times_ipa <-  as.data.frame(model12[[i]]) %>%
#   transmute(veg.Arm =  `(Intercept)`+ siteDOK + siteEDG + siteFAM + siteSHR + siteTUR, 
#             noveg.Arm = `(Intercept)` + vegAbsent + siteDOK + siteEDG + siteFAM + siteSHR + siteTUR,
#             veg.Nat = `(Intercept)` + ipaNatural + siteDOK + siteEDG + siteFAM + siteSHR + siteTUR
#             
#             vegAbsent 
#             + vegPresent + ipaNatural + vegPresent:ipaNatural,
#              Present.Restored = `(Intercept)` + vegPresent + ipaRestored + vegPresent:ipaRestored,
#              Present.Armored = `(Intercept)` + vegPresent,
#              Absent.Natural = `(Intercept)` + ipaNatural,
#              Absent.Restored =`(Intercept)` +  ipaRestored + ipaRestored:vegPresent)%>%
#   gather() %>%
#   group_by(key) %>%
#   #mean_hdi(.width=0.5) %>%
#   median_hdi(.width=0.5) %>%
#   mutate_if(is.double, round, digits = 2) %>%
#   mutate(median = exp(value)) %>%
#   mutate(.lower = exp(.lower)) %>%
#   mutate(.upper = exp(.upper)) %>%
#   separate(key, c("Vegetation", "Shoreline")) %>%
#   select(-c("value")) %>%
#   mutate(species = paste(as.character(NAME)))


#   i=4 
#   df = as.data.frame(df_list[[i]]) # to id each species
#   #  NAME <- unique(df$species)
#   NAME <- "Surf Smelt"   
#   smelt_model[[i]] <-  as.data.frame(model11[[i]]) %>%
#     transmute(Per1.SHR.Natural = `(Intercept)` + siteSHR +ipaNatural + `siteSHR:ipaNatural`,
#               Per1.SHR.Armored = `(Intercept)` + siteSHR ,
#               Per1.SHR.Restored = `(Intercept)` + siteSHR +ipaRestored + `siteSHR:ipaRestored`,
#               Per1.TUR.Natural = `(Intercept)` + siteTUR +ipaNatural + `siteTUR:ipaNatural`,
#               Per1.TUR.Armored = `(Intercept)` + siteTUR ,
#               Per1.TUR.Restored = `(Intercept)` + siteTUR +ipaRestored + `siteTUR:ipaRestored`,
#               Per1.COR.Natural = `(Intercept)` + ipaNatural ,
#               Per1.COR.Armored = `(Intercept)` ,
#               Per1.COR.Restored = `(Intercept)` +ipaRestored,
#               Per1.DOK.Natural = `(Intercept)` + siteDOK + ipaNatural + `siteDOK:ipaNatural` + periodperiod2,
#               Per2.SHR.Natural = `(Intercept)` + siteSHR +ipaNatural + `siteSHR:ipaNatural` + periodperiod2,
#               Per2.SHR.Armored = `(Intercept)` + siteSHR + periodperiod2,
#               Per2.SHR.Restored = `(Intercept)` + siteSHR +ipaRestored + `siteSHR:ipaRestored` + periodperiod2,
#               Per2.TUR.Natural = `(Intercept)` + siteTUR +ipaNatural + `siteTUR:ipaNatural` + periodperiod2,
#               Per2.TUR.Armored = `(Intercept)` + siteTUR + periodperiod2,
#               Per2.TUR.Restored = `(Intercept)` + siteTUR +ipaRestored + `siteTUR:ipaRestored` + periodperiod2,
#               Per2.COR.Natural = `(Intercept)` + ipaNatural + periodperiod2,
#               Per2.COR.Armored = `(Intercept)`+ periodperiod2,
#               Per2.COR.Restored = `(Intercept)` +ipaRestored + periodperiod2,
#               Per3.SHR.Natural = `(Intercept)` + siteSHR +ipaNatural + `siteSHR:ipaNatural` + periodperiod3,
#               Per3.SHR.Armored = `(Intercept)` + siteSHR + periodperiod3,
#               Per3.SHR.Restored = `(Intercept)` + siteSHR +ipaRestored + `siteSHR:ipaRestored` + periodperiod3,
#               Per3.TUR.Natural.Pres = `(Intercept)` + siteTUR +ipaNatural + `siteTUR:ipaNatural` + periodperiod3,
#               Per3.TUR.Armored.Pres = `(Intercept)` + siteTUR + periodperiod3,
#               Per3.TUR.Restored.Pres = `(Intercept)` + siteTUR +ipaRestored + `siteTUR:ipaRestored` + periodperiod3,
#               Per3.COR.Natural.Pres = `(Intercept)` + ipaNatural + periodperiod3,
#               Per3.COR.Armored.Pres = `(Intercept)` + periodperiod3,
#               Per3.COR.Restored.Pres = `(Intercept)` +ipaRestored + periodperiod3,
#               Per4.SHR.Natural.Pres = `(Intercept)` + siteSHR +ipaNatural + `siteSHR:ipaNatural` + periodperiod4,
#               Per4.SHR.Armored.Pres = `(Intercept)` + siteSHR + periodperiod4,
#               Per4.SHR.Restored.Pres = `(Intercept)` + siteSHR +ipaRestored + `siteSHR:ipaRestored` + periodperiod4,
#               Per4.TUR.Natural.Pres = `(Intercept)` + siteTUR +ipaNatural + `siteTUR:ipaNatural` + periodperiod4,
#               Per4.TUR.Armored.Pres = `(Intercept)` + siteTUR + periodperiod4,
#               Per4.TUR.Restored.Pres = `(Intercept)` + siteTUR +ipaRestored + `siteTUR:ipaRestored` + periodperiod4,
#               Per4.COR.Natural.Pres = `(Intercept)` + ipaNatural + periodperiod4,
#               Per4.COR.Armored.Pres = `(Intercept)` + periodperiod4,
#               Per4.COR.Restored.Pres = `(Intercept)` +ipaRestored + periodperiod4,
#               Per1.DOK.Natural.Abs = `(Intercept)` + siteDOK + ipaNatural + `siteDOK:ipaNatural` + vegAbsent,
#               Per1.DOK.Armored.Abs = `(Intercept)` + siteDOK + vegAbsent,
#               Per1.DOK.Restored.Abs = `(Intercept)` + siteDOK + ipaRestored + `siteDOK:ipaRestored` + vegAbsent,
#               Per1.EDG.Natural.Abs = `(Intercept)` + siteEDG + ipaNatural + `siteEDG:ipaNatural` + vegAbsent,
#               Per1.EDG.Armored.Abs = `(Intercept)` + siteEDG + vegAbsent,
#               Per1.EDG.Restored.Abs = `(Intercept)` + siteEDG +ipaRestored + `siteEDG:ipaRestored` + vegAbsent,
#               Per1.FAM.Natural.Abs = `(Intercept)` + siteFAM +ipaNatural + `siteFAM:ipaNatural` + vegAbsent,
#               Per1.FAM.Armored.Abs = `(Intercept)` + siteFAM + vegAbsent,
#               Per1.FAM.Restored.Abs = `(Intercept)` + siteFAM +ipaRestored + `siteFAM:ipaRestored` + vegAbsent,
#               Per1.DOK.Natural.Abs = `(Intercept)` + siteDOK + ipaNatural + `siteDOK:ipaNatural` + periodperiod2 + vegAbsent,
#               Per2.DOK.Armored.Abs = `(Intercept)` + siteDOK + periodperiod2 + vegAbsent,
#               Per2.DOK.Restored.Abs = `(Intercept)` + siteDOK + ipaRestored + `siteDOK:ipaRestored` + periodperiod2 + vegAbsent,
#               Per2.EDG.Natural.Abs = `(Intercept)` + siteEDG + ipaNatural + `siteEDG:ipaNatural` + periodperiod2 + vegAbsent,
#               Per2.EDG.Armored.Abs = `(Intercept)` + siteEDG + periodperiod2 + vegAbsent,
#               Per2.EDG.Restored.Abs = `(Intercept)` + siteEDG +ipaRestored + `siteEDG:ipaRestored` + periodperiod2 + vegAbsent,
#               Per2.FAM.Natural.Abs = `(Intercept)` + siteFAM +ipaNatural + `siteFAM:ipaNatural` + periodperiod2 + vegAbsent,
#               Per2.FAM.Armored.Abs = `(Intercept)` + siteFAM + periodperiod2 + vegAbsent,
#               Per2.FAM.Restored.Abs = `(Intercept)` + siteFAM +ipaRestored + `siteFAM:ipaRestored` + periodperiod2 + vegAbsent,
#               Per3.DOK.Natural.Abs = `(Intercept)` + siteDOK + ipaNatural + `siteDOK:ipaNatural` + periodperiod3 + vegAbsent,
#               Per3.DOK.Armored.Abs = `(Intercept)` + siteDOK + periodperiod3 + vegAbsent,
#               Per3.DOK.Restored.Abs = `(Intercept)` + siteDOK + ipaRestored + `siteDOK:ipaRestored` + periodperiod3 + vegAbsent,
#               Per3.EDG.Natural.Abs = `(Intercept)` + siteEDG + ipaNatural + `siteEDG:ipaNatural` + periodperiod3 + vegAbsent,
#               Per3.EDG.Armored.Abs = `(Intercept)` + siteEDG + periodperiod3 + vegAbsent,
#               Per3.EDG.Restored.Abs = `(Intercept)` + siteEDG +ipaRestored + `siteEDG:ipaRestored` + periodperiod3 + vegAbsent,
#               Per3.FAM.Natural.Abs = `(Intercept)` + siteFAM +ipaNatural + `siteFAM:ipaNatural` + periodperiod3 + vegAbsent,
#               Per3.FAM.Armored.Abs = `(Intercept)` + siteFAM + periodperiod3 + vegAbsent,
#               Per3.FAM.Restored.Abs = `(Intercept)` + siteFAM +ipaRestored + `siteFAM:ipaRestored` + periodperiod3 + vegAbsent,
#               Per4.DOK.Natural.Abs = `(Intercept)` + siteDOK + ipaNatural + `siteDOK:ipaNatural` + periodperiod4 + vegAbsent,
#               Per4.DOK.Armored.Abs = `(Intercept)` + siteDOK + periodperiod4 + vegAbsent,
#               Per4.DOK.Restored.Abs = `(Intercept)` + siteDOK + ipaRestored + `siteDOK:ipaRestored` + periodperiod4 + vegAbsent,
#               Per4.EDG.Natural.Abs = `(Intercept)` + siteEDG + ipaNatural + `siteEDG:ipaNatural` + periodperiod4 + vegAbsent,
#               Per4.EDG.Armored.Abs = `(Intercept)` + siteEDG + periodperiod4 + vegAbsent,
#               Per4.EDG.Restored.Abs = `(Intercept)` + siteEDG +ipaRestored + `siteEDG:ipaRestored` + periodperiod4 + vegAbsent,
#               Per4.FAM.Natural.Abs = `(Intercept)` + siteFAM +ipaNatural + `siteFAM:ipaNatural` + periodperiod4 + vegAbsent,
#               Per4.FAM.Armored.Abs = `(Intercept)` + siteFAM + periodperiod4 + vegAbsent,
#               Per4.FAM.Restored.Abs = `(Intercept)` + siteFAM +ipaRestored + `siteFAM:ipaRestored` + periodperiod4 + vegAbsent) %>%
#     gather()%>%
#     group_by(key) %>%
#     mean_hdi(.width=0.89) %>%
#     mutate_if(is.double, round, digits = 2) %>%
#     mutate(median = exp(value)) %>%
#     mutate(.lower = exp(.lower)) %>%
#     mutate(.upper = exp(.upper)) %>%
#     separate(key, c("period", "site", "shoreline", "veg")) %>%
#     select(-c("value")) %>%
#     mutate(species = paste(as.character(NAME))) 
# #}

```

#
# Figure code not used
# #site*ipa + veg + period
# cols =c(brewer.pal(4,"Paired"))
# smeltmodel <-  bind_rows(smelt_model[[4]])
# smelt_plot <- ggplot(smeltmodel, aes(median, shoreline)) +
#   geom_point(aes(color = shoreline)) +
#   geom_errorbarh(aes(xmin = .lower, xmax=.upper,color=shoreline),width=0.1,size =0.3) +
#   scale_color_manual(values=rev(cols), name=" "
#                      #wes_palette(n=3, name="GrandBudapest1")
#   )+
#   guides(color = guide_legend(reverse = TRUE))+
#   
#   facet_grid(site~period, scales = "free_x")+
#   # scale_x_continuous(limits = c(0,1))+
#   xlab("Estimated Median Fish Count") +
#   ylab("Shoreline Type")+
#   theme_classic() +
#   theme(panel.spacing.y = unit(0,"line"), legend.position = "left")
# smelt_plot

# cols =c(brewer.pal(3,"Paired"))
# ipa_plot <-  bind_rows(ipa) %>%
# #  filter(!species=="Shiner Perch") %>%
#   ggplot(aes(median, key)) +
#   geom_point(aes(color = key)) +
#   geom_errorbarh(aes(xmin = .lower, xmax=.upper,color=key),width=0.1,size =0.3) +
#   scale_color_manual(values=rev(cols), name=" "
#                      #wes_palette(n=3, name="GrandBudapest1")
#   )+
#   guides( color = guide_legend(reverse = TRUE))+
#   facet_wrap(~species, ncol =1, nrow =5, strip.position = "right")+
#   # scale_x_continuous(limits = c(0,1))+
#   xlab("Estimated Median Fish Count") +
#   ylab("Shoreline Type")+
#   theme_classic() +
#   theme(legend.position = "left",panel.spacing.y = unit(0,"line")) 
# ipa_plot
# 
# #  save_plot("shoreline.jpg", shoreline, base_height = 4, base_width = 7)
# 
# cols =rev(c("steelblue1", "steelblue2", "steelblue3", "steelblue", "steelblue4", "navyblue"))
# 
# site 
# my_color = c("snow3","snow3","black","black","black","snow3")
#     site.fig <- bind_rows(site) %>%
#                  mutate(key= factor(key, levels = rev(c("FAM", "TUR", "COR", "SHR", "DOK", "EDG")))) %>%
#                 # filter(!species=="Shiner Perch") %>%
#                       ggplot(aes(median,  key)) +
#                       geom_point(aes(color = key)) +
#                       geom_errorbarh(aes(xmin = .lower, xmax=.upper, color=key)) +
#                       facet_wrap(~species, ncol =1, nrow =5, strip.position = "right")+
#                       scale_color_manual(name=" ",values = my_color
# #                                         values=wes_palette(n=6, name="GrandBudapest1", type = "continuous")
#                                          )+
#                       xlab("Estimated Median Fish Count") +
#                       ylab(" ")+
#                       #ggtitle("Net")+
#                       guides( color = guide_legend(reverse = TRUE))+
#                       theme_classic() +
#                       theme(legend.position = "none",
#                             panel.spacing.y = unit(0,"line"))
#       site.fig
#      save_plot("site_May22.jpg", site.fig, base_height = 4, base_width = 6)
# 

# site
# trying w estimate_means/posterior means
     # em_3_chi <- estimate_means(model9[[1]], cis=0.8)
     # em_3_chu <- estimate_means(model9[[2]], cis=0.8)
     # em_6 <- estimate_means(model6[[3]], cis=0.8)
     # em_9 <- estimate_means(model9[[4]], cis=0.9)
     # # cols site/Mean/CI_low/CI_high
     # # add species column to estimate_means tables
     # em_3_chi$species <- "Chinook"
     # em_3_chu$species <- "Chum"
     # em_3 <- bind_rows(em_3_chi, em_3_chu)
     # 
     # my_color = c("snow3","snow3","darkolivegreen","darkolivegreen","darkolivegreen","snow3")
     # site2.fig <- em_3 %>% 
     #  # mutate(site = recode(site, COR="Site.3", DOK="Site.5", EDG="Site.6",FAM="Site.1", SHR="Site.4", TUR="Site.2")) %>%
     #   mutate(site = fct_relevel(site, "EDG", "DOK", "SHR", "COR", "TUR", "FAM")) %>%
     #   ggplot(aes(Mean,  site)) +
     #   geom_point(aes(color = site)) +
     #   geom_errorbarh(aes(xmin = CI_low, xmax=CI_high, color=site)) +
     #   facet_wrap(~species, ncol =1, nrow =5, strip.position = "right")+
     #   scale_color_manual(name=" ",values = my_color
     #                      #                                         values=wes_palette(n=6, name="GrandBudapest1", type = "continuous")
     #   )+
     #   xlab("Estimated Mean Fish Count") +
     #   ylab(" ")+
     #   #ggtitle("Net")+
     #   guides( color = guide_legend(reverse = TRUE))+
     #   theme_classic() +
     #   theme(legend.position = "left",
     #         panel.spacing.y = unit(0,"line"))
     # site2.fig
     # save_plot("site_May22.jpg", site.fig, base_height = 4, base_width = 6)
     # #
     
     
#site*ipa
# cols =c(brewer.pal(3,"Paired"))
#     F.siteipa <-  bind_rows(site_ipa[[4]])
#       F.siteipa_plot <- ggplot(F.siteipa, aes(median, shoreline)) +
#                 geom_point(aes(color = shoreline)) +
#                 geom_errorbarh(aes(xmin = .lower, xmax=.upper,color=shoreline),width=0.1,size =0.3) +
#                 scale_color_manual(values=rev(cols), name=" "
#                                      #wes_palette(n=3, name="GrandBudapest1")
#                                      )+
#                 guides(color = guide_legend(reverse = TRUE))+
# 
#                 facet_grid(site~species, scales = "free_x")+
#                # scale_x_continuous(limits = c(0,1))+
#                 xlab("Estimated Median Fish Count") +
#                 ylab("Shoreline Type")+
#                 theme_classic() +
#                 theme(panel.spacing.y = unit(0,"line"), legend.position = "left")
# # 
# # site*ipa using estimate_means
#       cols =c(brewer.pal(3,"Paired"))
#       F.siteipa2 <-  em_9 %>%
#      # mutate(site = recode(site, COR="Site.3", DOK="Site.5", EDG="Site.6",FAM="Site.1", SHR="Site.4", TUR="Site.2")) %>%
#       mutate(site = fct_relevel(site, "FAM", "TUR", "COR", "SHR", "DOK", "EDG")) %>%
#       mutate(species="Surf Smelt")
#       F.siteipa_plot2 <- ggplot(F.siteipa2, aes(Mean, ipa)) +
#         geom_point(aes(color = ipa)) +
#         geom_errorbarh(aes(xmin = CI_low, xmax=CI_high,color=ipa),width=0.1,size =0.3) +
#         scale_color_manual(values=rev(cols), name=" "
#                            #wes_palette(n=3, name="GrandBudapest1")
#         )+
#         guides(color = guide_legend(reverse = TRUE))+
#         
#         facet_grid(site~species, scales = "free_x")+
#         # scale_x_continuous(limits = c(0,1))+
#         xlab("Estimated Median Fish Count") +
#         ylab("Shoreline Type")+
#         theme_classic() +
#         theme(panel.spacing.y = unit(0,"line"), legend.position = "none")
#     F.siteipa_plot2
#       
#       
# #site + ipa using estimate_means
#     cols =c(brewer.pal(2,"Paired"))
#     siteipa2 <- em_6 %>%
#       #mutate(site = recode(site, COR="Site.3", DOK="Site.5", EDG="Site.6",FAM="Site.1", SHR="Site.4", TUR="Site.2")) %>%
#       mutate(site = fct_relevel(site, "FAM", "TUR", "COR", "SHR", "DOK", "EDG")) %>%
#       mutate(species="Herring")
#     siteipa_plot2 <- ggplot(siteipa2, aes(Mean, ipa)) +
#       geom_point(aes(color = ipa)) +
#       geom_errorbarh(aes(xmin = CI_low, xmax=CI_high,color=ipa),width=0.1,size =0.3) +
#       scale_color_manual(values=rev(cols), name=" "
#                          #wes_palette(n=3, name="GrandBudapest1")
#       )+
#       guides(color = guide_legend(reverse = TRUE))+
#       facet_grid(site~species, scales = "free_x")+
#       # scale_x_continuous(limits = c(0,1))+
#       xlab("Estimated Median Fish Count") +
#       ylab("Shoreline Type")+
#       theme_classic() +
#       theme(panel.spacing.y = unit(0,"line"), legend.position = "none")
#     siteipa_plot2
#     # 
#     
    

#site_ipa and site*ipa in 2-panel plot
# ggarrange(siteipa_plot2, F.siteipa_plot2, ncol=2, nrow=1,common.legend = FALSE, 
          labels = c("A", "B"))

# # veg
# #my_color = c("snow3","darkolivegreen","snow3","darkolivegreen","snow3","darkolivegreen","snow3","darkolivegreen")
# my_color = c("snow3","black","snow3","black","snow3","black","snow3","black")
#  
# veg_only <- bind_rows(veg)
# veg_only_plot <- veg_only %>%
#   mutate(Eelgrass = factor(key, levels = c("Absent", "Present"))) %>%
#   ggplot(aes(x=Eelgrass, y=median, fill = Eelgrass)) +
#   geom_bar(colour ="black",position="dodge", stat="identity")+
#      scale_fill_manual(values = my_color)+
#   geom_errorbar(aes(ymin = .lower, ymax=.upper, group = key), width = 0.01, position=position_dodge(.9)) +
#    facet_wrap(~species, scales = "free_y", ncol = 5) +
#   ylab("Mean Fish Count") +
#   xlab(" ")+
#   theme_classic() +
#   theme(strip.background =element_rect(fill="white"), axis.text.x = element_text(angle = 90)) 
  
# veg
#trying w estimate_means
my_color = c("darkolivegreen","snow3","darkolivegreen","snow3","darkolivegreen","snow3","darkolivegreen")

em_1_chi <- estimate_means(model1[[1]], cis=0.9)
em_1_chu <- estimate_means(model1[[2]], cis=0.9)
# cols site/Mean/CI_low/CI_high
# add species column 
em_1_chi$species <- "Chinook"
em_1_chu$species <- "Chum"
em_1 <- bind_rows(em_1_chi, em_1_chu)

veg_only_plot2 <- em_1 %>%
  ggplot(aes(x=veg, y=Mean, fill = veg)) +
  geom_bar(colour ="black",position="dodge", stat="identity")+
  scale_fill_manual(values = my_color)+
  geom_errorbar(aes(ymin = CI_low, ymax=CI_high, group = veg), width = 0.01, position=position_dodge(.9)) +
  facet_wrap(~species, scales = "free_y", ncol = 5) +
  ylab("Est. Median Fish Count") +
  xlab(" ")+
  theme_classic() +
  theme(strip.background =element_rect(fill="white"), axis.text.x = element_text(angle = 90))
veg_only_plot2


# # veg + ipa
# veg_ipa <- bind_rows(veg_plus_ipa)
# veg_plot <- veg_ipa %>%
#    ggplot(aes(Shoreline,median, fill=Vegetation)) +
#   geom_bar(aes(group = Vegetation,Shoreline), colour ="black",position="dodge", stat="identity") +
#   scale_fill_manual(values = my_color)+
#   geom_errorbar(aes(ymin = .lower, ymax=.upper, group = Vegetation), width = 0.01, position=position_dodge(.9)) +
#   # facet_wrap(~species, scales = "free_y", ncol = 5) +
#   ylab("Est. Median Fish Count") +
#   xlab(" ")+
#   theme_classic() +
#   theme(strip.background =element_rect(fill="white"), axis.text.x = element_text(angle = 90))
# veg_plot
# 
# 
# # veg*ipa
#  veg_times_ipa <- bind_rows(veg_times_ipa)
#   veg_ipa_plot <- veg_times_ipa %>%
#     ggplot(aes(Shoreline, median, fill=Vegetation)) +
#         geom_bar(aes(group = Vegetation, Shoreline), colour ="black",position="dodge", stat="identity") +
#     scale_fill_manual(values = my_color)+
#     geom_errorbar(aes(ymin = .lower, ymax=.upper, group = Vegetation), width = 0.01, position=position_dodge(.9)) +
#    # facet_wrap(~species, scales = "free_y", ncol = 5) +
#     ylab("Est. Median Fish Count") +
#     xlab(" ")+
#     theme_classic() +
#     theme(strip.background =element_rect(fill="white"), axis.text.x = element_text(angle = 90))
# veg_ipa_plot
# #save_plot("Eelgrass Fish .jpg", eg_plot, base_height = 4, base_width = 6)
# ```
# site*IPA interaction only (change depending on which species/model)
site.ipa.inxnplot.Chin <- emmip(model9[[1]], ipa ~ site, CIs=TRUE) +
  ggtitle("Chinook")
site.ipa.inxnplot.Chin

site.ipa.inxnplot.Chum <- emmip(model9[[1]], ipa ~ site, CIs=TRUE) +
  ggtitle("Chinook")
site.ipa.inxnplot.Chin

# other plots for paper

# posterior distributions for coefficients
mc_smelt <- mcmc_intervals(model9[[4]], prob_outer = 0.95) +
  ggplot2::labs(
    title = "Smelt", xlab = "Effect size" )
mc_smelt

# site*IPA interaction only (change depending on which species/model)
site.ipa.inxnplot.Chin <- emmip(model9[[1]], ipa ~ site) +
  ggtitle("Chinook")
sii.C <- site.ipa.inxnplot.Chin

site.ipa.inxnplot.chum <- emmip(model9[[2]], ipa ~ site) +
  ggtitle("Chum")
sii.cu <- site.ipa.inxnplot.chum

site.ipa.inxnplot.herr <- emmip(model9[[3]], ipa ~ site) +
  ggtitle("Herring")
sii.he <- site.ipa.inxnplot.herr

site.ipa.inxnplot.smelt <- emmip(model9[[4]], ipa ~ site, CIs=TRUE) +
  ggtitle("Smelt")
sii.sm <- site.ipa.inxnplot.smelt

ggarrange(sii.C, sii.cu, sii.he, sii.sm, ncol=2, nrow=2, common.legend = TRUE)


# veg*IPA interaction only (change depending on which species/model)
veg.ipa.inxnplot.Chin <- emmip(model12[[1]], ipa ~ veg, CIs=TRUE) +
  ggtitle("Chin")
veg.ipa.inxnplot.Chin


# site + IPA 
site.ipa.Herring <- emmip(model6[[3]], ipa ~ site, CIs=TRUE) +
  ggtitle("Herring")
site.ipa.Herring

# site*IPA interaction only (change depending on which species/model)
# bar plot
cols = c(brewer.pal(12, "Paired"))
# first pass model into a dataframe
site.ipa.inxnplot.Chin <- emmip(model9[[1]], ipa ~ site, CIs=TRUE,
                                plotit=FALSE)
p <- ggplot(site.ipa.inxnplot.Chin, aes(x=site,y=yvar, fill=ipa)) + 
  geom_bar(stat="identity",position="dodge", color="black")
p1 <- p + geom_errorbar(position=position_dodge(.9),width=.25, 
                        aes(ymax=UCL, ymin=LCL),alpha=0.3)
p1 <- p1 + scale_fill_manual(values=rev(cols))
p1 <- p1 + theme_classic()
p1 <- p1  + labs(x="Site", y="Log(Count)", fill="Shoreline",title="Chinook")


# first pass model into a dataframe
site.ipa.inxnplot.chum <- emmip(model9[[2]], ipa ~ site, CIs=TRUE,
                                plotit=FALSE)
q <- ggplot(site.ipa.inxnplot.chum, aes(x=site,y=yvar, fill=ipa)) + 
  geom_bar(stat="identity",position="dodge", color="black")
q2 <- q + geom_errorbar(position=position_dodge(.9),width=.25, 
                        aes(ymax=UCL, ymin=LCL),alpha=0.3)
q2 <- q2 + scale_fill_manual(values=rev(cols))
q2 <- q2 + theme_classic()
q2 <- q2  + labs(x="Site", y="Log(Count)", fill="Shoreline", title="Chum")


# first pass model into a dataframe
site.ipa.inxnplot.herr <- emmip(model9[[3]], ipa ~ site, CIs=TRUE,
                                plotit=FALSE)
r <- ggplot(site.ipa.inxnplot.herr, aes(x=site,y=yvar, fill=ipa)) + 
  geom_bar(stat="identity",position="dodge", color="black")
r1 <- r + geom_errorbar(position=position_dodge(.9),width=.25, 
                        aes(ymax=UCL, ymin=LCL),alpha=0.3)
r1 <- r1 + scale_fill_manual(values=rev(cols))
r1 <- r1 + theme_classic()
r1 <- r1  + labs(x="Site", y="Effect Size", fill="Shoreline", title="Herring")


# first pass model into a dataframe
site.ipa.inxnplot.smelt <- emmip(model9[[4]], ipa ~ site, CIs=TRUE,
                                 plotit=FALSE)
s <- ggplot(site.ipa.inxnplot.smelt, aes(x=site,y=yvar, fill=ipa)) + 
  geom_bar(stat="identity",position="dodge", color="black")
s1 <- s + geom_errorbar(position=position_dodge(.9),width=.25, 
                        aes(ymax=UCL, ymin=LCL),alpha=0.3)
s1 <- s1 + scale_fill_manual(values=rev(cols))
s1 <- s1 + theme_classic()
s1 <- s1  + labs(x="Site", y="Log(Abundance)", fill="Shoreline", title="Smelt")


# first pass model into a dataframe
veg.ipa.inxnplot.Chum <- emmip(model8[[2]], ipa ~ veg, CIs=TRUE,
                                 plotit=FALSE)
s <- ggplot(veg.ipa.inxnplot.Chum, aes(x=veg,y=yvar, fill=ipa)) + 
  geom_bar(stat="identity",position="dodge", color="black")
s1 <- s + geom_errorbar(position=position_dodge(.9),width=.25, 
                        aes(ymax=UCL, ymin=LCL),alpha=0.3)
s1 <- s1 + scale_fill_manual(values=rev(cols))
s1 <- s1 + theme_classic()
s1 <- s1  + labs(x="Eelgrass", y="Log(Abundance)", fill="Shoreline", title="Chum")

ggarrange(p1, q2, ncol=2, nrow=1, common.legend = TRUE)


# 
# # general plots
# #abundance on Y axis, year/month on X axis
# 
# abundance <- net %>%
#   # filter(!month %in% c("04","05")) %>% 
#   filter(species %in% c("Chinook", "Chum", "Herring", "Surf Smelt"))%>%
#   group_by(year, month, species) %>%
#   summarise(count=sum(species_count)) %>%
#   ggplot(aes(as.factor(month, year), count, fill=species)) +
#   geom_bar(position="dodge", stat="identity") +
#   facet_wrap(~species)
# 
# plot(abundance)
# 
# #histograms of fish sizes by species
# p <- ggplot(net_import2[which(net_import2$length_mm>0),], aes(x=species, y=length_mm, fill=species)) + 
#   geom_violin()
# p + coord_flip()+ stat_summary(fun.y=mean, geom="point", size=2, color="black")
# p 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
#  # calculate effect sizes
#  # prepare lists; change below to fit this model
#  # create one list for each effect, and one for the SDs
#   sd=list()
#   site_eff = list()
#   site_ipa_eff = list()
# 
#  # start the effect size calculation; iterate through species
#   for(i in 1: length(df_list)){
#     df = as.data.frame(df_list[[i]]) #to id each species
#    NAME <- c("Shoreline Effects")
# 
#  #get standard deviation of draws to calculate effect size; CHANGE BELOW TO FIT THIS MODEL
#  sd <- as.data.frame(model8[[i]])  %>%
#   # this transmute was for veg*ipa effects; not sure if it's needed here
#    transmute(Absent.Armored = `(Intercept)` ,
#              Present.Natural =  `(Intercept)`+ vegPresent + ipaNatural + vegPresent:ipaNatural,
#              Present.Restored = `(Intercept)` + vegPresent + ipaRestored + vegPresent:ipaRestored,
#              Present.Armored = `(Intercept)` + vegPresent,
#              Absent.Natural = `(Intercept)` + ipaNatural,
#              Absent.Restored =`(Intercept)` +  ipaRestored + ipaRestored:vegPresent)%>%
#    gather() %>%
#    group_by(key)  %>%
#    mutate(value = exp(value)) %>%
#    summarise(sd = sd(value)) %>%
#    mutate(species = paste(as.character(NAME)))
# 
#  #change the "width" to change your credible interval
#  site_ipa_eff <-  as.data.frame(model8[[i]])  %>%
#    transmute(Absent.Armored = `(Intercept)` ,
#              Present.Natural =  `(Intercept)`+ vegPresent + ipaNatural + vegPresent:ipaNatural,
#              Present.Restored = `(Intercept)` + vegPresent + ipaRestored + vegPresent:ipaRestored,
#              Present.Armored = `(Intercept)` + vegPresent,
#              Absent.Natural = `(Intercept)` + ipaNatural,
#              Absent.Restored =`(Intercept)` +  ipaRestored + ipaRestored:vegPresent)%>%
#    gather() %>%
#    group_by(key)  %>%
#    mean_hdi(.width=0.50) %>%
#    #median_hdi(.width=0.5) %>%
#    mutate_if(is.double, round, digits = 2) %>%
#    mutate(mean = exp(value)) %>%
#    mutate(.lower = exp(.lower)) %>%
#    mutate(.upper = exp(.upper)) %>%
#    #separate(key, c("Vegetation", "Shoreline")) %>%
#    select(-c("value")) %>%
#    mutate(species = paste(as.character(NAME)))
# 
#   veg[[i]] <-  as.data.frame(model1[[i]])  %>%
#     transmute(Present = `(Intercept)`,
#               Absent = `(Intercept)` + vegAbsent)%>%
#     gather() %>%
#     group_by(key)  %>%
#     mean_hdi(.width=0.50) %>%
#     #median_hdi(.width=0.5) %>%
#     mutate_if(is.double, round, digits = 2) %>%
#     mutate(mean = exp(value)) %>%
#     mutate(.lower = exp(.lower)) %>%
#    mutate(.upper = exp(.upper)) %>%
#    select(-c("value")) %>%
#    mutate(species = paste(as.character(NAME)))
# # #
# ipa[[i]] <-  as.data.frame(model1[[i]])  %>%
#    transmute(Armored = `(Intercept)` ,
#              Natural =  `(Intercept)`+ ipaNatural,
#              Restored = `(Intercept)` + ipaRestored)  %>%
#    gather() %>%
#    group_by(key)  %>%
#    mean_hdi(.width=0.5) %>%
#    #median_hdi(.width=0.5) %>%
#    mutate_if(is.double, round, digits = 2) %>%
#    mutate(mean = exp(value)) %>%
#    mutate(.lower = exp(.lower)) %>%
#    mutate(.upper = exp(.upper)) %>%
#    select(-c("value")) %>%
#    mutate(species = paste(as.character(NAME)))
# }
# 
# #bind rows to calculate effect size and plot
# df<- site_eff %>%
#   bind_rows()
# 
# #get "effect sizes" or numbers for differences between groups
# #above needs to be mean_hdi for these values
# #probably want to plot everything using the mean instead of median
# #
# # #tidy the sd info, join with means to calculate effect size.
# # #pooled sd is only for more than 1 covariate I think
#  sd_pooled <- bind_rows(sd) %>%
#    group_by(species) %>%
#    spread(Vegetation, sd) %>%
#    summarise(pooled = sqrt((Absent^2)+ (Present^2)/2))
# #
# effect <-  df %>%
#   select(species, Shoreline,Vegetation,median) %>%
#   spread(Vegetation, median) %>%
#   group_by(species, Shoreline) %>%
#   summarise(m = Present-Absent) %>% #top half of equation
#   left_join(sd_pooled) %>% #add in pooled std dev from draws
#   mutate(effect = m/pooled) %>%
#   select(-c(m,pooled))
# #
# # # write.csv(effect, "output/effect_size.csv")
# 
# 
# ## plot effect sizes using plot(model) which plots credible intervals
# 
# 
# 
# 
# 
# 
# #use month as replicate...
# df<- net %>%
#   select(-1) %>%
#   separate(month, c("year", "month", "day"), sep = "-")%>%
#   filter(species %in% c("Chinook", "Chum", "Surf Smelt", "Herring")) %>%
#   within(species_count [species == 'none'] <- '0') %>%
#   mutate(veg = recode(site, 'COR' = "Present",
#                       'TUR'="Present",
#                       'FAM'="Absent",
#                       'DOK'="Absent",
#                       'EDG'="Absent",
#                       'SHR'="Present")) %>%
#   mutate(species_count = as.numeric(species_count))%>% 
#   group_by(year,month,veg,ipa, species) %>%
#   summarise(total = sum(species_count))  %>%
#   ungroup( ) %>%
#   spread(species, total, fill = 0)   %>%
#   gather(5:8, key = species, value = count) %>%
# #  filter(!month %in% c("04","05")) %>%
#   filter(!month %in% c("07","08","09")) %>%
# 
#   filter(!count == 357)
# 
# col.filters <- unique(df$species) 
# 
lapply(seq_along(col.filters), function(x) {
  filter(df, species == col.filters[x])
}
) -> df_list
# 
# names(df_list) <- col.filters
# 
# # test for site effects on count
# model1 = list()
# for(i in 1: length(df_list)) {
#   model1[[i]] <- stan_glm(count ~ ipa*veg ,
#                           data = df_list[[i]],
#                           family = neg_binomial_2(link="log"), init_r=1,
#                           adapt_delta = 0.999
#                           # prior_intercept = normal(0,10), prior = normal(0,10)
#   )
# }
# 
# print(model1[1]) #Chinook summary
# print(model1[2]) #chum summary
# 
# 
# # create one list for each effect, and one for the SDs
# sd=list()
# veg_ipa = list()
# ipa=list()
# veg=list()
# 
# # start the effect size calculation; iterate through species
# for(i in 1: length(df_list)){
#   df = as.data.frame(df_list[[i]]) #to id each species
#   NAME <- unique(df$species)
#   
#   #get standard deviation of draws to calculate effect size; CHANGE BELOW TO FIT THIS MODEL    
#   sd[[i]] <- as.data.frame(model1[[i]])  %>%
#     # this transmute was for veg*ipa effects; not sure if it's needed here
#     # transmute(siteCOR = `(Intercept)` ,
#     Present.Natural =  `(Intercept)`+ ipaNatural,
#   Present.Restored = `(Intercept)` + ipaRestored,
#   Absent.Armored = `(Intercept)` + vegAbsent,
#   Absent.Natural = `(Intercept)` + vegAbsent + ipaNatural +ipaNatural:vegAbsent,
#   Absent.Restored =`(Intercept)` + vegAbsent + ipaRestored + ipaRestored:vegAbsent) %>%
#   gather() %>%
#   group_by(key)  %>%
#   mutate(value = exp(value)) %>%
#   summarise(sd = sd(value)) %>%
#   separate(key, c("Vegetation", "Shoreline")) %>%
#   mutate(species = paste(as.character(NAME))) 
# 
# #change the "width" to change your credible interval
# veg_ipa[[i]] <-  as.data.frame(model1[[i]])  %>%
#   transmute(Present.Armored = `(Intercept)` ,
#             Present.Natural =  `(Intercept)`+ ipaNatural,
#             Present.Restored = `(Intercept)` + ipaRestored,
#             Absent.Armored = `(Intercept)` + vegAbsent,
#             Absent.Natural = `(Intercept)` + vegAbsent + ipaNatural +ipaNatural:vegAbsent,
#             Absent.Restored =`(Intercept)` + vegAbsent + ipaRestored + ipaRestored:vegAbsent) %>%
#   gather() %>%
#   group_by(key)  %>%
#   mean_hdi(.width=0.50) %>%
#   #median_hdi(.width=0.5) %>%
#   mutate_if(is.double, round, digits = 2) %>%
#   mutate(mean = exp(value)) %>%
#   mutate(.lower = exp(.lower)) %>%
#   mutate(.upper = exp(.upper)) %>%
#   separate(key, c("Vegetation", "Shoreline")) %>%
#   select(-c("value")) %>%
#   mutate(species = paste(as.character(NAME)))
# 
# veg[[i]] <-  as.data.frame(model1[[i]])  %>%
#   transmute(Present = `(Intercept)`,
#             Absent = `(Intercept)` + vegAbsent)%>%
#   gather() %>%
#   group_by(key)  %>%
#   mean_hdi(.width=0.50) %>%
#   #median_hdi(.width=0.5) %>%
#   mutate_if(is.double, round, digits = 2) %>%
#   mutate(mean = exp(value)) %>%
#   mutate(.lower = exp(.lower)) %>%
#   mutate(.upper = exp(.upper)) %>%
#   select(-c("value")) %>%
#   mutate(species = paste(as.character(NAME)))
# 
# ipa[[i]] <-  as.data.frame(model1[[i]])  %>%
#   transmute(Armored = `(Intercept)` ,
#             Natural =  `(Intercept)`+ ipaNatural,
#             Restored = `(Intercept)` + ipaRestored)  %>%
#   gather() %>%
#   group_by(key)  %>%
#   mean_hdi(.width=0.5) %>%
#   #median_hdi(.width=0.5) %>%
#   mutate_if(is.double, round, digits = 2) %>%
#   mutate(mean = exp(value)) %>%
#   mutate(.lower = exp(.lower)) %>%
#   mutate(.upper = exp(.upper)) %>%
#   select(-c("value")) %>%
#   mutate(species = paste(as.character(NAME)))
# }
# 
# #bind rows to calculate effect size and plot 
# df<- veg_ipa %>% 
#   bind_rows() 
# 
# #get "effect sizes" or numbers for differences between groups
# #above needs to be mean_hdi for these values 
# #probably want to plot everything using the mean instead of median
# 
# #tidy the sd info, join with means to calculate effect size. 
# sd_pooled <- bind_rows(sd) %>%
#   group_by(species, Shoreline) %>%
#   spread(Vegetation, sd) %>%
#   summarise(pooled = sqrt((Absent^2)+ (Present^2)/2))
# 
# effect <-  df %>%
#   select(species, Shoreline,Vegetation,median) %>%
#   spread(Vegetation, median) %>%
#   group_by(species, Shoreline) %>%
#   summarise(m = Present-Absent) %>% #top half of equation
#   left_join(sd_pooled) %>% #add in pooled std dev from draws
#   mutate(effect = m/pooled) %>%
#   select(-c(m,pooled))
# 
# # write.csv(effect, "output/effect_size.csv")
# 
# 
# ```{r exploratory - FISH ABUNDANCE BY YEAR - PLOTS, eval = FALSE, echo=FALSE, warning= FALSE}
# #abundance on Y axis, year on X grouped by species
# 
# abundance <- net %>%
#   # filter(!month %in% c("04","05")) %>% 
#   filter(species %in% c("Chinook", "Chum", "Herring", "Surf Smelt"))%>%
#   separate(month, c("year", "month", "day"), sep = "-")%>%
#   group_by(year, month, site, species) %>%
#   summarise(count=sum(species_count)) %>%
#   ggplot(aes(as.factor(month), count, fill=species)) +
#   geom_bar(position="dodge", stat="identity") +
#   facet_wrap(~species)
# 
# plot(abundance)
# 
# 
# # multi-plot of effect sizes for a single model, 1 panel per species
# # can't figure out how to put these pre-designed plots into a grid
# chi_8 <- plot(model8[[1]])
# chu_8 <- plot(model8[[2]])
# her_8 <- plot(model8[[3]])
# ss_8 <- plot(model8[[4]])
# 
# plot_grid(
#   chi_8, chu_8, her_8, ss_8, 
#   labels=c("A","B","C","D"),
#   
# )
# plot(model8[[1]])
# plot(model8[[2]])
# plot(model8[[3]])
# plot(model8[[4]])
# 
# for(i in 1: length(df_list)){
#   plot(model8[[i]])
# }

# model comparison - failed effort by Tessa to use lists
loo_res <- list() #make list for loo output; all models and all species
model_list <- list(model1, model2, model3, model4, model5, model6, 
                   model7, model8, model9, model10, model11, model12) #put all model estimates in list
for(j in seq_along(model_list))  {
  for(i in 1:4)  {
    modelf <- model_list[[j]]
    loo_res[[j]][[i]] <- loo(modelf[[i]]) 
  }
}

########
# LENGTH MODELS; NOT USED
# #for length models, use different df
# lapply(seq_along(col.filters), function(x) {
#   filter(ldf, species == col.filters[x])
# }
# ) -> ldf_list
# 
# names(ldf_list) <- col.filters
# 
# 
# model13 = list()
# for(i in 1: length(ldf_list)) {
#   model13[[i]] <- stan_glm(total ~ ipa*perc_sm,
#                            data = ldf_list[[i]],
#                            family = neg_binomial_2(link="log"), init_r=1,
#                            adapt_delta = 0.999
#                            # prior_intercept = normal(0,10), prior = normal(0,10)
#   )
# }
# print(model13)
# 
# model14 = list()
# for(i in 1: length(ldf_list)) {
#   model14[[i]] <- stan_glm(total ~ ipa*perc_sm + veg,
#                            data = ldf_list[[i]],
#                            family = neg_binomial_2(link="log"), init_r=1,
#                            adapt_delta = 0.999
#                            # prior_intercept = normal(0,10), prior = normal(0,10)
#   )
# }
# print(model14)
# 
# model15 = list()
# for(i in 1: length(ldf_list)) {
#   model15[[i]] <- stan_glm(total ~ site + ipa*perc_sm,
#                            data = ldf_list[[i]],
#                            family = neg_binomial_2(link="log"), init_r=1,
#                            adapt_delta = 0.999
#                            # prior_intercept = normal(0,10), prior = normal(0,10)
#   )
# }
# print(model15)
# 
# model16 = list()
# for(i in 1: length(ldf_list)) {
#   model16[[i]] <- stan_glm(total ~ site + ipa*perc_sm + veg,
#                            data = ldf_list[[i]],
#                            family = neg_binomial_2(link="log"), init_r=1,
#                            adapt_delta = 0.999
#                            # prior_intercept = normal(0,10), prior = normal(0,10)
#   )
# }
# print(model16)
# 
# model17 = list()
# for(i in 1: length(ldf_list)) {
#   model17[[i]] <- stan_glm(total ~ ipa*perc_sm + veg + veg*ipa,
#                            data = ldf_list[[i]],
#                            family = neg_binomial_2(link="log"), init_r=1,
#                            adapt_delta = 0.999
#                            # prior_intercept = normal(0,10), prior = normal(0,10)
#   )
# }
# print(model17)
# 
# model18 = list()
# for(i in 1: length(ldf_list)) {
#   model18[[i]] <- stan_glm(total ~ site + ipa*perc_sm + site*ipa,
#                            data = ldf_list[[i]],
#                            family = neg_binomial_2(link="log"), init_r=1,
#                            adapt_delta = 0.999
#                            # prior_intercept = normal(0,10), prior = normal(0,10)
#   )
# }
# print(model18)
# 
# model19 = list()
# for(i in 1: length(ldf_list)) {
#   model19[[i]] <- stan_glm(total ~ site + ipa*perc_sm + veg + site*ipa + veg*ipa,
#                            data = ldf_list[[i]],
#                            family = neg_binomial_2(link="log"), init_r=1,
#                            adapt_delta = 0.999
#                            # prior_intercept = normal(0,10), prior = normal(0,10)
#   )
# }
# print(model19)

#UNUSED CODE
# length data
#use this if you are using googlesheets4
# net_import <- drive_get("https://docs.google.com/spreadsheets/d/1OhsndJNLAlHxT0TTHp7dmgbKgPZ7dgd-xNnE-aAxLdc/edit#gid=0") %>%
#   read_sheet( ) 
# 
# net_import2 <-net_import %>%
#   separate(date, into = c("year","month", "day"), sep = "-") %>%
#   mutate(month = str_pad(month, 2, side = c("left"), pad = "0")) %>%
#   mutate(day = str_pad(day, 2, side = c("left"), pad = "0"))   %>%
#   mutate(ipa = as.factor(ipa))   %>%
#   mutate(year = as.integer(year)) %>%
#   mutate(month = as.integer(month)) %>%
#   select(-c(transect_notes)) %>%
#   mutate(species = replace_na(species, "none"),length_mm = replace_na(length_mm, 0)) %>%
#   filter(species %in% c("Chinook", "Chum", "Surf Smelt", "Herring")) %>%
#   filter(length_mm > 0) %>% #only want to categorize positive catches
#   mutate(count = 1) %>% #for tally of all fish at each month/site/ipa/species
#   mutate(small = case_when(length_mm <100  ~ 1, T ~0)) %>% #assign each fish to small category (1) or not (0)
#   group_by(year,month,site,ipa, species) %>%
#   summarise(perc_sm = (sum(small)/sum(count))*100) %>% #percentage of all fish that are small
#   ungroup( ) %>%
#   mutate_if(is.character, as.factor) %>%
#   merge(all, all=TRUE) #so can be merged with df
# 
# #merge with df to pull in count and veg columns, for use in models that include length
# new_df <- full_join(net_import2, df) %>% 
#   filter(!is.na(perc_sm)) #make this df only include positive catches
#   
# ldf <- new_df %>%
#   mutate(species=as.factor(species))

#for snorkel data
#set up dataframe for count data using snorkel data
# df_snorkel<- snorkel %>%
#   select(-1) %>%
#   group_by(year,month,site,ipa, species) %>%
#   # how do I retain mean_length_mm while grouping/summarizing
#   summarise(total = sum(mean_species_count)) %>% 
#   ungroup( ) %>%
#   within(total [species == 'NA'] <- '0') %>%
#   mutate(total = ceiling(as.numeric(total))) %>%
#   mutate_if(is.character, as.factor) %>%
#   # mutate(log_size = log(mean_length_mm)) %>%
#   #  mutate(log_size = as.numeric(log_size)) %>%
#   #  within(log_size [mean_length_mm =="NA"] <- '0')
#   spread(ipa, total, fill = 0)  %>%
#   gather(5:7, key=ipa, value=total) %>%
#   spread(species, total, fill = 0)   %>%
#   gather(5:62, key = species, value = total) %>%
#   filter(species %in% c("Chinook", "Chum", "Surf Smelt", "Herring")) %>%
#   mutate(veg = recode(site, 'COR' = "Present",
#                       'TUR'="Present",
#                       'FAM'="Absent",
#                       'DOK'="Absent",
#                       'EDG'="Absent",
#                       'SHR'="Present")) 
# df<- df_snorkel


#for presence absence data 
# net_PA <- net %>% 
#   select(-1) %>%
#   separate(month, c("year", "month", "day"), sep = "-")%>%
#   filter(species %in% c("Chinook", "Chum", "Surf Smelt", "Herring")) %>%
#   within(species_count [species == 'none'] <- '0') %>% 
#   mutate(pres_abs = case_when(species_count < 0.9 ~ 0,
#                               TRUE ~ 1)) %>%
#   select(-c("tax_group", "org_type", "species_count", "mean_length_mm")) %>%
#  
#   #spread(species, pres_abs, fill = 0) %>% 
#   group_by(ipa, species) %>%
#   summarise(total = sum(pres_abs))  %>%
#   ungroup( ) %>%
#   select(1:6, 9:11, 26,45, 50,52,55,62) %>%
#   gather(7:15, key = species, value = pres_abs) %>%
#   filter(!species %in% c("Pink", "Coho", "Sand Lance", "Snake Prickleback")) %>% 
#   #filter(!month %in% c("04","05")) %>% 
#   #filter(!station == 2) %>%
#   mutate(year = as.factor(year))

# model 5 w Chum and Chinook

# bar plot
cols = c(brewer.pal(12, "Paired"))
# use estimate_means output
p <- ggplot(em_3_chu, aes(x=site,y=Mean, fill=ipa)) + 
  geom_bar(stat="identity",position="dodge", color="black")
p1 <- p + geom_errorbar(position=position_dodge(.9),width=.25, 
                        aes(ymax=10, ymin=CI_low),alpha=0.3)
p1 <- p1 + scale_fill_manual(values=rev(cols))
p1 <- p1 + theme_classic()
p1 <- p1  + labs(x="Site", y="Estimated Fish Count", fill="Shoreline",title="Chinook")
p1

#Kernel density plots
k_ch <- mcmc_dens_overlay(model9[[1]], pars=c("(Intercept)", "ipaNatural", "ipaRestored")) 
k_ch + geom_vline(aes(xintercept = 0), linetype = "dashed") 
k_ch + labs(x="Effect size")

k_c <- mcmc_dens_overlay(model9[[2]], pars=c("(Intercept)", "ipaNatural", "ipaRestored")) + geom_vline(aes(xintercept = 0), linetype = "dashed") 
ggarrange(k_ch, k_c, nrow=2)
