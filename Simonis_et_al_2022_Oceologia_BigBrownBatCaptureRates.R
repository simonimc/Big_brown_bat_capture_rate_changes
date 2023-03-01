#This R script file provides code for all statistical analyses from Simonis et al. 2022, Oceologia. 

#Data and its description can be found and downloaded on the Dryad Digital Repository at https://doi.org/10.5061/dryad.ngf1vhhvv 

#Due to the large amount of data, the following analyses were run remotely on a university computing core.
##install packages
#install.packages('emmeans')
#install.packages('brms')
#install.packages('rstan')
#install.packages('dplyr')
#install.packages('plyr')
#install.packages('tidyr')
#install.packages('tidyverse')

##turn on packages
library(emmeans)
library(brms)
library(rstan)
library(dplyr)
library(plyr)
library(tidyr)
library(tidyverse)

#remove environment to start b/c cluster does not automatically remove it each time
rm(list=ls())

EPFU_data<- read.csv('SIMONIS_et_al_BigBrownBatData_Dryad.csv', header = T, sep = ",")

head(EPFU_data)
str(EPFU_data)

#Get capture rates (capture count per year per site) for adults
EPFUcount<- EPFU_data %>% 
  select(state, site_mask, month, year, age, sex, repstat, years_Pd, disease_time_step, county_centroid_lat, county_centroid_lon, rMapCounty)  

EPFUcount<- count(EPFUcount, vars = c('state', 'site_mask', 'month', 'year', 'sex', 'age', 'repstat', 'years_Pd', 'disease_time_step', 'county_centroid_lat', 'county_centroid_lon', 'rMapCounty')) %>% 
  group_by(year, site_mask)

str(EPFUcount)
head(EPFUcount)

adult_EPFUcount<- EPFUcount %>%
  filter(age == 'adult') %>%
  droplevels()


#Adjust variables and levels where needed
adult_EPFUcount$disease_time_step<- factor(adult_EPFUcount$disease_time_step, levels = c('pre-invasion', 'invasion', 'epidemic', 'established'))


#Begin creating first initial models for all adults and for females only
#get priors for males and females
get_prior(freq~sex*disease_time_step*county_centroid_lat + (1|site_mask), data = adult_EPFUcount, family = poisson())

adult_prior<-set_prior("student_t(3, 0.7, 2.5)", class = "Intercept")


#Create initial model for adult counts as a function of disease time-step and sex
adult_count_mod1<- brm(freq ~ sex*disease_time_step*county_centroid_lat + (1|site_mask), data = adult_EPFUcount, family = poisson(), prior = adult_prior, chains = 4, warmup = 5000, iter = 10000, sample_prior = 'yes', control = list(max_treedepth = 15))


#run prior_summary() to check match
prior_summary(adult_count_mod1)
#still a match at Intercept and sd

summary(adult_count_mod1)

#Get conditional effects--bayes equivalent of frequentist fixed effects
cond_ef_adult_count1<- emmeans(adult_count_mod1, ~sex*disease_time_step*county_centroid_lat)

cond_ef_adult_count1

#try stan_code() so it spits out your betas and all the pre model parameters that brms eliminates the process of
stancode(adult_count_mod1)



#now do same steps but for females only and include repstat
#Get count datasets for females only
fem_EPFUcount<- EPFUcount %>%
  filter(sex == 'female', age == 'adult') %>%
  droplevels()

#Adjust variables and levels where needed
fem_EPFUcount$repstat<- factor(fem_EPFUcount$repstat, levels= c("non-reproductive", "pregnant", "lactating", "post-lactating"))

fem_EPFUcount$disease_time_step<- factor(fem_EPFUcount$disease_time_step, levels = c('pre-invasion', 'invasion', 'epidemic', 'established'))

#get priors for females
get_prior(freq~county_centroid_lat*repstat*disease_time_step + (1|site_mask), data = fem_EPFUcount, family = poisson())

fem_prior<-set_prior("student_t(3, 0.7, 2.5)", class = "Intercept")


#Create initial model for females 
fem_count_mod1<- brm(freq ~ county_centroid_lat*repstat*disease_time_step + (1|site_mask), data = fem_EPFUcount, family = poisson(), prior = fem_prior, chains = 4, warmup = 5000, iter = 10000, sample_prior = 'yes', control = list(max_treedepth = 15))

#run prior_summary() to check match
prior_summary(fem_count_mod1)
#still a match at Intercept and sd


summary(fem_count_mod1)


#Get conditional effects--which is really a post-hoc type thing for the bayes equivalent of fixed effects
cond_ef_fem_count1<- emmeans(fem_count_mod1, ~county_centroid_lat*disease_time_step*repstat)

cond_ef_fem_count1

#try stan_code() so it spits out your betas and all the pre model parameters that brms eliminates the process of
stancode(fem_count_mod1)



#Now create and test secondary models for all adults and females only, replacing latitude with north/south of interaction point
#first for all adults
#add latitudinal category variable from first model outputs in step 1 and latitudinal category weights var
adult_EPFUcount$lat_cat<- ifelse(adult_EPFUcount$county_centroid_lat > 39.3, 'north', 'south')

#get priors for males and females
get_prior(freq~sex*disease_time_step*lat_cat + (1|site_mask), data = adult_EPFUcount, family = poisson())

adult_prior<-set_prior("student_t(3, 0.7, 2.5)", class = "Intercept")


#Create second model for adult counts as a function of disease time-step, sex and latitudinal category
adult_count_mod2<- brm(freq ~ sex*disease_time_step*lat_cat + (1|site_mask), data = adult_EPFUcount, family = poisson(), prior = adult_prior, chains = 4, warmup = 5000, iter = 10000, sample_prior = 'yes', control = list(max_treedepth = 15))


#run prior_summary() to check match
prior_summary(adult_count_mod2)
#still a match at Intercept and sd

summary(adult_count_mod2)

#Get conditional effects--bayes equivilent of frequentist fixed effects
cond_ef_adult_count2<- emmeans(adult_count_mod2, ~sex*disease_time_step*lat_cat)

cond_ef_adult_count2

#try stan_code() so it spits out your betas and all the pre model parameters that brms eliminates the process of
stancode(adult_count_mod2)


#now for females only
#add latitudinal category variable 
fem_EPFUcount$lat_cat<- ifelse(fem_EPFUcount$county_centroid_lat > 39.3, 'north', 'south')


#get priors for females
get_prior(freq~lat_cat*repstat*disease_time_step + (1|site_mask), data = fem_EPFUcount, family = poisson())

fem_prior2<-set_prior("student_t(3, 0.7, 2.5)", class = "Intercept")


#Create second model for females 
fem_count_mod2<- brm(freq ~ lat_cat*repstat*disease_time_step + (1|site_mask), data = fem_EPFUcount, family = poisson(), prior = fem_prior2, chains = 4, warmup = 5000, iter = 10000, sample_prior = 'yes', control = list(max_treedepth = 15))


#run prior_summary() to check match
prior_summary(fem_count_mod2)
#still a match at Intercept and sd

summary(fem_count_mod2)

#Get conditional effects
cond_ef_fem_count2<- emmeans(fem_count_mod2, ~lat_cat*disease_time_step*repstat)

cond_ef_fem_count2

#try stan_code() so it spits out your betas and all the pre model parameters that brms eliminates the process of
stancode(fem_count_mod2)

#save environment for later
save.image(file='BigBrownBat_CaptureRate_Bayes_Cluster.RData')