#This R script file provides code for all statistical analyses from Simonis et al. 2024, Ecology and Evolution. 

#Data and its description can be found and downloaded on the Dryad Digital Repository at https://doi.org/10.5061/dryad.ngf1vhhvv 

#Due to the large amount of highly variable data, the following analyses should be run 
#remotely on a university computing core. It takes ~48 hours to run this script on a full node (24 cores)

##install packages
#install.packages('emmeans')
#install.packages('brms')
#install.packages('rstan')
#install.packages('dplyr')
#install.packages('lubridate')
#install.packages('tidyr')


##turn on packages
library(emmeans)
library(brms)
library(rstan)
library(dplyr)
library(lubridate)
library(tidyr)

#remove environment to start b/c cluster does not automatically remove it each time
rm(list=ls())

EPFU_data<- read.csv('SIMONIS_et_al_BigBrownBatData_Dryad.csv', header = T, sep = ",")

head(EPFU_data)
str(EPFU_data)

#make a dates column
EPFU_data$day<- as.character(EPFU_data$day)
EPFU_data$day<- ifelse(EPFU_data$day == '1', 
                       paste('0', EPFU_data$day, sep = ''), EPFU_data$day)
EPFU_data$day<- ifelse(EPFU_data$day == '2', 
                       paste('0', EPFU_data$day, sep = ''), EPFU_data$day)
EPFU_data$day<- ifelse(EPFU_data$day == '3', 
                       paste('0', EPFU_data$day, sep = ''), EPFU_data$day)
EPFU_data$day<- ifelse(EPFU_data$day == '4', 
                       paste('0', EPFU_data$day, sep = ''), EPFU_data$day)
EPFU_data$day<- ifelse(EPFU_data$day == '5', 
                       paste('0', EPFU_data$day, sep = ''), EPFU_data$day)
EPFU_data$day<- ifelse(EPFU_data$day == '6', 
                       paste('0', EPFU_data$day, sep = ''), EPFU_data$day)
EPFU_data$day<- ifelse(EPFU_data$day == '7', 
                       paste('0', EPFU_data$day, sep = ''), EPFU_data$day)
EPFU_data$day<- ifelse(EPFU_data$day == '8', 
                       paste('0', EPFU_data$day, sep = ''), EPFU_data$day)
EPFU_data$day<- ifelse(EPFU_data$day == '9', 
                       paste('0', EPFU_data$day, sep = ''), EPFU_data$day)
unique(EPFU_data$day)

EPFU_data$month<- as.character(EPFU_data$month)
EPFU_data$month<- ifelse(EPFU_data$month == '1', 
                         paste('0', EPFU_data$month, sep = ''), EPFU_data$month)
EPFU_data$month<- ifelse(EPFU_data$month == '2', 
                         paste('0', EPFU_data$month, sep = ''), EPFU_data$month)
EPFU_data$month<- ifelse(EPFU_data$month == '3', 
                         paste('0', EPFU_data$month, sep = ''), EPFU_data$month)
EPFU_data$month<- ifelse(EPFU_data$month == '4', 
                         paste('0', EPFU_data$month, sep = ''), EPFU_data$month)
EPFU_data$month<- ifelse(EPFU_data$month == '5', 
                         paste('0', EPFU_data$month, sep = ''), EPFU_data$month)
EPFU_data$month<- ifelse(EPFU_data$month == '6', 
                         paste('0', EPFU_data$month, sep = ''), EPFU_data$month)
EPFU_data$month<- ifelse(EPFU_data$month == '7', 
                         paste('0', EPFU_data$month, sep = ''), EPFU_data$month)
EPFU_data$month<- ifelse(EPFU_data$month == '8', 
                         paste('0', EPFU_data$month, sep = ''), EPFU_data$month)
EPFU_data$month<- ifelse(EPFU_data$month == '9', 
                         paste('0', EPFU_data$month, sep = ''), EPFU_data$month)

unique(EPFU_data$month)

EPFU_data$year<- as.character(EPFU_data$year)

EPFU_data$date<- paste(EPFU_data$year, EPFU_data$month, EPFU_data$day, sep = '')
EPFU_data$date<- ymd(EPFU_data$date)

#make a second copy of adult counts
EPFU_surveys<- EPFU_data

#summarize data by unique site_mask and date
EPFU_surveys<- EPFU_surveys %>%
  group_by(site_mask, date) %>%
  dplyr::summarize(.groups = 'drop') %>%
  droplevels()

#order dates by site, create groupings for consecutive nights, summarize the number of nights, and create vars for first and last night  
EPFU_surveys <- EPFU_surveys %>% 
  group_by(site_mask) %>% 
  arrange(date) %>% 
  dplyr::mutate(survey_group = cumsum(c(TRUE, diff(date) > 1))) %>% 
  group_by(survey_group, .add = TRUE) %>% 
  dplyr::summarise(num_nights = n(), start_date = first(date), end_date = last(date), .groups = 'keep')

#work to merge this and full adult count df
colnames(EPFU_data)[19]<- 'start_date'

temp<- left_join(x = EPFU_data, y = EPFU_surveys)

#visualize data by site_mask and start_date
temp<- temp %>%
  group_by(site_mask) %>%
  dplyr::arrange(start_date)

#fill in NAs
temp$end_date <- ave(temp$end_date, cumsum(!is.na(temp$end_date)), FUN = function(x) x[1])
temp$num_nights <- ave(temp$num_nights, cumsum(!is.na(temp$num_nights)), FUN = function(x) x[1])
temp$survey_group <- ave(temp$survey_group, cumsum(!is.na(temp$survey_group)), FUN = function(x) x[1])

#remove and replace proper survey start_date for a final adult_EPFUcount df
temp<- subset(temp, select = -c(survey_group, start_date))
EPFU_data<- left_join(x = temp, y = EPFU_surveys)
EPFU_data<- subset(EPFU_data, select = -survey_group)

#remove temporary df
rm(temp)

#re-fill in start dates where needed
EPFU_data<- EPFU_data %>%
  group_by(site_mask) %>%
  dplyr::arrange(end_date)

EPFU_data$start_date<- ave(EPFU_data$start_date, cumsum(!is.na(EPFU_data$start_date)), FUN = function(x) x[1])

#make survey dates a single column
EPFU_data$survey_dates<- paste(EPFU_data$start_date, ' to ', EPFU_data$end_date)

#transform into capture counts 
EPFUcount<- EPFU_data %>% 
  #select(state, site_mask, survey_dates, num_nights, year, age, sex, repstat, years_Pd, 
  #disease_time_step, county_centroid_lat, county_centroid_lon, rMapCounty)  %>%
  group_by(state, site_mask, survey_dates, num_nights, age, sex, repstat, disease_time_step, 
           county_centroid_lat) %>%
  summarise(freq = n(), .groups = 'keep')

str(EPFUcount)
head(EPFUcount)

adult_EPFUcount<- EPFUcount %>%
  filter(age == 'adult') %>%
  droplevels()

#Adjust variables and levels where needed
adult_EPFUcount$disease_time_step<- factor(adult_EPFUcount$disease_time_step, levels = c('pre-invasion', 'invasion', 'epidemic', 'established'))


#capture rate will be capture counts per site per number of nights per survey
#check distribution of response (capture count per number of net nights)
hist(adult_EPFUcount$freq/adult_EPFUcount$num_nights)

mean(adult_EPFUcount$freq/adult_EPFUcount$num_nights)

var(adult_EPFUcount$freq/adult_EPFUcount$num_nights)

#negative binomial-- should be good given the amount of variation in the data.
#Will use a gamma (since freq/num_nights is a non-integer), and select 
#weakly informative priors with a negative binomial shape for gamma prior (still weak/default)


#Create initial model for adult counts as a function of disease time-step and sex
#create/run model with weakly informative priors for a gamma distribution
adult_count_mod1<- brm(freq/num_nights ~ sex*disease_time_step*county_centroid_lat + (1|site_mask),
                       data = adult_EPFUcount, family = Gamma(link="log"), 
                       prior = c(prior(normal(0, 10), class = Intercept), #normal distribution with mean of 0 and SD of 10                                                                              
                                 prior(normal(0, 1), class = b), #population level effects
                                 prior(gamma(0.01, 0.01), class = shape)), #shape of distribution has most values close to 0 and tail of larger numbers 
                       chains = 4, warmup = 5000, iter = 10000, sample_prior = 'yes', control = list(max_treedepth = 15))



summary(adult_count_mod1)

#Get conditional effects--bayes equivalent of frequentist fixed effects
cond_ef_adult_count1<- emmeans(adult_count_mod1, ~sex*disease_time_step*county_centroid_lat)

cond_ef_adult_count1



#now do same steps but for females only and include repstat
#Get count datasets for females only
fem_EPFUcount<- EPFUcount %>%
  filter(sex == 'female', age == 'adult') %>%
  droplevels()

#Adjust variables and levels where needed
fem_EPFUcount$repstat<- factor(fem_EPFUcount$repstat, levels= c("non-reproductive", "pregnant", "lactating", "post-lactating"))

fem_EPFUcount$disease_time_step<- factor(fem_EPFUcount$disease_time_step, levels = c('pre-invasion', 'invasion', 'epidemic', 'established'))




#Create initial model for females 
#use the same weakly informed priors with a gamma shape repping neg binomial distribution
fem_count_mod1<- brm(freq/num_nights ~ repstat*disease_time_step*county_centroid_lat + (1|site_mask),
                     data = fem_EPFUcount, family = Gamma(link="log"), 
                     prior = c(prior(normal(0, 10), class = Intercept), #normal distribution with mean of 0 and SD of 10                                                                              
                               prior(normal(0, 1), class = b), #population level effects
                               prior(gamma(0.01, 0.01), class = shape)), #shape of distribution has most values close to 0 and tail of larger numbers 
                     chains = 4, warmup = 5000, iter = 10000, sample_prior = 'yes', control = list(max_treedepth = 15))


summary(fem_count_mod1)


#Get conditional effects--which is really a post-hoc type thing for the bayes equivalent of fixed effects
cond_ef_fem_count1<- emmeans(fem_count_mod1, ~county_centroid_lat*disease_time_step*repstat)

cond_ef_fem_count1



#Now create and test secondary models for all adults and females only, replacing latitude with north/south of interaction point
#first for all adults
#add latitudinal category variable from first model outputs in step 1 and latitudinal category weights var
adult_EPFUcount$lat_cat<- ifelse(adult_EPFUcount$county_centroid_lat > 39.5, 'north', 'south')



#Create second model for adult counts as a function of disease time-step, sex and latitudinal category
adult_count_mod2<- brm(freq/num_nights ~ sex*disease_time_step*lat_cat + (1|site_mask), 
                       data = adult_EPFUcount, family = Gamma(link="log"), 
                       prior = c(prior(normal(0, 10), class = Intercept), #normal distribution with mean of 0 and SD of 10                                                                              
                                 prior(normal(0, 1), class = b), #population level effects
                                 prior(gamma(0.01, 0.01), class = shape)), #shape of distribution has most values close to 0 and tail of larger numbers 
                       chains = 4, warmup = 5000, iter = 10000, sample_prior = 'yes', control = list(max_treedepth = 15))


summary(adult_count_mod2)

#Get conditional effects--bayes equivilent of frequentist fixed effects
cond_ef_adult_count2<- emmeans(adult_count_mod2, ~sex*disease_time_step*lat_cat)

cond_ef_adult_count2

#try stan_code() so it spits out your betas and all the pre model parameters that brms eliminates the process of
stancode(adult_count_mod2)


#now for females only
#add latitudinal category variable 
fem_EPFUcount$lat_cat<- ifelse(fem_EPFUcount$county_centroid_lat > 39.5, 'north', 'south')



#Create second model for females 
fem_count_mod2<- brm(freq ~ lat_cat*repstat*disease_time_step + (1|site_mask), 
                     data = fem_EPFUcount, family = Gamma(link="log"), 
                     prior = c(prior(normal(0, 10), class = Intercept), #normal distribution with mean of 0 and SD of 10                                                                              
                               prior(normal(0, 1), class = b), #population level effects
                               prior(gamma(0.01, 0.01), class = shape)), #shape of distribution has most values close to 0 and tail of larger numbers 
                     chains = 4, warmup = 5000, iter = 10000, sample_prior = 'yes', control = list(max_treedepth = 15))



summary(fem_count_mod2)

#Get conditional effects
cond_ef_fem_count2<- emmeans(fem_count_mod2, ~lat_cat*disease_time_step*repstat)

cond_ef_fem_count2

#try stan_code() so it spits out your betas and all the pre model parameters that brms eliminates the process of
stancode(fem_count_mod2)

#save environment for plotting figures later
save.image(file='BigBrownBat_CaptureRate_Bayes_Cluster.RData')