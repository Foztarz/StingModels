# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2023 04 26
#     MODIFIED:	James Foster              DATE: 2023 09 24
#
#  DESCRIPTION: Toy example of a binomial experiment spread over 4 trials
#               with consistent individual biases
#               
#       INPUTS: Starting parameters for the simulation
#               file to load.
#               
#      OUTPUTS: Saves model summary
#
#	   CHANGES: - 
#
#   REFERENCES: Bates D., Maechler M., Bolker B. & Walker S. (2015). 
#               Fitting Linear Mixed-Effects Models Using lme4.
#               Journal of Statistical Software, 67(1), 1-48. doi:10.18637/jss.v067.i01.
# 
#       USAGE:  
#TODO   ---------------------------------------------
#TODO   
#- Set up simulation  +
#- Test analysis      +
#- Individual slopes
#- Plotting
#- Nonlinear version?



# Load packages -----------------------------------------------------------
require(lme4)#package for modelling linear "mixed-effects" models
# mixed-effects models include a mixture of both population and group-level effects
# a population effect could be a change in probability proportional to some measured variable
# a group effect could a bias specific to some group, e.g. all trials by the same animal


# Seed random number generator --------------------------------------------
#simulation should be repeatable
set.seed(20230426)# date this script was started

# Starting parameters -----------------------------------------------------
n_indiv = 80 # total individuals
n_trials = 4 # trials per individual
n_treat = 2 # number of treatments (different individuals)

dt_sim = data.frame(indiv = sort(rep(x = 1:n_indiv, 
                                     times = n_trials)), #n_indiv each perform n_trials, sort by indiv number
                    treat = rep(x = sapply(X = 1:n_treat,
                                           FUN = rep,
                                           times = n_trials), #there are n treatments, shared equally across n_indiv
                                times = n_indiv/n_trials),
                    trial = rep(x = 1:n_trials,
                                times = n_indiv) #each indiv performs trials 1 to n
                    )

# Simulate curve -----------------------------------------------------------
zero_crossing = 2#  probability for treatment 1 would reach 0.5 on the 5th trial
pop_trial_slope = -1.20 #change in log odds with every trial (70% reduction ≈ 1-exp(-1.20))
pop_treat_slope = 1.10 #change in log odds with treatment 2 (tripled odds ≈ exp(1.10))
pop_treat_trial_slope = 0.30 #change in trial effect with treatment 2 (quite large)
pop_intercept = -zero_crossing*pop_trial_slope # set intercept to match zero crossing (inflection point of sigmoid)
ind_intercept_sd = 0.5 #standard deviation of individual intercepts (log odds)
ind_trial_sd = 0.3 #standard deviation of individual trial slopes
residual_sd = 0.05#small source of unnaccounted error
# ind_treat_sd = 0.1 #standard deviation of individual treatment slopes #indistinguishable from intercept
ind_intercepts = rnorm(n = 1:n_indiv,
                       mean = 0,
                       sd = ind_intercept_sd)
ind_trial_bias = rnorm(n = 1:n_indiv,
                       mean = 0,
                       sd = ind_trial_sd)
#simulate continuous trial sequence
xx = seq(from  = 0, to  = n_trials, length.out = 1e3)
#simulate log odds for treatment 1
yy_lodds = pop_trial_slope*xx + pop_intercept 
#convert to average choice probability
yy_prob = plogis(yy_lodds)
#convert to probability
#simulate log odds for treatment 2
yy_lodds2 = (pop_trial_slope + pop_treat_trial_slope)*xx + 
              pop_intercept +pop_treat_slope
#convert to average choice probability
yy_prob2 = plogis(yy_lodds2)

#plot function
plot(x = xx, 
     y = yy_prob,
     ylim = c(0,1),
     type = 'l',
     col = 2,
     lwd = 3,
     ylab = 'probability',
     xlab = 'trial')
lines(x = xx,
      y = yy_prob2,
      col = 7,
      lwd = 3)
abline(h = c(1,0.5,0),
       lty = c(1,3,1))
abline(v = c(zero_crossing, #original p(50%)
             -(pop_intercept + pop_treat_slope) / 
               (pop_trial_slope + pop_treat_trial_slope)), #treatment 2 p(50%)
       lty = 3,
       col = c(2,3))


# . Generate simulated data -----------------------------------------------

# #Make an empty vector
# response_y = rep(x = NA,
#                  length = n_indiv * n_trials)

dt_sim = within(dt_sim,
              {
                 #each data point is generated from the animal's intercept
                 #plus the effect of treatment
                 #plus the effect of successive trials
                 #plus the animal's response (slope) to changing treatment
                 #scaled by the effect of successive trials
                 #plus a small additional source or error
                 
       response_y = pop_intercept + ind_intercepts[indiv] +
                    (trial) * (pop_trial_slope + 
                                 (treat-1) * pop_treat_trial_slope) + 
                     ind_trial_bias[indiv] +
                    (treat-1) *  +pop_treat_slope
                    rnorm(n = n_indiv * n_trials,
                          mean = 0,
                          sd = residual_sd)
              }
              )
#convert to probability and simulate binomial
dt_sim = within(dt_sim,
              {
                #each response on the log-odds scale gives us
                #response probability on a "logit" scale
                 p_response = plogis(response_y)
                 #that response probabilty determines the sting event
                 #which we observe once for each trial
                 sting = rbinom(n = n_indiv * n_trials, #total trial number
                                size = 1, # one observation per trial
                                prob = p_response)
              }
              )

# . Inspect simulated data ------------------------------------------------

# View(dt_sim)
summary(dt_sim)

#plot the expected probabilities for each individual (after added noise)
plot(x = xx, 
     y = yy_prob,
     ylim = c(0,1),
     type = 'l',
     col = 2,
     lwd = 3,
     ylab = 'probability',
     xlab = 'trial')
lines(x = xx,
      y = yy_prob2,
      col = 7,
      lwd = 3)
with(subset(x = dt_sim,
            subset = treat == 1),
     points(x = trial,
            y = p_response,
            col = 2,
            pch = 3)
      )

with(subset(x = dt_sim,
            subset = treat == 2),
     points(x = trial,
            y = p_response,
            col = 7,
            pch = 4)
      )

abline(h = c(1,0.5,0),
       lty = c(1,3,1))
abline(v = c(zero_crossing, #original p(50%)
             -(pop_intercept + pop_treat_slope) / 
               (pop_trial_slope + pop_treat_trial_slope)), #treatment 2 p(50%)
       lty = 3,
       col = c(2,3))


# . aggregate simulated data to summarise ---------------------------------
dt_agg = aggregate(sting ~ trial*treat,
                   data = dt_sim,
                   FUN = mean)


#plot the expected probabilities for each individual (after added noise)
plot(x = xx, 
     y = yy_prob,
     ylim = c(0,1),
     type = 'l',
     col = 2,
     lwd = 3,
     ylab = 'probability',
     xlab = 'trial')
lines(x = xx,
      y = yy_prob2,
      col = 7,
      lwd = 3)
with(subset(x = dt_agg,
            subset = treat == 1),
     points(x = trial,
            y = sting,
            col = 2,
            pch = 3)
      )

with(subset(x = dt_agg,
            subset = treat == 2),
     points(x = trial,
            y = sting,
            col = 7,
            pch = 4)
      )


# Load packages required for mixed-effects modelling ----------------------
require(lme4)
# Fit a generalised linear model ------------------------------------------

#set some useful optimiser settings for models with many paramters
ctrl_opt = glmerControl(optimizer = 'bobyqa') #efficient optimiser for models with many parameters

#Maximal model with random intercepts for individuals
glmm.max = glmer(formula = sting~
                    treat * trial +
                    (1 + treat * trial | indiv),
                 data = dt_sim,
                 control = ctrl_opt,
                 family = binomial(link = 'logit')
)
#boundary (singular) fit: see ?isSingular
#Some random effects are too small to estimate properly, common warning
#Null model, with only random effects
glmm.null = glmer(formula = sting~
                     1 + (1|indiv),
                  data = dt_sim,
                  family = binomial(link = 'logit')
)
# Model comparison --------------------------------------------------------
#To prove any fixed effects, the mixed-effects model has to
#describe the data better than the random effects model alone.
extractAIC(glmm.max)[2]#2nd component is the AIC, 1st is the d.f.
extractAIC(glmm.null)[2]
#If AIC is lower for the maximal model, then the maximal model fits
# 329.1871 < 444.2998
#We can also perform a likelihood ratio test, confusingly called "anova"
anova(glmm.max, glmm.null, test = 'Chisq')
#This gives the same answer as the AIC, but the difference can be reported with
#a test statistic and p-value
##           npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
## glmm.null     2 444.30 451.84 -220.15   440.30                                        
## glmm.max     14 329.19 381.94 -150.59   301.19 139.11 12  < 2.2e-16 ***

#report:
#change in deviance = 139.11   
#change in degrees of freedom = 12
#p < 0.0001
#(often p < 2.2e-16; the smallest number the computer can think of)

# Check model reduction options -------------------------------------------
#If we are not sure that we need all of our parameters, we might consider model reduction.
#This is not always necessary, and without a good rationale it may be preferable
#to use the maximal model, which controls for all possible combinations of parameters.

#can we remove one or more fixed effect to improve the model fit?

#effect of stimulus-type interaction
no_int = update(object = glmm.max, 
                .~. - treat:trial - #remove the interaction
                  (1 + treat * trial | indiv) + #remove random effects that include it
                  (1 + treat + trial | indiv) #replace them with random effect that exclude it
)
anova(glmm.max,
      no_int)

# no_int: sting ~ treat + trial + (1 + treat + trial | indiv)
# glmm.max: sting ~ treat * trial + (1 + treat * trial | indiv)
#           npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)   
# no_int       9 320.70 354.62 -151.35   302.70          
# glmm.max    14 329.19 381.94 -150.59   301.19 1.5158  5     0.9112
#lower AIC following the parameter's removal indicates 
#that the model without an interaction is a better fit.
#However, the difference is not significant, so there is
#no significant evidence for or against an interaction.

#We simulated quite a substantial effect 
#(change in log odds of 0.30, an exp(0.30) ≈ 1.35 x increase in relative sting prob per trial),
#but it appears hard to detect in this size of dataset.


#We can also use this method to report the strength of our fixed effects
#we would report this as:
#likelihood ratio test, change in deviance  = 1.5158, d.f. = 5, p = 0.9112

#effect of treatment
no_treatment = update(object = no_int, 
                     .~. - treat - #remove the effect of treatment
                       (1 + treat + trial | indiv) + #remove random effects that include it
                       (1 + trial | indiv) #replace them with random effect that exclude it
)
anova(no_int,
      no_treatment)
##              npar   AIC   BIC   logLik deviance  Chisq Df Pr(>Chisq)    
## no_treatment    5 351.82 370.66 -170.91   341.82                          
## no_int          9 320.70 354.62 -151.35   302.70 39.114  4  6.598e-08 ***

#larger AIC following the parameter's removal indicates 
#that the model without an effect of treatment is a significantly poorer fit
#therefore, there is a significant effect of treatment
#we would report this as:
#likelihood ratio test, change in deviance  = 39.114, d.f. = 4, p < 0.001

#effect of trial
no_trial = update(object = no_treatment, 
                 .~. - trial - #remove the effect of type
                   (1 + trial | indiv) + #remove random effects that include it
                   (1 | indiv) #replace them with random effect that exclude it
)
anova(no_treatment,
      no_trial)
##              npar   AIC   BIC   logLik deviance  Chisq Df Pr(>Chisq)    
## no_trial        2 444.30 451.84 -220.15   440.30                       
## no_treatment    5 351.82 370.66 -170.91   341.82 98.483  3  < 2.2e-16 ***

#larger AIC following the parameter's removal indicates 
#that the model without an effect of trial is a significantly poorer fit
#therefore, there is a significant effect of trial
#we would report this as:
#likelihood ratio test, change in deviance  = 98.483, d.f. = 3, p < 0.001


# Discarded ---------------------------------------------------------------
# #Loop through experiment design and generate theoretical response
# for (ii in 1:length(response_y))
# {
#    #each data point is generated from the animal's intercept
#    #plus the effect of treatment
#    #plus the effect of successive trials
#    #plus the animal's response (slope) to changing treatment
#    #scaled by the effect of successive trials
#    #plus a small additional source or error
#    response_y[ii] = ind_intercepts[animal[ii]] +
#       (stype[ii] - 1) * type_mean +
#       (slope_animal[animal[ii]] +
#           (stype[ii] - 1) * type_slope * type_animal[animal[ii]]) * stimulus[ii] +
#       rnorm(n = 1,
#             mean = 0,
#             sd = residual_sd)
# }
# #convert to binomial
# correct_incorrect = rbinom(n = length(response_y),
#                            size = 1,
#                            prob = plogis(q = response_y - mean(response_y)) )
# 
# 
# # Combine this into a single data.frame format object,
# # with the experimental design and the measured response data
# sim_data = data.frame(
#    correct_incorrect = correct_incorrect,
#    animal = LETTERS[animal],
#    #animal names will now be capital letters
#    stimulus = stimulus,
#    type = factor(stype,
#                  labels = c('alpha',
#                             'beta'))
# )


