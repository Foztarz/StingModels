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
#- Set up simulation
#- Test analysis
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
n_indiv = 40 # total individuals
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
pop_trial_slope = -0.9 #change in log odds with every trial
pop_treat_slope = 1.10 #change in log odds with treatment 2 (tripled odds ≈ exp(1.1))
pop_treat_trial_slope = 0.1 #change in trial effect with treatment 2
pop_intercept = -zero_crossing*pop_trial_slope # set intercept to match zero crossing (inflection point of sigmoid)
ind_intercept_sd = 1.0 #standard deviation of individual intercepts (log odds)
ind_trial_sd = 0.5 #standard deviation of individual trial slopes
residual_sd = 0.10#small source of unnaccounted error
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

#Make an empty vector
response_y = rep(x = NA,
                 length = n_indiv * n_treat * n_trials)

dt_sim = within(dt_sim,
              {
                 #each data point is generated from the animal's intercept
                 #plus the effect of treatment
                 #plus the effect of successive trials
                 #plus the animal's response (slope) to changing treatment
                 #scaled by the effect of successive trials
                 #plus a small additional source or error
                 
                 #WIP SOMETHING IS WRONG HERE
                 response_y = pop_intercept + ind_intercepts[indiv] +
                    (treat - 1) * pop_treat_slope +
                    (trial - 1) * pop_trial_slope + ind_trial_bias[indiv] +
                    (trial - 1) * pop_treat_trial_slope +
                    rnorm(n = n_indiv * n_treat * n_trials,
                          mean = 0,
                          sd = residual_sd)
              }
              )
#Loop through experiment design and generate theoretical response
for (ii in 1:length(response_y))
{
   #each data point is generated from the animal's intercept
   #plus the effect of treatment
   #plus the effect of successive trials
   #plus the animal's response (slope) to changing treatment
   #scaled by the effect of successive trials
   #plus a small additional source or error
   response_y[ii] = ind_intercepts[animal[ii]] +
      (stype[ii] - 1) * type_mean +
      (slope_animal[animal[ii]] +
          (stype[ii] - 1) * type_slope * type_animal[animal[ii]]) * stimulus[ii] +
      rnorm(n = 1,
            mean = 0,
            sd = residual_sd)
}
#convert to binomial
correct_incorrect = rbinom(n = length(response_y),
                           size = 1,
                           prob = plogis(q = response_y - mean(response_y)) )


# Combine this into a single data.frame format object,
# with the experimental design and the measured response data
sim_data = data.frame(
   correct_incorrect = correct_incorrect,
   animal = LETTERS[animal],
   #animal names will now be capital letters
   stimulus = stimulus,
   type = factor(stype,
                 labels = c('alpha',
                            'beta'))
)
View(sim_data)


# Load packages required for mixed-effects modelling ----------------------
require(lme4)
# Fit a generalised linear model ------------------------------------------

#set some useful optimiser settings for models with many paramters
ctrl_opt = glmerControl(optimizer = 'bobyqa') #efficient optimiser for models with many parameters

#Maximal model with random intercepts for individuals
glmm.max = glmer(formula = correct_incorrect~
                    stimulus * type +
                    (1 + stimulus * type | animal),
                 data = dta,
                 control = ctrl_opt,
                 family = binomial(link = 'logit')
)
#boundary (singular) fit: see ?isSingular
#Some random effects are too small to estimate properly, common warning
#Null model, with only random effects
glmm.null = glmer(formula = correct_incorrect~
                     1 + (1|animal),
                  data = dta,
                  family = binomial(link = 'logit')
)
# Model comparison --------------------------------------------------------
#To prove any fixed effects, the mixed-effects model has to
#describe the data better than the random effects model alone.
extractAIC(glmm.max)[2]#2nd component is the AIC, 1st is the d.f.
extractAIC(glmm.null)[2]
#If AIC is lower for the maximal model, then the maximal model fits
# 685.9584 < 3587.245
#We can also perform a likelihood ratio test, confusingly called "anova"
anova(glmm.max, glmm.null, test = 'Chisq')
#This gives the same answer as the AIC, but the difference can be reported with
#a test statistic and p-value
##           npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
## glmm.null    2 3587.2 3599.1 -1791.62   3583.2                         
## glmm.max    14  686.0  769.1  -328.98    658.0 2925.3 12  < 2.2e-16 ***

#report:
#change in deviance = 2925.3 
#change in degrees of freedom = 12
# p < 2.2e-16 (the smallest number the computer can think of)
