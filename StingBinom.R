# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2023 04 26
#     MODIFIED:	James Foster              DATE: 2023 12 14
#
#  DESCRIPTION: Toy example of a binomial experiment spread over 4 trials
#               with consistent individual biases
#               
#       INPUTS: Starting parameters for the simulation
#               file to load.
#               
#      OUTPUTS: Saves model summary
#
#	   CHANGES: - post-hoc tests
#	            - plot predictions
#
#   REFERENCES: Bates D., Maechler M., Bolker B. & Walker S. (2015). 
#               Fitting Linear Mixed-Effects Models Using lme4.
#               Journal of Statistical Software, 67(1), 1-48. doi:10.18637/jss.v067.i01.
#               
#               Lenth R. V. (2021). emmeans: Estimated Marginal Means, 
#               AKA Least-Squares Means. R package version 1.5.4. 
#               https://cran.r-project.org/web/packages/emmeans/vignettes/basics.html
#
#               Šidák, Z. K. (1967). "Rectangular Confidence Regions for the 
#               Means of Multivariate Normal Distributions". J. Am. Statist. Assoc. 
#               62 (318): 626–633. doi:10.1080/01621459.1967.10482935.
#       USAGE:  
#TODO   ---------------------------------------------
#TODO   
#- Set up simulation  +
#- Test analysis      +
#- Individual slopes  +
#- Plotting +
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
zero_crossing = 1.5#  probability for treatment 1 would reach 0.5 between the 1st and 2nd trial
pop_trial_slope = -1.20 #change in log odds with every trial (70% reduction ≈ 1-exp(-1.20))
pop_treat_slope = 1.10 #change in log odds with treatment 2 (tripled odds ≈ exp(1.10))
pop_treat_trial_slope = 0.45 #change in trial effect with treatment 2 (quite large)
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
# 323.4593 < 439.8845
#We can also perform a likelihood ratio test, confusingly called "anova"
anova(glmm.max, glmm.null, test = 'Chisq')
#This gives the same answer as the AIC, but the difference can be reported with
#a test statistic and p-value
##           npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
## glmm.null     2 439.88 447.42 -217.94   435.88                                        
## glmm.max     14 323.46 376.22 -147.73   295.46 140.43 12  < 2.2e-16 ***

#report:
#change in deviance = 140.43   
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
# no_int       9 314.97 348.88 -148.49   296.97           
# glmm.max    14 323.46 376.22 -147.73   295.46 1.5103  5     0.9119
#lower AIC following the parameter's removal indicates 
#that the model without an interaction is a better fit.
#However, the difference is not significant, so there is
#no significant evidence for or against an interaction.

#We simulated quite a substantial effect 
#(change in log odds of 0.45, an exp(0.45) ≈ 1.57 x increase in relative sting prob per trial),
#but it appears hard to detect in this size of dataset.


#We can also use this method to report the strength of our fixed effects
#we would report this as:
#likelihood ratio test, change in deviance  = 1.5103, d.f. = 5, p = 0.9119

#effect of treatment
no_treatment = update(object = no_int, 
                     .~. - treat - #remove the effect of treatment
                       (1 + treat + trial | indiv) + #remove random effects that include it
                       (1 + trial | indiv) #replace them with random effect that exclude it
)
anova(no_int,
      no_treatment)
##              npar   AIC   BIC   logLik deviance  Chisq Df Pr(>Chisq)    
## no_treatment    5 359.78 378.62 -174.89   349.78                           
## no_int          9 314.97 348.88 -148.49   296.97 52.811  4  9.332e-11 ***

#larger AIC following the parameter's removal indicates 
#that the model without an effect of treatment is a significantly poorer fit
#therefore, there is a significant effect of treatment
#we would report this as:
#likelihood ratio test, change in deviance  = 52.811, d.f. = 4, p < 0.001

#effect of trial
no_trial = update(object = no_treatment, 
                 .~. - trial - #remove the effect of type
                   (1 + trial | indiv) + #remove random effects that include it
                   (1 | indiv) #replace them with random effect that exclude it
)
anova(no_treatment,
      no_trial)
##              npar   AIC   BIC   logLik deviance  Chisq Df Pr(>Chisq)    
## no_trial        2 439.01 446.55 -217.50   435.01                       
## no_treatment    5 359.78 378.62 -174.89   349.78 85.229  3  < 2.2e-16 ***

#larger AIC following the parameter's removal indicates 
#that the model without an effect of trial is a significantly poorer fit
#therefore, there is a significant effect of trial
#we would report this as:
#likelihood ratio test, change in deviance  = 85.229, d.f. = 3, p < 0.001


# Post-hoc comparisons ----------------------------------------------------
#Load the package for post-hoc comparisons
require(emmeans)
require(pbkrtest)# Used for calculating DF in mixed-effects model

# Calculate overall discrete effects
emm_intercepts = emmeans(glmm.max,#we'll use the maximal model, as model selection was not clearly against it
                         specs = list(pairwise~treat),
                         #compare the overall effects of different types
                         adjust = 'sidak')
summary(emm_intercepts)
## $`pairwise differences of treat`
## 1               estimate    SE  df z.ratio p.value
## treat1 - treat2    -3.16 0.661 Inf  -4.782  <.0001
#As for the likelihood ratio test, we see that there is a large effect of 
#treatment on probability of a response (as there should be).
#The estimated difference is that for treat1 the odds of a correct response are 
#16 log units lower than for beta for the same stimulus intensity.
#That is an odds ratio of exp(3.16) = 24 to 1 (treat2 vs treat1)

# Calculate post-hoc comparisons for continuous effects (slope)
emm_slopes_interact = emtrends(glmm.max,
                               specs = list(pairwise~treat),
                               var = 'trial',
                               #compare the effects of different types on response/stimulus slopes
                               adjust = 'sidak')
summary(emm_slopes_interact)
## `pairwise differences of treat`
## 1               estimate    SE  df z.ratio p.value
## treat1 - treat2   -0.366 0.555 Inf  -0.658  0.5103
#as for the model comparison, the model finds not significant difference in the trial-response relationship 
#between the two treatments. This is unlike what was specified in our input variables (treat_slope) 
#for the simulation. The estimated difference is that treat1 has an average 
#trial-response slope that is weaker by "-0.366". That is for every trial 
#the odds of a response decrease 1.44x faster
#for treatment 1 than treatment 2.


# Plot predictions --------------------------------------------------------

# . Extract predictions ---------------------------------------------------
#check all relevant variables for predictions
formula(glmm.max)
## sting ~ treat * trial + (1 + treat * trial | indiv)
newdta = with(dt_sim,
              expand.grid(trial = seq(from = min(trial),
                                         to = max(trial),
                                         length.out = length(trial)),
                          treat = unique(treat),
                          indiv = unique(indiv)
              ) 
)
#predictions (mean estimate)
prd = predict(glmm.max,
              newdata = newdta,
              type = 'response'
)
#bootstrap the confidence intervals (can take a while...)
pfun = function(x)
{
  predict(x,
          newdata = newdta, 
          type = "response")
}
# library(snow)#parallel processing requires "snow" on Windows
#Really benefits from some parallel processing
avail.cores = parallel::detectCores() - 1
clt = parallel::makeCluster(spec = avail.cores, 
                            type=if(Sys.info()[['sysname']] == 'Windows'){"SOCK" # 'sock' type on windows
                            }else{"FORK"}#'fork' type on mac & linux
)
parallel::clusterExport(clt,
                        list('glmm.max',
                             'dt_sim',
                             'newdta'
                        )
)
#this takes a long time to simulate
message('starting large simulation\n',
        'expect to wait at least 10 minutes...\n'
)
system.time({
  bt = bootMer(glmm.max,
               FUN = pfun,
               nsim = 100,#100 takes ≈15 minutes. Minimum of 20 to be able to calculate 95%CI. Increase number for greater detail.
               re.form = NULL,#NA for fixed effects, NULL to include random effects
               #fixed effects gives "confidence interval", population level effects
               #random effects gives "prediction interval", expected values for these individuals
               #on the logit scale, we may be able to improve precision by including individuals in our prediction
               parallel = ifelse(test = Sys.info()[['sysname']] == 'Windows',
                                 yes =  "snow",
                                 no =  "multicore"),
               ncpus = parallel::detectCores()-1, #leave one processor available for user
               cl = clt #the parallel cluster prepared above
  )
})
message('simulation finished!')
#now it has been used, close the cluster
parallel::stopCluster(clt)
#To get CI for only fixed effects, align and aggregate bootstrap estimates
pred_q = aggregate(pred ~ treat * trial, #aggregate predictions by fixed effect
                   FUN = quantile, #calculate as quantiles of bootstrap predictions
                   data = with(bt, 
                               data.frame(newdta, 
                                          pred = c(t(t)) #convert from rows to columns and align 
                               ) 
                   ),
                   probs = c(0,1) + c(1,-1)*0.05/2) #using alpha = 0.05, make two-tailed confidence intervals
# find mean prediction across 
mod_mean = aggregate(prd ~ treat * trial, 
                     data = cbind(newdta, prd), 
                     FUN = function(x){plogis(mean(qlogis(x)))})
#merge together
param_data = merge(merge(newdta, pred_q, all = TRUE), 
                   mod_mean)
#rename for plotting
mod_pred = within(param_data,
                  {
                    mod_mean = prd
                    CI_02.5 = pred[,'2.5%']
                    CI_97.5 = pred[,'97.5%']
                    rm(list = c('pred','prd'))
                  }
)

# . Plot predictions ------------------------------------------------------

#plot all data together
with(dt_agg,
     {
       plot(x = trial,
            y = sting,
            bg = adjustcolor(col = c('red4', 'orange2')[ 1 + treat %in% '2'],
                             alpha.f = 0.5), # 50% opacity
            col = 'black',
            pch = 21,
            ylim = c(0,1)) # dots
     }
)
#add model and shaded confidence intervals
#treatment 1
with(subset(mod_pred, treat %in% '1'),
     {
       polygon(x = c(sort(trial), sort(trial,decreasing = TRUE)),
               y = 
                 c(lowerCI = CI_02.5[order(trial)],
                   upperCI = CI_97.5[order(trial,decreasing = TRUE)]
                 ),
               col = adjustcolor(col = 'red',
                                 alpha.f = 0.2),
               border = NA
       )
       lines(sort(trial),
             mod_mean[order(trial)],#
             col = 'red',
             pch = 10
       )
     }
)
#stimulus type beta
with(subset(mod_pred, treat %in% '2'),
     {
       polygon(x = c(sort(trial), sort(trial,decreasing = TRUE)),
               y = 
                 c(lowerCI = CI_02.5[order(trial)],
                   upperCI = CI_97.5[order(trial,decreasing = TRUE)]
                 ),
               col = adjustcolor(col = 'orange2',
                                 alpha.f = 0.2),
               border = NA
       )
       lines(sort(trial),
             mod_mean[order(trial)],
             col = 'orange2',
             pch = 10
       )
     }
)


