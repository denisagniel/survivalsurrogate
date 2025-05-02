

#required packages, install first

library(tidyverse)
library(here)
library(glue)
library(simcausal)
library(devtools)
devtools::install_github("denisagniel/survivalsurrogate")
library(survivalsurrogate)
library(mlr3verse)
#remotes::install_github("mlr-org/mlr3extralearners@*release")
library(mlr3extralearners)
library(quantreg)
library(gbm)

#########################
## Set up the simulation parameters
###################################

#setting 1 is R=0; setting 2 is R=1, setting 3 is R in between
setting = 3
######################

n <- 1000
K <- 3
n_sim <- 1000
tt <- 5
t0 <- 4
yvars <- paste0('Y_', 0:tt)

################################
## Set up  learners to use
###################################################
lrn <- 'gbm'
lrnc <- glue('regr.{lrn}')
lrnb <- glue('classif.{lrn}')

#######################
## Set up some mean functions
###############################
if(setting == 1) {
	yfn <- function(s_tm1, g, x) plogis(-2 +0.5*s_tm1 - 1*g + 0.3*x)
	sfn <- function(s_tm1, g, x) -0.1*g + 0.5*x + 0.25*s_tm1
	gfn <- function(x) plogis(x)
}
if(setting == 2) {
		
	yfn <- function(s_tm1, g, x) plogis(-5 + 4.5*s_tm1 - 0.05*g*s_tm1 - 0.05*g+ 0.3*x)
	sfn <- function(s_tm1, g, x) -0.5*g +0.5*x + 0.25*s_tm1
	gfn <- function(x) plogis(x)
}
if(setting == 3) {
	yfn <- function(s_tm1, g, x) plogis(-5 + 4*s_tm1 - 0.1*g*s_tm1 - 1*g + 0.3*x)
	sfn <- function(s_tm1, g, x) -0.5*g + 0.5*x + 0.25*s_tm1
	gfn <- function(x) plogis(x)
}

rnorm_mix <- function(n, mu1, mu2, sd1 = 1, sd2 = 1, p_mix) {
  index <- rbinom(n, size = 1, prob = p_mix)
  n1 <- sum(index)
  n2 <- n - n1
  out <- rep(0, n)
  if (n1 > 0) 
    out[index == 1] <- rnorm(n1, mu1, sd1)
  if (n2 > 0) 
    out[index == 0] <- rnorm(n2, mu2, sd2)
  out
}
dnorm_mix <- function(x, mu1, mu2, sd1 = 1, sd2 = 1, p_mix) {
  # browser()
  p_mix*dnorm(x, mean = mu1, sd = sd1) + (1-p_mix)*dnorm(x, mean = mu2, sd = sd2)
}

#######################################
## what do the weighting functions look like?
#######################################
f_s0_g0 <- function(s, x) dnorm(s, sfn(0, 0, x))
f_s0_g1 <- function(s, x) dnorm(s, sfn(0, 1, x))
f_s0_both <- function(s, x) dnorm_mix(s, mu1 = sfn(0, 1, x), mu2 = sfn(0, 0, x), p_mix = gfn(x))

x.temp=rnorm(10000,0,1)
g.temp = gfn(x.temp)
s.temp = rnorm(10000, sfn(0, g.temp, x.temp))
sx_grid <- expand_grid(s = seq(quantile(s.temp,0.05), quantile(s.temp,0.95), length = 21),
                       x = seq(quantile(x.temp,0.05), quantile(x.temp,0.95), length = 21)) %>%

  mutate(mu_s0 = sfn(0, 0, x),
         mu_s1 = sfn(0, 1, x),
         f_00 = f_s0_g0(s, x),
         f_10 = f_s0_g1(s, x),
         e = gfn(x),
         f_b0 = f_s0_both(s, x),
         wt_00 = f_b0/f_00,
         wt_10 = f_b0/f_10)




######################
## Specify the DAG
######################
D <- DAG.empty() + 
  node("X", t=0, distr="rnorm") +
  node("G", t=0, distr="rbern", prob = gfn(X_0)) +
  node("Y", t=0, distr="rbern", prob=yfn(0, G_0, X_0), EFU = TRUE) +
  node("S", t=0, distr="rnorm", mean=sfn(0, G_0, X_0)) +
  node("Y", t=1:(t0+1), distr="rbern", prob=yfn(S[t-1], G_0, X_0), EFU=TRUE) +
  node("S", t=1:t0, distr="rnorm", mean=sfn(S[t-1], G_0, X_0)) 
if (tt > t0+1) {
  D <- D +
    node("Y", t=(t0+2):tt, distr="rbern", prob=yfn(S[t0], G_0, X_0), EFU=TRUE) 
}
D_0 <- set.DAG(D, vecfun = c('yfn', 'sfn', 'gfn'))

#######################
## Evaluate the true parameters
##################################
s_a1 <- c(
  node('G', t = 0, distr = 'rbern', prob = 1),
  node('S', t = 0, distr = 'rnorm_mix', 
       mu1 = sfn(0, 1, X_0), mu2 = sfn(0, 0, X_0), p_mix = gfn(X_0), sd1 = 1, sd2 = 1),
  node('S', t = 1:t0, distr = 'rnorm_mix', 
       mu1 = sfn(S[t-1], 1, X_0), mu2 = sfn(S[t-1], 0, X_0), p_mix = gfn(X_0), sd1 = 1, sd2 = 1)
)
s_a0 <- c(
  node('G', t = 0, distr = 'rbern', prob = 0),
  node('S', t = 0, distr = 'rnorm_mix', 
       mu1 = sfn(0, 1, X_0), mu2 = sfn(0, 0, X_0), p_mix = gfn(X_0), sd1 = 1, sd2 = 1),
  node('S', t = 1:t0, distr = 'rnorm_mix', 
       mu1 = sfn(S[t-1], 1, X_0), mu2 = sfn(S[t-1], 0, X_0), p_mix = gfn(X_0), sd1 = 1, sd2 = 1)
)

set_a1 <- c(
  node('G', t = 0, distr = 'rbern', prob = 1)
)
set_a0 <- c(
  node('G', t = 0, distr = 'rbern', prob = 0)
)

Delta_S_DAG <- D_0 + 
  action('s_a1', nodes = s_a1) +
  action('s_a0', nodes = s_a0) 
Delta_S_DAG <- set.targetE(Delta_S_DAG, outcome = 'Y', t = tt, param = 's_a1 - s_a0')
Delta_S <- eval.target(Delta_S_DAG, n = 1000000)$res

Delta_DAG <- D_0 + 
  action('set_a1', nodes = set_a1) +
  action('set_a0', nodes = set_a0) 
Delta_DAG <- set.targetE(Delta_DAG, outcome = 'Y', t = tt, param = 'set_a1 - set_a0')
Delta <- eval.target(Delta_DAG, n = 1000000)$res

R_0 <- 1 - Delta_S/Delta
Delta
Delta_S
R_0

#########################
## Make new function to give both boot and asymptotic
###################################

estimate_R = function (delta_if, delta_s_if, delta = NULL, delta_s = NULL, 
          se_type = "asymptotic", n_boot = 1000, alpha = 0.05) 
{
  if (is.null(delta)) 
    delta <- mean(delta_if)
  if (is.null(delta_s)) 
    delta_s <- mean(delta_s_if)
  n <- length(delta_if)
  sigma_sq <- delta^(-2) * mean(delta_if^2) + delta_s^2 * delta^(-4) * 
    mean(delta_s_if^2) - 2 * delta_s * delta^(-3) * mean(delta_if * 
                                                           delta_s_if)
  asymptotic_res <- tibble(estimand = c("Delta", "Delta_S", 
                                        "R"), estimate = c(delta, delta_s, 1 - delta_s/delta), 
                           se = c(sd(delta_if)/sqrt(n), sd(delta_s_if)/sqrt(n), 
                                  sqrt(sigma_sq)/sqrt(n)), ci_l = estimate - qnorm(1 - 
                                                                                     alpha/2) * se, ci_h = estimate + qnorm(1 - alpha/2) * 
                             se)
    gmat <- rBeta2009::rdirichlet(n_boot, rep(1, n)) * n
    boot_res_l <- list()
    for (b in 1:n_boot) {
      wt_b <- gmat[b, ]
      delta_b <- sum(delta_if * wt_b)/sum(wt_b) + delta - 
        mean(delta_if)
      delta_s_b <- sum(delta_s_if * wt_b)/sum(wt_b) + delta_s - 
        mean(delta_s_if)
      boot_res_l[[b]] <- tibble(estimand = c("Delta", "Delta_S", 
                                             "R"), boot_estimate = c(delta_b, delta_s_b, 1 - 
                                                                       delta_s_b/delta_b))
    }
    boot_res <- bind_rows(boot_res_l) %>% group_by(estimand) %>% 
      summarise(estimate = unique(case_when(estimand == 
                                              "Delta" ~ delta, estimand == "Delta_S" ~ delta_s, 
                                            estimand == "R" ~ 1 - delta_s/delta)), se_inf_fn = unique(asymptotic_res$se[asymptotic_res$estimand == 
                                                                                                                          unique(estimand)]), se = median(abs(boot_estimate - 
                                                                                                                            median(boot_estimate))), ci_l = quantile(boot_estimate, 
                                                                                                                                                                     alpha/2), ci_h = quantile(boot_estimate, 1 - 
                                                                                                                                                                                                 alpha/2)) 

return(list("asymptotic_res" = asymptotic_res, "boot_res" = boot_res))
}

#########################
## Set up the simulation parameters
###################################
yvars <- paste0('Y_', 0:tt)

##############
## Write a big function to run one simulation
#############################################
sim_fn <- function(this_run, set_seed = TRUE) {
  if (set_seed) set.seed(this_run)
  print(paste("Starting iteration = ", this_run))
  sim_ds <- sim(D_0, n = n) |>
    mutate(ff = sample(1:K, replace = TRUE, size = n)) |>
    mutate_at(vars(all_of(yvars)), ~ 1- .) |>
    as_tibble()
  
  #add censoring
  	cens = rexp(n, 0.1)
  	
  	#setting up the Y's and A's. 
  	sim_ds = mutate(sim_ds, "Y_0" = 1- (1*(is.na(Y_0) | Y_0 == 0)), "Y_1" = 1- (1*(is.na(Y_1) | Y_1 == 0)), "Y_2" = 1- (1*(is.na(Y_2) | Y_2 == 0)), "Y_3" = 1- (1*(is.na(Y_3) | Y_3 == 0)), "Y_4" = 1- (1*(is.na(Y_4) | Y_4 == 0)), "Y_5" = 1- (1*(is.na(Y_5) | Y_5 == 0)))
  	sim_ds = mutate(sim_ds, "A_0" = 1*(cens>0), "A_1" = 1*((cens>1 & Y_1 == 1) | (cens > 1 & Y_1 ==0 & Y_0 == 1)), "A_2" = 1*((cens>2 & Y_2 == 1) | (cens > 2 & Y_2 ==0 & Y_1 == 1)), "A_3" = 1*((cens>3 & Y_3 == 1) | (cens > 3 & Y_3 ==0 & Y_2 == 1)),"A_4" = 1*((cens>4 & Y_4 == 1) | (cens > 4 & Y_4 ==0 & Y_3 == 1)),"A_5" = 1*((cens>5 & Y_5 == 1) | (cens > 5 & Y_5 ==0 & Y_4 == 1)))
  	sim_ds = mutate(sim_ds, "Y_0" = Y_0*(cens>0), "Y_1" = Y_1*(cens>1), "Y_2" = Y_2*(cens>2), "Y_3" = Y_3*(cens>3), "Y_4" = Y_4*(cens>4), "Y_5" = Y_5*(cens>5))

  
  
  	tryCatch({
    p_deltahat <- plugin_delta(
    data = sim_ds,
    folds = 'ff',
    id = 'ID',
    x = 'X_0',
    g = 'G_0',
    a = paste0('A_', 0:tt),
    y = yvars,
    s = paste0('S_', 0:t0),
    binary_lrnr = lrn(lrnb, predict_type = 'prob'),
    cont_lrnr = lrn(lrnc),
    truncate_e = 0.005,
    verbose = FALSE
  )

  #
  p_deltahat_s <- plugin_delta_s(
    data = sim_ds,
    folds = 'ff',
    id = 'ID',
    x = 'X_0',
    g = 'G_0',
    a = paste0('A_', 0:tt),
    y = yvars,
    s = paste0('S_', 0:t0),
    binary_lrnr = lrn(lrnb, predict_type = 'prob'),
    cont_lrnr = lrn(lrnc),
    truncate_pi = 0.005,
    truncate_e = 0.005,
    verbose = FALSE
  )

tml_deltahat <- tmle_delta(data = sim_ds,
                             folds = 'ff',
                             id = 'ID',
                             x = 'X_0',
                             g = 'G_0',
                             a = paste0('A_', 0:tt),
                             y = yvars,
                             s = paste0('S_', 0:t0),
                             binary_lrnr = lrn(lrnb, predict_type = 'prob'),
                             cont_lrnr = lrn(lrnc),
                             truncate_e = 0.005,
                             verbose = FALSE)
  tml_deltahat_s <- tmle_delta_s(data = sim_ds,
                                 folds = 'ff',
                                 id = 'ID',
                                 x = 'X_0',
                                 g = 'G_0',
                                 a = paste0('A_', 0:tt),
                                 y = yvars,
                                 s = paste0('S_', 0:t0),
                                 binary_lrnr = lrn(lrnb, predict_type = 'prob'),
                                 cont_lrnr = lrn(lrnc),
                                 truncate_e = 0.005,
                                 verbose = FALSE)

  est.r.p = estimate_R(p_deltahat$if_data[[1]]$eif,
                       p_deltahat_s$if_data[[1]]$eif)
  est.r.t = estimate_R(p_deltahat$if_data[[1]]$eif,
                       p_deltahat_s$if_data[[1]]$eif, delta = tml_deltahat$tmle_est, delta_s = tml_deltahat_s$tmle_est)
  r_res.a <- bind_rows(est.r.p$asymptotic_res %>% mutate(method = 'plugin'), est.r.t$asymptotic_res %>% mutate(method = 'tmle')
  ) %>%
    mutate(true_val = case_when(estimand == 'Delta' ~ -Delta,
                                estimand == 'Delta_S' ~ -Delta_S,
                                estimand == 'R' ~ R_0),
           ci_covers = ci_h > true_val & ci_l < true_val)
  r_res.b <- bind_rows(est.r.p$boot_res %>% mutate(method = 'plugin'), est.r.t$boot_res %>% mutate(method = 'tmle')
  ) %>%
    mutate(true_val = case_when(estimand == 'Delta' ~ -Delta,
                                estimand == 'Delta_S' ~ -Delta_S,
                                estimand == 'R' ~ R_0),
           ci_covers = ci_h > true_val & ci_l < true_val)
  

   names(r_res.b) = c("estimand",  "estimate" , "se_inf_fn_boot", "se_boot","ci_l_boot", "ci_h_boot","method" , "true_val" ,"ci_covers_boot")
   r_res = inner_join(r_res.b,r_res.a, by = c("estimand","estimate","method","true_val"))
   r_res
   }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

sim_res <- map_dfr(1:n_sim, sim_fn)


sim_res_table =  sim_res %>%
  group_by(method, estimand) %>%
  summarise(estimate_m = round(mean(estimate),3),
  			bias = round(mean(estimate - true_val),3),
            empirical_se = round(sd(estimate),3),
            avg_est_se_a = round(mean(se),3),
            ci_coverage_a = round(mean(ci_covers),3),
            adj_ci_coverage_a = round(mean(ci_h > true_val + bias & ci_l< true_val + bias),3),
            rmse_a = round(sqrt(mean((estimate - true_val)^2)),3), 
  			  avg_estinf_se_boot = round(mean(se_inf_fn_boot),3),
  			avg_est_se_boot = round(mean(se_boot),3),
  			ci_coverage_boot = round(mean(ci_covers_boot),3))


sim_res_table

sim_res_table.small =  sim_res %>%
  group_by(method, estimand) %>%
  summarise(bias = round(mean(estimate - true_val),3),
            ci_coverage = round(mean(ci_covers),3))
sim_res_table.small.side= cbind(sim_res_table.small[1:3,-1],sim_res_table.small[4:6,-c(1:2)])
            
sim_res_table.small.side
 
sim_res_table.small.side[,1] = c("$\\Delta(t)$", "$\\Deltastt$","$R_{\\bS}(t,t_0)$")


latex.table(as.matrix(sim_res_table.small.side), paste("results_table_",n,"_",setting,sep=""), row.names = T, col.names = F, caption = "", dcolumn = T)

