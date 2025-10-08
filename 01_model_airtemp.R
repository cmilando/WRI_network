#' ============================================================================
#' ////////////////////////////////////////////////////////////////////////////
#' 
#' Bootstrap the linear model so the confidence intervals are more robust
#' 
#' ////////////////////////////////////////////////////////////////////////////
#' ============================================================================
library(boot)

library(future)
library(future.apply)
plan(multisession)

N_boot <- 500  # the number of bootstrap iterations
sz     <- 0.5   # the fraction of monitors to be used each bootstrap

#' ============================================================================
#' ////////////////////////////////////////////////////////////////////////////
# use all the data
lm1 <- glm(air_temp ~ daymet + green + albedo, data = df)


#' ============================================================================
#' ////////////////////////////////////////////////////////////////////////////

# do bootstrapped cross-validation to get more robust coefficients
# by only using a fraction of each monitors each time

boot.fn = function (i, sz) {
  
  # get indices
  monitor_id_subset = sample(1:N, size = round(sz * N),
                 replace = F)
  
  # get new glm for this subset of monitors
  lm_boot <- glm(air_temp ~ daymet + green + albedo, 
                 data = df %>% filter(monitor_id %in% monitor_id_subset))
  
  return (coef(lm_boot))
}

boot.l <- future_lapply(1:N_boot, FUN = boot.fn, sz = sz, 
                        future.seed=TRUE)

boot.out <- do.call(rbind, boot.l)

# get bootstrapped confidence intervals and compare to full model
boot.lb <- apply(boot.out, 2, \(x) quantile(x, 0.025))
boot.ub <- apply(boot.out, 2, \(x) quantile(x, 0.975))

out <- rbind(boot.lb,  boot.ub)
row.names(out) = c('eCI 2.5 %', 'eCI 97.5 %')
t(out)

confint(lm1)
