# ---------------------------------------------------------- #
# -------------------- Global functions -------------------- #
# ---------------------------------------------------------- #

formatting_fun <- function(idx,or,ci.lb,ci.ub){
  if(!is.na(or[idx])){
    return(paste0(sprintf(or[idx], fmt="%.2f")," (",sprintf(ci.lb[idx], fmt="%.2f"),",",sprintf(ci.ub[idx], fmt="%.2f"),")"))
  } else {
    return("N/A")
  }
}

spline_predictions <- function(mod,dat_pred,varnames){
  # mod = clogit object
  # dat_pred = values at which to make predictions
  # varnames = names of the variables included in predictions
  
  # get predicted values
  vars <- mod$assign
  vars_num <- as.integer(unlist(vars[varnames]))
  beta <- mod$coefficients
  vars_beta <- beta[vars_num]
  lp <- dat_pred %*% vars_beta
  # Get standard errors
  v <- mod$var
  v_vars <- v[vars_num,vars_num]
  var_lp <- dat_pred %*% v_vars %*% t(dat_pred)
  se_lp <- sqrt(diag(var_lp))
  # Put into data frame for plotting
  dat <- data.frame(pred=as.numeric(exp(lp)))
  dat$ci.lb <- as.numeric(exp(lp - 1.96*se_lp))
  dat$ci.ub <- as.numeric(exp(lp + 1.96*se_lp))
  
  # return
  return(dat)
}
