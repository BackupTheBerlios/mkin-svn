mkinfit <- function(mkinmod, observed, 
  parms.ini = rep(0.1, length(mkinmod$parms)),
  state.ini = c(100, rep(0, length(mkinmod$diffs) - 1)), 
  fixed_parms = rep(FALSE, length(mkinmod$parms)),
  fixed_initials = c(FALSE, rep(TRUE, length(mkinmod$diffs) - 1)), 
  plot = FALSE, quiet = FALSE,
  err = NULL, weight = "none", scaleVar = FALSE,
  ...)
{
  # Get names of observed variables
  obs_vars = unique(as.character(observed$name))

  # Name the parameters if they are not named yet
  if(is.null(names(parms.ini))) names(parms.ini) <- mkinmod$parms
  # Create a function calculating the differentials specified by the model
  mkindiff <- function(t, state, parms) {
    diffs <- vector()
    for (box in names(mkinmod$diffs))
    {
      diffname <- paste("d", box, sep="_")      
      diffs[diffname] <- with(as.list(c(state, parms)),
        eval(parse(text=mkinmod$diffs[[box]])))
    }
    return(list(c(diffs)))
  } 

  # Name the inital parameter values if they are not named yet
  if(is.null(names(state.ini))) names(state.ini) <- names(mkinmod$diffs)

  # TODO: Collect parameters to be optimised
  parms.optim <- parms.ini[!fixed_parms]
  parms.fixed <- parms.ini[fixed_parms]

  state.ini.optim <- state.ini[!fixed_initials]
  state.ini.optim.boxnames <- names(state.ini.optim)
  names(state.ini.optim) <- paste(names(state.ini.optim), "0", sep="_")
  state.ini.fixed <- state.ini[fixed_initials]

  cost.old <- 1e100
  calls <- 0
  out_predicted <- NA
  # Define the model cost function
  cost <- function(P)
  {
    assign("calls", calls+1, inherits=TRUE)
    if(length(state.ini.optim) > 0) {
      odeini <- c(P[1:length(state.ini.optim)], state.ini.fixed)
      names(odeini) <- c(state.ini.optim.boxnames, names(state.ini.fixed))
    } else odeini <- state.ini.fixed

    odeparms <- c(P[(length(state.ini.optim) + 1):length(P)], parms.fixed)

    # Get more timepoints if plotting is desired
    if(plot) {
      outtimes = unique(sort(c(observed$time, 
        seq(min(observed$time), max(observed$time), length.out=100))))
    } else outtimes = unique(observed$time)

    # Solve the ode
    out <- ode(
      y = odeini,
      times = outtimes,
      func = mkindiff, 
      parms = odeparms)
     
  
    # Output transformation for models with ghost compartments like SFORB
    out_transformed <- data.frame(time = out[,"time"])
    for (var in names(mkinmod$map)) {
      if(length(mkinmod$map[[var]]) == 1) {
        out_transformed[var] <- out[, var]
      } else {
        out_transformed[var] <- rowSums(out[, mkinmod$map[[var]]])
      }
    }    
    assign("out_predicted", subset(out_transformed, time %in% observed$time), inherits=TRUE)

    mC <- modCost(out_transformed, observed, y = "value",
      err = err, weight = weight, scaleVar = scaleVar)

    # Report and/or plot if the model is improved
    if (mC$model < cost.old) {
      if(!quiet) cat("Model cost at call ", calls, ": ", mC$model, "\n")

      # Plot the data and the current model output if requested
      if(plot) {
        plot(0, type="n", 
          xlim = range(observed$time), ylim = range(observed$value, na.rm=TRUE),
          xlab = "Time", ylab = "Observed")
        col_obs <- pch_obs <- 1:length(obs_vars)
        names(col_obs) <- names(pch_obs) <- obs_vars
        for (obs_var in obs_vars) {
          points(subset(observed, name == obs_var, c(time, value)), 
            pch = pch_obs[obs_var], col = col_obs[obs_var])
        }
        matlines(out_transformed$time, out_transformed[-1])
        legend("topright", inset=c(0.05, 0.05), legend=obs_vars, 
          col=col_obs, pch=pch_obs, lty=1:length(pch_obs))
      }
    
      assign("cost.old", mC$model, inherits=TRUE)
    }
    return(mC)
  }
  fit <- modFit(cost, c(state.ini.optim, parms.optim), ...)

  fit$observed <- mkin_long_to_wide(observed)
  fit$predicted <- out_predicted

  # Calculate chi2 error levels
  means <- aggregate(value ~ time + name, data = observed, mean, na.rm=TRUE)
  predicted_long <- mkin_wide_to_long(out_predicted, time = "time")
  errdata <- merge(means, predicted_long, by = c("time", "name"), suffixes = c("_mean", "_pred"))
  errdata <- errdata[order(errdata$time, errdata$name), ]
  fit$errmin.overall <- mkinerrmin(errdata, length(parms.ini))
  fit$errmin <- vector()
  for (obs_var in obs_vars)
  {
    errdata.var <- subset(errdata, name == obs_var)
    n.parms.optim <- length(grep(paste("k", obs_var, sep="_"), names(parms.optim)))
    n.initials.optim <- length(grep(paste(obs_var, 0, sep="_"), names(state.ini.optim)))
    n.optim <- n.parms.optim + n.initials.optim
    fit$errmin[obs_var] <- mkinerrmin(errdata.var, n.optim)
  }

  # Calculate chi2
  class(fit) <- c("mkinfit", "modFit") 
  return(fit)
}
