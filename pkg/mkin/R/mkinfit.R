mkinfit <- function(mkinmod, observed, 
  parms.ini = rep(0.1, length(mkinmod$parms)),
  fixed_parms = rep(FALSE, length(mkinmod$parms)),
  state.ini = c(100, rep(0, length(mkinmod$diffs) - 1)), 
  fixed_initials = c(FALSE, rep(TRUE, length(mkinmod$diffs) - 1)), 
  ...)
{
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
  if(is.null(names(state.ini))) names(state.ini) <- paste(names(mkinmod$diffs), "0", sep="_")

  # TODO: Collect parameters to be optimised

  # Define the model cost function
  cost <- function(P)
  {
    inistates <- c(P[[1]], rep(0, length(mkinmod$diffs) - 1))
    names(inistates) = names(mkinmod$diffs)

    odeparms <- P[2:length(P)]
    names(odeparms) <- mkinmod$parms
    # Solve the ODE
    out <- ode(
      y = inistates,
      times = unique(observed$time),
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
    
    return(modCost(out_transformed, observed, y = "value"))
  }
  modFit(cost, c(state.ini[1], parms.ini))
}
