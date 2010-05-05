mkingen <- function(spec = list(parent = list(type = "SFO", input = NA, output = NA)))
{
  parms <- list() 
  for (varname in names(mkinlist))
  {
    # New (sub)compartments needed for the model type
    new <- switch(mkinlist[[varname]]$type,
      SFO = varname,
      SFORB = paste(varname, c("free", "bound"), sep="_")
    )
  }
  model <- list(diff = diff, parms = parms)
  class(model) <- mkingen 
}
