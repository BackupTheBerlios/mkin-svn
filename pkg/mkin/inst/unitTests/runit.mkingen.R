test.mkingen <- function(
{  
  SFO.diff <- c(
    parent = "d_parent <- - k_parent_sink * parent"
  )
  SFO.parms <- c("parent_0", "k_parent_sink")
  SFO <- list(diff = SFO.diff, parms = SFO.parms)
  class(SFO) <- "mkingen"
  checkIdentical(SFO, mkingen(spec = list(parent = list(type = "SFO", input = NA, output = NA))))
}
