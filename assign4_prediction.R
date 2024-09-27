### ONLY past your predict.y function here

predict.y <- function(x)
{
  load("fit_params.Rdata")
  x <- data.matrix(x)
  f.x <- t(beta_min_p_pos)%*%(x)
  return(f.x)
}
