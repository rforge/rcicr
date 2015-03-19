#' Simulate pixel intensity range for noise 
#' 
#' @export
#' @param nrep Number of replications
#' @param img_size Size of noise pattern in pixels (one value equal for width and height)
#' @return Matrix with range of noise intensities for each replication 
simulateNoiseIntensities <- function(nrep=1000, img_size=512) {
  
  results <- matlab::zeros(nrep, 2)
  s <- generateNoisePattern(img_size=512)
  
  pb <- tcltk::tkProgressBar(title="Running simulations", min=0, max=nrep, initial=0)
  for (i in 1:nrep) {
    tcltk::setTkProgressBar(pb, i)
    
    params <- (runif(4096) * 2) - 1
    
    noise <- generateNoiseImage(params, s) 
    results[i,] <- range(noise)
  }
  close(pb)
  boxplot(results)
  return(results)
}
