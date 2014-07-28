#' Simulate pixel intensity range for noise 
#' 
#' @export
#' @param nrep Number of replications
#' @param img_size Size of noise pattern in pixels (one value equal for width and height)
#' @param distribution Which distribution for parameters to use (normal/uniform)
#' @return Matrix with range of noise itensities for each replication 
simulateNoiseIntensities <- function(nrep=1000, img_size=512, distribution='uniform') {
  
  results <- zeros(nrep, 2)
  s <- generateNoisePattern(img_size=512)
  
  pb <- tkProgressBar(title="Running simulations", min=0, max=nrep, initial=0)
  for (i in 1:nrep) {
    setTkProgressBar(pb, i)
    
    if (distribution == 'normal') {
      params <- rnorm(4096)        
    } 
    
    if (distribution == 'uniform') {
      params <- (runif(4096) * 2) - 1
    }
    
    noise <- generateNoiseImage(params, s) 
    results[i,] <- range(noise)
  }
  close(pb)
  boxplot(results)
  return(results)
}
