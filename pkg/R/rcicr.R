# Script by Ron Dotsch, based on Matlab code by Oliver Langner and Python code by Ron Dotsch
# r.dotsch@psych.ru.nl

#' Generate single sinusoid patch
#'
#' @export
#' @param img_size Integer specifying size of sinusoid patch in number of pixels
#' @param cycles Integer specifying number of cycles sinusoid should span
#' @param angle Value specifying the angle (rotation) of the sinusoid
#' @param phase Value specifying phase of sinusoid
#' @param contrast Value between -1.0 and 1.0 specifying contrast of sinusoid
#' @return The sinusoid image with size \code{img_size}.
#' @examples
#' generateSinusoid(512, 2, 90, pi/2, 1.0)
generateSinusoid <- function(img_size, cycles, angle, phase, contrast) {
  
  # Generates an image matrix containing a sinusoid, angle (in degrees) of 0 will give vertical, 90 horizontally oriented sinusoid  
  angle <- aspace::as_radians(angle)
  sinepatch = matlab::repmat(matlab::linspace(0, cycles, img_size), img_size, 1)
  sinusoid <- (sinepatch * cos(angle) + t(sinepatch) * sin(angle)) * 2 * pi
  sinusoid <- contrast * sin(sinusoid + phase)
  return(sinusoid)
}


#' Generate sinusoid noise pattern
#' 
#' @export
#' @param img_size Integer specifying size of the noise pattern in number of pixels
#' @return List with two elements: the 3D sinusoid matrix with size \code{img_size}, and and indexing
#' matrix with the same size to easily change contrasts.
#' @examples
#' generateNoisePattern(256)
generateNoisePattern <- function(img_size=512) {
  # Settings of sinusoids
  scales <- c(1, 2, 4, 8, 16)
  orientations <- c(0, 30, 60, 90, 120, 150)
  phases <- c(0, pi/2)
  
  # Size of sinusoids per scale
  mg <- matlab::meshgrid(1:img_size, 1:img_size,1:length(scales))
  x <- mg$x
  y <- mg$y
  rm(mg)
  sinSize = x / y
  
  # Number of sinsoids needed
  nrSin = length(scales) * length(orientations) * length(phases)
  
  # Pre allocate memory
  sinusoids = matlab::zeros(c(img_size, img_size, nrSin))
  sinIdx = matlab::zeros(c(img_size, img_size, nrSin))
  
  # counters
  co = 0 # sinusoid layer counter
  idx = 0 # contrast index counter
  
  for (scale in scales) {
    for (orientation in orientations) {
      for (phase in phases) {
        # Generate single sinusoid
        size <- sinSize[scale, img_size]
        s <- generateSinusoid(size, 2, orientation, phase, 1)
        
        # Repeat to fill scale
        sinusoids[,,co] <- matlab::repmat(s, scale, scale)
        
        # Create index matrix
        for (col in 1:scale) {
          for (row in 1:scale) {
            
            # insert absolute index for later contrast weighting
            sinIdx[(size * (row-1) + 1) : (size * row), (size * (col-1) + 1) : (size * col), co] = idx
            
            # Update contrast counter
            idx = idx + 1
            
          }
        } 
        
        # Update layer counter
        co = co + 1
        
      }
    }
  }
  return(list(sinusoids=sinusoids, sinIdx=sinIdx))
}


#' Generate single noise image based on parameter vector
#' 
#' @export
#' @param params Vector with 4096 values specifying the contrast of each sinusoid in noise
#' @param s 3D sinusoid matrix (generated using \code{generateNoisePattern()})
#' @return The noise pattern as pixel matrix
#' @examples
#' params <- rnorm(4096) # generates 4096 normally distributed random values
#' s <- generateNoisePattern(img_size=256)
#' noise <- generateNoiseImage(params, s)
generateNoiseImage <- function(params, s) {
  noise <- apply(s$sinusoids * array(params[s$sinIdx], dim(s$sinusoids)), 1:2, mean)
  return(noise)
}

#' Generate classification image based on set of stimuli (matrix: trials, parameters), responses (vector), and sinusoid
#' 
#' @export
#' @param stimuli Matrix with one row per trial, each row containing the 4096 parameters for the original stimulus
#' @param responses Vector containing the response to each trial (1 if participant selected original , -1 if participant selected inverted;
#' this can be changed into a scale)
#' @param s 3D sinusoid matrix (generated using \code{generateNoisePattern()})
#' @return The classification image as pixel matrix
generateCI <- function(stimuli, responses, s) {
  weighted <- responses * stimuli
  params <- colMeans(weighted)
  return(generateNoiseImage(params, s))
}

#' Determines optimal scaling constant for a list of ci's
#' 
#' @export
#' @param cis List of cis, each of which are a list containing the pixel matrices of at least the noise pattern (\code{$ci}) and if the noise patterns need to be written to jpegs, als the base image (\code{$base})
#' @param saveasjpegs Boolean, when set to true, the autoscaled noise patterns will be combined with their respective base images and saved as jpegs (using the key of the list as name)
#' @return List of scaled noise patterns and determind scaling factor
autoscale <- function(cis, saveasjpegs=TRUE) {
  # Get range of each ci
  ranges <- matlab::zeros(length(names(cis)), 2)
  for (ciname in names(cis)) {
    ranges[which(ciname==names(cis)), ] <- range(cis[[ciname]]$ci)
  }  

  # Determine the lowest possible scaling factor constant
  if (abs(min(ranges[,1])) > max(ranges[,2])) {
    constant <- abs(min(ranges[,1]))
  }  else {
    constant <- max(ranges[,2])
  }

  print(paste0("Using scaling factor constant:", constant))
  
  # Scale all noise patterns
  for (ciname in names(cis)) {
    cis[[ciname]]$scaled <-  (cis[[ciname]]$ci + constant) / (2*constant)
    
    # Combine and save to jpeg if necessary
    if (saveasjpegs) {
      ci <- (cis[[ciname]]$scaled + cis[[ciname]]$base) / 2
      jpeg::writeJPEG(ci, paste0(ciname, '_autoscaled.jpg'), quality=1.0)
    }
  
  }

  cis[['autoscaling.constant']] <- constant
  return(cis)
}