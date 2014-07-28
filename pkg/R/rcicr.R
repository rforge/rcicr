# Script by Ron Dotsch, based on Matlab code by Oliver Langner and Python code by Ron Dotsch
# r.dotsch@psych.ru.nl

#' @import matlab
library(matlab)

#' @import aspace
library(aspace)

#' @import biOps
library(biOps)

#' @import tcltk
library(tcltk)

#' @import reshape
library(reshape)

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
  angle <- as_radians(angle)
  sinepatch = repmat(linspace(0, cycles, img_size), img_size, 1)
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
#' generateNoisePattern(512)
#' generateNoisePattern(256)
#' generateNoisePattern(128)
generateNoisePattern <- function(img_size=512) {
  # Settings of sinusoids
  scales <- c(1, 2, 4, 8, 16)
  orientations <- c(0, 30, 60, 90, 120, 150)
  phases <- c(0, pi/2)
  
  # Size of sinusoids per scale
  mg <- meshgrid(1:img_size, 1:img_size,1:length(scales))
  x <- mg$x
  y <- mg$y
  rm(mg)
  sinSize = x / y
  
  # Number of sinsoids needed
  nrSin = length(scales) * length(orientations) * length(phases)
  
  # Pre allocate memory
  sinusoids = zeros(c(img_size, img_size, nrSin))
  sinIdx = zeros(c(img_size, img_size, nrSin))
  
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
        sinusoids[,,co] <- repmat(s, scale, scale)
        
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
#' @return The noise pattern as image
#' @examples
#' params <- rnorm(4096) # generates 4096 normally distributed random values
#' s <- generateNoisePattern(img_size=512)
#' noise <- generateNoiseImage(params, s)
#' generateNoisePattern(256)
#' generateNoisePattern(128)
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
#' @return The classification image
generateCI <- function(stimuli, responses, s) {
  weighted <- responses * stimuli
  params <- colMeans(weighted)
  return(generateNoiseImage(params, s))
}

