#' Generates 2IFC stimuli 
#' 
#' Generate stimuli for 2 images forced choice reverse correlation task. 
#' 
#' Will save the stimuli as
#' jpeg's to a folder, including .Rdata file needed for analysis of data after data collection. This
#' .Rdata file contains the parameters that were used to generate each stimulus.
#' 
#' @export
#' @param base_face_files List containing base face file names (jpegs) used as base images for stimuli
#' @param n_trials Number specifying how many trials the task will have (function will generate two images for each trial per base image: original and inverted/negative noise)
#' @param img_size Number specifying the number of pixels that the stimulus image will span horizontally and vertically (will be square, so only one integer needed)
#' @param stimulus_path Path to save stimuli and .Rdata file to
#' @param label Label to prepend to each file for your convenience
#' @param use_same_parameters Boolean specifying whether for each base image, the same set of parameters is used, or unique set is created for each base image
#' @param seed Integer seeding the random number generator (for reproducibility)
#' @param maximize_baseimage_contrast Boolean specifying wheter the pixel values of the base image should be rescaled to maximize its contrast. 
#' @return Nothing, everything is saved to files. 
generateStimuli2IFC <- function(base_face_files, n_trials=770, img_size=512, stimulus_path='./stimuli', label='rcic', use_same_parameters=TRUE, seed=1, maximize_baseimage_contrast=TRUE) {
  
  # Initalize #
  s <- generateNoisePattern(img_size)
  dir.create(stimulus_path, recursive=T)
  set.seed(seed)
  
  stimuli_params <- list()
  base_faces <- list()
  
  for (base_face in names(base_face_files)) {
    # Read base face
    img <- jpeg::readJPEG(base_face_files[[base_face]])    
    
    # Change base face to grey scale if necessary
    if (length(dim(img)) == 3) {
      img <- apply(img, c(1, 2), mean)
    }
    
    # Adjust size of base face
    #base_faces[[base_face]] <- biOps::imgMedianShrink(img, x_scale=img_size/ncol(img), y_scale=img_size/nrow(img))
    
    # If necessary, rescale to maximize contrast
    if (maximize_baseimage_contrast) {
      img <- (img - min(img)) / (max(img) - min(img))
    }
    
    # Save base image to list
    base_faces[[base_face]] <- img
  }
  
  # Generate parameters #
  if (use_same_parameters) {
    
    # Generate stimuli parameters, one set for all base faces
    params <- matlab::zeros(n_trials, 4096)
    for (trial in 1:n_trials) {  
      params[trial,] <- (runif(4096) * 2) - 1
    }    
    
    # Assign to each base face the same set
    for (base_face in names(base_faces)) {
      stimuli_params[[base_face]] <- params
    }
    
    rm(params)
  } else {
    for (base_face in names(base_faces)) {
      # Generate stimuli parameters, unique to each base face
      stimuli_params[[base_face]] <- matlab::zeros(n_trials, 4096)  
      for (trial in 1:n_trials) { 
        stimuli_params[[base_face]][trial,] <- (runif(4096) * 2) - 1
      }    
    }
    
  }
  
  # Generate stimuli #
  
  pb <- tcltk::tkProgressBar(title="Generating stimuli for all base faces", label=paste0("trials:", n_trials, " base faces:", length(base_faces)), min=0, max=n_trials, initial=0)
  stimuli <- matlab::zeros(img_size, img_size, n_trials)
  
  for (trial in 1:n_trials) {
    tcltk::setTkProgressBar(pb, trial)
    
    if (use_same_parameters) {
      # compute noise pattern, can be used for all base faces
      stimuli[,,trial] <- generateNoiseImage(stimuli_params[[base_face]][trial,], s) 
    }
    
    for (base_face in names(base_faces)) {
      if (!use_same_parameters) {
        # compute noise pattern unique to this base face
        stimuli[,,trial] <- generateNoiseImage(stimuli_params[[base_face]][trial,], s)        
      }
      
      # Scale noise (based on simulations, most values fall within this range [-0.3, 0.3], test
      # for yourself with simulateNoiseIntensities())
      stimulus <- ((stimuli[,,trial] + 0.3) / 0.6)
            
      # add base face
      combined <- (stimulus + base_faces[[base_face]]) / 2
      
      # write to file
      jpeg::writeJPEG(combined, paste(stimulus_path, paste(label, base_face, seed, sprintf("%05d_ori.jpg", trial), sep="_"), sep='/'), quality = 1.0)
      
      # compute inverted stimulus
      stimulus <- ((-stimuli[,,trial] + 0.3) / 0.6)
      
      # add base face
      stimulus <- (stimulus + base_faces[[base_face]]) / 2
      
      # write to file
      jpeg::writeJPEG(stimulus, paste(stimulus_path, paste(label, base_face, seed, sprintf("%05d_inv.jpg", trial), sep="_"), sep='/'), quality = 1.0)
    }
  }
  
  close(pb)  
  
  # Save all to image file (IMPORTANT, this file is necessary to analyze your data later and create classification images)
  save(base_face_files, base_faces, img_size, label, n_trials, s, seed, stimuli_params, stimulus_path, trial, use_same_parameters, file=paste(stimulus_path, paste(label, "seed", seed, "time", format(Sys.time(), format="%b_%d_%Y_%H_%M.Rdata"), sep="_"), sep='/'), envir=environment())
  
  
}

#' Generates 2IFC classification image 
#' 
#' Generate classification image for 2 images forced choice reverse correlation task. 
#' 
#' This funcions saves the classification image as jpeg to a folder and returns the CI. Your choice of scaling
#' matters. The default is \code{'matched'}, and will match the range of the intensity of the pixels to
#' the range of the base image pixels. This scaling is non linear and depends on the range of both base image
#' and noise pattern. It is truly suboptimal, because it shifts the 0 point of the noise (that is, pixels that would
#' have not changed base image at all before scaling may change the base image after scaling and vice versa). It is
#' however the quick and dirty way to see how the CI noise affects the base image.
#' 
#' For more control, use \code{'constant'} scaling, where the scaling is independent of 
#' the base image and noise range, but where the choice of constant is arbitrary (provided by the user with t
#' the \code{constant} parameter). The noise is then scale as follows: \code{scaled <- (ci + constant) / (2*constant)}.
#' Note that pixels can take intensity values between 0 and 1 If your scaled noise exceeds those values,
#' a warning will be given. You should pick a higher constant (but do so consistently for different classification images
#' that you want to compare). The higher the constant, the less visible the noise will be in the resulting image.
#' 
#' When creating multiple classification images a good strategy is to find the lowest constant that works for all 
#' classification images. This can be automatized using the \code{autoscale} function.
#' 
#' @export
#' @param stimuli Vector with stimulus numbers (should be numeric) that were presented in the order of the response vector. Stimulus numbers must match those in file name of the generated stimuli
#' @param responses Vector specifying the responses in the same order of the stimuli vector, coded 1 for original stimulus selected and -1 for inverted stimulus selected.
#' @param baseimage String specifying which base image was used. Not the file name, but the key used in the list of base images at time of generating the stimuli.
#' @param rdata String pointing to .RData file that was created when stimuli were generated. This file contains the contrast parameters of all generated stimuli.
#' @param saveasjpeg Boolean stating whether to additionally save the CI as jpeg image
#' @param filename Optional string to specify a file name for the jpeg image
#' @param antiCI Optional boolean specifying whether antiCI instead of CI should be computed
#' @param scaling Optional string specifying scaling method: \code{none}, \code{constant}, or \code{matched} (default)
#' @param constant Optional number specifying the value used as constant scaling factor for the noise (only works for \code{scaling='constant'})
#' @return List of pixel matrix of classification noise only, scaled classification noise only, base image only and combined 
generateCI2IFC <- function(stimuli, responses, baseimage, rdata, saveasjpeg=TRUE, filename='', antiCI=FALSE, scaling='constant', constant=0.1) {
  
  # Load parameter file (created when generating stimuli)
  load(rdata)
  
  # Check whether critical variables have been loaded
  if (!exists('s', envir=environment(), inherits=FALSE)) {
    stop('File specified in rdata argument did not contain s variable.')
  }
  
  if (!exists('base_faces', envir=environment(), inherits=FALSE)) {
    stop('File specified in rdata argument did not contain base_faces variable.')
  }
  
  if (!exists('stimuli_params', envir=environment(), inherits=FALSE)) {
    stop('File specified in rdata argument did not contain stimuli_params variable.')
  }
  
  # Get base image
  base <- base_faces[[baseimage]]
  
  # Get parameters of actually presented stimuli (this will work with non-consecutive stims as well)
  params <- stimuli_params[[baseimage]][stimuli,]
  
  # Compute classification image
  if (antiCI) {
    params = -params
  } 
  ci <- generateCI(params, responses, s)
  
  # Scale 
  if (scaling == 'none') {
    scaled <- ci
  } else if (scaling == 'constant') {
    scaled <- (ci + constant) / (2*constant)
    if (max(scaled) > 1.0 | min(scaled) < 0) {
      warning('Chosen constant value for constant scaling made noise of classification image exceed possible intensity range of pixels (<0 or >1). Choose a lower value, or clipping will occur.')
    } 
  } else if (scaling == 'matched') {
    scaled <- min(base) + ((max(base) - min(base)) * (ci - min(ci)) / (max(ci) - min(ci)))
  } else {
    warning(paste0('Scaling method \'', scaling, '\' not found. Using none.'))
    scaled <- ci
  }
  
  # Combine with base image
  combined <- (scaled + base) / 2
  
  # Save to file
  if (saveasjpeg) {
    
    if (filename == '') {
      filename <- paste0(baseimage, '.jpg')
    }
    
    if (antiCI) {
      filename <- paste0('antici_', filename)
    } else {
      filename <- paste0('ci_', filename)
    }
    
    jpeg::writeJPEG(combined, filename)
    
  }
  
  # Return list
  return(list(ci=ci, scaled=scaled, base=base, combined=combined))
}



#' Generates multiple 2IFC classification images by participant or condition 
#' 
#' Generate classification image for 2 images forced choice reverse correlation task. 
#' 
#' This funcions saves the classification images by participant or condition as jpeg to a folder and returns the CIs.
#' 
#' @export
#' @param data Data frame 
#' @param by String specifying column name that specifies the smallest unit (participant, condition) to subset the data on and calculate CIs for
#' @param stimuli String specifying column name in data frame that contains the stimulus numbers of the presented stimuli
#' @param responses String specifying column name in data frame that contains the responses coded 1 for original stimulus selected and -1 for inverted stimulus selected.
#' @param baseimage String specifying which base image was used. Not the file name, but the key used in the list of base images at time of generating the stimuli.
#' @param rdata String pointing to .RData file that was created when stimuli were generated. This file contains the contrast parameters of all generated stimuli.
#' @param saveasjpeg Boolean stating whether to additionally save the CI as jpeg image
#' @param antiCI Optional boolean specifying whether antiCI instead of CI should be computed
#' @param scaling Optional string specifying scaling method: \code{none}, \code{constant},  \code{matched} or \code{autoscale} (default)
#' @param constant Optional number specifying the value used as constant scaling factor for the noise (only works for \code{scaling='constant'})
#' @return List of classification image data structures (which are themselves lists of pixel matrix of classification noise only, scaled classification noise only, base image only and combined) 
batchGenerateCI2IFC <- function(data, by, stimuli, responses, baseimage, rdata, saveasjpeg=TRUE, antiCI=FALSE, scaling='autoscale', constant=0.1) {
 
  if (scaling == 'autoscale') {
    doAutoscale <- TRUE
    scaling <- 'none'
  } else {
    doAutoscale <- FALSE
  }
  
  pb <- tcltk::tkProgressBar(title="Computing classification images",  min=0, max=length(unique(data[,by])), initial=0)
  cis <- list()
  counter <- 0
  
  for (unit in unique(data[,by])) {
    
    # Update progress bar
    counter <- counter + 1
    tcltk::setTkProgressBar(pb, counter, label=unit)
    
    # Get subset of data 
    unitdata <- data[data[,by] == unit, ]

    # Specify filename for CI jpeg
    filename <- paste0(baseimage, '_', by, '_', unitdata[1,by])
    
    # Compute CI with appropriate settings for this subset (Optimize later so rdata file is loaded only once)
    cis[[filename]] <- generateCI2IFC(unitdata[,stimuli], unitdata[,responses], baseimage, rdata, saveasjpeg, paste0(filename, '.jpg'), antiCI, scaling, constant)
  }
  
  if (doAutoscale) {
    tcltk::setTkProgressBar(pb, counter, label="Autoscaling...")
    cis <- autoscale(cis, saveasjpegs=saveasjpeg)
  }
  
  close(pb)
  return(cis)

}

# Suppress checking notes for variables loaded at runtime from .RData files
if(getRversion() >= "2.15.1")  utils::globalVariables(c("s", "base_faces", "stimuli_params"))

