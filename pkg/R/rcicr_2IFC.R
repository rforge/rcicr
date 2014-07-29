#' Generates 2IFC stimuli 
#' 
#' Generate stimuli for 2 images forced choice reverse correlation task. 
#' 
#' Will save the stimuli as
#' jpeg's to a folder, including .Rdata file needed for analysis of data after data collection. This
#' .Rdata file contains the parameters that were used to generate each stimulus.
#' 
#' @export
#' @param base_face_file List containing base face file names (jpegs) used as base images for stimuli
#' @param n_trials Number specifying how many trials the task will have (function will generate two images for each trial per base image: original and inverted/negative noise)
#' @param img_size Number specifying the number of pixels that the stimulus image will span horizontally and vertically (will be square, so only one integer needed)
#' @param stimulus_path Path to save stimuli and .Rdata file to
#' @param label Label to prepend to each file for your convenience
#' @param use_same_parameters Boolean specifying whether for each base image, the same set of parameters is used, or unique set is created for each base image
#' @param seed Integer seeding the random number generator (for reproducibility)
#' @param distribution String specifying whether \code{uniform} or \code{normal} distribution should be used for contrasts.
#' @return Nothing, everything is saved to files. 
generateStimuli2IFC <- function(base_face_files, n_trials=770, img_size=512, stimulus_path='./stimuli', label='rcic', use_same_parameters=TRUE, seed=1, distribution='uniform') {
  
  # Initalize #
  s <- generateNoisePattern(img_size)
  dir.create(stimulus_path, recursive=T)
  set.seed(seed)
  
  stimuli_params <- list()
  base_faces <- list()
  
  for (base_face in names(base_face_files)) {
    # Read base face
    img <- biOps::readJpeg(base_face_files[[base_face]])    
    
    # Change base face to grey scale if necessary
    if (biOps::imageType(img) != "grey") {
      img <- biOps::imgRGB2Grey(img)
    } 
    
    # Adjust size of base face
    base_faces[[base_face]] <- biOps::imgMedianShrink(img, x_scale=img_size/ncol(img), y_scale=img_size/nrow(img))
    
  }
  
  # Generate parameters #
  if (use_same_parameters) {
    
    # Generate stimuli parameters, one set for all base faces
    params <- zeros(n_trials, 4096)
    for (trial in 1:n_trials) {  
      if (distribution == 'normal') {
        params[trial,] <- rnorm(4096)        
      } 
      if (distribution == 'uniform') {
        params[trial,] <- (runif(4096) * 2) - 1
      }
    }    
    
    # Assign to each base face the same set
    for (base_face in names(base_faces)) {
      stimuli_params[[base_face]] <- params
    }
    
    rm(params)
  } else {
    for (base_face in names(base_faces)) {
      # Generate stimuli parameters, unique to each base face
      stimuli_params[[base_face]] <- zeros(n_trials, 4096)  
      for (trial in 1:n_trials) { 
        if (distribution == 'normal') {
          stimuli_params[[base_face]][trial,] <- rnorm(4096)        
        } 
        if (distribution == 'uniform') {
          stimuli_params[[base_face]][trial,] <- (runif(4096) * 2) - 1
        }
      }    
    }
    
  }
  
  # Generate stimuli #
  
  pb <- tcltk::tkProgressBar(title="Generating stimuli for all base faces", label=paste0("trials:", n_trials, " base faces:", length(base_faces)), min=0, max=n_trials, initial=0)
  stimuli <- zeros(img_size, img_size, n_trials)
  
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
      
      stimulus <- stimuli[,,trial]
      
      # apply distribution dependent scaling
      if (distribution == 'normal') {
        # centered around 0, bring to [0, 255]
        stimulus <- (stimulus + 0.5) * 255
      }
      
      if (distribution == 'uniform') {
        # values are based on simulations, most values will fall withinin this range: [-0.3, 0.3]
        # test for yourself with simulateNoiseIntensities() function
        stimulus <- ((stimulus + 0.3) / (0.6)) * 255
      }
      
      stimulus <- biOps::imagedata(stimulus)
      
      # add base face
      stimulus <- biOps::imgAverage(list(stimulus, base_faces[[base_face]]))
      
      # write to file
      biOps::writeJpeg(paste(stimulus_path, paste(label, base_face, seed, sprintf("%05d_ori.jpg", trial), sep="_"), sep='/'), stimulus)
      
      # compute inverted stimulus
      stimulus <- -stimuli[,,trial]
      
      # apply distribution dependent scaling
      if (distribution == 'normal') {
        # centered around 0, bring to [0, 255]
        stimulus <- (stimulus + 0.5) * 255
      }
      
      if (distribution == 'uniform') {
        # values are based on simulations, most values will fall withinin this range: [-0.3, 0.3]
        # test for yourself with simulateNoiseIntensities() function
        stimulus <- ((stimulus + 0.3) / (0.6)) * 255
      }
      
      stimulus <- biOps::imagedata(stimulus)      
      
      # add base face
      stimulus <- biOps::imgAverage(list(stimulus, base_faces[[base_face]]))
      
      # write to file
      biOps::writeJpeg(paste(stimulus_path, paste(label, base_face, seed, sprintf("%05d_inv.jpg", trial), sep="_"), sep='/'), stimulus)
    }
  }
  
  close(pb)  
  
  # Save all to image file (IMPORTANT, this file is necessary to analyze your data later and create classification images)
  save(base_face_files, base_faces, distribution, img_size, label, n_trials, s, seed, stimuli_params, stimulus_path, trial, use_same_parameters, file=paste(stimulus_path, paste(label, "seed", seed, "time", format(Sys.time(), format="%b_%d_%Y_%H_%M.Rdata"), sep="_"), sep='/'), envir=environment())
  
  
}