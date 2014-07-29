\name{rcicr-package}
\alias{rcicr-package}
\alias{rcicr}
\docType{package}
\title{
Reverse correlation image classification
}
\description{
Toolbox with functions to generate stimuli and analyze data of reverse correlation image classification experiments.
}
\details{
\tabular{ll}{
Package: \tab rcicr\cr
Type: \tab Package\cr
Version: \tab 0.2.1\cr
Date: \tab 2014-07-26\cr
License: \tab GPL-2\cr
}
A short tutorial: Load the package with \code{library(rcicr)}. Then create stimuli with: \code{generateStimuli2IFC(base_face_files, n_trials = 770)}and analyze your data by computing classification images with: \code{generateCI2IFC(stimuli, responses, baseimage, rdata)}. Use the \code{autoscale} function to optimize the scaling of the noise relative to the base image in such a way that different ci's autoscaled simultaneously are scaled in identical ways, thus making them easy to compare. See the help pages for those functions for more information using \code{?generateStimuli2IFC} and \code{?generateCI2IFC}. A longer tutorial will follow soon.

Currently, the package is still in alpha stage. Much may still change. It only supports 2 Image Forced Choice tasks, although the underlying functions can be used for other versions of the reverse correlation task. It also only supports sinusoid noise. In the future, it will support Gaussian white noise, as well as additional variants of the task.

If you use this package for your experiments, please cite it as: 

Dotsch, R. (2014). rcicr Reverse Correlation Image Classification Toolbox for R. http://r-forge.r-project.org/projects/rcicr/

}
\author{
Ron Dotsch <rdotsch@gmail.com> (http://ron.dotsch.org/)
Maintainer: Ron Dotsch <rdotsch@gmail.com>
}
\references{
Dotsch, R., & Todorov, A. (2012). Reverse correlating social face perception. Social Psychological and Personality Science, 3 (5), 562-571.

Dotsch, R., Wigboldus, D. H. J., Langner, O., & Van Knippenberg, A. (2008). Ethnic out-group faces are biased in the prejudiced mind. Psychological Science, 19, 978-980.
}
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
~~ simple examples of the most important functions ~~
}
