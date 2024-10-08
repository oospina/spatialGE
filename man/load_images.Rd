% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{load_images}
\alias{load_images}
\title{load_images: Place tissue images within STlist}
\usage{
load_images(x = NULL, images = NULL)
}
\arguments{
\item{x}{an STlist}

\item{images}{a string indicating a folder to load images from}
}
\value{
an STlist with images
}
\description{
Loads the images from tissues to the appropriate STlist slot.
}
\details{
This function looks for \code{.PNG} or \code{.JPG} files within a folder matching the
sample names in an existing STlist. Then, loads the images to the STlist which
can be used for plotting along with other spatialGE plots.
}
