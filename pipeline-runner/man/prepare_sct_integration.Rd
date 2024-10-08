% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qc-6-integrate_scdata-seuratv4.R
\name{prepare_sct_integration}
\alias{prepare_sct_integration}
\title{Prepare for integration after SCTransform}
\usage{
prepare_sct_integration(data.split, reduction, normalization, k.filter, npcs)
}
\arguments{
\item{data.split}{list of Seurat objects}

\item{reduction}{reduction method}

\item{normalization}{normalization method}

\item{k.filter}{number of neighbors (k) to use when filtering anchors}

\item{npcs}{number of PCs}
}
\value{
data.anchors to use for integration
}
\description{
This function runs the steps required to prepare the list of Seurat object normalized with
SCTransform for integration, and finds the integration anchors.
For further details see the documentation for
\code{\link[Seurat:SelectIntegrationFeatures]{Seurat::SelectIntegrationFeatures()}},
\code{\link[Seurat:PrepSCTIntegration]{Seurat::PrepSCTIntegration()}},
and \href{https://satijalab.org/seurat/articles/sctransform_v2_vignette.html#perform-integration-using-pearson-residuals-1}{sctransform_v2 vignette}.
}
