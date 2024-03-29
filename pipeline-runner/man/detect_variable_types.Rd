% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qc-7-embed_and_cluster.R
\name{detect_variable_types}
\alias{detect_variable_types}
\title{Detect cell-level metadata variable types}
\usage{
detect_variable_types(cl_metadata)
}
\arguments{
\item{cl_metadata}{data.table}
}
\value{
list of cell-level metadata types
}
\description{
detect cell level metadata types
\itemize{
\item group (CLMPerSample, like metadataCategorical)
\item "cellset type" (CLM, like cellSets)
\item excludes continuous/high cardinality (to avoid infinite cellsets)
}
}
