% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gem2s-6-construct_qc_config.R
\name{construct_qc_config}
\alias{construct_qc_config}
\title{Constructs default QC configuration}
\usage{
construct_qc_config(scdata_list, unfiltered_samples, technology)
}
\arguments{
\item{scdata_list}{list of Seurat objects}

\item{unfiltered_samples}{character vector of unfiltered sample ids}
}
\value{
list of QC configuration parameters
}
\description{
This function returns the default parameters used during QC as a nested list.
It is sent to the API, which in turn saves it as a jsonb object in the
PostgreSQL database.
}
\details{
For each pipeline step, it customizes and adds the QC parameters for all samples
in the running experiment.
}
