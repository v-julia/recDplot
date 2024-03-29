% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rec_plots.R
\name{plot_PDCP_with_control}
\alias{plot_PDCP_with_control}
\title{Pairwise nucleotide Distance Correspondence Plot with control}
\usage{
plot_PDCP_with_control(dna_object, st1, e1, st2, e2)
}
\arguments{
\item{dna_object}{list of DNA sequences produced by \code{read.dna} function of \code{ape} package (\code{as.character = TRUE} mode)}

\item{st1}{start position of genome region 1}

\item{e1}{end position of genome region 1}

\item{st2}{start position of genome region 2}

\item{e2}{end position of genome region 2}
}
\value{
list with ggplot object with PDC and control plot for two regions, dataframe with distances in two regions of alignment and in concatenated odd and even positions of alignment
}
\description{
This function plots PDC plot and its control on the same figure. Each dot corresponds to a pair of nucleotide distances between
the same pair of genomes in two genomic regions - with alignment positions \code{st1-e1} and \code{st2-e2} (see axis).
Returns list with ggplot and dataframe with pairwise distances calculated for \code{st1-e1} and \code{st2-e2} regions
}
\examples{
\dontrun{
alignment = read.dna(path/to/file, format="fasta", as.character=TRUE)
alignment[alignment=='-'] <- NA
plot_PDCP_with_control(alignment, 1, 500, 600, 1000)}
}
