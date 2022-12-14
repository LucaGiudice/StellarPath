% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funs_wrappers.R
\name{enrichment_analysis}
\alias{enrichment_analysis}
\title{Enrichment analysis of the results of classification}
\usage{
enrichment_analysis(out_l, extended = FALSE)
}
\arguments{
\item{out_l}{The output list obtained from the classification of the user-defined data with the function: prepare_data, analyse_training and classification}

\item{extended}{Default FALSE, the overlap between pathways is computed with only the significant features falling into that pathways}
}
\value{
For each unique omic (e.g. gene expression), it provides a dataframe (SDR) which is a result of the function
find_SDR, for example $TRAINING$gex$SDR. For each omic type (e.g. transcriptomics which can be combination of genes, lncRNAs and miRNAs), it adds
a list of two elements in the $TRAINING$pathway_analysis$transcriptomics section of the data list: top_pathways_l and top_pathways_df
top_pathways_df details the pathways/sets which have a significant PSN. The list is the result of the function find_PSN.
- EnrichmentMap similarity network between the pathways/sets of targets that have been found significant with their PSNs
}
\description{
It takes in input the resulting list of the classification function. It finds the run of cross validation in which the model classified
the samples with the best performance. It finds the similarity networks that separate the two sample's groups.
It adds annotations to the significant gene pathways/non-coding target sets/PSNs:
- ID and database source for canonical gene pathways
- Regulation pathway targeted by the non-coding significant elements
}
\details{
It creates a similarity network between the pathways/sets of targets that have been found significant with their PSNs.
For each pair of pathway/target set, there is a weighted edge which measures how much is strong their overlap and so their similarity.
This network is save in an object called EnrichmentMap because it is analog to the idea presented with Cytoscape.
}
