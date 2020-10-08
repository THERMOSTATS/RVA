#' This is data to be included in package
#' @name Sample_summary_statistics_table
#' @docType data
#' @format An example summary statistics table as dataframe, row names are gene ID
#' \describe{
#'   \item{logFC}{log2 fold change from comparison}
#'   \item{AveExpr}{Average expression for this gene}
#'   \item{P.Value}{p value}
#'   \item{adj.P.Val}{adjusted p value or FDR}
#'   ...
#' }
#'
"Sample_summary_statistics_table"

#' This is data to be included in package
#' @name Sample_summary_statistics_table1
#' @docType data
#' @format Second example summary statistics table as dataframe, row names are gene ID
#' \describe{
#'   \item{logFC}{log2 fold change from comparison}
#'   \item{AveExpr}{Average expression for this gene}
#'   \item{P.Value}{p value}
#'   \item{adj.P.Val}{adjusted p value or FDR}
#'   ...
#' }
#'
"Sample_summary_statistics_table1"

#' This is data to be included in package
#' @name Sample_summary_statistics_table2
#' @docType data
#' @format Third example summary statistics table as dataframe, row names are gene ID
#' \describe{
#'   \item{logFC}{log2 fold change from comparison}
#'   \item{AveExpr}{Average expression for this gene}
#'   \item{P.Value}{p value}
#'   \item{adj.P.Val}{adjusted p value or FDR}
#'   ...
#' }
#'
"Sample_summary_statistics_table2"

#' This is data to be included in package
#' @name Sample_disease_gene_set
#' @docType data
#' @format An example disease gene set from summary statistics table as dataframe, row names are gene ID
#' the summary statistics can be calculated from disease vs healthy, which is this example.
#' \describe{
#'   \item{logFC}{log2 fold change from comparison}
#'   \item{AveExpr}{Average expression for this gene}
#'   \item{P.Value}{p value}
#'   \item{adj.P.Val}{adjusted p value or FDR}
#'   ...
#' }
#'
"Sample_disease_gene_set"

#' This is data to be included in package
#' @name c2BroadSets
#' @docType data
#' @format GeneSetCollection
#' \describe{
#'   \item{Genesetcollection}{GeneSetCollection from BroadCollection}
#' }
#'
"c2BroadSets"

#' This is data to be included in package
#' @name count_table
#' @docType data
#' @format An example count table where row names are gene ID, each column is a sample
#' \describe{
#'   \item{counttable}{count table}
#'   ...
#' }
"count_table"

#' This is data to be included in package
#' @name sample_count_cpm
#' @docType data
#' @format An example cpm table where row names are gene ID, each column is a sample
#' \describe{
#'   \item{counttable}{count cpm table}
#'   ...
#' }
"sample_count_cpm"

#' This is data to be included in package
#' @name sample_annotation
#' @docType data
#' @format Sample annotation document
#' \describe{
#'   \item{sample_id}{sample name}
#'   \item{tissue}{tissue for comparison}
#'   \item{subject_id}{subject id}
#'   \item{day}{time points}
#'   ...
#' }
"sample_annotation"

#' This is data to be included in package
#' @name wpA2020
#' @docType data
#' @format Rwikipathway data downloaded version 2020
#' \describe{
#'   \item{name}{pathway name}
#'   \item{version}{version}
#'   \item{wpid}{pathway id}
#'   \item{org}{host name}
#'   ...
#' }
"wpA2020"
