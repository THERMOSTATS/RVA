#' @title Validate Pathways DB
#' @description To ensure selected db name is correct.
#' @param pathway.db The databse to be used for encrichment analysis. Can be one of the following, "rWikiPathways", "KEGG", "REACTOME", "Hallmark","rWikiPathways_aug_2020"
#' @importFrom stringr str_c
#' @references Xingpeng Li, Tatiana Gelaf Romer & Siddhartha Pachhai RVA - RNAseq Visualization Automation tool.
validate.pathways.db <- function(pathway.db){
  valid.db <- c("rWikiPathways", "KEGG", "REACTOME", "Hallmark","rWikiPathways_aug_2020")
  if (!pathway.db %in% valid.db){
    stop(paste0("The DB name you have entered is cannot be used please select one of the follwing",
                str_c(valid.db, collapse = " , " )))
  }
  if (pathway.db == "rWikiPathways_aug_2020"){
    warning(paste0("\n\n Using rWikiPathways from August 2020, please use",
                   " pathway.db = rWikiPathways for latest version \n\n" ))
  }
  cat(paste0("\n\n Currently using the ", pathway.db ," database for enrichment \n\n"))

}


#' @title Validate pval flag
#' @description To ensure p value flags are the same accross datasets.
#' @param data A list of summary statistics table (data.frame) from limma or DEseq2, where rownames are gene id.
#' @param value P value flag.
#' @importFrom purrr map
#' @references Xingpeng Li, Tatiana Gelaf Romer & Siddhartha Pachhai RVA - RNAseq Visualization Automation tool.
validate.pvalflag <- function(data, value) {
  if(!inherits(data, 'list')) {
    options=colnames(data)
    if(!(value %in% options)) {
      stop(paste0('\n \n The supplied p-value flag was not found in dataset. If a signle data frame was passed make sure that p.value.flag is one your columns. If a list was passed make sure that p.value.flag is present in all data frames. \n'))
    }
  }else if(inherits(data, 'list')) {
    map(data, validate.pvalflag, value=value)
  }

}


#' @title Validate Single Table is not list
#' @description Makes sure the summary table being input is of the right class and format.
#' @param data summary statistics table (data.frame) from limma or DEseq2, where rownames are gene id.
#' @importFrom purrr map
#' @references Xingpeng Li, Tatiana Gelaf Romer & Siddhartha Pachhai RVA - RNAseq Visualization Automation tool.
validate.single.table.isnotlist <- function(data) {

  if (inherits(data,"list") & length(data) == 1){
    stop(paste0("There is only 1 summary table being passed in as a list, please convert to data frame"))
  }

  if (inherits(data,"list")){
    map(data,validate.single.table.isnotlist)
  }
}






#' @title Validate Comp Names
#' @description This function ensures that when a list of data frames are used as input the
#'        the number of comp names are the same as the number of data frames.
#' @param comp.names a character vector contain the comparison names corresponding to the same order to the \code{dat.list}. default = NULL.
#' @param data summary statistics table (data.frame) from limma or DEseq2, where rownames are gene id.
#' @references Xingpeng Li, Tatiana Gelaf Romer & Siddhartha Pachhai RVA - RNAseq Visualization Automation tool.
validate.comp.names <- function(comp.names, data) {
  if (!inherits(data,"list")){
    if (inherits(comp.names, "list") | inherits(comp.names, "vector")){
      stop(paste0("If input is a data frame please input comp.names as character or NULL"))
    }
  }
  else {
    if(length(data) != length(comp.names)){
      stop(paste0("Please make sure the provided summary statistics",
                " list dat.list has the same length as comp.names"))
    }

    if(length(unique(comp.names)) != length(comp.names)){
        stop(paste0("Number of unique comp names is not equal to the number of comp names provided"))

    }
  }
}




#' @title Validate Geneset
#' @description This function ensures that the input geneset to check.cutoff
#'        is formatted properly and in a usable form.
#' @param data summary statistics table or a list contain multiple summary statistics tables from limma or DEseq2, where each row is a gene.
#' @param geneset a summary statistic table contain the genes want to be highlighted, the gene name format (in row names)
#' needs to be consistent to the main summary statistics table). For example, this summary statistics
#'  table coulb be the output summary statistics table from Disease vs Healthy comparison (Only contain
#'  the subsetted significant genes want to be highlighted).
#' @param highlight.1 genes want to be highlighted, in the format of a vector consists of gene names. The gene name format
#'  needs to be consistent to the main summary statistics table.
#' @param highlight.2 genes want to be highlighted, in the format of a vector consists of gene names. The gene name format
#'  needs to be consistent to the main summary statistics table.

#'
#' @importFrom purrr set_names map2 map
#'
#' @return A character value indicating if the geneset was passed as a
#'         dataframe (`df`) or two vectors (`vec`), if a list is input
#'         the number of returned values equal the length of the list
#'
#'
#' @details The function ensures that only a dataframe or vectors are supplied,
#'         that at least one or the other is supplied, and that their formatting
#'         is correct if supplied. It also checks if any of the genes overlap
#'         with the genes in the datanames.
#' @references Xingpeng Li, Tatiana Gelaf Romer & Siddhartha Pachhai RVA - RNAseq Visualization Automation tool.
validate.geneset <- function(data,
                             geneset,
                             highlight.1,
                             highlight.2) {


        #added recursively check list
        if(inherits(data, "list")){

          cat ("Checking gene sets for listof data frames")

          map(data,validate.geneset, geneset=geneset, highlight.1 = highlight.1, highlight.2 = highlight.2)


        }else{


          if (class(highlight.1) == "list"){
            highlight.1 = unlist(highlight.1, use.names=FALSE)
          }

          if (class(highlight.2) == "list"){
            highlight.2 = unlist(highlight.2, use.names=FALSE)
          }

          if(is.null(geneset) & is.null(highlight.1) & is.null(highlight.2)) {
                  geneset.type <- "No.highlight"
                  return(geneset.type)
          }
          if(!is.null(geneset) & (!is.null(highlight.1) | !is.null(highlight.2))) {
                  stop(paste0("Please supply only ONE of: geneset parameter OR",
                              " highlight.1 and/or highlight.2 parameter."))
          }

          if(is.data.frame(geneset)) {
                  geneset.genes <- row.names(geneset)
                  if(is.null(geneset.genes)) {
                          stop(paste0("Input gene set must have rows labeled",
                                      " with gene names."))
                  }
                geneset.type <- "df"

          }

          else if(is.character(geneset)) {
                geneset.genes <- geneset

          }

          else if (is.null(geneset)) {

                if(is.null(highlight.1) | is.null(highlight.2)) {
                        stop(paste0("You must supply either a geneset parameter",
                                    " or a highlight.1 AND a highlight.2",
                                    " parameter."))
                }
                if(!is.character(highlight.1) | !is.character(highlight.2)) {
                        stop(paste0("Invalid data type supplied as gene sets to",
                                    " highlight.1 and highlight.2. Please",
                                    " supply two character vectors."))
                }

                geneset.genes <- c(highlight.1, highlight.2)
                geneset.type <- "vec"

          }
          else {
                stop(paste0("Invalid data type supplied as gene set. Please",
                            " supply a character vector or dataframe."))
          }

        validate.genes.present(rownames(data), geneset.genes)
        return(geneset.type)

        }
}


#' @title Validate genes present
#' @description Checks how many of the gene id's in the dataset are there in the geneset.
#' @param data.genes The gene id's.
#' @param geneset a summary statistic table contain the genes want to be highlighted, the gene name format (in row names)
#' needs to be consistent to the main summary statistics table). For example, this summary statistics
#'  table coulb be the output summary statistics table from Disease vs Healthy comparison (Only contain
#'  the subsetted significant genes want to be highlighted).
#' @references Xingpeng Li, Tatiana Gelaf Romer & Siddhartha Pachhai RVA - RNAseq Visualization Automation tool.
validate.genes.present <- function(data.genes, geneset) {

        missing.genes <- geneset[!(geneset %in% data.genes)]
        if(length(missing.genes) == length(geneset)) {
                stop(paste0("None of the genes in the gene set exist in dataset."))
        }
        if(length(missing.genes) > 0) {
                warning(paste0("\n\n",length(missing.genes)," genes from the gene set are missing",
                                " in the dataset and will be",
                                " ignored: \n\n"))
        }
}

#' @title Validate Foldchange
#' @description This function ensures the fold change minimum, maximum, and step
#'         are valid.
#' @inheritParams plot.cutoff
#'
#' @details Specifically it checks that the FCmax is greater than the FCmin,
#'         that at least 1 FCstep can fit within the FCmax and FCmin, that
#'         FCmax and FCmin values are non-negative, and that FCstep is positive.
#' @references Xingpeng Li, Tatiana Gelaf Romer & Siddhartha Pachhai RVA - RNAseq Visualization Automation tool.
validate.FC <- function(FCmin,
                        FCmax,
                        FCstep) {
        if(FCmax < FCmin) {
                stop("Choose an FCmax greater than FCmin.")
        }
        if(FCmax < 0 | FCmin < 0) {
                stop("FCmin and FCmax must be non-negative.")
        }
        if(FCstep <= 0) {
                stop("FCstep must be positive.")
        }
        if(abs(FCstep) >= abs(FCmax - FCmin)) {
                stop("FCstep does not fit within the FC bounds.")
        }
}

#' @title Validate Pvalues
#' @description This function ensures the fold change minimum, maximum, and step
#'         are valid.
#' @inheritParams plot.cutoff
#'
#' @details Specifically it checks that the pvalues are between 0-1, and that
#'         at least 1 `p.step` fits within the `p.min` and `p.max` bounds and
#'         is positive.
#' @references Xingpeng Li, Tatiana Gelaf Romer & Siddhartha Pachhai RVA - RNAseq Visualization Automation tool.
validate.pvals <- function(p.min,
                           p.max,
                           p.step) {
        validate.pval.range(p.min, "p.min")
        validate.pval.range(p.max, "p.max")

        if(p.step <= 0 | p.step > (p.max-p.min)) {
                stop("Invalid p.step value. p.step values must be positive and",
                     " less than p-value range to plot.")
        }
}

validate.signif <- function(signif.vals) {
        for(p.val in signif.vals) {
                validate.pval.range(p.val, "signif.val")
        }
}

#' @title Validate P-value Range
#'
#' @description Error-handling for invalid p-value.
#'
#' @param pval The pvalue
#' @param name The name of the value to include in the error.
#' @references Xingpeng Li, Tatiana Gelaf Romer & Siddhartha Pachhai RVA - RNAseq Visualization Automation tool.
validate.pval.range <- function(pval,
                                name) {
        if(pval < 0 | pval > 1) {
                stop(paste0("Invalid ", name, " value entered;",
                            " p-values must be between 0 and 1."))
        }
}

#' @title Validate Summary Statistics File
#'
#' @description Check for required column names and types.
#'
#' @inheritParams plot.cutoff
#' @param datin the summary statistics file.
#' @references Xingpeng Li, Tatiana Gelaf Romer & Siddhartha Pachhai RVA - RNAseq Visualization Automation tool.
validate.stats <- function(datin,
                           name = 1,
                           ...) {

        validate.stats.cols(datin, name = name, list(...))
        validate.col.types(datin, name = name, list(...))
}

#' @title Check Summary Statistics Required Columns
#' @description Required columns are `FCflag` and `FDRflag`
#'
#' @inheritParams validate.stats
#' @references Xingpeng Li, Tatiana Gelaf Romer & Siddhartha Pachhai RVA - RNAseq Visualization Automation tool.
validate.stats.cols <- function(datin,
                                name = 1,
                                req.cols) {

        existing.cols <- colnames(datin)

        missing.cols <- req.cols[!(req.cols %in% existing.cols)]

        if(length(missing.cols) > 0) {
                stop(paste0("The following required columns are missing from",
                            " summary statistics file ", name, ": ",
                            paste0(missing.cols, collapse = ",")))
        }
}

#' @title Check Summary Statistics Required Column Types
#' @description `FCflag` and `FDRflag` must be numeric.
#'
#' @inheritParams validate.stats
#' @references Xingpeng Li, Tatiana Gelaf Romer & Siddhartha Pachhai RVA - RNAseq Visualization Automation tool.
validate.col.types <- function(datin,
                               name = 1,
                               flags) {
        for(flag in flags) {
                validate.numeric(datin, flag, name = name)
        }
}


#' @title Validate Flag Value Is Valid
#' @description Enures that the `value` is one of `options` and throws an error
#'         otherwise.
#' @param value The user-input value for the parameter
#' @param options A vector of valid values for `value`
#' @param name The name of the parameter to be displayed in the error
#' @references Xingpeng Li, Tatiana Gelaf Romer & Siddhartha Pachhai RVA - RNAseq Visualization Automation tool.
validate.flag <- function(value, name, options) {
        if(!(value %in% options)) {
                stop(paste0("Invalid ", name, " value \"", value, "\". ",
                            name, " values must be one of [\"",
                            str_c(options, collapse = "\", \""), "\"]."))
        }

}

validate.geneid.flag <- function(value, name) {
        options <-c("ACCNUM", "ALIAS", "ENSEMBL", "ENSEMBLPROT, ENSEMBLTRANS",
                    "ENTREZID", "ENZYME", "EVIDENCE", "EVIDENCEALL", "GENENAME",
                    "GO", "GOALL", "IPI", "MAP", "OMIM, ONTOLOGY",
                    "ONTOLOGYALL", "PATH", "PFAM", "PMID", "PROSITE", "REFSEQ",
                    "SYMBOL", "UCSCKG", "UNIGENE", "UNIPROT")
        validate.flag(value, name, options)
}

#' @title Validate Numeric Column
#' @description Ensures that a column in a dataframe which must be numeric is
#'         numeric and throws an error otherwise.
#'
#' @param datin The data in question.
#' @param col The column to validate as numeric.
#'
#' @details This specifically checks if any of the values in the column can
#'         be coerced as numeric.
#' @references Xingpeng Li, Tatiana Gelaf Romer & Siddhartha Pachhai RVA - RNAseq Visualization Automation tool.
validate.numeric <- function(datin,
                             col,
                             name = 1) {
        vals <- unlist(datin[,col])
        non.na.vals <- vals[!is.na(vals)]

        vals.num <- as.numeric(vals)
        if(all(is.na(vals.num))) {
                stop(paste0("Non-numeric data detected in ", col, "column of",
                            " dataset ", name))
        }
}

#' @title Validate Annotation Table
#' @description Ensure that an annotation has all of the required columns.
#'
#' @param data The input count data.
#' @param annot The annotation dataframe.
#' @param annot.flags The vector of annotation flags passed by the user.
#' @param fill The fill value indicated by the user.
#' @param baseline.flag The baseline.flag passed by the user.
#' @param baseline.val The baseline value passed by the user.
#'
#' @details The function will check the following:
#'         \itemize{
#'                 \item The `annot.flags` values are columns in `annot`
#'                 \item If `fill` = "cfb": validate the `baseline.flag` and
#'                         `baseline.val` parameters.
#'                 \item `sample.id` is a column in `annot`.
#'         }
#' @references Xingpeng Li, Tatiana Gelaf Romer & Siddhartha Pachhai RVA - RNAseq Visualization Automation tool.
validate.annot <- function(data, annot, annot.flags, sample.id,
                           fill = "CPM",
                           baseline.flag = NULL,
                           baseline.val = NULL) {
        cols <- colnames(annot)
        if(!is.null(sample.id) && !(sample.id %in% cols)) {
                stop(paste0("Required column [", sample.id,
                            "] is missing from your annotation input data.")
                )
        }
        missing.flags <- annot.flags[!(annot.flags %in% cols)]
        if(length(missing.flags) > 0) {
                stop(paste0("Annotation data missing flag column(s) [",
                            str_c(missing.flags, collapse = ", "),
                            "].")
                )
        }

        if(!all(colnames(data) == annot[,sample.id])) {
          stop("Please make sure the column names of count data matched \'sample.id\' in sample annotation file, they need to be in the same length and order.")
        }

        if(fill == "CFB") {
                validate.baseline(annot, baseline.val, baseline.flag)
        }
}

#' @title Validate Baseline Values
#' @description Ensures that user-input `baseline.val` and `baseline.flag`
#'         parameters are valid with respect to the `annot` dataframe.
#' @inheritParams validate.annot
#'
#' @details Specifically, validates that `baseline.flag` value(s) are columns
#'         in `annot`, and that `baseline.val` value(s) occur at least once in
#'         their respective `baseline.flag` columns.
#' @references Xingpeng Li, Tatiana Gelaf Romer & Siddhartha Pachhai RVA - RNAseq Visualization Automation tool.
validate.baseline <- function(annot, baseline.val, baseline.flag) {
        cols <- colnames(annot)

        missing.baseline <- baseline.flag[!(baseline.flag %in% cols)]
        if(length(missing.baseline) > 0) {
                stop(paste0("Baseline flag column(s) [",
                            str_c(missing.baseline, collapse = ", "),
                            "] are missing from you sample annotation",
                            " file.")
                )
        }

        if(length(baseline.val) != length(baseline.flag)) {
                stop(paste0("Unequal lengths for baseline.flag and",
                            " baseline.val input. Each flag must have a",
                            " corresponding baseline value.")
                )
        }

        for(i in seq_along(baseline.flag)) {
                flag <- baseline.flag[[i]]
                val <- baseline.val[[i]]
                data <- annot[[flag]]
                if(all(data != val)) {
                        stop(paste0("Baseline value ", val, " for flag ",
                                    "\"", flag, "\" is not one of the",
                                    " values in column \"", flag, "\"",
                                    " of your annotation data.")
                        )
                }
        }
}

#' @title Validate Data Input
#' @description Ensures that the data input has the required formatting.
#'
#' @param data The wide-format dataframe with input data.
#'
#' @details Specifically, checks if `data` has rownaems and that all other
#'        columns can be coerced to numeric.
#' @references Xingpeng Li, Tatiana Gelaf Romer & Siddhartha Pachhai RVA - RNAseq Visualization Automation tool.
validate.data <- function(data) {
        if(is.null(rownames(data))) {
                stop("Data must have rownames labeling genes.")
        }
        for(data.col in colnames(data)) {
                validate.numeric(data, data.col)
        }
}

#' @title Validate Data in the Context of Annotation
#' @description Ensures that the annotation file matches the data file with
#'         respect to sample IDs. Throws warnings if there are discrepencies.
#' @references Xingpeng Li, Tatiana Gelaf Romer & Siddhartha Pachhai RVA - RNAseq Visualization Automation tool.
validate.data.annot <- function(data, annot, sample.id) {

        annot.samples <- annot[[sample.id]]
        data.samples <- colnames(data)

        extra.samples <- data.samples[!(data.samples %in% annot.samples)]
        if(length(extra.samples) > 0) {
                warning(paste0("The following samples in your data input data",
                               " are missing from your annotation input: [",
                               str_c(extra.samples, collapse = ", "),
                               "]. These samples will be omitted from analysis.")
                )
        }

        missing.samples <- annot.samples[!(annot.samples %in% data.samples)]
        if(length(missing.samples) > 0) {
                warning(paste0("The following samples in your annotation input",
                               " are missing from your data input data: [",
                               str_c(missing.samples, collapse = ", "),
                               "]. These samples will be omitted from analysis.")
                        )
        }
}
