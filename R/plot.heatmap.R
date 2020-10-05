#' @title Plot a CPM Heatmap
#'
#' @description An alias for `plot.heatmap.expr(annot, cpm, fill = "CPM", ...)`.
#'
#' @inheritParams plot.heatmap.expr
#'
#' @export
plot.heatmap.cpm <- function(cpm,
                             annot,
                             title = "RVA CPM Heatmap",
                             ...) {
  plot.heatmap.expr(cpm, annot, fill = "CPM", title = title, ...)
}

#' @title Plot a CFB Heatmap
#'
#' @description An alias for `plot.heatmap.expr(annot, cpm, fill = "CFB", ...)`.
#'
#' @inheritParams plot.heatmap.expr
#'
#' @export
plot.heatmap.cfb <- function(cpm,
                             annot,
                             title = "RVA CFB Heatmap",
                             ...) {
  plot.heatmap.expr(cpm, annot, fill = "CFB", title = title, ...)
}

#' @title Plot Heatmap From Raw CPM
#'
#' @description Create a heatmap with either CFB or CPM averaged across
#'         individual samples.
#'
#' @param data A wide-format dataframe with geneid rownames, sample column
#'         names, and fill data matching `input.type`.
#' @param annot A long-format dataframe with any pertinent treatment data about
#'         the samples. The only required column is one titled the `sample.id`
#'         value with values matching the column names of sample IDs in `data`.
#'         Additional columns can contain information such as treament
#'         compounds, dates of sample collection, or dosage quantities.
#' @param sample.id The column name to specify sample ID.
#' @param annot.flags A vector of column names corresponding to column names
#'         in `annot` which will be used to define the x-axis for the heatmap.
#'         Default = `c("day", "dose")`.
#' @param ct.table.id.type The gene id format in `data` should be one of:
#'         ACCNUM, ALIAS, ENSEMBL, ENSEMBLPROT, ENSEMBLTRANS,
#'         ENTREZID, ENZYME, EVIDENCE, EVIDENCEALL, GENENAME, GO, GOALL, IPI,
#'         MAP, OMIM, ONTOLOGY, ONTOLOGYALL, PATH, PFAM, PMID, PROSITE, REFSEQ,
#'         SYMBOL, UCSCKG, UNIGENE, UNIPROT.
#' @param gene.id.type The gene id format of `gene.names`, should be one of:
#'         ACCNUM, ALIAS, ENSEMBL, ENSEMBLPROT, ENSEMBLTRANS,
#'         ENTREZID, ENZYME, EVIDENCE, EVIDENCEALL, GENENAME, GO, GOALL, IPI,
#'         MAP, OMIM, ONTOLOGY, ONTOLOGYALL, PATH, PFAM, PMID, PROSITE, REFSEQ,
#'         SYMBOL, UCSCKG, UNIGENE, UNIPROT.
#' @param gene.names A character vector or list of ensembl IDs for which to
#'         display gene information. If `NULL`, all genes will be included.
#'         Default = `NULL`.
#' @param gene.count The number of genes to include, where genes are selected
#'         based on ranking by values in `fill`. Default = 10.
#' @param title A title for the heatmap. Default = "RVA Heatmap".
#' @param fill One of `c("CPM", "CFB")` to fill the heatmap cells with.
#'         Default = "CFB".
#' @param baseline.flag A character vector of column names. If `fill = "CFB"`,
#'         these columns in `annot` contain the values to compare
#'         across. Ignored if `fill = "CPM"`. Default = "timepoint".
#' @param baseline.val A character vector of values. This vector must be the
#'         same length as `baseline.flag`, and the value at each index must
#'         represent a value from the column given by the corresponding index
#'         in `baseline.flag`. The samples corresponding to these values will
#'         be used as a baseline when calculating CFB. Ignored if
#'         `fill = "CPM"`. Default = "Week 0".
#' @param plot.save.to The address to save the heatmap plot.
#' @param input.type One of `count` or `cpm` indicating what the input data type
#'         is. If `count`, the CPM of the input data will be calculated using
#'         [edgeR::cpm()]. Default = `count`.
#'
#' @return The function returns a list with 2 items:
#' \item{df.sub}{"A data frame of change from baselines values (fill = CFB in this example) for each gene id that is divided by a combination of treatment group and time point}
#' \item{gp}{A Heatmap object from ComplexHeatmap which can be plotted}
#'
#'
#'
#' @importFrom dplyr left_join mutate filter select_at vars group_by_at ungroup group_by select summarize arrange_at
#' @import edgeR
#' @importFrom edgeR cpm
#' @import tidyr
#' @importFrom tidyr unite
#' @importFrom ggplot2 facet_wrap ggsave ggplot
#' @importFrom stringr word str_flatten
#' @importFrom tibble as_tibble
#' @import org.Hs.eg.db
#' @import clusterProfiler
#' @importFrom ComplexHeatmap Heatmap columnAnnotation
#' @importFrom circlize colorRamp2
#' @importFrom data.table as.data.table dcast.data.table melt.data.table dcast melt
#' @importFrom grid grid.text gpar
#'
#' @export plot.heatmap.expr
#'
#'
#' @references Xingpeng Li,Tatiana Gelaf Romer & Aliyah Olaniyan, RVA - RNAseq Visualization Automation tool.
#'
#' @details The function takes raw CPM data and returns both a list containing a data frame
#'          with values based on the fill parameter and a heatmap plot.
#'
#'
#' @examples
#' plot.heatmap.expr(data = example_count_table,
#' annot = sample_annotation_table,
#' plot.save.to = "~/address_to_folder/heatmap_plot.png")
#'
#' #Save figures using ComplexHeatmap
#' library(ComplexHeatmap)
#' png("heatmap_plot_name.png", width = 500, height = 500)
#' draw(hm.expr$gp)
#' dev.off()
#'


plot.heatmap.expr <- function(data = count,
                              annot = meta,
                              sample.id = "sample_id",
                              annot.flags = c("day", "Treatment", "tissue"),
                              ct.table.id.type = "ENSEMBL",
                              gene.id.type = "SYMBOL",
                              gene.names = NULL,
                              gene.count = 10,
                              title = "RVA Heatmap",
                              fill = "CFB",
                              baseline.flag = "day",
                              baseline.val = "0",
                              plot.save.to = NULL,
                              input.type = "count") {

  options(warn=-1)
  suppressWarnings({
  suppressMessages({

  validate.geneid.flag(ct.table.id.type, "ct.table.id.type")
  validate.geneid.flag(gene.id.type, "gene.id.type")
  validate.flag(fill, "fill", c("CFB", "CPM"))
  validate.flag(input.type, "input.type", c("count", "cpm"))
  validate.annot(data, annot, annot.flags, sample.id, fill, baseline.flag, baseline.val)
  validate.data(data)
  validate.data.annot(data, annot, sample.id)

  user.title = title

  if(input.type == "count") {
      data <- cpm(data, log = TRUE) %>%
          as_tibble(rownames = "geneid")
  } else {
      data <- data %>%
          as_tibble(rownames = "geneid")
  }

  data <- reformat.ensembl(data, ct.table.id.type)

  if(!is.null(gene.names)) {
    gene.names <- transform.geneid(gene.names,
                                   from = gene.id.type,
                                   to = ct.table.id.type)
    data <- data %>% filter(geneid %in% gene.names)
  }

  long <- data %>%
    as.data.table %>%
    melt(id.vars = "geneid",
         variable.name = "sample_id",
         value.name = "cpm") %>%
    merge(annot, by = "sample_id") %>%
    select_at(vars(c("geneid", "subject_id", annot.flags, "cpm")))

  long$cpm <- as.numeric(long$cpm)

  if(fill == "CFB") {
    long <- calc.cfb(long, annot, baseline.flag, baseline.val)
    fill.var <- "cfb"
  } else if (fill == "CPM") {
    long <- long %>%
      filter(!is.na(cpm))
    fill.var <- "cpm"
  }

  long <- long %>%
    as.data.table %>%
    group_by_at(vars(c("geneid", annot.flags))) %>%
    dplyr::summarize(!!fill.var := mean(UQ(as.name(fill.var)))) %>%
    ungroup()

  if(is.null(gene.names)) {
    gene.names <- long %>%
      group_by(geneid)
    if(fill == "CFB") {
      gene.names <- gene.names %>%
        dplyr::summarize(cfb = max(abs(cfb)))
    }
    if(fill == "CPM") {
      gene.names <- gene.names %>%
        dplyr::summarize(cpm = max(cpm))
    }
    gene.names <- gene.names %>%
      ungroup() %>%
      arrange_at(fill.var, dplyr::desc) %>%
      .$geneid %>%
      .[1:gene.count]

    long <- long %>%
      filter(geneid %in% gene.names)
  }

  annot <- long %>%
    select_at(vars(annot.flags)) %>%
    unique %>%
    as.data.frame

  wide.df <- long %>%
    as.data.table %>%
    dcast(geneid ~ ..., value.var = fill.var)

  wide <- wide.df %>%
    dplyr::select(-geneid) %>%
    as.matrix

  #new code - if cpm then do z-score transofrm
  if (fill == "CPM"){
    wide = t(scale(t(wide)))
  }
  #change anoot flags into factors
  annot[,annot.flags] = data.frame(lapply(annot[,annot.flags], as.factor))

  #new code - change color scale to adjust for z - score
  range_cpm =c(wide)

  colors <- switch(fill,
                   CFB = colorRamp2(c(-2, 0, 2), c("blue", "grey", "red")),
                   CPM = get.cpm.colors(range_cpm))# use to be long$cpm inside
  title <- switch(fill,
                  CFB = "CFB (log2CPM)",# originally "log2(CFB)",
                  CPM = "z-score (log2CPM)") # use to be log2 (CPM)
  colnames(wide) <- NULL
  gene.display <- transform.geneid(wide.df$geneid,
                                     from = ct.table.id.type,
                                     to = "SYMBOL")
  ctid <- gene.display
  colnames(ctid)[1] <- "FROM"
  ctid <- ctid[!duplicated(ctid$FROM), ] #if there are multiple symbols, pick the first one
  gene.display <- ctid %>% left_join(ctid)
  colnames(gene.display)[1] = ct.table.id.type

  rownames(wide) <- gene.display[,2]
  set.seed(200) #fix the color scheme
  gp <- Heatmap(wide, col = colors, name = "Heatmap",
                na_col = "black",
                cluster_columns = F,
                heatmap_legend_param = list(title = title),
                column_title = user.title,
                row_title = "Genes",
                bottom_annotation = columnAnnotation(df=annot))
  gp
  if(is.null(plot.save.to)){
    print("Plot file name not specified, a plot in Heatmap object will be output to the first object of the return list!")
  }else{
    png(filename = plot.save.to)
    print(gp)
    dev.off()
  }
  #wide.df
  return(list(gp = gp,
              df.sub = wide.df))
  })
  })
}

#' @title Get CPM Colors
#' @description This function creates the color gradient for the cpm data.
#' @param data The CPM dataset.
#' @importFrom stats median
#'
get.cpm.colors <- function(data) {
  max.val <- max(data)
  med.val <- median(data)
  min.val <- min(data)

  colorRamp2(c(min.val, med.val, max.val), c("blue", "grey", "red"))
}

#' @title Reformat Ensembl GeneIDs
#' @description This is the function to exclude the version number from the input ensembl type gene ids.
#' @param logcpm The input count table transformed into log counts per million.
#' @param ct.table.id.type The gene id format in `logcpm` should be one of: ACCNUM, ALIAS, ENSEMBL, ENSEMBLPROT,
#'  ENSEMBLTRANS, ENTREZID, ENZYME, EVIDENCE, EVIDENCEALL, GENENAME, GO, GOALL, IPI, MAP, OMIM,
#'   ONTOLOGY, ONTOLOGYALL, PATH, PFAM, PMID, PROSITE, REFSEQ, SYMBOL, UCSCKG, UNIGENE, UNIPROT.
#' @importFrom stringr word
#' @importFrom dplyr %>% mutate
#'
reformat.ensembl <- function(logcpm, ct.table.id.type){
  if(ct.table.id.type == "ENSEMBL"){
    logcpm <- logcpm %>% mutate(geneid = word(geneid, sep = '\\.'))}
  return(logcpm)
}

#' @title Transform GeneIDs
#' @description This is the function to transform the input gene id type to another gene id type.
#' @param gene.names Genes,in the format of character vector, to be transformed.
#' @param from The gene id format of `gene.names`, should be one of: ACCNUM, ALIAS, ENSEMBL, ENSEMBLPROT,
#'  ENSEMBLTRANS, ENTREZID, ENZYME, EVIDENCE, EVIDENCEALL, GENENAME, GO, GOALL, IPI, MAP, OMIM,
#'   ONTOLOGY, ONTOLOGYALL, PATH, PFAM, PMID, PROSITE, REFSEQ, SYMBOL, UCSCKG, UNIGENE, UNIPROT.
#' @param to The new gene id format should be one of: ACCNUM, ALIAS, ENSEMBL, ENSEMBLPROT,
#'  ENSEMBLTRANS, ENTREZID, ENZYME, EVIDENCE, EVIDENCEALL, GENENAME, GO, GOALL, IPI, MAP, OMIM,
#'   ONTOLOGY, ONTOLOGYALL, PATH, PFAM, PMID, PROSITE, REFSEQ, SYMBOL, UCSCKG, UNIGENE, UNIPROT.
#' @importFrom dplyr %>%
#' @import org.Hs.eg.db
#' @importFrom clusterProfiler bitr
#'
transform.geneid <- function(gene.names,
                             from = gene.id.type,
                             to = ct.table.id.type){
  out <- gene.names %>%
    clusterProfiler::bitr(fromType = from, toType = to,OrgDb = org.Hs.eg.db)
  return(out)
}

#' @title Calculate CFB
#' @description This function calculates the change from baseline.
#' @param data Dataframe with subject id, annotation flag, gene id and cpm value
#'             (from count tables) columns.
#' @param annot A long-format dataframe with any pertinent treatment data about
#'         the samples. The only required column is one titled the `sample.id`
#'         value with values matching the column names of sample IDs in `data`.
#'         Additional columns can contain information such as treament
#'         compounds, dates of sample collection, or dosage quantities.
#' @param baseline.flag A character vector of column names. These columns in `annot`
#'         contain the values to compare across.
#' @param baseline.val A character vector of values. This vector must be the
#'         same length as `baseline.flag`, and the value at each index must
#'         represent a value from the column given by the corresponding index
#'         in `baseline.flag`.
#' @importFrom dplyr %>% mutate filter select_at
#' @importFrom stringr str_flatten
#' @importFrom tidyr unite separate
#' @importFrom data.table dcast melt
#' @importFrom rlang UQ
#'
calc.cfb <- function(data, annot, baseline.flag, baseline.val) {
  cast.formula <- paste0("... ~ ", paste0(baseline.flag, collapse = "+"))
  relevant.vars <- annot %>%
    unite(col = "var", baseline.flag, remove = F) %>%
    .$var %>%
    unique
  baseline.var <- str_flatten(baseline.val, collapse = "_")
  relevant.vars <- relevant.vars[relevant.vars != baseline.var]

  data %>%
    dcast(cast.formula, value.var = "cpm", fun.aggregate = mean) %>%
    melt(value.name = "cpm", variable.name = "variable",
         measure.vars = relevant.vars) %>%
    filter(!is.na(UQ(as.name(baseline.var))) & !is.na(variable)) %>%
    separate(variable, into = baseline.flag, sep = "_") %>%
    dplyr::mutate(cpm = as.numeric(cpm),
                  !!baseline.var := as.numeric(UQ(as.name(baseline.var)))) %>%
    dplyr::mutate(cfb = cpm-UQ(as.name(baseline.var))) %>%
    filter(!is.na(cfb)) %>%
    select_at(vars(-c(baseline.var, "cpm")))
}
