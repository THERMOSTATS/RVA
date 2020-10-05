#' @title Plot gene expression
#' @description This is the function to process the gene count table to show gene expression variations over time or across groups.
#'
#' @param data Count table in the format of dataframe with gene id as row.names.
#' @param anno Annotation table that provides design information.
#' @param gene.names Genes to be visualized, in the format of character vector.
#' @param ct.table.id.type The gene id format in `data` should be one of: ACCNUM, ALIAS, ENSEMBL, ENSEMBLPROT,
#'  ENSEMBLTRANS, ENTREZID, ENZYME, EVIDENCE, EVIDENCEALL, GENENAME, GO, GOALL, IPI, MAP, OMIM,
#'   ONTOLOGY, ONTOLOGYALL, PATH, PFAM, PMID, PROSITE, REFSEQ, SYMBOL, UCSCKG, UNIGENE, UNIPROT.
#' @param gene.id.type The gene id format of `gene.names`, should be one of: ACCNUM, ALIAS, ENSEMBL, ENSEMBLPROT,
#'  ENSEMBLTRANS, ENTREZID, ENZYME, EVIDENCE, EVIDENCEALL, GENENAME, GO, GOALL, IPI, MAP, OMIM,
#'   ONTOLOGY, ONTOLOGYALL, PATH, PFAM, PMID, PROSITE, REFSEQ, SYMBOL, UCSCKG, UNIGENE, UNIPROT.
#' @param treatment The column name to specify treatment groups.
#' @param sample.id The column name to specify sample IDs.
#' @param time The column name to specify different time points.
#' @param log.option Logical option, whether to log2 transform the CPM as y-axis. Default = True.
#' @param plot.save.to The address to save the plot from simplified cutoff combination with FDR of 0.01, 0.05, 0.1, and 0.2.
#' @param input.type One of `count` or `cpm` indicating what the input data type is. If `count`, the CPM of the input data will be
#'  calculated using [edgeR::cpm()]. Default = `count`.
#'
#' @return The function returns a ggplot object.
#'
#' @import data.table
#' @importFrom dplyr left_join mutate n_distinct filter select
#' @importFrom tidyr pivot_longer
#' @importFrom stringr word
#' @importFrom tibble as_tibble
#' @importFrom edgeR cpm
#' @import org.Hs.eg.db
#' @importFrom ggplot2 facet_wrap ggsave ggplot aes_ geom_boxplot geom_point scale_color_brewer theme_classic
#' @importFrom clusterProfiler bitr
#'
#' @export plot.gene
#'
#' @references Xingpeng Li,Tatiana Gelaf Romer & Aliyah Olaniyan, RVA - RNAseq Visualization Automation tool.
#'
#' @details The function takes the gene counts and returns a ggplot that shows gene expression variation over time or group.
#'
#'
#' @examples
#' plot <- plot.gene(data = example_count_table,
#' anno = sample_annotation_table,
#' plot.save.to = "~/address_to_folder/gene_plot.png")
#'
#' #Save figures using ggplot2
#' library(ggplot2)
#' ggsave(plot,
#' "gene_plot_name.png",
#' device = "png",
#' width = 100,
#' height = 100,
#' dpi = 200,
#' limitsize = FALSE)
#'

plot.gene <- function(data = dat,
                      anno = meta,
                      gene.names = c("AAAS", "A2ML1", "AADACL3"),
                      ct.table.id.type = "ENSEMBL",
                      gene.id.type = "SYMBOL",
                      treatment = "Treatment",
                      sample.id = "sample_id",
                      time = "day",
                      log.option = T,
                      plot.save.to = NULL,
                      input.type = "count") {

        options(warn=-1)
        suppressWarnings({
        suppressMessages({

        validate.geneid.flag(ct.table.id.type, "ct.table.id.type")
        validate.geneid.flag(gene.id.type, "gene.id.type")
        validate.flag(input.type, "input.type", c("count", "cpm"))
        validate.annot(data, anno, c(time, treatment, sample.id), sample.id)
        validate.data(data)
        validate.data.annot(data, anno, sample.id)


        if(input.type == "count") {
                logcpm <- cpm(data, log = TRUE) %>%
                        as_tibble(rownames = "geneid")
        } else {
                logcpm <- data %>%
                        as_tibble(rownames = "geneid")
        }

        logcpm <- reformat.ensembl(logcpm, ct.table.id.type)

        gene.names <- transform.geneid(gene.names,
                                       from = gene.id.type,
                                       to = ct.table.id.type) %>% pull(ct.table.id.type)

        validate.genes.present(logcpm$geneid, gene.names)

        #clean meta data
        anno <- anno[,c(sample.id, treatment, time)]

        ##merge with meta data and clean up for ggplot
        pd.dat <- logcpm %>%
                filter(geneid %in% gene.names) %>%
                pivot_longer(-geneid,
                             values_to = "exprs",
                             names_to = sample.id) %>%
                left_join(anno, by = sample.id) %>%
                mutate(exprs = as.numeric(exprs)) %>%
                mutate(!!time := as.factor(UQ(as.name(time))))

        #convert the gene names back to gene symbol for the plot
        lib <- transform.geneid(pd.dat$geneid,
                                from = ct.table.id.type,
                                to = gene.id.type)

        pd.dat <- pd.dat %>% left_join(lib, by = c("geneid" = ct.table.id.type))
        pd.dat$new.geneid <- pd.dat[[gene.id.type]]



        #ggplot
        pgene <- ggplot(pd.dat, aes_(x=as.name(time),
                                     y=~exprs,
                                     color=as.name(treatment))) +
                geom_boxplot(position=position_dodge(0.8)) +
                geom_point(position=position_dodge(0.8)) +
                scale_color_brewer(palette="Dark2") +
                theme_classic()

        #check if there is only 1 gene input
        if(n_distinct(pd.dat$geneid) > 1){
                pgene <- pgene + facet_wrap(~new.geneid)
        }
        #check if y-axis log transformed
        if(log.option){
                pgene <- pgene + ylab("Gene expression (log2CPM)")
        }else{
                pgene <- pgene + ylab("Gene expression (CPM)")
        }


        if(is.null(plot.save.to)){
                print("Plot file name not specified, a plot in ggplot object will be output to the second object of the return list!")
        }else{
                ggsave(filename = plot.save.to,
                       plot = pgene,
                       dpi = 300,
                       units = "in",
                       device='png')
        }
        pgene
        pd.dat <- pd.dat %>% dplyr::select(-new.geneid)
        return(list(pgene, pd.dat))

        })
        })
}

#' @title Reformat Ensembl GeneIDs
#' @description This is the function to exclude the version number from the input ensembl type gene ids.
#' @param logcpm The input count table transformed into log counts per million.
#' @param ct.table.id.type The gene id format in `logcpm` should be one of: ACCNUM, ALIAS, ENSEMBL, ENSEMBLPROT,
#'  ENSEMBLTRANS, ENTREZID, ENZYME, EVIDENCE, EVIDENCEALL, GENENAME, GO, GOALL, IPI, MAP, OMIM,
#'   ONTOLOGY, ONTOLOGYALL, PATH, PFAM, PMID, PROSITE, REFSEQ, SYMBOL, UCSCKG, UNIGENE, UNIPROT.
#' @importFrom stringr word
#' @importFrom dplyr %>% mutate

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
