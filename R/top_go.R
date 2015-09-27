#' Get binary factor vector for genes
#'
#' @param genes data.frame of gene metadata (data$voom or data$dgelist)
#' @param interaction.tbl gene interaction table constructed from genes classified into
#' an interaction mode
#' @param filter.mode Boolean FALSE when all genes in the interaction table are used
#' to construct the gene vector. TRUE when the gene vector should be constructed only
#' from genes classified into a particular mode
#' @return factor vector where genes of interest are represented by 1 and the opposite
#' are represented by 0
#' @export getGeneVector
#' @examples
#' # construct gene vector from genes exhibiting high stabilization mode
#' \dontrun{
#' this_gene_vector <- getGeneVector(data$dgelist$genes, edgeR_interaction.tbl,
#' filter.mode = "High.Stab")
#' }
#' 
getGeneVector <- function(genes, interaction.tbl, filter.mode = FALSE) {
  
    # all DEGs or just genes with Mode matching filter.mode 
    filterInteractionTable <- function(interaction.tbl, filter.mode) {
      if (filter.mode == FALSE)
        interaction.tbl
    
      else {
        dplyr::filter(interaction.tbl, Mode == filter.mode)
      }
    }
  
  interaction.tbl <- filterInteractionTable(interaction.tbl, filter.mode)

  makeGeneVector <- function(genes, interaction.tbl) {
    # gene universe
    gene_vector <- rep(-1, length(genes$gene_symbol))
    
    # differentially expressed
    de_indices <- which(genes$gene_symbol %in% interaction.tbl$gene_symbol)
    nde_indices <- which(!(genes$gene_symbol %in% interaction.tbl$gene_symbol))
    # genes of interest
    gene_vector[de_indices] <- "1"
    gene_vector[nde_indices] <- "0"
    gene_vector <- factor(gene_vector)
    names(gene_vector) <- genes$gene_symbol
    return(gene_vector)
  }
  makeGeneVector(genes, interaction.tbl)
}

#' Generate gene ontology table and data with topGO from interaction data.frame
#' using Fisher's Exact Test
#'
#' @import topGO
#' @param data EList or DGEList with massaged counts and all 0 rows removed 
#'        (see dge2voom_q())
#' @param interaction.tbl data.frame where rows are genes and 
#'        (see interactionGeneTable())
#' @param ontology character vector of length n = 1 indicating gene ontology
#' @param filter.mode Boolean FALSE when all genes in the interaction table are used
#' to construct the gene vector. TRUE when the gene vector should be constructed only
#' from genes classified into a particular mode
#' @return list where list[[1]] is a topGO table and list[[2]] is a topGO data object
#' @export getGOTable
#' @examples
#' # not really
getGOTable <- function(data, interaction.tbl, ontology, filter.mode = FALSE) {
  if (class(interaction.tbl) == "Interaction" & missing(v)) {
    interaction.tbl <- interaction.tbl@interaction.tbl
    v <- interaction.tbl@voom
  }

#  if (!(is.data.frame(interaction.tbl)) & is.list(interaction.tbl)) {
#    interaction.tbl <- interaction.tbl$interaction.tbl
#  }

  if (is.data.frame(interaction.tbl)) {

    # every gene NI
    if (nrow(interaction.tbl) == 0) {
      print("no genes with interactions classified")
      return(list(go.results = 0,
                  go.data = 0))
    }
    # > 0 genes classified
    if (nrow(interaction.tbl) > 0) {
      # gene vector where 1 are genes classified into a mode
      if (filter.mode == FALSE) {
        interaction.tbl <- interaction.tbl
      }
      # gene vector where 1 are genes classified in filter.mode
      if (filter.mode != FALSE) {
        stopifnot(is.character(filter.mode))
        interaction.tbl <- dplyr::filter(interaction.tbl, Mode == filter.mode)
      }

      int_genes <- getGeneVector(data$genes, interaction.tbl)

      # > 0 genes found -------------------------------------------------------
      if (is.factor(int_genes) & nlevels(int_genes) == 2) {

        go_data <- new("topGOdata", ontology = ontology, allGenes = int_genes,
                       annot = annFUN.org, mapping = "org.Hs.eg.db",
                       ID = "symbol")
        sig_genes <- topGO::sigGenes(go_data)
        go_stats <- topGO::termStat(go_data)

        # Fisher's Exact Test ---------------------------------------------------------
        test_stat <- new("classicCount", testStatistic = GOFisherTest,
                         name = "Fisher test")
        fishers_result <- topGO::getSigGroups(go_data, test_stat)

        all_results <- topGO::GenTable(go_data, classicFisher = fishers_result,
                                orderBy = "classicFisher", 
                                ranksOf = "classicFisher", topNodes = 20)
        
        return(list(go.results = all_results,
                    sig.genes = sig_genes,
                    term.stat = go_stats,
                    fishers.result = fishers_result))
      }
      else {
        return(list(go.results = NULL,
                    go.data = NULL))
      }
    }
  }
else {
  return(list(go.results = 0,
              go.data = 0))
    }
}

#' save topGO table as .csv and topGO DAG as .png
#'
#' @param go.res list where go.res[[1]] is a topGO table and go.res[[2]] is the topGO
#' data object
#' @return saved .csv file with topGO table in signal combination output directory
#' and saved .png file with topGO graph in output dir
#' @export saveTopGoRes

saveTopGoRes <- function(go.res) {
  signal_combo <- paste(T1, T2, sep = "_")
  this_dir <- file.path(root_dir, "signal_combinations", signal_combo,
                        "gene_ontology")

  if(!(file.exists(this_dir))) {
    dir.create(this_dir)
  }

  table_path <- file.path(this_dir, 'topGO_table.csv')
  write.csv(go.res[[1]], table_path)

  # go graph
  graph_path <- file.path(this_dir, 'topGO_graph.png')
  png(graph_path, height = 1200, width = 1200)
  topGO::showSigOfNodes(go.res[[2]], score(result_ks_elim))
  dev.off()
}

#' get topGO tables and significant genes for a particular ontology with gene
#' vectors for every interaction mode present in the interaction table
#'
#' @param v Elist constructed with dge2voom_q()
#' @param interaction.tbl interaction data.frame generated from
#' interaction mode classification
#' @param ontology character vector of length 1 indicating ontology ("MF", "CC", "BP")
#' @return list where first element is a topGO table and the second element is a
#' character vector of significant genes
#' @export getGOTermsForModes
#' @examples
#' \dontrun{
#' v <- dge2voom_q(counts, genes)
#' interactionModeClassification()
#' mf_terms <- getGOTermsForModes(v, edgeR_interaction.tbl, ontology = "MF")
#' }
getGOTermsForModes <- function(data, interaction.tbl, ontology) {
  go_terms <- list()
  ADDITIVE_VECTORS <- c("Sym.Right", "Sym.Left",
                        "Step.Up", "Step.Down")
  MODES <- MODES[MODES != "NI" & MODES != "A" & MODES != "UC"
                 & !(MODES %in% ADDITIVE_VECTORS)]
  
  # gene set enrichment for each interaction mode
  for (i in 1:length(MODES)) {
    print(paste0("Getting GO terms for ", MODES[i]))
    this_table <- getGOTable(data, interaction.tbl, ontology, filter.mode = MODES[i])
    go_terms[[i]] <- this_table
  }
  names(go_terms) <- MODES
  return(go_terms)
}

#' get GO tables for all ontologies, given an interaction table data.frame constructed
#' from geneInteractionTable()
#'
#' @param v EList constructed with voom()
#' @param interaction.tbl data.frame constructed with geneInteractionTable()
#' @return nested list where each element is a list containing GO tables
#' and significant genesgenerated with topGO using each mode as a gene vector
#' @export getGOTermsForAllModes
#' @examples
#' \dontrun{
#' v <- dge2voom_q(counts, genes)
#' interactionModeClassification()
#' all_go_terms <- getGOTermsForModesAll(v, edgeR_interaction.tbl)
#' }
getGOTermsForAllModes <- function(data, interaction.tbl) {
  all_go_terms <- list(BP = getGOTermsForModes(data, interaction.tbl, "BP"),
                       MF = getGOTermsForModes(data, interaction.tbl, "MF"),
                       CC = getGOTermsForModes(data, interaction.tbl, "CC"))
  names(all_go_terms) <- c("BP", "MF", "CC")
  return(all_go_terms)
}





