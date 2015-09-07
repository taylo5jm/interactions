#' Check to see if a drug combination exists in eld_v4

#' @param T1 Signal X can be a treatment_name in treatment table
#' @param T2 Signal Y ""
#' @param DRUG_COMBINATION Boolean TRUE if T1 and T2 are both drugs
#' @return A character vector of length 1 indicating signal combination is present
#' @export isValidCombination

isValidCombination <- function(T1, T2, DRUG_COMBINATION) {

 if (DRUG_COMBINATION == TRUE) {
   # male/drug x/drug y/drug xy
   all_males <- dbGetQuery(CON, paste0("SELECT * FROM treatment
                                       NATURAL JOIN hsexperiment
                                       NATURAL JOIN donor_pair
                                       NATURAL JOIN donor
                                       NATURAL JOIN cell
                                       NATURAL JOIN sample
                                       NATURAL JOIN assay2sample
                                       NATURAL JOIN rnaseq_assay
                                       WHERE donor_sex = 'M'"))
   all_males_tbl <- table(all_males$treatment_name) %>%
     as.data.frame %>% arrange(Freq)
   all_m_drug_combos <- all_males_tbl[grep("/", all_males_tbl$Var1),]
   is_available <- grep(paste(T1, T2, sep = " / "), all_m_drug_combos$Var1)
   if (length(is_available) == 1) {
     T12 <- all_m_drug_combos$Var1[is_available]
     return(T12)
   }
   else {
     is_available <- grep(paste(T2, T1, sep = " / "), all_m_drug_combos$Var1)
     if (length(is_available) == 1) {
       T12 <- all_m_drug_combos$Var1[is_available]
       return(T12)
     }
     else {
       return("")
       #warning("Drug Combination in males is not available")
     }
   }
 }

 if (DRUG_COMBINATION == FALSE) {
   return("Error")
 }
}

#' Get all drug combinations present in eld_v4
#'
#' @param CON MySQLConnection object
#' @return list where first element are the names of all drug combinations in the database
#' @examples
#' \dontrun{
#' # CON <- dbConnect (#)
#' all_combos <- allDrugCombinations(CON)

allDrugCombinations <- function(CON) {
  this_combos <- dbGetQuery(CON, "SELECT * FROM treatment
                                 NATURAL JOIN hsexperiment
                                 NATURAL JOIN donor_pair
                                 NATURAL JOIN donor
                                 NATURAL JOIN cell
                                 NATURAL JOIN sample
                                 NATURAL JOIN assay2sample
                                 NATURAL JOIN rnaseq_assay
                                 WHERE treatment_name != 'None'")
  all_drug_combos <- this_combos[grep("/", this_combos$treatment_name),]
  all_drug_combos_tbl <- dplyr::count(all_drug_combos, treatment_name, donor_sex)
  return(list(all.drugs = all_drug_combos,
              all.drugs.tbl = all_drug_combos_tbl))
}

#' Retrieve libraries treated with a particular condition and filter by various parameters
#'
#' @param CON MySQLConnection object
#' @param condition treatment_name in treatment table. If searching for a drug
#'  combination, condition should be a character vector of length n > 1.
#' @param conc character vector of length 1 indicating concentration of drug
#' delivered.
#' @param sex character vector of length 1 indicating sex of donors to be retrieved
#' @param cell.type character vector of length 1 indicating cell type of donors
#' @return data.frame with each row representing a library that satisfies all conditions
#' specified in parameters
#' @examples
#' # not yet
getCondition <- function(CON, condition, conc = "", sex = "", cell.type = "") {
  if (length(condition) > 1) {
    condition <- paste(condition, sep = " / ")
  }
  this_combos <- dbGetQuery(CON, paste0("SELECT * FROM treatment
                                        NATURAL JOIN hsexperiment
                                        NATURAL JOIN donor_pair
                                        NATURAL JOIN donor
                                        NATURAL JOIN cell
                                        NATURAL JOIN sample
                                        NATURAL JOIN assay2sample
                                        NATURAL JOIN rnaseq_assay
                                        WHERE treatment_name RLIKE ", shQuote(condition)))
  this_combos <- filterMetadata(this_combos, conc, sex, cell.type)
  return(this_combos)
}

#' Filter data.frame of libraries retrieved by attributes such as concentration,
#' cell type, and donor sex
#'
#' @param conc character vector indicating concentration at which the drug was a
#' dministered (currently supports {"Low", "Mid", "High", "Super", "Cmax"})
#' @param sex character vector indicating donor sex ("M" or "F")
#' @param cell.type character vector of length 1 indicating cell type
#' @return filtered data.frame
#' @export filterMetadata

filterMetadata <- function(this_combos, conc, sex, cell.type) {
  if (conc != "") {
    this_combos <- this_combos[grep(
      conc, this_combos$alternate_treatment_name,
      ignore.case = TRUE),]
  }
  # filter sex
  if (sex != "") {
    this_combos <- dplyr::filter(this_combos,
                                 donor_sex == str_to_upper(sex))
  }
  # filter cell type
  if (cell.type != "") {
    this_combos <- dplyr::filter(this_combos,
                                 cell_name == str_to_upper(cell.type))
  }
  return(this_combos)
}

#' Get libraries for a specified signal interaction
#'
#' @param CON MySQLConnection
#' @param T1 Signal X, character vector of length 1 corresponding to treatment name in the treatment table
#' @param T2 Signal Y, character vector of length 1 corresponding to treatment name in the treatment table
#' @param T12 Signal X + Y, character vector of length 1 corresponding to treatment name for signal
#' combination in treatment table
#' @param cell.type Character vector of length 1 indicating cell type of interest
#' @return data.frame for libraries of interest where rows are libraries and columns are variables associated
#' with those libraries
#' @export getInteraction

getInteraction <- function(CON, CONDITIONS, cell.type, sex) {
   metadata <- dbGetQuery(CON, paste0("SELECT * FROM treatment
                                      NATURAL JOIN hsexperiment
                                      NATURAL JOIN donor_pair
                                      NATURAL JOIN donor
                                      NATURAL JOIN cell
                                      NATURAL JOIN sample
                                      NATURAL JOIN assay2sample
                                      NATURAL JOIN rnaseq_assay
                                      WHERE
                                      treatment_name RLIKE '", CONDITIONS$signal.xy, "'
                                      OR treatment_name LIKE '", CONDITIONS$signal.x, "'
                                      OR treatment_name LIKE '", CONDITIONS$signal.y, "'",
                                      "AND donor_sex = ", shQuote(sex)))
   return(metadata)
}

#' Generate a table of all treatment/concentration conditions for a given signal X and signal Y
#'
#' @param metadata data.frame of metadata retrieved from getInteraction()
#' @return base::table with number of libraries retrieved for each
#' treatment/concentration condition
#' @examples
#' \dontrun{
#' treatmentConcTable(metadata)
#' }
treatmentConcTable <- function(metadata) {
 treatment_concs <- strsplit(metadata$alternate_treatment_name, "_") %>%
   lapply(function(x) (paste(x[2], x[3], sep = ":"))) %>% unlist
 treatment_concs_tbl <- table(treatment_concs)
 return(treatment_concs_tbl)
}

#' Set conditions for interaction analysis
#' 
#' @param null Character vector indicating signal theta, often no stimuli or none
#' @param signal.x Character vector indicating signal x
#' @param signal.y Character vector indicating signal y
#' @param signal.xy Character vector indicating signal x + y. If missing, 
#'    function will assume signal.xy = signal x + signal y. (ex. "Sitagliptin / Metformin")
#' @return CONDITIONS list of length 4 where the elements are character vectors represnting
#'    signal theta, signal x, signal y, and signal x + y
#' @export setConditions
#' @examples
#' setConditions("None", "Sitagliptin", "Metformin")
#' setConditions("F:None", "M:None", "F:Sitagliptin", "M:Sitagliptin")

setConditions <- function(null, signal.x, signal.y, signal.xy) {
  if (missing(signal.xy)) {
    signal.xy <- paste(signal.x, signal.y, sep = " / ") 
  }  
  CONDITIONS <<- list(
    signal.theta = null,
    signal.x = signal.x,
    signal.y = signal.y,
    signal.xy = signal.xy
  )
}

#' return signal combination indices in treatment/concentration table
#'
#' @param T1 character vector of length 1 indicating Signal X
#' @param T2 character vector of length 1 indicating Signal Y
#' @param treatment.conc.tbl table (1d vector)
#' @return vector indices in table where Signal X/Signal Y combination
#' is present
#'
#' @examples
#' \dontrun{
#' comboIndices("Atorvastatin", "Pioglitazone", treatment_conc_tbl)
#' }
comboIndices <- function(T1, T2, treatment.conc.tbl, conc) {
 # treatment combo names in table
 t1_combo_name <- strsplit(T1, c())[[1]][1:5]
 t2_combo_name <- strsplit(T2, c())[[1]][1:5]
 # any concentration
 if (missing(conc)) {
   combo_name <- paste(paste0(t1_combo_name, collapse = ""),
                       paste0(t2_combo_name, collapse = ""),
                       sep = ":")
 }
 # filter for concentration
 if (is.character(conc)) {
   combo_name <- paste(paste0(t1_combo_name, collapse = ""),
                       paste0(t2_combo_name, collapse = ""),
                       conc,
                       sep = ":")
 }
 combo_tbl_index <- grep(combo_name, names(treatment.conc.tbl))
 return(combo_tbl_index)
}

#' return boolean indicating if drug combination was found in the database
#' (HS)
#' @param combo.tbl.index numeric vector indicating vector indices of signal combination
#' in library frequency distribution table
#' @return boolean indicating if the combination exists in the database
#' @export comboExists

comboExists = function(combo.tbl.index) {
 if (length(combo_tbl_index > 0)) {
   return(TRUE)
 }
 else {
   return(FALSE)
 }
}

#' Get vector of unique batches retrieved from a particular database
#' @param metadata data.frame with metadata retrieved for 4 conditions
#' @return vector of unique batches for a given metadata \code{data.frame}
#' @export uniqueBatches

uniqueBatches <- function(metadata) {
 if (is.character(metadata$alternate_treatment_name)) {
   this_batches <- strsplit(metadata$alternate_treatment_name, ":") %>%
     lapply(function(x) (x[1])) %>% unlist %>% unique
   this_batches <- this_batches[!is.na(this_batches)]
   return(this_batches)
 }
 else {
   warning("Error when attempting to retrieve unique batches")
 }
}

#' Get vehicles for libraries of interest
#'
#' @param CON MySQLConnection object
#' @param cell.type Character vector indicating which cell type to retrieve
#' @param metadata data.frame with metadata retrieved for 4 conditions
#' @return data.frame with libraries from \code{metadata} and respective vehicles
#' @export getVehicles
getVehicles <- function (CON, cell.type, metadata) {
 evenDonors <- function(metadata) {
   donors <- unique(metadata$donor_id)
   # libs for signal x, signal y, signal x + y have been retrieved
   if (try(nrow(metadata) / length(donors) == 3)) {
     return(donors)
   }
   else {
     warning("One or more donors are not mapped to > 1 conditions")
     return("")
   }

 }
 donors <- evenDonors(metadata)

 getVehicle <- function(donor, CON, cell.type) {
   n_batches <- uniqueBatches(metadata)
   if (length(n_batches) == 1) {
     dbGetQuery(CON, paste0("SELECT * FROM treatment
                            NATURAL JOIN hsexperiment
                            NATURAL JOIN donor_pair
                            NATURAL JOIN donor
                            NATURAL JOIN cell
                            NATURAL JOIN sample
                            NATURAL JOIN assay2sample
                            NATURAL JOIN rnaseq_assay
                            WHERE
                            treatment_name LIKE 'None'
                            AND treatment_label RLIKE 'ELD0107'
                            AND donor_id = ", donor,
                            " AND cell_name = ", shQuote(cell.type),
                            " LIMIT 1"))
   }
   else {
     dbGetQuery(CON, paste0("SELECT * FROM treatment
                            NATURAL JOIN hsexperiment
                            NATURAL JOIN donor_pair
                            NATURAL JOIN donor
                            NATURAL JOIN cell
                            NATURAL JOIN sample
                            NATURAL JOIN assay2sample
                            NATURAL JOIN rnaseq_assay
                            WHERE
                            treatment_name LIKE 'None'
                            AND donor_id = ", shQuote(donor),
                            "AND cell_name = ", shQuote(cell.type),
                            " LIMIT 1"))
   }
 }
 vehicles <- lapply(donors, getVehicle, CON, cell.type)
 vehicles <- data.frame(matrix(unlist(vehicles), nrow = length(vehicles), byrow = TRUE))
 colnames(vehicles) <- colnames(metadata)
 metadata <- rbind(metadata, vehicles, deparse.level = 0) %>%
   filter(cell_name == cell.type, donor_sex == "M")
 return(metadata)
}

#' Get matrix of expected (raw) counts from gene_rsem_count with RNA Seq assay-id 
#' (HS)
#' @param CON MySQLConnection object (should be connected to database with
#' table "gene_rsem_count" containing field "rnaseq_assay_id")
#' @param metadata data.frame where rows are libraries and columns are
#' variables associated with libraries
#' @param level character vector of length n = 1 indicating if "gene" or
#' "isoform" expected counts should be used. NOTE: ONLY SUPPORTS "gene" argument
#' for now
#' @return matrix where each row is a gene and each column is a library
#' @export getExpectedCounts

getExpectedCounts <- function(CON, metadata, level = "gene") {
 getCounts <- function(rnaseq_assay_id, level) {
   dbGetQuery(CON, paste0('SELECT * from ', paste0(level, '_rsem_count '),
                          'WHERE rnaseq_assay_id = ', rnaseq_assay_id))
 }
 counts <- sapply(metadata$rnaseq_assay_id, getCounts, level) %>%
   as.matrix()
 return(counts)
}

#' Get data.frame with gene metadata
#'
#' @param CON MySQLConnection object
#' @param counts matrix of counts where rows are genes and columns are libraries
#' @param level character vector indicating whether gene or isoform counts should be retrieved. NOTE:
#' only "gene" is supported as of now.
#' @return data.frame with gene metadata
#' @export getGeneMetadata

getGeneMetadata <- function(CON, counts, level = "gene") {
 gene_id <- counts[,2] %>% unlist()
 dbGetQuery(CON, paste0('SELECT * from hg19_', level)) %>%
   filter(hg19_gene_id %in% gene_id)
}

#' Strip counts data.frame of any extra variables
#'
#' @param counts matrix of counts where each row is a gene and each column is a library
#' @return counts matrix of where each row is a gene and each column is a library
#' @export asCountsMatrix
asCountsMatrix <- function(counts) {
 # # just expected counts
 counts <- counts[4,1:ncol(counts)] %>% sapply(unlist)
 return(counts)
}


checkMetadata <- function(metadata, all.concs = TRUE) {
  if (is.numeric(metadata)) {
    tryAllConcs <- function(CON, T1, T2, T12, cell_type, conc) {
      conc <- c("Low", "Mid", "High", "Super", "Cmax")
      this_list <- list()
      for (i in 1:length(this_list)) {
        this_list[[i]] <- getAllMetadata(CON, T1, T2, T12, cell_type, conc = conc[i])
      }
      return(this_list)
    }
    tryAllConcs(CON, T1, T2, T12, cell_type, conc)
    return("Metadata not retrieved for treatment/concentration condition")
  }
  else {
    return(metadata)
  }
}


#' Set conditions for interaction analysis
#' 
#' @param null Character vector indicating signal theta, often no stimuli or none
#' @param signal.x Character vector indicating signal x
#' @param signal.y Character vector indicating signal y
#' @param signal.xy Character vector indicating signal x + y. If missing, 
#'    function will assume signal.xy = signal x + signal y. (ex. "Sitagliptin / Metformin")
#' @return CONDITIONS list of length 4 where the elements are character vectors represnting
#'    signal theta, signal x, signal y, and signal x + y
#' @export setConditions
#' @examples
#' setConditions("None", "Sitagliptin", "Metformin")
#' setConditions("F:None", "M:None", "F:Sitagliptin", "M:Sitagliptin")

setConditions <- function(null, signal.x, signal.y, signal.xy) {
  if (missing(signal.xy)) {
    signal.xy <- paste(signal.x, signal.y, sep = " / ") 
  }  
  CONDITIONS <<- list(
    null = null,
    signal.x = signal.x,
    signal.y = signal.y,
    signal.xy = signal.xy
  )
}


