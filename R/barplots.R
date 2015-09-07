#' Color bars in geom_bar plot given an interaction class and ggplot2
#'
#' @param int_class character vector of length n = 1 indicating interaction class ("Positive", "Negative", "NI",
#' "Null", "Ambiguous)
#' @param gg_plot ggplot2 object
#' @return ggplot2 object with bars colored with respect to interaction class
#'

colorBarplot <- function(int_class, gg_plot) {
  if (int_class == "Positive") {
    return(gg_plot + theme(panel.background = element_rect(fill = "#99CCFF")))
  }
  if (int_class == "Negative") {
    return(gg_plot + theme(panel.background = element_rect(fill = "#FF6666")))
  }
}


geom_bar_constructor <- function(logfc, int_mode, gene) {
  stopifnot(is.matrix(logfc), is.character(int_mode), is.character(gene))
  gg_plot <- qplot(x = COND, y = as.numeric(logfc),
                   geom = "bar", stat = "identity", xlab = "Condition",
                   ylab = "logFC", main = paste0(gene,
                                                    " Interaction Profile - ",
                                                 int_mode))
  colorBarplot(findClass(int_mode), gg_plot)
}


# barplots for all modes
expedite_mode_gg_bar <- function(modes, int_genes) {
  lapply(unique(modes), mode_geom_bar, int_genes)
}

#' Make barplot of outcome vector distribution with base::barplot
#'
#' @param pw_s character vector
#' @param remove.ni Boolean TRUE when null outcome vector should be removed from the frequency distribution table
#' @param geom_bar Boolean TRUE when ggplot2 with geom_bar should be used instead of base R graphics
#' @param log character vector of length 1 indicating which variable (if any) to plot log of
#' @return base::barplot of outcome vector frequency distribution
#'


qplotOutcomeVectorDist <- function(pw_s, remove.ni = TRUE, geom_bar = TRUE, log = "") {
  if (!(is.character(pw_s))) {
    pw_s <- apply(pw_s, 1, stringCoerce)
  }

  if (is.character(pw_s)) {
    outcome_dist <- ovTable(pw_s)
    if (remove.ni == TRUE) {
      outcome_dist <- filter(outcome_dist, Outcome.Vector != "0,0,0,0,0,0")
    }
  }

#   if (is.table(pw_s) | is.data.frame(pw_s)) {
#     outcome_dist <- pw_s
#     colnames(outcome_dist) <- c("Outcome.Vector", "Freq")
#     if (remove.ni == TRUE) {
#       outcome_dist <- filter(outcome_dist, Outcome.Vector != "0,0,0,0,0,0")
#     }
#   }
  if (geom_bar == FALSE) {
    if (log == "") {
      print(barplot(outcome_dist$Freq))
    }
    else {
      print(barplot(outcome_dist$Log.Freq))
    }
  }

  if (geom_bar == TRUE) {
    if (log == TRUE) {
      this_log <- "y"
      this_main <- "Log Frequency Distribution of Outcome Vectors"
    }
    else {
      this_log <- ""
      this_main <- "Frequency Distribution of Outcome Vectors"
    }
      p <- qplot(x = 1:nrow(outcome_dist), y = Freq, data = outcome_dist,
                 geom = "bar", stat = "identity", xlab = "Outcome Vector", ylab = "Frequency",
                 main = this_main, log = this_log)
    return(p)
  }
}


expedite_mode_geom_bar3 <- function(interaction_table, cis, file_type = "png") {
  stopifnot(nrow(interaction_table) < 1000)
  setwd(root_dir)
  signal_combo <- paste(T1, T2, sep = "_")
  this_dir <- paste(root_dir, "signal_combinations", signal_combo, cell_type,
                    "mode_barplots", sep = "/")

  if(!file.exists(file.path(this_dir))) {
    dir.create(file.path(this_dir))
  }

  print('Producing interaction profile barplots in
        ~/mode_barplots/limma_classification. Y-axis scale may vary
        drastically by plot.')
  if (getwd() == root_dir) {
    setwd(this_dir)
  }

  for (i in seq(1, nrow(interaction_table), 12)) {
    interaction_table <- interaction_table[i:(i+11),]
    par(mfrow = c(3, 4))
    pdf('mode_barplots.pdf')

    for (j in 1:12) {
      p <- geom_bar_constructor3(as.numeric(interaction_table[j,9:12]),
                                 cis[[1]][j,1:4],
                                 cis[[2]][j,1:4],
                                 as.character(interaction_table[j,7]),
                                 interaction_table[j,2])
      print(p)
      dev.off()
    }
  }
}


geom_bar_constructor <- function(log.fc, ci.l, ci.r, int.mode, gene) {
  #stopifnot(is.numeric(log.fc), is.character(int.mode), is.character(gene))
  gg_plot <- qplot(x = COND, y = as.numeric(log.fc),
                   geom = "bar", stat = "identity", xlab = "Condition",
                   ylab = "logFC", main = paste0(gene,
                                                 " Interaction Profile - ",
                                                 int.mode)) +
    geom_errorbar(aes(y = log.fc, ymin = ci.l, ymax = ci.r))
  colorBarplot(findClass(int.mode), gg_plot)
}


