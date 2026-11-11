#### SpaMTP Lipid Simplifier Functions #################################################################################################################################################################################

#' Uses common lipid nomenclature to simplify lipid annotations
#'
#' Based on rgoslin, this function has the ability to simplify lipid annotations
#' into more general terms, groups and classes.
#' This function can input a data frame containing a column storing the
#' annotation and refine only the entries which are annotated as known lipid
#' structures.
#'
#' @param data data.frame containing a column with lipid annotations
#' @param annotationColumn Character string matching the column name containing
#'   the lipid annotations in the provided df (default = "annotations").
#' @param database Character string defining the database matching the lipid
#'   annotations. Possible entries include
#'   c('Shorthand2020','Goslin','FattyAcids','LipidMaps','SwissLipids','HMDB')
#'   (default = "HMDB").
#' @param verbose Logical; show a progress bar (default = TRUE).
#' @param lipidInfo Character string defining what level of information to return for each lipid. Options are either "all" or "simple". "Simple" returns a smaller, simplified list of annotations for each lipid. Please visit https://bioconductor.org/packages/release/bioc/vignettes/rgoslin/inst/doc/introduction.html to see all possible outputs (default = "simple")).
#'
#' @return Data.frame containing additional columns with simplified lipid names
#' @export
#'
#' @examples
#' lipids <- data.frame(annotations = c("PC 16:0/18:1", "PE 18:0/20:4"))
#' refineLipids(lipids, database = "Shorthand2020", verbose = FALSE)
refineLipids <- function(
    data, annotationColumn = "annotations", database = "HMDB",
    lipidInfo = "simple", verbose = TRUE
) {
    if (!(lipidInfo %in% c("simple", "all"))) {
        warning("Invalid input for 'lipidInfo'!... Must be either 'all' or 'simple'. Will use default == 'simple'")
        lipidInfo <- "simple"
    }
    simple_columns <- c(
        "Lipid.Maps.Category", "Lipid.Maps.Main.Class",
        "Species.Name", "Species.Name.Simple"
    )
    pb <- .progressBar(nrow(data), verbose)
    rows_list <- lapply(seq_len(nrow(data)), function(idx) {
        annotations <- strsplit(data[[annotationColumn]][idx], "; ")[[1]]
        lipid.df <- suppressMessages(rgoslin::parseLipidNames(annotations, grammar = database))
        add_information <- if (lipidInfo == "all") {
            c(colnames(lipid.df), "Species.Name.Simple")
        } else {
            simple_columns
        }
        .progressTick(pb, idx)
        if (identical(unique(lipid.df$Grammar), "NOT_PARSEABLE")) {
            return(as.data.frame(matrix(NA,
                ncol = length(add_information), nrow = 1,
                dimnames = list(NULL, add_information)
            )))
        }
        lipid.df$Species.Name.Simple <- ifelse(
            is.na(lipid.df$Lipid.Maps.Main.Class) | is.na(lipid.df$Total.C) | is.na(lipid.df$Total.DB),
            NA,
            paste0(
                lipid.df$Lipid.Maps.Main.Class,
                "(",
                lipid.df$Total.C,
                ":",
                lipid.df$Total.DB,
                ")"
            )
        )
        lipid.df <- lipid.df[add_information]
        data.frame(lapply(
            lipid.df,
            function(col) paste(unique(col), collapse = "; ")
        ))
    })
    .progressDone(pb)

    final_df <- dplyr::bind_rows(rows_list)
    add_information <- colnames(final_df)
    data[add_information] <- lapply(add_information, function(col) {
        .tidyLipidField(final_df[[col]], species = col == "Species.Name.Simple")
    })
    data
}

# Remove empty and NA entries left by collapsing parsed names with "; ".
.tidyLipidField <- function(x, species = FALSE) {
    x <- stringr::str_replace_all(x, "\\bNA\\b", "")
    x <- gsub("(^|\\W)\\;\\s+|\\s+\\;($|\\W)", "\\1\\2", x)
    x <- gsub("\\;\\s+$", "", x)
    if (species) {
        x <- gsub("[; ]", "", x)
        x <- gsub("\\)", "); ", x)
        x <- sub(";\\s*$", "", x)
    }
    stringr::str_trim(x)
}
