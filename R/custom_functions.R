ExtinctionCoefficients <- function(aa_sequence){
  ## Split sequence into pieces
  aa_character <- strsplit(toupper(aa_sequence),"",fixed=T)[[1]]

  ## Calculate Trp(W), Tyr(Y) and Cystine(C) counts
  Numb_Trp <- length(which(aa_character == "W"))
  Numb_Tyr <- length(which(aa_character == "Y"))
  Numb_Cystine <- floor(length(which(aa_character == "C"))/2)

  ## Define Extinction Coefficients of Tyr, Trp and Cystine
  Ext_Tyr <- 1490
  Ext_Trp <- 5500
  Ext_Cystine <- 125

  ## Compute the molecular weight
  mw_aa <- Peptides::mw(aa_sequence)


  ## Calculate Extinction Coefficients of sequence given all of Cystine are paired
  E_Prot_paired <- Numb_Tyr * Ext_Tyr + Numb_Trp * Ext_Trp + Numb_Cystine * Ext_Cystine
  Absorb_paired <- round(E_Prot_paired/mw_aa, digits = 3)


  ## Calculate Extinction Coefficients of sequence given all of Cystine are unpaired
  E_Prot_unpaired <- Numb_Tyr * Ext_Tyr + Numb_Trp * Ext_Trp
  Absorb_unpaired <- round(E_Prot_unpaired/mw_aa, digits = 3)

  results <- c(E_Prot_paired, Absorb_paired, E_Prot_unpaired, Absorb_unpaired)
  names(results) <- c("Ext_Coef_paired", "Absorb_paired", "Ext_Coef_unpaired", "Absorb_unpaired")
  return(results)

}

### detect upload sequence type
detect_sequence_type <- function(file_path) {
  sequences <- readLines(file_path)
  seq_lines <- grep("^>", sequences, invert = TRUE, value = TRUE)
  seq_sample <- paste(seq_lines, collapse = "")
  if (grepl("[^ATCGNatcgn]", seq_sample, perl = TRUE)) {
    return("Protein")
  } else {
    return("DNA")
  }
}
