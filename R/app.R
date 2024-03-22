source("R/custom_functions.R")


launchApp <- function(pythonPath = NULL, outDirPath = NULL, host = '127.0.0.1', port = 3839) {
  library(shiny)
  library(reticulate)
  library(DT)
  library(Peptides)
  library(stringr)
  library(tableHTML)

  if (is.null(pythonPath)) {
    pythonPath <- "~/miniconda3/bin/python"  # Please adjust the default paths as necessary
  }
  if (is.null(outDirPath)) {
    outDirPath <- "~/share/igblast"  # Please adjust the default paths as necessary
  }
  reticulate::use_python(pythonPath)

  ui <- shiny::fluidPage(
    titlePanel("Antibody Analysis Tool"),
    tabsetPanel(
      tabPanel("Sequence Analysis",
               #titlePanel("Antibody Sequence and IgBLAST Analysis"),
               fluidRow(
                 column(12,
                        wellPanel(
                          fileInput("fastaFile", "Choose a FASTA File", accept = c(".fasta", ".fa")),
                          hr(),
                          textAreaInput("fastaText", "Or paste FASTA content here:", rows = 3),
                          hr(),
                          fluidRow(
                            column(3,
                                   selectInput("removeDuplicates", "RemoveDuplicate",
                                               choices = list("No" = "no", "Yes" = "yes"),
                                               selected = "yes")),
                            column(3,
                                   selectInput("analysisType", "SequenceType:",
                                               choices = c("Protein" = "igblastp", "DNA" = "igblastn"),
                                               selected = "")),
                            column(3,
                                   selectInput("species", "Species",
                                               choices = c("Mouse" = "mouse", "Human" = "human",
                                                           "Rat" = "rat", "Rhesus monkey" = "rhesus_monkey",
                                                           "Rabbit" = "rabbit"))),
                            column(3,
                                   selectInput("cdrScheme", "CDR regions:",
                                               choices = list("IMGT" = "imgt", "Chothia" = "chothia", "Kabat" = "kabat"),
                                               selected = "imgt"))
                          ),
                          actionButton("analyzeBtn", "Analyze"),
                          downloadButton("downloadData", "Download CSV"),
                          downloadButton("downloadIgBlastResultsV", "Download IgBLAST V Gene Results"),
                          conditionalPanel(
                            condition = "input.analysisType == 'igblastp'",
                            downloadButton("downloadIgBlastResultsJ", "Download IgBLAST J Gene Results")
                          )
                        )
                 ),
                 column(12,
                        dataTableOutput("resultsTable")
                 )
               )
      ),
      tabPanel("Physicochemical Property",
               checkboxGroupInput("propertySelection", "Select properties to analyze:",
                                  choices = list("Length" = "Length",
                                                 "Molecular weight" = "Molecular weight",
                                                 "Net charge" = "Net charge",
                                                 "Hydrophobicity index" = "Hydrophobicity index",
                                                 "Instability index" = "Instability index",
                                                 "Isoelectric point" = "Isoelectric point",
                                                 "Aliphatic index" = "Aliphatic index",
                                                 "Extinction coefficients" = "Extinction coefficients"),
                                  selected = c("Length", "Molecular weight")),
               actionButton("analyzePropertiesBtn", "Analyze"),
               dataTableOutput("propertyResultsTable"),
               downloadButton("downloadMergedResults", "Download Results")

      ),
      tabPanel("PTM",
               actionButton("analyzePTMBtn", "Analyze"),
               uiOutput("ptmResultsTable")
      )


    )
  )


  server <- function(input, output, session) {
    df_seq <- reactiveVal()
    propertiesResult <- reactiveVal()
    analysisResults <- reactiveVal()
    ptmResults <- reactiveVal()
    igBlastResultsFileV <- reactiveVal()
    igBlastResultsFileJ <- reactiveVal()
    detectedType <- reactiveVal()

    observe({

      req(input$fastaFile)
      fasta_path <- input$fastaFile$datapath
      detected_type <- detect_sequence_type(fasta_path)
      detectedType(detected_type)

      if (detected_type != input$analysisType) {
        updateSelectInput(session, "analysisType", selected = detected_type)
      }
      print(detectedType())
    })




    observeEvent(input$analyzeBtn, {
      #req(input$fastaFile)
      req(detectedType())
      scheme <- input$cdrScheme
      cdr_definition <- input$cdrScheme

      sequenceType <- if (!is.null(input$analysisType)) input$analysisType else detectedType()


      # Perform the initial sequence analysis
      withProgress(message = 'Analysis in progress...', value = 0, {
        setProgress(value = 0.1, message = "Preparing analysis...")
        Sys.sleep(0.5)  # Simulating delay

        fastaPath <- if (!is.null(input$fastaFile)) {
          input$fastaFile$datapath
        } else if (input$fastaText != "") {
          tempFile <- tempfile(fileext = ".fasta")
          writeLines(unlist(strsplit(input$fastaText, "\r\n|\r|\n")), tempFile)
          tempFile
        } else {
          stop("No input provided")
        }

        if (!file.exists(fastaPath) || file.size(fastaPath) == 0) {
          stop("The FASTA file is missing or empty.")
        }

        inputData <- if (!is.null(input$fastaFile)) {
          list(type = "file", content = readLines(input$fastaFile$datapath))
        } else if (input$fastaText != "") {
          list(type = "text", content = strsplit(input$fastaText, "\r\n|\r|\n")[[1]])
        } else {
          stop("No input provided")
        }

        setProgress(value = 0.3, message = "Performing IgBLAST analysis...")

        # IgBLAST analysis logic
        OUTDIR <- normalizePath(outDirPath)
        Sys.setenv(IGDATA = OUTDIR)

        SPECIES <- input$species
        SEQTYPE <- "Ig"
        #fastaPath <- input$fastaFile$datapath
        outputFileV <- tempfile(fileext = ".tsv")
        outputFileJ <- tempfile(fileext = ".tsv")


        if (detectedType() == "DNA" | input$analysisType == "igblastn") {
          cmd <- paste(
            "igblastn",
            "-germline_db_V", shQuote(paste0(OUTDIR, "/database/imgt_", SPECIES, "_ig_v")),
            "-germline_db_D", shQuote(paste0(OUTDIR, "/database/imgt_", SPECIES, "_ig_d")),
            "-germline_db_J", shQuote(paste0(OUTDIR, "/database/imgt_", SPECIES, "_ig_j")),
            "-auxiliary_data", shQuote(paste0(OUTDIR, "/optional_file/", SPECIES, "_gl.aux")),
            "-domain_system imgt -ig_seqtype", SEQTYPE, "-organism", SPECIES,
            "-outfmt 19 -query", shQuote(fastaPath),
            "-out", shQuote(outputFileV)
          )

          system(cmd, intern = TRUE)
          igBlastResultsFileV(outputFileV)

          igBlastResults <- read.table(igBlastResultsFileV(), sep = "\t", header = T)
          if (input$removeDuplicates == "yes"){
            igBlastResults <- igBlastResults[which(!duplicated(igBlastResults$sequence_alignment_aa)),]
          }else{
            igBlastResults <- igBlastResults
          }
          df_seq(data.frame(
            seq_id = igBlastResults$sequence_id,
            sequence = igBlastResults$sequence_alignment_aa,
            chain = igBlastResults$locus,
            fwr1 = igBlastResults$fwr1_aa,
            cdr1 = igBlastResults$cdr1_aa,
            fwr2 = igBlastResults$fwr2_aa,
            cdr2 = igBlastResults$cdr2_aa,
            fwr3 = igBlastResults$fwr3_aa,
            cdr3 = igBlastResults$cdr3_aa,
            fwr4 = igBlastResults$fwr4_aa,
            v_call = igBlastResults$v_call,
            d_call = igBlastResults$d_call,
            j_call = igBlastResults$j_call))
          df <- df_seq()[,which(colnames(df_seq()) != "sequence")]
          analysisResults(df)

        } else if (detectedType() == "Protein" | input$analysisType == "igblastp") {  # igblastp
          # For V gene
          cmdV <- paste(
            "igblastp",
            "-germline_db_V", shQuote(paste0(OUTDIR, "/database/imgt_aa_", SPECIES, "_ig_v")),
            "-query", shQuote(fastaPath),
            "-outfmt", "7",
            "-organism", SPECIES,
            "-out", shQuote(outputFileV)
          )

          system(cmdV, intern = TRUE)
          igBlastResultsFileV(outputFileV)
          igBlastResultsV <- read.table(igBlastResultsFileV(), sep = "\t", header = F, fill = T)

          # For J gene
          cmdJ <- paste(
            "igblastp",
            "-germline_db_V", shQuote(paste0(OUTDIR, "/database/imgt_aa_", SPECIES, "_ig_j")),
            "-query", shQuote(fastaPath),
            "-outfmt", "7",
            "-organism", SPECIES,
            "-out", shQuote(outputFileJ)
          )
          system(cmdJ, intern = TRUE)
          igBlastResultsFileJ(outputFileJ)
          igBlastResultsJ <- read.table(igBlastResultsFileJ(), sep = "\t", header = F, fill = T)

          setProgress(value = 0.3, message = "Performing VDJ Region analysis...")
          fastaContent <- inputData$content
          removeDuplicatesFlag <- input$removeDuplicates == "yes"

          py_run_string(paste0("
from Bio import SeqIO
from abnumber import Chain
from io import StringIO
import pandas as pd

def process_fasta(fasta_content, remove_duplicates, scheme, cdr_definition):
    data = []  # Data without sequence
    data_seq = []  # Data with sequence
    observed_sequences = set()
    fasta_io = StringIO('\\n'.join(fasta_content))
    for record in SeqIO.parse(fasta_io, 'fasta'):
        if remove_duplicates and record.seq in observed_sequences:
            continue
        observed_sequences.add(record.seq)
        chain = Chain(str(record.seq), scheme=scheme, cdr_definition=cdr_definition)

        # Append data with sequence
        data_seq.append({
            'seq_id': record.id,
            'chain': chain.chain_type,
            'sequence': str(record.seq),  # Add the sequence here
            'fwr1': chain.fr1_seq,
            'cdr1': chain.cdr1_seq,
            'fwr2': chain.fr2_seq,
            'cdr2': chain.cdr2_seq,
            'fwr3': chain.fr3_seq,
            'cdr3': chain.cdr3_seq,
            'fwr4': chain.fr4_seq
        })

    df_seq = pd.DataFrame(data_seq)
    return df_seq  # Return both dataframes
"))

          # removeDuplicatesFlag to Python function
          py$fasta_content <- fastaContent
          py$remove_duplicates <- removeDuplicatesFlag
          py$scheme <- scheme
          py$cdr_definition <- cdr_definition
          df_seq(py$process_fasta(py$fasta_content, py$remove_duplicates, py$scheme, py$cdr_definition))

          v.df <- igBlastResultsV[match(df_seq()$seq_id, igBlastResultsV$V2),]
          j.df <- igBlastResultsJ[match(df_seq()$seq_id, igBlastResultsJ$V2),]
          current_df <- df_seq()
          current_df$v_call <- v.df$V3
          current_df$j_call <- j.df$V3
          df_seq(current_df)
          df <- df_seq()[,which(colnames(df_seq()) != "sequence")]

          analysisResults(df)

        }
        setProgress(value = 1, message = "Analysis complete!")
      })
    })

    observeEvent(input$analyzePropertiesBtn, {

      req(df_seq())
      print("Analyzing properties...")
      selectedProperties <- input$propertySelection


      # Initialize an empty data frame to store results
      results <- data.frame(seq_id = df_seq()$seq_id)

      # Dynamically add selected properties to the results data frame
      if ("Length" %in% selectedProperties) {
        results$Length <- sapply(df_seq()$sequence, Peptides::lengthpep)
      }
      if ("Molecular weight" %in% selectedProperties) {
        results$`Molecular weight` <- round(sapply(df_seq()$sequence, Peptides::mw), 3)
      }
      if ("Net charge" %in% selectedProperties) {
        results$`Net charge` <- round(sapply(df_seq()$sequence, Peptides::charge), 3)
      }
      if ("Hydrophobicity index" %in% selectedProperties) {
        results$`Hydrophobicity index` <- round(sapply(df_seq()$sequence, Peptides::hydrophobicity), 3)
      }
      if ("Instability index" %in% selectedProperties) {
        results$`Instability index` <- round(sapply(df_seq()$sequence, Peptides::instaIndex), 3)
      }
      if ("Isoelectric point" %in% selectedProperties) {
        results$`Isoelectric point` <- round(sapply(df_seq()$sequence, function(seq) Peptides::pI(seq, pKscale="Bjellqvist")), 3)
      }
      if ("Aliphatic index" %in% selectedProperties) {
        results$`Aliphatic index` <- round(sapply(df_seq()$sequence, Peptides::aIndex), 3)
      }
      if ("Extinction coefficients" %in% selectedProperties) {
        extinction_values <- lapply(df_seq()$sequence, ExtinctionCoefficients)
        extinction_df <- do.call(rbind, extinction_values)
        colnames(extinction_df) <- c("Ext_Coef_paired", "Absorb_paired", "Ext_Coef_unpaired", "Absorb_unpaired")
        results <- cbind(results, extinction_df)
      }

      propertiesResult(results)  # Store the results in a reactive value for rendering
    })

    observeEvent(input$analyzePTMBtn, {
      req(df_seq())
      print("Analyzing properties...")
      data <- data.frame(
        seq_id = df_seq()$seq_id,
        fwr1 = df_seq()$fwr1,
        cdr1 = df_seq()$cdr1,
        fwr2 = df_seq()$fwr2,
        cdr2 = df_seq()$cdr2,
        fwr3 = df_seq()$fwr3,
        cdr3 = df_seq()$cdr3,
        fwr4 = df_seq()$fwr4
      )

      data.re <- apply(data[,-1], 2, function(x) stringr::str_replace_all(x, "RGD", '<span style="color:red; font-weight:bold">RGD</span>'))

      data.re <- apply(data.re, 2, function(x) stringr::str_replace_all(x, "NG", '<span style="background-color:red; font-weight:bold">NG</span>'))
      data.re <- apply(data.re, 2, function(x) stringr::str_replace_all(x, "DG", '<span style="background-color:red; font-weight:bold">DG</span>'))
      data.re <- apply(data.re, 2, function(x) stringr::str_replace_all(x, "NSG", '<span style="background-color:red; font-weight:bold">NSG</span>'))
      data.re <- apply(data.re, 2, function(x) stringr::str_replace_all(x, "N([^P])S", '<span style="background-color:red; font-weight:bold">N\\1S</span>'))
      data.re <- apply(data.re, 2, function(x) stringr::str_replace_all(x, "N([^P])T", '<span style="background-color:red; font-weight:bold">N\\1T</span>'))

      data.re <- apply(data.re, 2, function(x) stringr::str_replace_all(x, "DP", '<span style="background-color:orange; font-weight:bold">DP</span>'))
      data.re <- apply(data.re, 2, function(x) stringr::str_replace_all(x, "DS", '<span style="background-color:orange; font-weight:bold">DS</span>'))
      data.re <- apply(data.re, 2, function(x) stringr::str_replace_all(x, "NS", '<span style="background-color:orange; font-weight:bold">NS</span>'))


      data.re <- apply(data.re, 2, function(x) stringr::str_replace_all(x, "DDD", '<span style="color:blue; font-weight:bold">DDD</span>'))
      data.re <- apply(data.re, 2, function(x) stringr::str_replace_all(x, "EEE", '<span style="color:blue; font-weight:bold">EEE</span>'))
      data.re <- apply(data.re, 2, function(x) stringr::str_replace_all(x, "WWW", '<span style="color:blue; font-weight:bold">WWW</span>'))
      data.re <- apply(data.re, 2, function(x) stringr::str_replace_all(x, "YYY", '<span style="color:blue; font-weight:bold">WWW</span>'))
      data.re <- apply(data.re, 2, function(x) stringr::str_replace_all(x, "FFF", '<span style="color:blue; font-weight:bold">WWW</span>'))
      data.re <- apply(data.re, 2, function(x) stringr::str_replace_all(x, "RRR", '<span style="color:blue; font-weight:bold">WWW</span>'))

      data.re <- apply(data.re, 2, function(x) stringr::str_replace_all(x, "(?<![C])C(?![C])", '<span style="background-color:yellow">C</span>'))

      data.renew <- cbind(data["seq_id"], data.re)

      ptmResults(data.renew)
    })



    # Render the initial analysis results table
    output$resultsTable <- DT::renderDT({
      req(analysisResults())
      datatable(analysisResults(),
                options = list(
                  pageLength = 10,
                  scrollX = TRUE,
                  columnDefs = list(
                    list(targets = which(colnames(analysisResults()) == "cdr1"), className = 'cdr1'),
                    list(targets = which(colnames(analysisResults()) == "cdr2"), className = 'cdr2'),
                    list(targets = which(colnames(analysisResults()) == "cdr3"), className = 'cdr3')
                  ),
                  initComplete = JS(
                    "function(settings, json) {",
                    "$(this.api().table().node()).css('font-family', 'Courier New');",
                    "$(this.api().table().node()).css('border-collapse', 'collapse');",
                    "$(this.api().table().node()).find('th').css('border', '1px solid black');",
                    "$(this.api().table().node()).find('td').css('border', '1px solid black');",
                    "$(this.api().table().node()).find('.cdr1').css({'color': 'green', 'font-weight': 'bold'});",
                    "$(this.api().table().node()).find('.cdr2').css({'color': 'blue', 'font-weight': 'bold'});",
                    "$(this.api().table().node()).find('.cdr3').css({'color': 'red', 'font-weight': 'bold'});",
                    "}"
                  )
                ))
    })



    output$propertyResultsTable <- DT::renderDT({
      req(propertiesResult())
      datatable(propertiesResult(),
                options = list(pageLength = 10, scrollX = TRUE))
    })


    output$ptmResultsTable <- shiny::renderUI({
      req(ptmResults())
      htmlTable <- tableHTML(ptmResults(), escape = FALSE, rownames = FALSE, collapse = 'separate')
      tags$div(style = 'font-family: "Courier New";', HTML(htmlTable))
    })


    # Download handler for the initial analysis results
    output$downloadData <- downloadHandler(
      filename = function() {
        paste("antibody_sequences", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        req(analysisResults())
        write.csv(analysisResults(), file, row.names = FALSE)
      }
    )


    output$downloadIgBlastResultsV <- downloadHandler(
      filename = function() {
        paste0("igblast_results_v_", Sys.Date(), ".txt")
      },
      content = function(file) {
        req(input$analysisType)
        fastaPath <- if (!is.null(input$fastaFile)) {
          input$fastaFile$datapath
        } else if (input$fastaText != "") {
          tempFile <- tempfile(fileext = ".fasta")
          writeLines(unlist(strsplit(input$fastaText, "\r\n|\r|\n")), tempFile)
          tempFile
        } else {
          stop("No input provided")
        }

        outputFileV <- tempfile(fileext = ".txt")
        OUTDIR <- normalizePath(outDirPath)
        Sys.setenv(IGDATA = OUTDIR)
        SPECIES <- input$species
        SEQTYPE <- "Ig"

        if (detectedType() == "DNA" | input$analysisType == "igblastn") {
          cmd <- paste(
            "igblastn",
            "-germline_db_V", shQuote(paste0(OUTDIR, "/database/imgt_", input$species, "_ig_v")),
            "-germline_db_D", shQuote(paste0(OUTDIR, "/database/imgt_", input$species, "_ig_d")),
            "-germline_db_J", shQuote(paste0(OUTDIR, "/database/imgt_", input$species, "_ig_j")),
            "-auxiliary_data", shQuote(paste0(OUTDIR, "/optional_file/", input$species, "_gl.aux")),
            "-domain_system imgt -ig_seqtype Ig",
            "-num_alignments_V 15",
            "-organism", input$species,
            "-outfmt 4 -query", shQuote(fastaPath),
            "-out", shQuote(outputFileV)
          )
        } else if (detectedType() == "Protein" | input$analysisType == "igblastp") {
          cmd <- paste(
            "igblastp",
            "-germline_db_V", shQuote(paste0(OUTDIR, "/database/imgt_aa_", input$species, "_ig_v")),
            "-query", shQuote(fastaPath),
            "-num_alignments_V 15",
            "-outfmt", "4",
            "-organism", input$species,
            "-out", shQuote(outputFileV)
          )
        }

        system(cmd, intern = TRUE)

        file.copy(outputFileV, file)
      }
    )

    # output$downloadIgBlastResultsJ <- downloadHandler(
    #   filename = function() {
    #     paste0("igblast_results_j_", Sys.Date(), ".tsv")
    #   },
    #   content = function(file) {
    #     req(igBlastResultsFileJ())
    #     if (file.exists(igBlastResultsFileJ())) {
    #       file.copy(igBlastResultsFileJ(), file)
    #     } else {
    #       stop("J gene results file does not exist.")
    #     }
    #   }
    # )

    # Download handlers for IgBLAST J gene results - with distinction between DNA and Protein analysis
    output$downloadIgBlastResultsJ <- downloadHandler(
      filename = function() {
        paste0("igblast_results_j_", Sys.Date(), ".txt")
      },
      content = function(file) {
        req(input$analysisType)
        fastaPath <- if (!is.null(input$fastaFile)) {
          input$fastaFile$datapath
        } else if (input$fastaText != "") {
          tempFile <- tempfile(fileext = ".fasta")
          writeLines(unlist(strsplit(input$fastaText, "\r\n|\r|\n")), tempFile)
          tempFile
        } else {
          stop("No input provided")
        }

        outputFileJ <- tempfile(fileext = ".txt")
        OUTDIR <- normalizePath(outDirPath)
        Sys.setenv(IGDATA = OUTDIR)
        SPECIES <- input$species
        SEQTYPE <- "Ig"

        if (detectedType() == "Protein" | input$analysisType == "igblastp") {
          cmd <- paste(
            "igblastp",
            "-germline_db_V", shQuote(paste0(OUTDIR, "/database/imgt_aa_", input$species, "_ig_j")),
            "-query", shQuote(fastaPath),
            "-num_alignments_V 15",
            "-outfmt", "4",
            "-organism", input$species,
            "-out", shQuote(outputFileJ)
          )
        }

        system(cmd, intern = TRUE)

        file.copy(outputFileJ, file)
      }
    )


    output$downloadMergedResults <- downloadHandler(
      filename = function() {
        paste("merged_results", Sys.Date(), ".tsv", sep = "_")
      },
      content = function(file) {
        req(df_seq(), propertiesResult())
        merged_df <- cbind(df_seq(), propertiesResult())
        write.table(merged_df, file, sep = "\t", row.names = FALSE, quote = FALSE)
      }
    )
  }

  shinyApp(ui = ui, server = server, options = list(host = host, port = port))
}
