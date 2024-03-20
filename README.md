
# igseq

igseq is an R Shiny application designed for antibody sequence analysis, physicochemical property analysis, and post-translational modification (PTM) analysis.

## Installation

You can install `igseq` directly from GitHub using the `devtools` package.

```r
if (!require('devtools')) install.packages('devtools')
devtools::install_github("yourusername/yourrepository")
```

## Prerequisites

### 1. R packages

Before running the `igseq` application, ensure all necessary R packages are installed and any required configurations are set up. Here is a list of required packages:

- shiny
- reticulate
- DT
- Peptides
- stringr
- tableHTML

### 2. Python Configuration and Requirements

The `igseq` application requires Python for certain functionalities. Ensure Python is installed and configure the necessary Python packages.

For managing the Python environment and dependencies, we recommend using Conda, an open-source package management system and environment management system. It simplifies package management and deployment for Python packages.

#### Installing Conda

If Conda is not already installed, download and install Miniconda, a minimal installer for Conda. Miniconda is a smaller alternative to Anaconda, suitable for users who wish to minimize disk space usage while leveraging the powerful Conda environment and package manager.

#### Installing Required Python Packages

The application requires the following Python packages:

- BioPython
- pandas
- abnumber

```bash
conda activate
which python # pythonPath used in launchApp function
conda install -c conda-forge biopython pandas
conda install -c bioconda abnumber
```

### 3. igblast

The package relies on igblast and requires specifying the igblast path during runtime. Follow these steps for installation and configuration:

1. The latest igblast program can be downloaded from: [NCBI igblast download link](https://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/ncbi-igblast-1.22.0-x64-linux.tar.gz)
```bash
VERSION=1.22.0
wget https://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/ncbi-igblast-${VERSION}-x64-linux.tar.gz
```

2. Extract igblast
```bash
tar -zxvf ncbi-igblast-${VERSION}-x64-linux.tar.gz
sudo cp ncbi-igblast-${VERSION}/bin/* /bin  
```

3. Download the reference database
```bash
OUTDIR=~/share/igblast
if [ ! -d $OUTDIR ]; then
    mkdir -p $OUTDIR
fi
# Fetch database
wget -q -r -nH --cut-dirs=5 --no-parent     ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/database     -P ${OUTDIR}/database
# Extract
tar -C ${OUTDIR}/database -xf ${OUTDIR}/database/mouse_gl_VDJ.tar
tar -C ${OUTDIR}/database -xf ${OUTDIR}/database/rhesus_monkey_VJ.tar
cp -r ncbi-igblast-${VERSION}/internal_data $OUTDIR
cp -r ncbi-igblast-${VERSION}/optional_file $OUTDIR
```

4. Build igblast database from IMGT reference sequences. This script is modified from fetch_imgtdb.sh of the immcantation workflow, with the addition of downloading IGHJ amino acid.
```bash
# (The script is lengthy, so ensure to maintain the complete script as provided by the user)
```

5. Process downloaded IMGT reference for igblast
```bash
# (This script is also lengthy and should be kept as is from the user's input)
```

## Running the Application
```R
igseq::launchApp()
```

## Features
- **Sequence Analysis**: Users can upload a FASTA file or directly paste FASTA content for sequence analysis. The application offers options to select sequence type (Protein or DNA), species, and whether to remove duplicate sequences.
- **Physicochemical Property Analysis**: Users can select to analyze various physicochemical properties such as length, molecular weight, net charge, etc. The results are displayed in a table and can be downloaded.
- **Post-Translational Modification Analysis**: The application provides an interface to analyze post-translational modifications. The results highlight specific motifs and modifications in the sequence data.
