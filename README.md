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

Before running the `igseq` application, you need to ensure that all the necessary R packages are installed and any required configurations are set up. Here is a list of required packages:

- shiny
- reticulate
- DT
- Peptides
- stringr
- tableHTML

### 2. Python Configuration and Requirements
The `igseq` application requires Python for certain functionalities. To ensure that these features work correctly, you need to have Python installed and configure the necessary Python packages.

For managing the Python environment and dependencies, we recommend using Conda, an open-source package management system and environment management system. It simplifies package management and deployment for Python packages, especially when dealing with packages that have numerous dependencies.

#### Installing Conda

If you don't already have Conda installed, you can download and install Miniconda, a minimal installer for Conda. Miniconda is a smaller alternative to Anaconda, suitable for users who wish to minimize disk space usage while still leveraging the powerful Conda environment and package manager.


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

该包运行依靠igblast，并在运行时声明igblast的路径，请按照以下步骤进行安装和配置

  1. latest igblast program can be downloaded from: https://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/ncbi-igblast-1.22.0-x64-linux.tar.gz
  ```bash
  VERSION=1.22.0
  wget https://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/ncbi-igblast-${VERSION}-x64-linux.tar.gz
  ```
  
  2. extract igblast
  ```bash
  tar -zxvf ncbi-igblast-${VERSION}-x64-linux.tar.gz
  sudo cp ncbi-igblast-${VERSION}/bin/* /bin  
  ```
  
  3. download reference database
  ```bash
  OUTDIR=~/share/igblast
  if [ ! -d $OUTDIR ]; then
    mkdir -p $OUTDIR
  fi
  # Fetch database
  wget -q -r -nH --cut-dirs=5 --no-parent \
      ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/database \
      -P ${OUTDIR}/database
  # Extract
  tar -C ${OUTDIR}/database -xf ${OUTDIR}/database/mouse_gl_VDJ.tar
  tar -C ${OUTDIR}/database -xf ${OUTDIR}/database/rhesus_monkey_VJ.tar
  cp -r ncbi-igblast-${VERSION}/internal_data $OUTDIR
  cp -r ncbi-igblast-${VERSION}/optional_file $OUTDIR

  ```
  4. build igblast database from IMGT reference sequences. this part of script is edit from fetch_imgtdb.sh of immcantation workflow (https://bitbucket.org/kleinstein/immcantation/src/master/scripts/fetch_imgtdb.sh), just add IGHJ amino acid was also downloaded.
  ```bash
OUTDIR_IMGT=~/share/germlines/imgt
REPERTOIRE="imgt"
if [ ! -d $OUTDIR_IMGT ]; then
  mkdir -p $OUTDIR_IMGT
fi
SPECIES_QUERY=("human:Homo+sapiens"
               "mouse:Mus"
               "rat:Rattus+norvegicus"
               "rabbit:Oryctolagus+cuniculus"
               "rhesus_monkey:Macaca+mulatta")
COUNT=0
for SPECIES in ${SPECIES_QUERY[@]}
do
    KEY=${SPECIES%%:*}
    VALUE=${SPECIES#*:}
    #REPLACE_VALUE=${SPECIES_REPLACE[$COUNT]#*:}
	echo "Downloading ${KEY} repertoires into ${OUTDIR_IMGT}"
	
	# Download VDJ
	echo "|- VDJ regions"
    FILE_PATH="${OUTDIR_IMGT}/${KEY}/vdj"
    FILE_PATH_AA="${OUTDIR_IMGT}/${KEY}/vdj_aa"
    FILE_PATH_LV="${OUTDIR_IMGT}/${KEY}/leader_vexon"
    mkdir -p $FILE_PATH $FILE_PATH_AA $FILE_PATH_LV
    
    # VDJ Ig
    echo "|---- Ig"
    for CHAIN in IGHV IGHD IGHJ IGKV IGKJ IGLV IGLJ
    do
        URL="https://www.imgt.org/genedb/GENElect?query=7.1+${CHAIN}&species=${VALUE}"
        FILE_NAME="${FILE_PATH}/${REPERTOIRE}_${KEY}_${CHAIN}.fasta"
        TMP_FILE="${FILE_NAME}.tmp"
        #echo $URL
        wget $URL -O $TMP_FILE -q
        awk '/<pre>/{i++}/<\/pre>/{j++}{if(j==2){exit}}{if(i==2 && j==1 && $0!~"^<pre>"){print}}' $TMP_FILE > $FILE_NAME
        # Make sed command work also for mac, see: https://stackoverflow.com/a/44864004
        sed -i.bak "$REPLACE_VALUE" $FILE_NAME && rm $FILE_NAME.bak
        rm $TMP_FILE
    done
    
    # Artificial spliced leader and V exon for Ig
    for CHAIN in IGHV IGKV IGLV
    do
        URL="https://www.imgt.org/genedb/GENElect?query=8.1+${CHAIN}&species=${VALUE}&IMGTlabel=L-PART1+V-EXON"
        FILE_NAME="${FILE_PATH_LV}/${REPERTOIRE}_lv_${KEY}_${CHAIN}.fasta"
        TMP_FILE="${FILE_NAME}.tmp"
        #echo $URL
        wget $URL -O $TMP_FILE -q
        awk '/<pre>/{i++}/<\/pre>/{j++}{if(j==2){exit}}{if(i==2 && j==1 && $0!~"^<pre>"){print}}' $TMP_FILE > $FILE_NAME
        # Make sed command work also for mac, see: https://stackoverflow.com/a/44864004
        sed -i.bak "$REPLACE_VALUE" $FILE_NAME && rm $FILE_NAME.bak
        rm $TMP_FILE
    done
        
    # V and J amino acid for Ig
    for CHAIN in IGHV IGHJ IGKV IGLV
    do
        URL="https://www.imgt.org/genedb/GENElect?query=7.3+${CHAIN}&species=${VALUE}"
        FILE_NAME="${FILE_PATH_AA}/${REPERTOIRE}_aa_${KEY}_${CHAIN}.fasta"
        TMP_FILE="${FILE_NAME}.tmp"
        #echo $URL
        wget $URL -O $TMP_FILE -q
        awk '/<pre>/{i++}/<\/pre>/{j++}{if(j==2){exit}}{if(i==2 && j==1 && $0!~"^<pre>"){print}}' $TMP_FILE > $FILE_NAME
        # Make sed command work also for mac, see: https://stackoverflow.com/a/44864004
        sed -i.bak "$REPLACE_VALUE" $FILE_NAME && rm $FILE_NAME.bak
        rm $TMP_FILE
    done
    # Download leaders
    echo "|- Spliced leader regions"
    FILE_PATH="${OUTDIR_IMGT}/${KEY}/leader"
    mkdir -p $FILE_PATH

    # Spliced leader Ig
    echo "|---- Ig"
    for CHAIN in IGH IGK IGL
    do
        URL="https://www.imgt.org/genedb/GENElect?query=8.1+${CHAIN}V&species=${VALUE}&IMGTlabel=L-PART1+L-PART2"
        FILE_NAME="${FILE_PATH}/${REPERTOIRE}_${KEY}_${CHAIN}L.fasta"
        TMP_FILE="${FILE_NAME}.tmp"
        #echo $URL
        wget $URL -O $TMP_FILE -q
        awk '/<pre>/{i++}/<\/pre>/{j++}{if(j==2){exit}}{if(i==2 && j==1 && $0!~"^<pre>"){print}}' $TMP_FILE > $FILE_NAME
        sed -i.bak "$REPLACE_VALUE" $FILE_NAME && rm $FILE_NAME.bak
        rm $TMP_FILE
    done
    # Download constant regions
    echo "|- Spliced constant regions"
    FILE_PATH="${OUTDIR_IMGT}/${KEY}/constant/"
    mkdir -p $FILE_PATH

    # Constant Ig
    echo "|---- Ig"
    for CHAIN in IGHC IGKC IGLC
    do
        # IMGT does not have artificially spliced IGKC / IGLC for multiple species
        if [ "$CHAIN" == "IGHC" ]; then
            QUERY=14.1
        else
            QUERY=7.5
        fi

        URL="https://www.imgt.org/genedb/GENElect?query=${QUERY}+${CHAIN}&species=${VALUE}"
        FILE_NAME="${FILE_PATH}/${REPERTOIRE}_${KEY}_${CHAIN}.fasta"
        TMP_FILE="${FILE_NAME}.tmp"
        #echo $URL
        wget $URL -O $TMP_FILE -q
        awk '/<pre>/{i++}/<\/pre>/{j++}{if(j==2){exit}}{if(i==2 && j==1 && $0!~"^<pre>"){print}}' $TMP_FILE > $FILE_NAME
        sed -i.bak "$REPLACE_VALUE" $FILE_NAME && rm $FILE_NAME.bak
        rm $TMP_FILE
    done
    echo ""
    ((COUNT++))
done
  ```
  
  5. process downloaded IMGT reference for igblast
  ```bash
OUTDIR_FA=${OUTDIR}/fasta
TMPDIR=${OUTDIR}/tmpdir
GERMDIR=$OUTDIR_IMGT
if [ ! -d $OUTDIR_FA ]; then
  mkdir -p $OUTDIR_FA
fi
if [ ! -d $TMPDIR ]; then
  mkdir -p $TMPDIR
fi

for SPECIES_Q in ${SPECIES_QUERY[@]}
do
    SPECIES=${SPECIES_Q%%:*}
    for CHAIN in IG
    do
        # VDJ nucleotides
        for SEGMENT in V D J
        do
	        # VDJ nucleotides
	        F=$(echo imgt_${SPECIES}_${CHAIN}_${SEGMENT}.fasta | tr '[:upper:]' '[:lower:]')
	        cat ${GERMDIR}/${SPECIES}/vdj/imgt_${SPECIES}_${CHAIN}?${SEGMENT}.fasta > ${TMPDIR}/${F}
	    done
		# C nucleotides
	    F=$(echo imgt_${SPECIES}_${CHAIN}_c.fasta | tr '[:upper:]' '[:lower:]')
	    cat ${GERMDIR}/${SPECIES}/constant/imgt_${SPECIES}_${CHAIN}?C.fasta > ${TMPDIR}/${F}
	    # V amino acids
        F=$(echo imgt_aa_${SPECIES}_${CHAIN}_v.fasta | tr '[:upper:]' '[:lower:]')
        cat ${GERMDIR}/${SPECIES}/vdj_aa/imgt_aa_${SPECIES}_${CHAIN}?V.fasta > ${TMPDIR}/${F}
	done
done

# J amino acids
for SPECIES_Q in ${SPECIES_QUERY[@]}
do
	SPECIES=${SPECIES_Q%%:*}
	CHAIN=IG
    F=$(echo imgt_aa_${SPECIES}_${CHAIN}_j.fasta | tr '[:upper:]' '[:lower:]')
    cat ${GERMDIR}/${SPECIES}/vdj_aa/imgt_aa_${SPECIES}_${CHAIN}?J.fasta > ${TMPDIR}/${F}
done


# Parse each created fasta file to create igblast database
cd ${TMPDIR}
NT_FILES=$(ls *.fasta | grep -E "imgt_(human|mouse|rhesus_monkey|rat|rabbit).+\.fasta")
for F in ${NT_FILES}; 
do
	edit_imgt_file.pl ${F} | awk '/^>/{if(seen[$0]++) next} {print}' > ${OUTDIR}/fasta/${F}
	makeblastdb -parse_seqids -dbtype nucl -in ${OUTDIR}/fasta/${F} \
        -out ${OUTDIR}/database/${F%%.*}
done

AA_FILES=$(ls *.fasta | grep -E "imgt_aa_(human|mouse|rhesus_monkey|rat|rabbit).+\.fasta")
for F in ${AA_FILES}; do
	edit_imgt_file.pl ${F} | awk '/^>/{if(seen[$0]++) next} {print}' > ${OUTDIR}/fasta/${F}
	makeblastdb -parse_seqids -dbtype prot -in ${OUTDIR}/fasta/${F} \
        -out ${OUTDIR}/database/${F%%.*}
done
cd -
# Remove temporary fasta files
#cd -; rm -rf $TMPDIR
  ```
  



## Running the Application
```R
igseq::launchApp()
```


## Features
Sequence Analysis: Users can upload a FASTA file or directly paste FASTA content for sequence analysis. The application offers options to select sequence type (Protein or DNA), species, and whether to remove duplicate sequences.
Physicochemical Property Analysis: Users can select to analyze various physicochemical properties such as length, molecular weight, net charge, etc. The results are displayed in a table and can be downloaded.
Post-Translational Modification Analysis: The application provides an interface to analyze post-translational modifications. The results highlight specific motifs and modifications in the sequence data.
