# SciX Genies

This is a Github repository for the UNSW SciX program 'Exploring the human genome' for learning Bioinformatics.

This tutorial will follow these steps:

1. Install R and Rstudio
2. Download the files from this github repository
3. Move the SciX-main.zip to Desktop and unzip the file
3. Open the .Rmd file in RStudio
4. Perform analysis by following instructions in the Rmarkdown file

# Overview
Transcriptomics is a field of biology interested in the structure and function of ribonucleic acid (RNA). RNA is genetic information used by our cells to instruct and regulate protein production. Proteins catalyze chemical reactions as enzymes, function as signalling molecules like neurotransmitters or cytokines, and provide the structures which make up cells, tissues and organs. 

By sequencing and measuring the abundance of RNA, scientists are able to infer cellular activity and the proteins being produced within tissues and organs. By comparing RNA expression levels between individuals with disease and healthy controls we can then identify the genes contributing to a disease and provide evidence for diagnostic biomarkers and drug targets. 

In this study, RNA sequencing data from the brains of individuals with neurodegenerative disease and healthy controls was checked for quality and pre-processed using peer reviewed and commonly used bioinformatics tools (see below). The data analysed in this study include:

# Parkinson's disease (PD)
Link to data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106608 \
Link to paper: NA \
Number of patients: 14 \
Number of controls: 18

# Alzheimer's disease (AD)
Link to data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159699 \
Link to original paper: https://www.nature.com/articles/s41588-020-0696-0 \
Number of patients: 12 \
Number of controls: Old:10 Young:8 

# Huntingtons disease (HD)
Link to data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129473 \
Link to original paper: https://bmcmedgenomics.biomedcentral.com/articles/10.1186/s12920-019-0581-9 \
Number of patients: 8 \
Number of controls: 5

# Preprocessing method
Trimmomatic v.0.39: Remove sequencing adapters and low quality reads \
STAR v.2.7.9a: Align sequencing reads to the GRch38 Homo sapiens reference genome \
Stringtie v.2.1.7: Annotate genes against GENCODE v.29 and quantify read counts 

# Install R and Rstudio
We will be using the R programming language and the Rstudio IDE for this analysis.

To install R please follow this link to download and install the correct version for your Mac, Windows or Linux computer.
https://cran.csiro.au/

Rstudio can be downloaded and installed at this link:
https://posit.co/download/rstudio-desktop/

# Download files from github
At the top of this page is a green button **CODE**. Click the button and then select **Download ZIP**.
This will download a SciX-main.zip file to your Downloads folder on your computer. This file contains everything you will need to perform this study. 
Unzip the file by double clicking and move the unzipped directory to your Desktop.

# Open the .Rmd in RStudio
Choose a disease and open the respective folder e.g AD, PD, HD \
There are three files labelled 1 --> 3 from metadata to gene ontology. The full script is found in 3.gene_ontology, however we will be learning the theory for metadata, differential expression and gene ontology seperately. \
Double click on the SciX.Genies.Rmd file to open in RStudio. Carefully follow the instructions.\
Enjoy!


Lachlan Gray \
Edited by: Helen King \
Edited by: Ben Hanrahan \
Edited by: Andrew Li







