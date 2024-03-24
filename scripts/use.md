# IPHOMED
Integrated Proteomics of HOst-MicrobiomE-Diet

# Installation
Clone the repository and enter it:
```text
git clone https://github.com/rvaldes/IPHOMED.git
cd IPHOMED
```

Create a miniconda environment with the included conda environment file.

```text
conda env create -n iphomed -f iphomed.yml
```

Modify the *config.yml* file to include all the databases required for the analysis, including:
* *Bowtie* databases (Human and mouse databases)
  * Pre-built databases are available in https://benlangmead.github.io/aws-indexes/bowtie
* *Kraken/Bracken* database
  * Pre-built databases are available in https://ccb.jhu.edu/software/kraken2/index.shtml?t=downloads
* Host Protein fasta files (Human and mouse proteomes, fasta format)
  * Reference proteomes are availabe in https://www.uniprot.org/proteomes/
* cRAP protein sequences (fasta format)
  * Availabe in https://www.thegpm.org/crap/
* **IPHOMED** installation folder
* *Metamorpheus* (CMD.dll) path.
  * Locate it in your **IPHOMED** conda environment folder. 

Memory, threads and queues of all rules can be modified in the Snakefile to optmize the resources of your server.

# Running IPHOMED

Copy the *Snakefile* and *config.yml* in your working directory.

Modify the *config.yml* file according to the details of your analysis:

* *PAIRED*: True/False. Single/paired end squencing configuration.
* *HOST*: HUMAN/MOUSE. Define the host in your experiment.
* *SUB_DEPTH*: subsampling depth for shotgun metagenomics analysis

The working directory must have the following structure:
```text
|- Snakefile
|- config.yml
|- samples.txt
|- FASTQs/
   |- [Sample]_R[12]_001.fastq.gz
|- mzML/
   |- *mzML
   |- ExperimentalDesign.tsv
```

Create the *samples.txt*, including the sample names [Sample] of all FASTQ files in FASTQs folder (shotgun metagenomics data)

Create the ExperimentalDesign.tsv, including the information of the proteomics data.

| FileName    | Condition | Biorep | Fraction | Techrep
| ----------- | ----------- | ----------- | ----------- | ----------- |
| SampleA.mzML | Group_A | 1 | 1 | 1
| SampleB.mzML | Group_A | 2 | 1 | 1
| SampleC.mzML | Group_B | 1 | 1 | 1
| SampleD.mzML | Group_B | 2 | 1 | 1

Activate the **IPHOMED** conda environment:

```text
conda activate iphomed
```

Finally, run the command (Adjust the parameters according to your cluster)
```text
snakemake
```