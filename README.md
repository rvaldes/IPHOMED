

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="images/logo-DM.png" width=80%>
  <img src="images/logo-LM.png" width=80%>
</picture>

# Installation
Clone the repository and enter it:
```text
git clone https://github.com/rvaldes/IPHOMED.git
cd IPHOMED
```

Create a conda environment with the included conda environment file.

```text
conda env create -n iphomed -f iphomed.yml
```

Modify the *config.yml* file to include all the databases required for the analysis, including:
* *Bowtie* databases (Human and mouse databases)
  * Pre-built databases are available at https://benlangmead.github.io/aws-indexes/bowtie
* *Kraken/Bracken* microbiome database
  * Pre-built databases are available at https://ccb.jhu.edu/software/kraken2/index.shtml?t=downloads
  * It is *highly recommended* to include human and mouse genomes in the database.
* Host proteomes (Human and mouse fasta files)
  * Reference proteomes are available at https://www.uniprot.org/proteomes/
* cRAP protein sequences (fasta format)
  * Available at https://www.thegpm.org/crap/
* Dietary proteome (refined IPHOMED database)
  * Available at https://drive.google.com/file/d/1fL0PpcpqG3ZobJfY-0om2o9_8IJ8RNIH/view?usp=drive_link
* **IPHOMED** installation folder
* *Metamorpheus* (CMD.dll) path.
  * Locate it in your **IPHOMED** conda environment folder. 

Memory, threads and queues of all rules can be modified in the Snakefile to optimize the resources of your server.

# Running IPHOMED

Copy the *Snakefile* and *config.yml* in your working directory.

Modify the *config.yml* file according to the details of your analysis:

* PAIRED: True/False. Single/paired-end sequencing configuration.
* HOST: HUMAN/MOUSE. Define the host in your experiment.
* SUB_DEPTH: subsampling depth for shotgun metagenomics analysis

The working directory must have the following structure:
```text
|- Snakefile
|- config.yml
|- samples.txt
|- FASTQs/
   |- *_R[12]_001.fastq.gz
|- mzML/
   |- *mzML
   |- ExperimentalDesign.tsv
```

Create the *samples.txt*, including the sample names "[Sample]_R1_001.fastq.gz" of all FASTQ files in FASTQs folder (shotgun metagenomics data)


|     |
| --- |
| SampleA |
| SampleB |
| ... |


Create the *ExperimentalDesign.tsv* (in the mzML folder), including the information on the proteomics data. [Example](https://github.com/smith-chem-wisc/MetaMorpheus/wiki/Experimental-Design):

| FileName    | Condition | Biorep | Fraction | Techrep
| ----------- | ----------- | ----------- | ----------- | ----------- |
| SampleA.mzML | Group_A | 1 | 1 | 1
| SampleB.mzML | Group_A | 2 | 1 | 1
| SampleC.mzML | Group_B | 1 | 1 | 1
| SampleD.mzML | Group_B | 2 | 1 | 1
| ... | ... | ... | ... | ...

Activate the **IPHOMED** conda environment:

```text
conda activate iphomed
```

Finally, run the command (adjust the parameters according to your cluster):
```text
snakemake
```
The output *iphomed* folder will contain the following files:
* *host.proteins.tsv*: host proteins detected by **IPHOMED**.
* *bacteria.proteins.tsv*: bacterial proteins detected by **IPHOMED**.
* *diet.proteins.tsv*: dietary proteins detected by **IPHOMED** without quality filtering.
* *diet.proteins.filtered.tsv*: filtered dietary proteins detected after **IPHOMED** refinement.
* *diet.proteins.filtered.bySample.tsv*: filtered dietary proteins detected in each sample independently.

