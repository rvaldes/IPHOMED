
### Integrated Proteomics of HOst-MicrobiomE-Diet (IPHOMED)
### Rafael Valdés-Mas (Elinav Lab)

import json
from os.path import join, basename, dirname
from os import getcwd
from subprocess import check_output
from glob import glob

# Globals ------------------------------------------------

configfile: 'config.yml'

# params
PAIRED = config['PAIRED']
HOST = config['HOST']
SUB_DEPTH = config['SUB_DEPTH']

# databases
KRAKEN_DB = config['KRAKEN_DB']
BRACKEN_DB = config['BRACKEN_DB']
BOW_MOUSE_DB = config['BOW_M_DB']
BOW_HUMAN_DB = config['BOW_H_DB']
UNIPROT_H_DB = config['UNIPROT_H_DB']
UNIPROT_M_DB = config['UNIPROT_M_DB']

if HOST == 'HUMAN':
	BOW_DB = BOW_HUMAN_DB
	UNIPROT_DB = UNIPROT_H_DB
else:
	BOW_DB = BOW_MOUSE_DB
	UNIPROT_DB = UNIPROT_M_DB

file = open("samples.txt", "r")
SAMPLES = [line.rstrip("\n") for line in file]


# Rules -------------------------------------------------

rule all:
        input:
                all_bracken = join('bracken', 'all.combined'),
                all_diamond = expand(join('diamond', '{sample}', 'blastout'), sample = SAMPLES),
                all_dimaond_iphomed = expand(join('diamond_iphomed', '{sample}', 'final.blastout'), sample = SAMPLES)


if PAIRED:

	rule counts_pre_paired:
        input:
                r1 = 'FASTQs/{sample}_R1.gz'
        output:
                counts = join('preprocess', '{sample}', '{sample}.pre.counts')
        log:
                join('preprocess', '{sample}', 'pre.counts.log')
        benchmark:
                join('preprocess', '{sample}', 'pre.counts.benchmark.tsv')
        threads:
                1
        resources:
                mem = 100,
                queue = 'new-short'
        shell:
                "zcat {input.r1} | wc -l | awk '{{print ($1/4)}}' > {output.counts}"

	rule preprocess_paired:
	        input:
	                r1 = 'FASTQs/{sample}_R1.gz',
	                r2 = 'FASTQs/{sample}_R2.gz'
	        params:
	                bowtieDB = BOW_DB
	        output:
	                r1 = join('preprocess', '{sample}', '{sample}_R1.fastq.gz'),
                    r2 = join('preprocess', '{sample}', '{sample}_R2.fastq.gz'),
                    unpaired = join('preprocess', '{sample}', '{sample}_unpaired.fastq.gz'),
                    singleton = join('preprocess', '{sample}', '{sample}_singleton.fastq.gz')
	        log:
	                join('preprocess', '{sample}', 'preprocess.log')
	        benchmark:
	                join('preprocess', '{sample}', 'preprocess.benchmark.tsv')
	        threads:
	                5
	        resources:
	                mem = 5000,
	                queue = 'new-short'
	        shell:
	                "fastp --in1 {input.r1} --in2 {input.r2} --detect_adapter_for_pe --stdout --length_required 50 --unpaired1 {output.unpaired} --unpaired2 {output.unpaired} | bowtie2 -x {params.bowtieDB} -p 5 --interleaved - --very-fast | samtools fastq -f 4 -N -c 6 -1 {output.r1} -2 {output.r2} -s {output.singleton} -"

	rule counts_paired:
	        input:
	                r1 = rules.preprocess.output.r1
	        output:
	                counts = join('preprocess', '{sample}', '{sample}.counts')
	        log:
	                join('preprocess', '{sample}', 'counts.log')
	        benchmark:
	                join('preprocess', '{sample}', 'counts.benchmark.tsv')
	        threads:
	                1
	        resources:
	                mem = 100,
	                queue = 'new-short'
	        shell:
	                "zcat {input.r1} | wc -l | awk '{{print ($1/4)}}' > {output.counts}"


	rule subsample_paired:
	        input:
	                r1 = rules.preprocess_paired.output.r1,
	                r2 = rules.preprocess_paired.output.r2
	        params:
	        		sub_depth = SUB_DEPTH
	        output:
	                r1 = join('subsampling', '{sample}', '{sample}_R1.fastq.gz'),
	                r2 = join('subsampling', '{sample}', '{sample}_R2.fastq.gz')
	        log:
	                join('subsampling', '{sample}', 'subsampling.log')
	        benchmark:
	                join('subsampling', '{sample}', 'subsampling.benchmark.tsv')
	        threads:
	                5
	        resources:
	                mem = 5000,
	                queue = 'new-short'
	        shell:
	                """
	                seqtk sample -s 100 {input.r1} {params.sub_depth} | pigz -p 5 > {output.r1}
	                seqtk sample -s 100 {input.r2} {params.sub_depth} | pigz -p 5 > {output.r2}
	                """

	rule kraken_paired:
	        input:
	                r1 = rules.subsample_paired.output.r1,
	                r2 = rules.subsample_paired.output.r2
	        params:
	                kdb = KRAKEN_DB
	        output:
	                out = join('kraken2', '{sample}', '{sample}.out'),
	                report = join('kraken2', '{sample}', '{sample}.report')
	        log:
	                join('kraken2', '{sample}', 'kraken2.log')
	        benchmark:
	                join('kraken2', '{sample}', 'kraken.benchmark.tsv')
	        threads:
	                1
	        resources:
	                mem = 400000,
	                queue = 'new-short'
	        shell:
	                "kraken --db {params.kdb} --threads 1 --output {output.out} --report {output.report} --use-names --paired {input.r1} {input.r2}"

else:

	rule counts_pre_single:
	        input:
	                r1 = 'FASTQs/{sample}.gz'
	        output:
	                counts = join('preprocess', '{sample}', '{sample}.pre.counts')
	        log:
	                join('preprocess', '{sample}', 'pre.counts.log')
	        benchmark:
	                join('preprocess', '{sample}', 'pre.counts.benchmark.tsv')
	        threads:
	                1
	        resources:
	                mem = 100,
	                queue = 'new-short'
	        shell:
	                "zcat {input.r1} | wc -l | awk '{{print ($1/4)}}' > {output.counts}"

	rule preprocess_single:
	        input:
	                r1 = 'FASTQs/{sample}.gz'
	        params:
	                bowtieDB = BOW_DB
	        output:
	                r1 = join('preprocess', '{sample}', '{sample}.fastq.gz')
	        log:
	                join('preprocess', '{sample}', 'preprocess.log')
	        benchmark:
	                join('preprocess', '{sample}', 'preprocess.benchmark.tsv')
	        threads:
	                5
	        resources:
	                mem = 5000,
	                queue = 'new-short'
	        shell:
	                "fastp --in1 {input.r1}  --stdout --length_required 50  | bowtie2 -x {params.bowtieDB} -p 5 - --very-fast | samtools fastq -f 4 -N -c 6 -0 {output.r1} -"

	rule counts_single:
	        input:
	                r1 = rules.preprocess_single.output.r1
	        output:
	                counts = join('preprocess', '{sample}', '{sample}.counts')
	        log:
	                join('preprocess', '{sample}', 'counts.log')
	        benchmark:
	                join('preprocess', '{sample}', 'counts.benchmark.tsv')
	        threads:
	                1
	        resources:
	                mem = 100,
	                queue = 'new-short'
	        shell:
	                "zcat {input.r1} | wc -l | awk '{{print ($1/4)}}' > {output.counts}"


	rule subsample_single:
	        input:
	                r1 = rules.preprocess_single.output.r1
	        params:
	        		sub_depth = SUB_DEPTH
	        output:
	                r1 = join('subsampling', '{sample}', '{sample}.fastq.gz')
	        log:
	                join('subsampling', '{sample}', 'subsampling.log')
	        benchmark:
	                join('subsampling', '{sample}', 'subsampling.benchmark.tsv')
	        threads:
	                5
	        resources:
	                mem = 5000,
	                queue = 'new-short'
	        shell:
	                "seqtk sample -s 100 {input.r1} {params.sub_depth} | pigz -p 5 > {output.r1}"


	rule kraken_single:
	        input:
	                r1 = rules.subsample_single.output.r1
	        params:
	                kdb = KRAKEN_DB
	        output:
	                out = join('kraken', '{sample}', '{sample}.out'),
	                report = join('kraken', '{sample}', '{sample}.report')
	        log:
	                join('kraken', '{sample}', 'kraken.log')
	        benchmark:
	                join('kraken', '{sample}', 'kraken.benchmark.tsv')
	        threads:
	                1
	        resources:
	                mem = 400000,
	                queue = 'new-short'
	        shell:
	                "kraken2 --db {params.kdb} --threads 1 --output {output.out} --report {output.report} --use-names {input.r1} "


rule bracken:
        input:
                report = rules.kraken2.output.report
        params:
                bdb = BRACKEN_DB
        output:
                out = join('bracken', '{sample}', 'bracken')
        log:
                join('bracken', '{sample}', 'bracken.log')
        benchmark:
                join('bracken', '{sample}', 'bracken.benchmark.tsv')
        threads:
                1
        resources:
                mem = 2000,
                queue = 'new-short'
        shell:
                "bracken -d {params.bdb} -i {input.report} -o {output.out} -r 100 -l S -t 0"


rule species_95:
        input:
                bracken = expand(join('bracken', '{sample}', 'bracken'), sample = SAMPLES)
        output:
                all = join('iphomed', 'species.list')
        log:
                join('iphomed', 'combine.log')
        benchmark:
                join('iphomed', 'combine.benchmark.tsv')
        threads:
                1
        resources:
                mem = 1000,
                queue = 'new-short'
        shell:
                """
                taxonkit list --ids 2 -I ""| taxonkit filter -E species -o bacteria-species.txt
                taxonkit lineage -L -nr bacteria-species.txt -o bacteria-species.detailed.txt
                wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt -O assembly_summary.refseq.txt
                sed -i s'/^#//' assembly_summary.refseq.txt
                wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt -O assembly_summary.genbank.txt
                sed -i s'/^#//' assembly_summary.genbank.txt
                Rscript scripts/species_95.R
                """


if PAIRED:

	rule diamond_initial_paired:
	        input:
	                r1 = rules.preprocess.output.r1,
	                r2 = rules.preprocess.output.r2,
	                db = rules**
	        output:
	                diamond = join('diamond', '{sample}', 'blastout')
	        log:
	                join('diamond', '{sample}', 'diamond.log')
	        benchmark:
	                join('diamond', '{sample}', 'diamond.benchmark.tsv')
	        threads:
	                1
	        resources:
	                mem = 50000,
	                queue = 'new-short'
	        shell:
	                "cat {input.r1} {input.r2} | diamond blastx -o {output.diamond} --db {input.db} -e 0.001 --threads 1 -k 1 -c 1"


else:

	rule diamond_initial_single:
	        input:
	                r1 = rules.preprocess.output.r1,
	                db = rules**
	        output:
	                diamond = join('diamond', '{sample}', 'blastout')
	        log:
	                join('diamond', '{sample}', 'diamond.log')
	        benchmark:
	                join('diamond', '{sample}', 'diamond.benchmark.tsv')
	        threads:
	                1
	        resources:
	                mem = 50000,
	                queue = 'new-short'
	        shell:
	                "diamond blastx -q {input.r1} -o {output.diamond} --db {input.db} -e 0.001 --threads 1 -k 1 -c 1"


rule diamond_iphomed:
        input:
                r1 = rules.subsample_1.output.r1,
                db = rules**
        output:
                diamond = join('diamond_iphomed', '{sample}', '{sample}.blastout')
        log:
                join('diamond_iphomed', '{sample}', 'diamond.log')
        benchmark:
                join('diamond_iphomed', '{sample}', 'diamond.benchmark.tsv')
        threads:
                5
        resources:
                mem = 5000,
                queue = 'elinav'
        shell:
                "diamond blastx -q {input.r1} -o {output.diamond} --db {input.db} -e 0.001 --threads 1 -k 1 -c 1"