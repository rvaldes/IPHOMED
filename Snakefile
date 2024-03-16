
### Integrated Proteomics of HOst-MicrobiomE-Diet (IPHOMED)
### Rafael Valdes-Mas (Elinav Lab)

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
UNIPROT_CRAP_DB = config['UNIPROT_CRAP_DB']

# IPHOMED directory
IPHOMED = config['IPHOMED']

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
		counts_pre = expand(join('preprocess', '{sample}', '{sample}.pre.counts'), sample = SAMPLES),
		counts_post = expand(join('preprocess', '{sample}', '{sample}.counts'), sample = SAMPLES),
		all_diamond_iphomed = expand(join('iphomed_diamond', '{sample}', '{sample}.blastout'), sample = SAMPLES)
		
if PAIRED:
	rule counts_pre:
		input:
			r1 = 'FASTQs/{sample}_R1.fastq.gz'
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

	rule preprocess:
		input:
			r1 = 'FASTQs/{sample}_R1.fastq.gz',
			r2 = 'FASTQs/{sample}_R2.fastq.gz'
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

	rule counts:
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

	rule subsample:
		input:
			r1 = rules.preprocess.output.r1,
			r2 = rules.preprocess.output.r2
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

	rule kraken:
		input:
			r1 = rules.subsample.output.r1,
			r2 = rules.subsample.output.r2
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
			mem = 320000,
			queue = 'elinav'
		shell:
			"kraken2 --db {params.kdb} --threads 1 --output {output.out} --report {output.report} --use-names --paired {input.r1} {input.r2}"

else:
	rule counts_pre:
		input:
			r1 = 'FASTQs/{sample}.fastq.gz'
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

	rule preprocess:
		input:
			r1 = 'FASTQs/{sample}.fastq.gz'
		params:
			bowtieDB = BOW_DB
		output:
			r1 = join('preprocess', '{sample}', '{sample}_R1.fastq.gz')
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

	rule counts:
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


	rule subsample:
		input:
			r1 = rules.preprocess.output.r1
		params:
			sub_depth = SUB_DEPTH
		output:
			r1 = join('subsampling', '{sample}', '{sample}_R1.fastq.gz')
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


	rule kraken:
		input:
				r1 = rules.subsample.output.r1
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
				mem = 320000,
				queue = 'elinav'
		shell:
				"kraken2 --db {params.kdb} --threads 1 --output {output.out} --report {output.report} --use-names {input.r1} "

rule bracken:
	input:
		report = rules.kraken.output.report
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
	params:
		iphomed_dir = IPHOMED
	output:
		species_taxIDs = join('species_95', 'bacterial-species.txt'),
		species_taxIDs_detailed = join('species_95', 'bacterial-species.detailed.txt'),
		genbank_assemblies = join('species_95', 'assembly_summary.genbank.txt'),
		genomes_95 = join('species_95', 'genomes_95.info.tsv'),
		species_abundances = join('species_95', 'bacterial_abundances.95.tsv'),
		assemblies_95 = join('species_95', 'assemblies_95.tsv'),
		assemblies_download = join('species_95', 'download_genomes.sh')
	log:
		join('species_95', 'species95.log')
	benchmark:
		join('species_95', 'species95.benchmark.tsv')
	threads:
		1
	resources:
		mem = 5000,
		queue = 'new-short'
	shell:
		"""
		taxonkit list --ids 2 -I \"\" | taxonkit filter -E species -o {output.species_taxIDs}
		taxonkit lineage -L -nr {output.species_taxIDs} -o {output.species_taxIDs_detailed}
		wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt -O {output.genbank_assemblies}
		sed -i s'/^#//' {output.genbank_assemblies}
		Rscript {params.iphomed_dir}/scripts/species_95.R {output.species_taxIDs_detailed} {output.genbank_assemblies} {output.species_abundances} {output.genomes_95} {output.assemblies_95}
		awk '{{print \"datasets  download genome accession \", $1,\" --filename \", $1, \".zip"}}'  {output.assemblies_95} | sed -e s'/ .zip/.zip/' > {output.assemblies_download}
		"""

rule download_genomes:
	input:
		script = rules.species_95.output.assemblies_download
	output:
		fasta = join('genome_sequences', 'genomes.fasta'),
		assembly_info = join('genome_sequences', 'assembly_info.txt')
	log:
		join('genome_sequences', 'download_genomes.log')
	benchmark:
		join('genome_sequences', 'download_genomes.benchmark.tsv')
	threads:
		1
	resources:
		mem = 1000,
		queue = 'new-short'
	shell:
		"""
		sh {input.script}
		ls *zip | awk '{{print \"7za -y x \"$1}}'  > unzip.sh
		sh unzip.sh
		rm *zip unzip.sh
		cat ncbi_dataset/*/*/*fna > {output.fasta}
		awk '{{print FILENAME\"\t\"$0}}' ncbi_dataset/*/*/*fna | grep \">\" | awk '{{print $1\"\t\"$2}}' | sed -e s'/[\._0-9a-zA-Z\/]*\///'| sed -e s'/>//' | sed -e s'/\_genomic.fna//' > {output.assembly_info}
		rm -rf ncbi_dataset
		"""

if PAIRED:
	rule aligment2genomes:
		input:
			reference = rules.download_genomes.output.fasta,
			assembly_info = rules.download_genomes.output.assembly_info
		params:
			iphomed_dir = IPHOMED
		output:
			bam = join('alignment2genome', 'alignment2genome.bam'),
			depth = join('alignment2genome', 'alignment2genome.depth'),
			summary_depth = join('alignment2genome', 'alignment2genome.summary')
		log:
			join('alignment2genome', 'aligment2genomes.log')
		benchmark:
			join('alignment2genome', 'aligment2genomes.benchmark.tsv')
		threads:
			20
		resources:
			mem = 1000,
			queue = 'new-short'
		shell:
			"""
			minimap2 -t 20 -ax sr {input.reference} <(cat preprocess/*/*R1*gz) <(cat preprocess/*/*R2*gz) | samtools view -Sb - | samtools sort - > {output.bam}
			samtools depth -a {output.bam} > {output.depth}
			python {params.iphomed_dir}/scripts/genome_coverage.py {input.assembly_info} {output.depth} > {output.summary_depth}
			"""

else:
	rule aligment2genomes:
		input:
			reference = rules.download_genomes.output.fasta,
			assembly_info = rules.download_genomes.output.assembly_info
		params:
			iphomed_dir = IPHOMED
		output:
			bam = join('alignment2genome', 'alignment2genome.bam'),
			depth = join('alignment2genome', 'alignment2genome.depth'),
			summary_depth = join('alignment2genome', 'alignment2genome.summary')
		log:
			join('alignment2genome', 'alignment2genomes.log')
		benchmark:
			join('alignment2genome', 'alignment2genomes.benchmark.tsv')
		threads:
			20
		resources:
			mem = 1000,
			queue = 'new-short'
		shell:
			"""
			minimap2 -t 20 -ax sr {input.reference} <(cat preprocess/*/*gz) | samtools view -Sb - | samtools sort - > {output.bam}
			samtools depth -a {output.bam} > {output.depth}
			python {params.iphomed_dir}/scripts/genome_coverage.py {input.assembly_info} {output.depth} > {output.summary_depth}
			"""

rule download_protein_sequences:
	input:
		depth = rules.aligment2genomes.output.summary_depth
	params:
		iphomed_dir = IPHOMED
	output:
		covered_genomes = join('protein_sequences', 'covered_genomes.list'),
		assemblies_download = join('protein_sequences', 'download_proteins.sh'),
		protein_sequences = join('protein_sequences', 'proteins.fasta'),
		diamond_db = join('protein_sequences', 'proteins_diamond.dmnd')
	log:
		join('protein_sequences', 'download_proteins.log')
	benchmark:
		join('protein_sequences', 'download_proteins.benchmark.tsv')
	threads:
		1
	resources:
		mem = 2000,
		queue = 'new-short'
	shell:
		"""
		Rscript {params.iphomed_dir}/scripts/covered_genomes.R {input.depth} {output.covered_genomes}
		awk '{{print \"datasets  download genome accession \", $0,\"--include protein --filename \", $1, \".zip\"}}' {output.covered_genomes} | sed -e s'/ .zip/.zip/' > {output.assemblies_download}
		sh {output.assemblies_download}
		ls *zip | awk '{{print \"7za -y x \"$0}}' > unzip.sh
		sh unzip.sh
		rm *zip unzip.sh
		cat ncbi_dataset/data/*/protein.faa > {output.protein_sequences}
		diamond makedb --in {output.protein_sequences} -d protein_sequences/proteins_diamond
		rm -rf ncbi_dataset
		"""

rule aligment2proteinsequence:
	input:
		database = rules.download_protein_sequences.output.diamond_db,
		fasta = rules.download_protein_sequences.output.protein_sequences
	params:
		iphomed_dir = IPHOMED
	output:
		diamond = join('alignment2protein', 'alignment2protein.diamond'),
		lengths = join('alignment2protein', 'proteins.fasta.lengths'),
		summary_depth = join('alignment2protein', 'alignment2protein.summary'),
		summary_depth_filtered = join('alignment2protein', 'alignment2protein.summary.filtered')
	log:
		join('alignment2protein', 'alignment2protein.log')
	benchmark:
		join('alignment2protein', 'alignment2protein.benchmark.tsv')
	threads:
		10
	resources:
		mem = 5000,
		queue = 'new-short'
	shell:
		"""
		diamond blastx -q <(zcat preprocess/*/*_R*gz) -o {output.diamond} --db {input.database}  -e 0.001 --threads 10 -k 1 -c 1
		bioawk  -c fastx '{{ print $name, length($seq) }}' < {input.fasta} > {output.lengths}
		cat {output.diamond} | python {params.iphomed_dir}/scripts/protein_coverage.py {output.lengths} > {output.summary_depth}
		awk '$5>0' {output.summary_depth} | awk '{{print $1}}' > {output.summary_depth_filtered}
		"""

rule iphomed_bacteria_database:
	input:
		fasta = rules.download_protein_sequences.output.protein_sequences,
		summary_depth_filtered = rules.aligment2proteinsequence.output.summary_depth_filtered
	output:
		fasta = join('iphomed_database', 'iphomed.bacteria.fasta'),
		diamond_db = join('iphomed_database', 'iphomed.bacteria.dmnd')
	log:
		join('iphomed_database', 'iphomed_bacteria_database.log')
	benchmark:
		join('iphomed_database', 'iphomed_bacteria_database.benchmark.tsv')
	threads:
		1
	resources:
		mem = 1000,
		queue = 'new-short'
	shell:
		"""
		seqkit grep -f {input.summary_depth_filtered} {input.fasta} > {output.fasta}
		diamond makedb --in {output.fasta} -d iphomed_database/iphomed.bacteria
		"""

if PAIRED:
	rule diamond_iphomed_bacteria:
		input:
			r1 = rules.subsample.output.r1,
			r2 = rules.subsample.output.r2,
			db = rules.iphomed_bacteria_database.output.diamond_db
		output:
			diamond = join('iphomed_bacteria_diamond', '{sample}', '{sample}.blastout')
		log:
			join('iphomed_bacteria_diamond', '{sample}', 'iphomed_bacteria_diamond.log')
		benchmark:
			join('iphomed_bacteria_diamond', '{sample}', 'iphomed_bacteria_diamond.benchmark.tsv')
		threads:
			10
		resources:
			mem = 5000,
			queue = 'new-short'
		shell:
			"""
			zcat {input.r1} {input.r2} | diamond blastx -o {output.diamond} --db {input.db} -e 0.001 --threads 10 -k 1 -c 1
			"""

else:
	rule diamond_iphomed_bacteria:
		input:
			r1 = rules.subsample.output.r1,
			db = rules.iphomed_bacteria_database.output.diamond_db
		output:
			diamond = join('iphomed_bacteria_diamond', '{sample}', '{sample}.blastout')
		log:
			join('iphomed_bacteria_diamond', '{sample}', 'iphomed_bacteria_diamond.log')
		benchmark:
			join('iphomed_bacteria_diamond', '{sample}', 'iphomed_bacteria_diamond.benchmark.tsv')
		threads:
			10
		resources:
			mem = 5000,
			queue = 'new-short'
		shell:
			"""
			diamond blastx -q {input.r1} -o {output.diamond} --db {input.db} -e 0.001 --threads 10 -k 1 -c 1
			"""

rule iphomed_proteomics_calibration:
	input:
		iphomed_host = UNIPROT_DB,
		mzml = 'mzML/ExperimentalDesign.tsv'
	output:
		calibrated = join('iphomed', 'Calibrated', 'Task1CalibrationTask/AutoGeneratedManuscriptProse.txt')
	log:
		join('iphomed', 'Calibrated', 'iphomed_proteomics_calibration.log')
	benchmark:
		join('iphomed', 'Calibrated', 'iphomed_proteomics_calibration.benchmark.tsv')
	threads:
		10
	resources:
		mem = 30000,
		queue = 'new-long'
	shell:
		"""
		dotnet /home/labs/elinav/rvaldes/miniconda3/envs/iphomed/lib/dotnet/tools/metamorpheus/CMD.dll -d {input.iphomed_host} -s iphomed/Calibrated/Task1Task1SearchTask/ -t {params.iphomed_dir}/TaskCalibrationTaskconfig.toml -o Search
		"""


rule iphomed_proteomics_search:
	input:
		iphomed_bacteria = rules.iphomed_bacteria_database.output.fasta,
		iphomed_diet = '{params.iphomed_dir}/iphomed.diet.fasta',
		iphomed_host = UNIPROT_DB,
		iphomed_crap = UNIPROT_CRAP_DB,
		calibrated = rules.iphomed_proteomics_calibration.output.calibrated
	output:
		search = join('iphomed', 'Search', 'Task1SearchTask/AutoGeneratedManuscriptProse.txt')
	log:
		join('iphomed', 'Search', 'iphomed_proteomics_search.log')
	benchmark:
		join('iphomed', 'Search', 'iphomed_proteomics_search.benchmark.tsv')
	threads:
		10
	resources:
		mem = 30000,
		queue = 'new-long'
	shell:
		"""
		dotnet /home/labs/elinav/rvaldes/miniconda3/envs/iphomed/lib/dotnet/tools/metamorpheus/CMD.dll -d {input.iphomed_host} {input.iphomed_bacteria} {input.iphomed_diet} {input.iphomed_crap} -s iphomed/Calibrated/Task1Task1SearchTask/ -t {params.iphomed_dir}/TaskSearchTaskconfig.toml -o iphomed/Search
		"""