#######M######I#####V#####E######G######E######C######I#######R######D#########
# Name:                 shave_markdup.smk
# Author:               Loïc TALIGNANI
# Affiliation:          IRD_UMR224_MIVEGEC
# Aim:                  Snakefile for SHort-read Alignment pipeline for VEctor 
# Date:                 2022.10.05
# Run:                  snakemake --snakefile workflow/rules/shave.smk --cores X --use-conda
# Latest modification:  2023.01.16
# Done:                 Added RealignerTargetCreator, IndelRealigner, 
#                       UnifiedGenotyper, compression

###############################################################################
# PUBLICATIONS #

###############################################################################
# CONFIGURATION #
configfile: "config/config.yaml"

from snakemake.utils import min_version
import shutil

min_version("5.18.0")

###############################################################################
# WILDCARDS #
SAMPLE, = glob_wildcards("resources/reads/{sample}_R1.fastq.gz")

###############################################################################
# RESOURCES #
OS = config["os"]                      # Operating system
CPUS = config["resources"]["cpus"]     # Threads
MEM_GB = config["resources"]["mem_gb"] # Memory (RAM) in Gb
TMPDIR = config["resources"]["tmpdir"] # Temporary directory

###############################################################################
# ENVIRONMENTS #
FASTQC = config["conda"][OS]["fastqc"]              # FastQC
FASTQSCREEN = config["conda"][OS]["fastq-screen"]   # Fastq-Screen
MULTIQC = config["conda"][OS]["multiqc"]            # MultiQC 1.14
CUTADAPT = config["conda"][OS]["cutadapt"]          # Cutadapt
SICKLETRIM = config["conda"][OS]["sickle-trim"]     # Sickle-trim
BOWTIE2 = config["conda"][OS]["bowtie2"]            # Bowtie2
BWA = config["conda"][OS]["bwa"]                    # Bwa
SAMTOOLS = config["conda"][OS]["samtools"]          # SamTools
BEDTOOLS = config["conda"][OS]["bedtools"]          # BedTools
BCFTOOLS = config["conda"][OS]["bcftools"]          # BcfTools
GAWK = config["conda"][OS]["gawk"]                  # Gawk
GATK = config["conda"][OS]["gatk"]                  # GATK 3.6
GATK4 = config["conda"][OS]["gatk4"]                # GATK 4.3.0
PICARD = config["conda"][OS]["picard"]              # Picard 2.27.5
QUALIMAP = config["conda"][OS]["qualimap"]          # Qualimap 2.2.0
TRIMMOMATIC = config["conda"][OS]["trimmomatic"]    # Trimmomatic 0.39

###############################################################################
# PARAMETERS #

BWAPATH = config["bwa"]["path"]                     # BWA path to indexes
BT2PATH = config["bowtie2"]["path"]                 # Bowtie2 path to indexes
SENSITIVITY = config["bowtie2"]["sensitivity"]      # Bowtie2 sensitivity preset

REFPATH = config["path"]                            # Path to genomes references
REFERENCE = config["reference"]                     # Genome reference sequence, in fasta format
INDEX = config["index"]                             # Genome reference index, in .fai format
DICTIONARY = config["dictionary"]                  # Genome reference dictionary, made w/ picard CreateSequenceDictionary, in .dict format

ALLELES = config["alleles"]["alleles_target"]       # Alleles against which to genotype (VCF format) 
MINCOV = config["consensus"]["mincov"]              # Minimum coverage, mask lower regions with 'N'
MINAF = config["consensus"]["minaf"]                # Minimum allele frequency allowed
IUPAC = config["consensus"]["iupac"]                # Output variants in the form of IUPAC ambiguity codes

###############################################################################
# FUNCTIONS AND COMMANDS #

# Output Field Separator = <tab>. if the third column is empty, print col 1, col 2 -1 (transform coordinates from 1-based to 0-based for IGV) and col 2.  
# else, print the three column and transform the second one in 0-based coordinates.
# look for "robust quoting for snakemake for more info".
#AWK_CMD_INTERVALS = r"""'BEGIN { OFS = "\t" } { if( $3 == "") { print $1, $2-1, $2 } else { print $1, $2-1, $3}}'"""

############################# O N S U C C E S S ##############################
onsuccess:
    shutil.rmtree(".snakemake")

################################## A L L #####################################
rule all:
    input:
        bam = expand("results/02_Mapping/{sample}_bwa_sorted-mark-dup.bam", sample=SAMPLE),
        metrics=expand("results/02_Mapping/{sample}_bwa_sorted-mark-dup_metrics.txt", sample=SAMPLE),

################################ R U L E S ####################################
rule convert_and_sort:
    message:
        "Samtools view conversion of sample sam in bam format and sorting by coordinates"
    conda:
        SAMTOOLS
    input:
        "results/02_Mapping/{sample}_bwa-mapped.sam",
    output:
        bam = "results/02_Mapping/{sample}_bwa_sorted.bam",
    threads:
        cpus = CPUS
    shell:
        "samtools view --threads {threads.cpus} -bS {input} | "
        "samtools sort -o {output} -T {sample}.temp &> {log} "

###############################################################################
rule mark_duplicates_spark:
    conda:
        GATK4
    input:
        "results/02_Mapping/{sample}_bwa_sorted.bam",
    output:
        bam = temp("results/02_Mapping/{sample}_bwa_sorted-mark-dup.bam"),
        metrics="results/02_Mapping/{sample}_bwa_sorted-mark-dup_metrics.txt",
    benchmark:
        "benchmarks/markduplicatesspark/{sample}_bwa.tsv"
    log:
        "results/11_Reports/markduplicatesspark/{sample}_bwa_sorted-mark-dup.log",
    params:
        extra="--remove-sequencing-duplicates",  # optional
        #java_opts=,  # optional
        #spark_runner="",  # optional, local by default
        #spark_v1.19.1="",  # optional
        #spark_extra="", # optional
    resources:
        mem_mb=16000,
    threads: CPUS
    shell:
        "gatk MarkDuplicatesSpark -I {input} -O {output.bam} -M {output.metrics} {params.extra} > {log} 2>&1"
    #wrapper:
    #    "v1.22.0/bio/gatk/markduplicatesspark"