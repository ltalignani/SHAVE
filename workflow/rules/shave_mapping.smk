#######M######I#####V#####E######G######E######C######I#######R######D#########
# Name:                 shave.smk
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
BOWTIE2 = config["conda"][OS]["bowtie2"]            # Bowtie2
BWA = config["conda"][OS]["bwa"]                    # Bwa

###############################################################################
# PARAMETERS #
TRIMMER = config["trimmer"]                         # Trimmers ('sickle' or 'trimmomatic')
ALIGNER = config["aligner"]                         # Aligners ('bwa' or 'bowtie2')
MARKDUP = config["markdup"]                         # Mark Duplicate Program ('picard' or 'samtools')

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
        mapped = expand("results/02_Mapping/{sample}_bwa-mapped.sam", sample= SAMPLE),

################################ R U L E S ####################################
rule bwa_mapping:
    # Aim: reads mapping against reference sequence
    # Use: bwa mem -t [THREADS] -x [REFERENCE] [FWD_R1.fq] [REV_R2.fq] 1> [MAPPED.sam]
    message:
        "BWA-MEM mapping {wildcards.sample} sample reads against reference genome sequence"
    conda:
        BWA
    resources:
        cpus = CPUS
    params:
        ref = REFPATH+REFERENCE,
        extra = r"'@RG\tID:{sample}\tSM:{sample}\tCN:SC\tPL:ILLUMINA'" # Manage ReadGroup
    input:
        fwdreads = "results/01_Trimmimg/trimmomatic/{sample}_trimmomatic_R1.fastq.gz",
        revreads = "results/01_Trimmimg/trimmomatic/{sample}_trimmomatic_R2.fastq.gz"
    output:
        mapped = "results/02_Mapping/{sample}_bwa-mapped.sam"
    benchmark:
        "benchmarks/bwa/{sample}.tsv"
    log:
        "results/11_Reports/bwa/{sample}.log"
    shell:
        "bwa mem -M -T 0 -t {resources.cpus} -v 1 -R {params.extra} {params.ref} {input.fwdreads} {input.revreads} 1> {output.mapped} 2> {log}"