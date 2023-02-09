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

SAMTOOLS = config["conda"][OS]["samtools"]          # SamTools
PICARD = config["conda"][OS]["picard"]              # Picard 2.27.5

###############################################################################
# PARAMETERS #

CONFIG = config["fastq-screen"]["config"]           # Fastq-screen --conf
MAPPER = config["fastq-screen"]["aligner"]          # Fastq-screen --aligner
SUBSET = config["fastq-screen"]["subset"]           # Fastq-screen --subset

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
        multiqc = "results/00_Quality_Control/MULTIQC/multiqc_report.html",
        check = expand("results/00_Quality_Control/validatesamfile/{sample}_bwa_md_realigned_fixed_ValidateSam.txt", sample=SAMPLE),
        flagstat = expand("results/00_Quality_Control/realigned/{sample}_bwa_md_realigned_fixed_bam.flagstat.txt", sample=SAMPLE),
        idxstats = expand("results/00_Quality_Control/realigned/{sample}_bwa_md_realigned_fixed.idxstats.txt", sample=SAMPLE),
        stats = expand("results/00_Quality_Control/realigned/{sample}_bwa_md_realigned_fixed_stats.txt", sample=SAMPLE)

################################ R U L E S ####################################
rule validate_sam:
    # Aim: Basic check for bam file validity, as interpreted by the Broad Institute.
    # Use: picard.jar ValidateSamFile \
    #      -I input.bam \
    #      -MODE SUMMARY
    message:
        "Picard ValidateSamFile"
    conda:
        PICARD
    input:
        bam = "results/04_Polishing/realigned/{sample}_bwa_md_realigned_fixed.bam",
    output:
        check = "results/00_Quality_Control/validatesamfile/{sample}_bwa_md_realigned_fixed_ValidateSam.txt"
    threads: CPUS
    log:
        "results/11_Reports/validatesamfiles/{sample}_bwa_realigned_fixed_validate_bam.log"
    shell:
        """
        picard ValidateSamFile -I {input.bam} -O {output.check} -M SUMMARY > {log} 2>&1 || true
        """

###############################################################################
rule samtools_stats:
    # Aim: Collects statistics from BAM files
    # Use: samtools stats -r ref.fa input.bam # -r: Reference sequence (required for GC-depth and mismatches-per-cycle calculation).
    message:
        "SamTools stats"
    conda:
        SAMTOOLS
    resources:
       cpus = CPUS
    input:
        bam = "results/04_Polishing/realigned/{sample}_bwa_md_realigned_fixed.bam",
        ref = REFPATH+REFERENCE,
    output:
        stats = "results/00_Quality_Control/realigned/{sample}_bwa_md_realigned_fixed_stats.txt"
    log:
        "results/11_Reports/samtools/{sample}_bwa_realigned_fixed_stats.log"
    shell:
        """
        samtools stats --threads {resources.cpus} -r {input.ref} {input.bam} 1> {output.stats} 2> {log}
        """

###############################################################################
rule samtools_idxstats:
    # Aim: samtools idxstats – reports alignment summary statistics
    #       Retrieve and print stats in the index file corresponding to the input file. Before calling idxstats, the input BAM file should be indexed by samtools index.
    #       If run on a SAM or CRAM file or an unindexed BAM file, this command will still produce the same summary statistics, but does so by reading through the entire file. 
    #       This is far slower than using the BAM indices.
    #       The output is TAB-delimited with each line consisting of reference sequence name, sequence length, # mapped read-segments and # unmapped read-segments. 
    #       It is written to stdout. Note this may count reads multiple times if they are mapped more than once or in multiple fragments.
    input:
        bam="results/04_Polishing/realigned/{sample}_bwa_md_realigned_fixed.bam",
        idx="results/04_Polishing/realigned/{sample}_bwa_md_realigned_fixed.bai",
    output:
        "results/00_Quality_Control/realigned/{sample}_bwa_md_realigned_fixed.idxstats.txt",
    log:
        "results/11_Reports/samtools/idxstats/{sample}_bwa_realigned_fixed_idxstats.log",
    params:
        extra="",  # optional params string
    wrapper:
         "v1.23.1/bio/samtools/idxstats"

###############################################################################
rule samtools_flagstat:
    conda:
        SAMTOOLS
    input:
        bam="results/04_Polishing/realigned/{sample}_bwa_md_realigned_fixed.bam",
    output:
        flagstat="results/00_Quality_Control/realigned/{sample}_bwa_md_realigned_fixed_bam.flagstat.txt",
    log:
        "results/11_Reports/samtools/flagstat/{sample}_bwa_realigned_fixed_bam.log",
    params:
        extra="",  # optional params string
    shell:
        "samtools flagstat {input.bam} > {output.flagstat} &> {log}"
    #wrapper:
    #    "v1.23.1/bio/samtools/flagstat"

###############################################################################
rule multiqc:
    # Aim: aggregates bioinformatics analyses results into a single report
    # Use: multiqc [OPTIONS] --output [MULTIQC/] [FASTQC/] [FASTQSCREEN/] [qulimaps/]
    priority: 42
    message:
        "MultiQC reports aggregating"
    input:
        expand("results/00_Quality_Control/realigned/{sample}_bwa_md_realigned_fixed_bam.flagstat.txt", sample=SAMPLE),
        expand("results/00_Quality_Control/realigned/{sample}_bwa_md_realigned_fixed.idxstats.txt", sample=SAMPLE),
        expand("results/00_Quality_Control/realigned/{sample}_bwa_md_realigned_fixed_stats.txt", sample=SAMPLE),
        expand("results/02_Mapping/{sample}_bwa_sorted-mark-dup_metrics.txt", sample=SAMPLE),
        "results/00_Quality_Control/fastqc/",
        "results/00_Quality_Control/fastq-screen/",
    output:
        "results/00_Quality_Control/MULTIQC/multiqc_report.html"
    params:
        extra=""  # Optional: extra parameters for multiqc.
    log:
        "results/11_Reports/multiqc/multiqc.log"
    wrapper:
        "v1.23.1/bio/multiqc"





















