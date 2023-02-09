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
BEDTOOLS = config["conda"][OS]["bedtools"]          # BedTools
BCFTOOLS = config["conda"][OS]["bcftools"]          # BcfTools
GAWK = config["conda"][OS]["gawk"]                  # Gawk
GATK = config["conda"][OS]["gatk"]                  # GATK 3.6
GATK4 = config["conda"][OS]["gatk4"]                # GATK 4.3.0
PICARD = config["conda"][OS]["picard"]              # Picard 2.27.5
QUALIMAP = config["conda"][OS]["qualimap"]          # Qualimap 2.2.0

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


################################## A L L #####################################
rule all:
    input:
        index_post_realign = expand("results/04_Polishing/realigned/{sample}_bwa_md_realigned_fixed.bai", sample=SAMPLE),
        fixmateinformation = expand("results/04_Polishing/realigned/{sample}_bwa_md_realigned_fixed.bam", sample=SAMPLE),
        #qualimap = expand("results/00_Quality_Control/qualimap/{sample}_bwa/qualimapReport.html", sample=SAMPLE),

################################ R U L E S ####################################
rule fixmateinformation:
    # Aim: This tool ensures that all mate-pair information is in sync between each read and its mate pair.
    #      If no #OUTPUT file is supplied then the output is written to a temporary file and then copied over 
    #      the #INPUT file (with the original placed in a .old file.)
    # Use: picard.jar FixMateInformation \
    #      -I input.bam \
    #      -O fixed_mate.bam \
    #      --ADD_MATE_CIGAR true
    message:
        "Picard FixMateInformation"
    conda:
        PICARD
    input:
        realigned = "results/04_Polishing/realigned/{sample}_bwa_md_realigned.bam",
    output:
        fixed = "results/04_Polishing/realigned/{sample}_bwa_md_realigned_fixed.bam"
    threads: CPUS
    log:
        "results/11_Reports/fixmateinformation/{sample}_bwa_realigned_fixed.log"
    shell:
        """
        picard FixMateInformation -I {input.realigned} -O {output.fixed} --ADD_MATE_CIGAR true &> {log}
        """

###############################################################################
rule samtools_index_post_realign:
    # Aim: indexing marked as duplicate BAM file fir samtools_idx_stats rule
    # Use: samtools index -@ [THREADS] -b [MARKDUP.bam] [INDEX.bai] # -b: Generate BAI-format index for BAM files (default)
    message:
        "SamTools indexing realigned fixed BAM file for Picard ValidateSamFile"
    conda:
        SAMTOOLS
    resources:
       cpus = CPUS
    input:
        fixedbam = "results/04_Polishing/realigned/{sample}_bwa_md_realigned_fixed.bam"
    output:
        index = "results/04_Polishing/realigned/{sample}_bwa_md_realigned_fixed.bai"
    log:
        "results/11_Reports/samtools/{sample}_bwa_realigned_fixed_indexed.log"
    threads: CPUS
    shell:
        """
        samtools index -@ {resources.cpus} -b {input.fixedbam} {output.index} &> {log}
        """
###############################################################################
rule qualimap:
    # Aim: Qualimap is a platform-independent application written in Java and R that provides both a Graphical User Interface (GUI) and a 
    # command-line interface to facilitate the quality control of alignment sequencing data. Shortly, Qualimap:
    #       1. Examines sequencing alignment data according to the features of the mapped reads and their genomic properties
    #       2. Povides an overall view of the data that helps to to the detect biases in the sequencing and/or mapping of the data and eases decision-making for further analysis.
    #
    # Use: qualimap bamqc -bam {input.bam} \ 
    #       -c \                        Paint chromosome limits inside charts
    #       -nt {threads} \             
    #       -outdir {output.report} \   Output directory for HTML report (default value is report.html)
    #       -outformat PDF \            Format of the ouput report (PDF or HTML, default is HTML)
    #       -sd                         Activate this option to skip duplicate alignments from the analysis                   
    conda:
        QUALIMAP
    input:
        bam = rules.fixmateinformation.output.fixed,                                 #"results/04_Polishing/realigned/{sample}_bwa_md_realigned_fixed.bam",
        index = rules.samtools_index_post_realign.output.index                      # not used in the command, but it's here so snakemake knows to run the rule after the indexing
    output:
        protected("results/00_Quality_Control/qualimap/{sample}_bwa/qualimapReport.html"),
        protected("results/00_Quality_Control/qualimap/{sample}_bwa/raw_data_qualimapReport/genome_fraction_coverage.txt"),
        protected("results/00_Quality_Control/qualimap/{sample}_bwa/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt"),
        protected("results/00_Quality_Control/qualimap/{sample}_bwa/genome_results.txt"),
        protected("results/00_Quality_Control/qualimap/{sample}_bwa/raw_data_qualimapReport/coverage_histogram.txt")
    params:
        outdir = "results/00_Quality_Control/qualimap/{sample}/"
    threads: CPUS
    resources: 
        mem_mb = 8000,
    log:
        "results/11_Reports/qualimap/logs/{sample}_bwa_qualimap.log" 
    shell:
        """
        qualimap bamqc -bam {input.bam} -nt {threads} --java-mem-size={resources.mem_mb}G -outdir {params.outdir} 
        """
