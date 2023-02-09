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
onsuccess:
    shutil.rmtree(".snakemake")

################################## A L L #####################################
rule all:
    input:
        multiqc = "results/00_Quality_Control/MULTIQC/multiqc_report.html",
        qualimap = expand("results/00_Quality_Control/qualimap/{sample}_{aligner}/qualimapReport.html", sample=SAMPLE, aligner=ALIGNER),
        check = expand("results/00_Quality_Control/validatesamfile/{sample}_{aligner}_md_realigned_fixed_ValidateSam.txt", sample=SAMPLE, aligner=ALIGNER),
        flagstat = expand("results/00_Quality_Control/realigned/{sample}_{aligner}_md_realigned_fixed_bam.flagstat.txt", sample=ALIGNER, aligner=ALIGNER),
        idxstats = expand("results/00_Quality_Control/realigned/{sample}_{aligner}_md_realigned_fixed.idxstats.txt", sample=SAMPLE, aligner=ALIGNER),
        stats = expand("results/00_Quality_Control/realigned/{sample}_{aligner}_md_realigned_fixed_stats.txt", sample=SAMPLE, aligner=ALIGNER),        
        igv_output = expand("results/04_Polishing/{sample}_{aligner}_realignertargetcreator.bed", sample=SAMPLE, aligner=ALIGNER),
        index_post_realign = expand("results/04_Polishing/realigned/{sample}_{aligner}_md_realigned_fixed.bai", sample=SAMPLE, aligner=ALIGNER),
        fixmateinformation = expand("results/04_Polishing/realigned/{sample}_{aligner}_md_realigned_fixed.bam", sample=SAMPLE, aligner=ALIGNER),

################################ R U L E S ####################################
rule SetNmMdAndUqTags:
    # Aim: This tool takes in a coordinate-sorted SAM or BAM and calculates the NM, MD, and UQ tags by comparing with the reference.
    # Use: picard.jar SetNmMdAndUqTags \
    #       R=reference_sequence.fasta
    #       I=sorted.bam \
    #       O=fixed.bam
    message:
        "Picard SetNmMdAndUqTags"
    conda:
        PICARD
    input:
        bam = "results/02_Mapping/{sample}_{aligner}_sorted-mark-dup.bam",
        ref = REFPATH+REFERENCE,                     # "resources/genomes/GCA_018104305.1_AalbF3_genomic.fasta",
    output:
        fix = "results/02_Mapping/{sample}_{aligner}_sorted-mark-dup-fx.bam"
    threads: CPUS
    log:
        "results/11_Reports/SetNmMdAndUqTags/{sample}_{aligner}_sorted-mark-dup-fx.log"
    benchmark:
        "benchmarks/setnmmdanduqtags/{sample}_{aligner}.tsv"
    shell:
        """
        picard SetNmMdAndUqTags R={input.ref} I={input.bam} O={output.fix} > {log} 2>&1 || true
        """

###############################################################################
rule samtools_index_markdup:
    # Aim: indexing marked as duplicate BAM file
    # Use: samtools index -@ [THREADS] -b [MARKDUP.bam] [INDEX.bai]
    message:
        "SamTools indexing marked as duplicate BAM file"
    conda:
        SAMTOOLS
    resources:
       cpus = CPUS
    input:
        markdup = "results/02_Mapping/{sample}_{aligner}_sorted-mark-dup-fx.bam",
    output:
        index = temp("results/02_Mapping/{sample}_{aligner}_sorted-mark-dup-fx.bai"),
    log:
        "results/11_Reports/samtools/{sample}_{aligner}_sorted-mark-dup-index.log"
    shell:
        "samtools index "     # Samtools index, tools for alignments in the SAM format with command to index alignment
        "-@ {resources.cpus} " # --threads: Number of additional threads to use (default: 1)
        "-b "                  # -b: Generate BAI-format index for BAM files (default)
        "{input.markdup} "     # Markdup bam input
        "{output.index} "      # Markdup index output
        "&> {log}"             # Log redirection

###############################################################################
rule realignertargetcreator:
    # Aim:      RealignerTargetCreator identify what regions need to be realigned.
    #           Local realignment around indels. Takes a coordinate-sorted and indexed BAM and a VCF of known indels and creates a target intervals file.
    # Use:      gatk3 -T RealignerTargetCreator \
    #           -R human_g1k_v37_decoy.fasta \
    #           -L 10:96000000-97000000 \
    #           -known INDEL_chr10_1Mb_b37_1000G_phase3_v4_20130502.vcf \
    #           -I 7156_snippet.bam \
    #           -o 7156_realignertargetcreator.intervals
    message:
        "RealignerTargetCreator creates a target intervals file for indel realignment"
    input:
        bam="results/02_Mapping/{sample}_{aligner}_sorted-mark-dup-fx.bam",
        bai="results/02_Mapping/{sample}_{aligner}_sorted-mark-dup-fx.bai",
        ref=REFPATH+REFERENCE,                                                           #"resources/genomes/GCA_018104305.1_AalbF3_genomic.fasta",
        fai="resources/genomes/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa.fai",
        dict="resources/genomes/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.dict",
    output:
        intervals=temp("results/04_Polishing/{sample}_{aligner}.intervals"),
    benchmark:
        "benchmarks/realignertargetcreator/{sample}_{aligner}.tsv"
    log:
        "results/11_Reports/realignertargetcreator/{sample}_{aligner}.log",
    params:
        extra="--defaultBaseQualities 20 --filter_reads_with_N_cigar",  # optional
    resources:
        mem_mb=16000,
    threads: CPUS
    wrapper:
        "v1.22.0/bio/gatk3/realignertargetcreator"

###############################################################################
rule awk_intervals_for_IGV:
    # Aim: View intervals on IGV
    # Use: awk -F '[:-]' 'BEGIN { OFS = "\t" } { if( $3 == "") { print $1, $2-1, $2 } else { print $1, $2-1, $3}}' \ 
    #      realignertargetcreator.intervals > realignertargetcreator.bed
    message:
        "Awk IGV intervals visualization"
    conda:
        GAWK
    input:
        intervals="results/04_Polishing/{sample}_{aligner}.intervals"
    params:
        cmd = r"""'BEGIN { OFS = "\t" } { if( $3 == "") { print $1, $2-1, $2 } else { print $1, $2-1, $3}}'"""
    output:
        bed = "results/04_Polishing/{sample}_{aligner}_realignertargetcreator.bed"
    log:
        "results/11_Reports/awk/{sample}_{aligner}_intervals_for_IGV.log"
    threads: 
        CPUS
    shell:
        "awk -F '[:-]' "                # Awk, a program that you can use to select particular records in a file and perform operations upon them
        "{params.cmd} "                 # {AWK_CMD_INTERVALS:q} :q : is asking snakemake to quote the awk command for me. 
        "{input.intervals} "            # Intervals input
        "1> {output.bed} "              # BedGraph output
        "2> {log} "                     # Log redirection

###############################################################################
rule indelrealigner:
    # Aim:  Mappers cannot “see” indels near ends of reads because mismatches are “cheaper” than a gap in this context.
    #       IndelRealigner takes a coordinate-sorted and indexed BAM and a target intervals file generated by RealignerTargetCreator. 
    #       IndelRealigner then performs local realignment on reads coincident with the target intervals using consenses 
    #       from indels present in the original alignment.
    # Use:  java -Xmx16G -Djava.io.tmpdir=/tmp -jar GenomeAnalysisTK.jar \
    #       -T IndelRealigner \
    #       -R human_g1k_v37_decoy.fasta \
    #       -targetIntervals realignertargetcreator.intervals \
    #       -known INDEL_chr10_1Mb_b37_1000G_phase3_v4_20130502.vcf \ 
    #       -I 7156_snippet.bam \
    #       -o 7156_snippet_indelrealigner.bam     
    message:
        "Indel realignment"
    input:
        bam="results/02_Mapping/{sample}_{aligner}_sorted-mark-dup-fx.bam",
        bai="results/02_Mapping/{sample}_{aligner}_sorted-mark-dup-fx.bai",
        ref=REFPATH+REFERENCE,                                                            #"resources/genomes/GCA_018104305.1_AalbF3_genomic.fasta",
        fai="resources/genomes/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa.fai",
        dict="resources/genomes/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.dict",
        target_intervals="results/04_Polishing/{sample}_{aligner}.intervals"
    output:
        bam= temp("results/04_Polishing/realigned/{sample}_{aligner}_md_realigned.bam"),
        bai= temp("results/04_Polishing/realigned/{sample}_{aligner}_md_realigned.bai"),
    benchmark:
        "benchmarks/indelrealigner/{sample}_{aligner}.tsv",
    log:
        "results/11_Reports/indelrealigner/{sample}_{aligner}.log",
    params:
        extra="--defaultBaseQualities 20 --filter_reads_with_N_cigar",  # optional
    threads: CPUS
    resources:
        mem_mb=16000,
    wrapper:
        "v1.22.0/bio/gatk3/indelrealigner"

###############################################################################
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
        realigned = "results/04_Polishing/realigned/{sample}_{aligner}_md_realigned.bam",
    output:
        fixed = "results/04_Polishing/realigned/{sample}_{aligner}_md_realigned_fixed.bam"
    threads: CPUS
    log:
        "results/11_Reports/fixmateinformation/{sample}_{aligner}_realigned_fixed.log"
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
        fixedbam = "results/04_Polishing/realigned/{sample}_{aligner}_md_realigned_fixed.bam"
    output:
        index = "results/04_Polishing/realigned/{sample}_{aligner}_md_realigned_fixed.bai"
    log:
        "results/11_Reports/samtools/{sample}_{aligner}_realigned_fixed_indexed.log"
    threads: 
        CPUS
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
        bam = "results/04_Polishing/realigned/{sample}_{aligner}_md_realigned_fixed.bam",
    output:
        protected("results/00_Quality_Control/qualimap/{sample}_{aligner}/qualimapReport.html"),
        protected("results/00_Quality_Control/qualimap/{sample}_{aligner}/raw_data_qualimapReport/genome_fraction_coverage.txt"),
        protected("results/00_Quality_Control/qualimap/{sample}_{aligner}/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt"),
        protected("results/00_Quality_Control/qualimap/{sample}_{aligner}/genome_results.txt"),
        protected("results/00_Quality_Control/qualimap/{sample}_{aligner}/raw_data_qualimapReport/coverage_histogram.txt")
    threads: CPUS,
    resources:
        mem_gb=8,
    log:
        stderr="results/11_Reports/qualimap/logs/{sample}_{aligner}_qualimap.stderr",
        stdout="results/11_Reports/qualimap/logs/{sample}_{aligner}_qualimap.stdout" #cp -r results/00_Quality_Control/qualimap/wildcards.sample}_{wildcards.aligner}/wildcards.sample}_{wildcards.aligner}_qualimapReport.html results/00_Quality_Control/qualimap/
    shell:
        """
        qualimap bamqc -bam {input.bam} -c -nt {threads} --java-mem-size={resources.mem_gb}G -outdir results/00_Quality_Control/qualimap/{wildcards.sample}_{wildcards.aligner} -sd > {log.stdout} 2> {log.stderr}
        """

###############################################################################
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
        bam = "results/04_Polishing/realigned/{sample}_{aligner}_md_realigned_fixed.bam",
    output:
        check = "results/00_Quality_Control/validatesamfile/{sample}_{aligner}_md_realigned_fixed_ValidateSam.txt"
    threads: CPUS
    log:
        "results/11_Reports/validatesamfiles/{sample}_{aligner}_realigned_fixed_validate_bam.log"
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
        bam = "results/04_Polishing/realigned/{sample}_{aligner}_md_realigned_fixed.bam",
        ref = REFPATH+REFERENCE,
    output:
        stats = "results/00_Quality_Control/realigned/{sample}_{aligner}_md_realigned_fixed_stats.txt"
    log:
        "results/11_Reports/samtools/{sample}_{aligner}_realigned_fixed_stats.log"
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
        bam="results/04_Polishing/realigned/{sample}_{aligner}_md_realigned_fixed.bam",
        idx="results/04_Polishing/realigned/{sample}_{aligner}_md_realigned_fixed.bai",
    output:
        "results/00_Quality_Control/realigned/{sample}_{aligner}_md_realigned_fixed.idxstats.txt",
    log:
        "results/11_Reports/samtools/idxstats/{sample}_{aligner}_realigned_fixed_idxstats.log",
    params:
        extra="",  # optional params string
    wrapper:
         "v1.22.0/bio/samtools/idxstats"

###############################################################################
rule samtools_flagstat:
    input:
        expand("results/04_Polishing/realigned/{sample}_{aligner}_md_realigned_fixed.bam", sample=SAMPLE, aligner=ALIGNER),
    output:
        "results/00_Quality_Control/realigned/{sample}_{aligner}_md_realigned_fixed_bam.flagstat.txt",
    log:
        "results/11_Reports/samtools/flagstat/{sample}_{aligner}_realigned_fixed_bam.log",
    params:
        extra="",  # optional params string
    wrapper:
        "v1.22.0/bio/samtools/flagstat"

###############################################################################
rule multiqc:
    # Aim: aggregates bioinformatics analyses results into a single report
    # Use: multiqc [OPTIONS] --output [MULTIQC/] [FASTQC/] [FASTQSCREEN/] [qulimaps/]
    priority: 42
    message:
        "MultiQC reports aggregating"
    input:
        expand("results/00_Quality_Control/realigned/{sample}_{aligner}_md_realigned_fixed_bam.flagstat.txt", sample=SAMPLE, aligner=ALIGNER),
        expand("results/00_Quality_Control/realigned/{sample}_{aligner}_md_realigned_fixed.idxstats.txt", sample=SAMPLE, aligner=ALIGNER),
        expand("results/00_Quality_Control/realigned/{sample}_{aligner}_md_realigned_fixed_stats.txt", sample=SAMPLE, aligner=ALIGNER),
        expand("results/02_Mapping/{sample}_{aligner}_sorted-mark-dup_metrics.txt", sample=SAMPLE, aligner=ALIGNER),
        "results/00_Quality_Control/fastqc/", #expand("results/00_Quality_Control/fastqc/{sample}_fastqc.html", sample=SAMPLE),
        "results/00_Quality_Control/fastq-screen/",
    output:
        "results/00_Quality_Control/MULTIQC/multiqc_report.html"
    params:
        extra=""  # Optional: extra parameters for multiqc.
    log:
        "results/11_Reports/multiqc/multiqc.log"
    wrapper:
        "v1.22.0/bio/multiqc"





















