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
# CUTADAPT:
LENGTHc = config["cutadapt"]["length"]                      # Cutadapt --minimum-length
TRUSEQ_FOR = config["cutadapt"]["kits"]["truseq_for"]       # Cutadapt --adapter Illumina TruSeq
TRUSEQ_REV = config["cutadapt"]["kits"]["truseq_rev"]       # Cutadapt --adapter Illumina TruSeq
NEXTERA = config["cutadapt"]["kits"]["nextera"]             # Cutadapt --adapter Illumina Nextera
SMALL = config["cutadapt"]["kits"]["small"]                 # Cutadapt --adapter Illumina Small

# SICKLE-TRIM:
COMMAND = config["sickle-trim"]["command"]          # Sickle-trim command
ENCODING = config["sickle-trim"]["encoding"]        # Sickle-trim --qual-type
QUALITY = config["sickle-trim"]["quality"]          # Sickle-trim --qual-threshold
LENGTH = config["sickle-trim"]["length"]            # Sickle-trim --length-treshold

# TRIMMOMATIC:
TRUSEQ2_PE: config['trimmomatic']["adapters"]["truseq2-pe"]     # Truseq2-PE adapters fasta
TRUSEQ2_SE: config['trimmomatic']["adapters"]["truseq2-se"]     # Truseq2-SE adapters fasta
TRUSEQ3_PE: config['trimmomatic']["adapters"]["truseq3-pe"]     # Truseq2-PE adapters fasta
TRUSEQ3_PE2: config['trimmomatic']["adapters"]["truseq3-pe-2"]   # Truseq3-PE2 adapters fasta
TRUSEQ3_SE: config['trimmomatic']["adapters"]["truseq3-se"]     # Truseq2-PE adapters fasta

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
        fastqc          = "results/00_Quality_Control/fastqc/", #expand("results/00_Quality_Control/fastqc/{sample}_fastqc.html", sample=SAMPLE),
        trimmed_fastqc  = "results/00_Quality_Control/trimmed_fastqc/",
        forward_reads   = expand("results/01_Trimming/sickle/{sample}_sickle-trimmed_R1.fastq.gz", sample=SAMPLE),
        reverse_reads   = expand("results/01_Trimming/sickle/{sample}_sickle-trimmed_R2.fastq.gz", sample=SAMPLE),
        #forwardUnpaired = expand("results/01_Trimmimg/trimmomatic/{sample}_trimmomatic_unpaired_R1.fastq.gz", sample=SAMPLE),
        #reverseUnpaired = expand("results/01_Trimmimg/trimmomatic/{sample}_trimmomatic_unpaired_R2.fastq.gz", sample=SAMPLE),
        fastqscreen     = "results/00_Quality_Control/fastq-screen/",
        multiqc         = "results/00_Quality_Control/MULTIQC/multiqc_report.html",

################################ R U L E S ####################################
rule fastqc_quality_control:
    # Aim: reads sequence files and produces a quality control report
    # Use: fastqc [OPTIONS] --output [DIR/] [SAMPLE_1.fastq] ... [SAMPLE_n.fastq]
    priority: 1
    message:
        "FastQC reads quality controling"
    conda:
        FASTQC
    resources:
        cpus = CPUS
    input:
        fastq = "resources/reads/"
    output:
        fastqc = directory("results/00_Quality_Control/fastqc/")
    log:
        "results/11_Reports/quality/fastqc.log"
    shell:
        "mkdir -p {output.fastqc} " # (*) this directory must exist as the program will not create it
        "2> /dev/null && "          # in silence and then...
        "fastqc "                    # FastQC, a high throughput sequence QC analysis tool
        "--quiet "                    # -q: Supress all progress messages on stdout and only report errors
        "--threads {resources.cpus} " # -t: Specifies files number which can be processed simultaneously
        "--outdir {output.fastqc} "   # -o: Create all output files in the specified output directory (*)
        "{input.fastq}/*.fastq.gz "   # Input file.fastq
        "&> {log}"                    # Log redirection

###############################################################################
rule fastqscreen_contamination_checking:
    # Aim: screen if the composition of the library matches with what you expect
    # Use fastq_screen [OPTIONS] --outdir [DIR/] [SAMPLE_1.fastq] ... [SAMPLE_n.fastq]
    message:
        "Fastq-Screen reads contamination checking"
    conda:
        FASTQSCREEN
    resources:
        cpus = CPUS
    params:
        config = CONFIG,
        mapper = MAPPER,
        subset = SUBSET
    input:
        fastq = "resources/reads/"
    output:
        fastqscreen = directory("results/00_Quality_Control/fastq-screen/")
    log:
        "results/11_Reports/quality/fastq-screen.log"
    shell:
        "fastq_screen "                  # FastqScreen, what did you expect ?
        "-q "                            # --quiet: Only show log warning
        "--threads {resources.cpus} "    # --threads: Specifies across how many threads bowtie will be allowed to run
        "--conf {params.config} "        # path to configuration file
        "--aligner {params.mapper} "     # -a: choose aligner 'bowtie', 'bowtie2', 'bwa'
        "--subset {params.subset} "      # Don't use the whole sequence file, but create a subset of specified size
        "--outdir {output.fastqscreen} " # Output directory
        "{input.fastq}/*.fastq.gz "      # Input file.fastq
        "&> {log}"                       # Log redirection

###############################################################################
rule cutadapt_adapters_removing:
    # Aim: finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequence from your high-throughput sequencing reads
    # Use: cutadapt [OPTIONS] -a/-A [ADAPTER] -o [OUT-FWD.fastq.gz] -p [OUT-REV.fastq.gz] [IN-FWD.fastq.gz] [IN-REV.fastq.gz]
    # Rmq: multiple adapter sequences can be given using further -a options, but only the best-matching adapter will be removed
    message:
        "Cutadapt adapters removing for {wildcards.sample} sample"
    conda:
        CUTADAPT
    resources:
        cpus = CPUS
    params:
        length = LENGTHc,
        truseq_for = TRUSEQ_FOR,
        truseq_rev = TRUSEQ_REV,
        nextera = NEXTERA,
        small = SMALL
    input:
        fwdreads = "resources/reads/{sample}_R1.fastq.gz",
        revreads = "resources/reads/{sample}_R2.fastq.gz"
    output:
        fwdreads = temp("results/01_Trimming/cutadapt/{sample}_cutadapt-removed_R1.fastq.gz"),
        revreads = temp("results/01_Trimming/cutadapt/{sample}_cutadapt-removed_R2.fastq.gz")
    log:
        "results/11_Reports/cutadapt/{sample}.log"
    shell:
       "cutadapt "                           # Cutadapt, finds and removes unwanted sequence from your HT-seq reads
        "-j {resources.cpus} "               # -j: Number of CPU cores to use. Use 0 to auto-detect (default: 1)
        "--trim-n "                          # --trim-n: Trim N's on ends of reads
        "--minimum-length {params.length} "  # -m: Discard reads shorter than length
        "--adapter {params.truseq_for} "     # -a: Sequence of an adapter ligated to the 3' end of the first read
        "-A {params.truseq_rev} "            # -A: 3' adapter to be removed from second read in a pair
        "--adapter {params.nextera} "        # -a: Sequence of an adapter ligated to the 3' end of the first read
        "-A {params.nextera} "               # -A: 3' adapter to be removed from second read in a pair
        "--adapter {params.small} "          # -a: Sequence of an adapter ligated to the 3' end of the first read
        "-A {params.small} "                 # -A: 3' adapter to be removed from second read in a pair
        "--output {output.fwdreads} "        # -o: Write trimmed reads to FILE
        "--paired-output {output.revreads} " # -p: Write second read in a pair to FILE
        "{input.fwdreads} "                  # Input forward reads R1.fastq
        "{input.revreads} "                  # Input reverse reads R2.fastq
        "&> {log}"                           # Log redirection

###############################################################################
rule sickle_trim_quality:
    # Aim: windowed adaptive trimming tool for FASTQ files using quality
    # Use: sickle [COMMAND] [OPTIONS]
    message:
        "Sickle-trim low quality sequences trimming for {wildcards.sample} sample"
    conda:
        SICKLETRIM
    params:
        command = COMMAND,
        encoding = ENCODING,
        quality = QUALITY,
        length = LENGTH
    input:
        fwdreads = "results/01_Trimming/cutadapt/{sample}_cutadapt-removed_R1.fastq.gz",
        revreads = "results/01_Trimming/cutadapt/{sample}_cutadapt-removed_R2.fastq.gz"
    output:
        fwdreads = "results/01_Trimming/sickle/{sample}_sickle-trimmed_R1.fastq.gz",
        revreads = "results/01_Trimming/sickle/{sample}_sickle-trimmed_R2.fastq.gz",
        single = temp("results/01_Trimming/sickle/{sample}_sickle-trimmed_SE.fastq.gz")
    log:
        "results/11_Reports/sickle-trim/{sample}.log"
    shell:
       "sickle "                 # Sickle, a windowed adaptive trimming tool for FASTQ files using quality
        "{params.command} "      # Paired-end or single-end sequence trimming
        "-t {params.encoding} "  # --qual-type: Type of quality values, solexa ; illumina ; sanger ; CASAVA, < 1.3 ; 1.3 to 1.7 ; >= 1.8
        "-q {params.quality} "   # --qual-threshold: Threshold for trimming based on average quality in a window (default: 20)
        "-l {params.length} "    # --length-threshold: Threshold to keep a read based on length after trimming (default: 20)
        "-f {input.fwdreads} "   # --pe-file1: Input paired-end forward fastq file
        "-r {input.revreads} "   # --pe-file2: Input paired-end reverse fastq file
        "-g "                    # --gzip-output: Output gzipped files
        "-o {output.fwdreads} "  # --output-pe1: Output trimmed forward fastq file
        "-p {output.revreads} "  # --output-pe2: Output trimmed reverse fastq file (must use -s option)
        "-s {output.single} "    # --output-single: Output trimmed singles fastq file
        "&> {log}"               # Log redirection

###############################################################################
rule trimmomatic:
    # Aim : Trimmomatic: a flexible read trimming tool for Illumina NGS data.
    message:
        "Trimming reads for {wildcards.sample}"
    conda:
        TRIMMOMATIC
    input:
        r1="resources/reads/{sample}_R1.fastq.gz",
        r2="resources/reads/{sample}_R2.fastq.gz",
        adapters = config["trimmomatic"]["adapters"]["truseq2-pe"]
    output:
        forward_reads   = "results/01_Trimmimg/trimmomatic/{sample}_trimmomatic_R1.fastq.gz",
        reverse_reads   = "results/01_Trimmimg/trimmomatic/{sample}_trimmomatic_R2.fastq.gz",
        forwardUnpaired = "results/01_Trimmimg/trimmomatic/{sample}_trimmomatic_unpaired_R1.fastq.gz",
        reverseUnpaired = "results/01_Trimmimg/trimmomatic/{sample}_trimmomatic_unpaired_R2.fastq.gz"
    log:
        "results/11_Reports/trimmomatic/{sample}.log"
    params:
        seedMisMatches =            str(config['trimmomatic']['seedMisMatches']),
        palindromeClipTreshold =    str(config['trimmomatic']['palindromeClipTreshold']),
        simpleClipThreshhold =      str(config['trimmomatic']['simpleClipThreshold']),
        LeadMinTrimQual =           str(config['trimmomatic']['LeadMinTrimQual']),
        TrailMinTrimQual =          str(config['trimmomatic']['TrailMinTrimQual']),
        windowSize =                str(config['trimmomatic']['windowSize']),
        avgMinQual =                str(config['trimmomatic']['avgMinQual']),
        minReadLen =                str(config['trimmomatic']['minReadLength']),
        phred = 		            str(config["trimmomatic"]["phred"])
    threads: CPUS
    shell:
        """
        trimmomatic PE \
        -threads {threads} \
        {params.phred} \
        {input.r1} {input.r2} \
        {output.forward_reads} {output.forwardUnpaired} {output.reverse_reads} {output.reverseUnpaired} \
        ILLUMINACLIP:{input.adapters}:{params.seedMisMatches}:{params.palindromeClipTreshold}:{params.simpleClipThreshhold} \
        LEADING:20 \
        TRAILING:3 \
        SLIDINGWINDOW:5:20 \
        AVGQUAL:20 \
        MINLEN:50 &>{log}
        """

###############################################################################
rule trimmed_fastqc:
    # Aim : reads sequence files and produces a quality control report after trimming
    # Use: fastqc [OPTIONS] --output [DIR/] [SAMPLE_1.fastq] ... [SAMPLE_n.fastq]
    message:
        "Quality check after trimming"
    conda:
        FASTQC
    threads:
        CPUS
    input:
        forward_reads = expand("results/01_Trimming/sickle/{sample}_sickle-trimmed_R1.fastq.gz", sample=SAMPLE),
        reverse_reads = expand("results/01_Trimming/sickle/{sample}_sickle-trimmed_R2.fastq.gz", sample=SAMPLE),
    output:
        fastqc = directory("results/00_Quality_Control/trimmed_fastqc/")
    log:
        "results/11_Reports/trimmed_fastqc/trimmed.fastqc.log"
    shell:
        "mkdir -p {output.fastqc} "     # (*) this directory must exist as the program will not create it
        "2> /dev/null && "              # in silence and then...
        "fastqc "                       # FastQC, a high throughput sequence QC analysis tool
        "--quiet "                      # -q: Supress all progress messages on stdout and only report errors
        "--threads {threads} "          # -t: Specifies files number which can be processed simultaneously
        "--outdir {output.fastqc} "     # -o: Create all output files in the specified output directory (*)
        "{input.forward_reads} {input.reverse_reads} "     # Input file.fastq
        "&> {log}"                      # Log redirection

###############################################################################
rule multiqc:
    # Aim: aggregates bioinformatics analyses results into a single report
    # Use: multiqc [OPTIONS] --output [MULTIQC/] [FASTQC/] [FASTQSCREEN/] [qulimaps/]
    priority: 42
    message:
        "MultiQC reports aggregating"
    input:
        "results/00_Quality_Control/fastqc/",
        "results/00_Quality_Control/fastq-screen/",
        "results/00_Quality_Control/trimmed_fastqc/",
    output:
        "results/00_Quality_Control/MULTIQC/multiqc_report.html"
    params:
        extra=""  # Optional: extra parameters for multiqc.
    log:
        "results/11_Reports/multiqc/multiqc.log"
    wrapper:
        "v1.21.2/bio/multiqc"
