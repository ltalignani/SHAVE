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
        multiqc = "results/00_Quality_Control/MULTIQC/multiqc_report.html",
        fastqc = "results/00_Quality_Control/fastqc/", #expand("results/00_Quality_Control/fastqc/{sample}_fastqc.html", sample=SAMPLE),
        trimmed_fastqc = "results/00_Quality_Control/trimmed_fastqc/",
        qualimap = expand("results/00_Quality_Control/qualimap/{sample}_{aligner}/qualimapReport.html", sample=SAMPLE, aligner=ALIGNER),
        fastqscreen = "results/00_Quality_Control/fastq-screen/",
        indexvcf = expand("results/05_Variants/{sample}_{aligner}.vcf.gz.tbi", sample=SAMPLE, aligner=ALIGNER),
        vcf_filtered ="results/05_Variants/merged_filtered/merged_hardfiltered.vcf.gz",
        table = "results/05_Variants/merged_raw/merged.table",
        combinegvcfs = "results/05_Variants/merged_raw/merged.vcf.gz",
        bgzip_vcfs = expand("results/05_Variants/{sample}_{aligner}.vcf.gz", sample=SAMPLE, aligner=ALIGNER),
        vcf = expand("results/05_Variants/{sample}_{aligner}.vcf", sample=SAMPLE, aligner=ALIGNER),
        check = expand("results/00_Quality_Control/validatesamfile/{sample}_{aligner}_md_realigned_fixed_ValidateSam.txt", sample=SAMPLE, aligner=ALIGNER),
        flagstat = expand("results/00_Quality_Control/realigned/{sample}_{aligner}_md_realigned_fixed_bam.flagstat.txt", sample=ALIGNER, aligner=ALIGNER),
        idxstats = expand("results/00_Quality_Control/realigned/{sample}_{aligner}_md_realigned_fixed.idxstats.txt", sample=SAMPLE, aligner=ALIGNER),
        stats = expand("results/00_Quality_Control/realigned/{sample}_{aligner}_md_realigned_fixed_stats.txt", sample=SAMPLE, aligner=ALIGNER),        
        igv_output = expand("results/04_Polishing/{sample}_{aligner}_realignertargetcreator.bed", sample=SAMPLE, aligner=ALIGNER),
        index_post_realign = expand("results/04_Polishing/realigned/{sample}_{aligner}_md_realigned_fixed.bai", sample=SAMPLE, aligner=ALIGNER),
        fixmateinformation = expand("results/04_Polishing/realigned/{sample}_{aligner}_md_realigned_fixed.bam", sample=SAMPLE, aligner=ALIGNER),

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
        forward_reads   = temp("results/01_Trimmimg/trimmomatic/{sample}_trimmomatic_R1.fastq.gz"),
        reverse_reads   = temp("results/01_Trimmimg/trimmomatic/{sample}_trimmomatic_R2.fastq.gz"),
        forwardUnpaired = temp("results/01_Trimmimg/trimmomatic/{sample}_trimmomatic_unpaired_R1.fastq.gz"),
        reverseUnpaired = temp("results/01_Trimmimg/trimmomatic/{sample}_trimmomatic_unpaired_R2.fastq.gz")
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
        forward_reads = expand("results/01_Trimmimg/trimmomatic/{sample}_trimmomatic_R1.fastq.gz", sample=SAMPLE),
        reverse_reads = expand("results/01_Trimmimg/trimmomatic/{sample}_trimmomatic_R2.fastq.gz", sample=SAMPLE),
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
rule bowtie2_mapping:
    # Aim: reads mapping against reference sequence
    # Use: bowtie2 -p [THREADS] -x [REFERENCE] -1 [FWD_R1.fq] -2 [REV_R2.fq] -S [MAPPED.sam]
    message:
        "Bowtie2 mapping {wildcards.sample} sample reads against reference genome sequence"
    conda:
        BOWTIE2
    resources:
        cpus = CPUS
    params:
        bt2path = BT2PATH,
        reference = REFERENCE,
        sensitivity = SENSITIVITY
    input:
        fwdreads = "results/01_Trimmimg/trimmomatic/{sample}_trimmomatic_R1.fastq.gz", 
        revreads = "results/01_Trimmimg/trimmomatic/{sample}_trimmomatic_R2.fastq.gz"
    output:
        mapped = temp("results/02_Mapping/{sample}_bowtie2-mapped.sam")
    benchmark:
        "benchmarks/bowtie2/{sample}.tsv"
    log:
        "results/11_Reports/bowtie2/{sample}.log"
    shell:
        "bowtie2 "                    # Bowtie2, an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences.
        "--threads {resources.cpus} "  # -p: Number of alignment threads to launch (default: 1)
        "--reorder "                   # Keep the original read order (if multi-processor option -p is used)
        "-x {params.bt2path}{params.reference} " # -x: Reference index filename prefix (minus trailing .X.bt2) [Bowtie-1 indexes are not compatible]
        "{params.sensitivity} "        # Preset (default: "--sensitive", same as [-D 15 -R 2 -N 0 -L 22 -i S,1,1.15])
        "-q "                          # -q: Query input files are FASTQ .fq/.fastq (default)
        "-1 {input.fwdreads} "         # Forward input reads
        "-2 {input.revreads} "         # Reverse input reads
        "1> {output.mapped} "          # -S: File for SAM output (default: stdout)
        "2> {log}"                     # Log redirection


##############################################################################
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
        extra = r"-R '@RG\tID:{sample}\tSM:{sample}\tCN:SC\tPL:ILLUMINA'" # Manage ReadGroup
    input:
        fwdreads = "results/01_Trimmimg/trimmomatic/{sample}_trimmomatic_R1.fastq.gz",
        revreads = "results/01_Trimmimg/trimmomatic/{sample}_trimmomatic_R2.fastq.gz"
    output:
        mapped = temp("results/02_Mapping/{sample}_bwa-mapped.sam")
    benchmark:
        "benchmarks/bwa/{sample}.tsv"
    log:
        "results/11_Reports/bwa/{sample}.log"
    shell:
        "bwa mem "                                                  # BWA-MEM algorithm, performs local alignment.
        "-M "                                                       # Mark shorter split hits as secondary (for Picard compatibility). 
        "-T 0 "                                                     # Don’t output alignment with score lower than INT. This option only affects output.
        "-t {resources.cpus} "                                      # -t: Number of threads (default: 12)
        "-v 1 "                                                     # -v: Verbosity level: 1=error, 2=warning, 3=message, 4+=debugging
        "{params.extra} "                                           # -R: Complete read group header line.
        "{params.ref} "                                             # Reference genome
        "{input.fwdreads} "                                         # Forward input reads
        "{input.revreads} "                                         # Reverse input reads
        "1> {output.mapped} "                                       # SAM output
        "2> {log}"                                                  # Log redirection

###############################################################################
rule samtools_view:
    # Aim: Convert or filter SAM/BAM.
    # Use : samtools view -bo aln.bam aln.sam
    message:
        "Samtools view conversion of sample sam in bam format"
    input:
        "results/02_Mapping/{sample}_{aligner}-mapped.sam",
    output:
        bam = temp("results/02_Mapping/{sample}_{aligner}.bam"),
    log:
        "results/11_Reports/samtools_view/{sample}_{aligner}.log",
    params:
        extra="",  # optional params string
        region="",  # optional region string
    threads: CPUS
    wrapper:
        "v1.22.0/bio/samtools/view"

###############################################################################
rule samtools_sort:
    # Aim: sorting
    # Use: samtools sort -@ [THREADS] -m [MEM] -T [TMPDIR] -O BAM -o [SORTED.bam] [FIXMATE.bam]
    message:
        "samtools sort {wildcards.sample} sample reads ({wildcards.aligner})"
    conda:
        SAMTOOLS
    resources:
       cpus = CPUS,
       mem_gb = MEM_GB
    params:
        tmpdir = TMPDIR
    input:
        bam = "results/02_Mapping/{sample}_{aligner}.bam"
    output:
        sorted = temp("results/02_Mapping/{sample}_{aligner}_sorted.bam")
    log:
        "results/11_Reports/samtools/{sample}_{aligner}_sorted.log"
    shell:
        "samtools sort "               # Samtools sort, tools for alignments in the SAM format with command to sort alignment file
        "--threads {resources.cpus} "  # -@: Number of additional threads to use (default: 1)
        "-m {resources.mem_gb}G "      # -m: Set maximum memory per thread, suffix K/M/G recognized (default: 768M)
        "-T {params.tmpdir} "          # -T: Write temporary files to PREFIX.nnnn.bam
        "--output-fmt BAM "            # -O: Specify output format: SAM, BAM, CRAM (here, BAM format)
        "-o {output.sorted} "          # Sorted bam output
        "{input.bam} "                 # bam input
        "&> {log}"                     # Log redirection

###############################################################################
rule mark_duplicates_spark:
    input:
        "results/02_Mapping/{sample}_{aligner}_sorted.bam",
    output:
        bam = temp("results/02_Mapping/{sample}_{aligner}_sorted-mark-dup.bam"),
        metrics="results/02_Mapping/{sample}_{aligner}_sorted-mark-dup_metrics.txt",
    benchmark:
        "benchmarks/markduplicatesspark/{sample}_{aligner}.tsv"
    log:
        "results/11_Reports/samtools/{sample}_{aligner}_sorted-mark-dup.log",
    params:
        extra="--remove-sequencing-duplicates",  # optional
        #java_opts=,  # optional
        #spark_runner="",  # optional, local by default
        #spark_v1.19.1="",  # optional
        #spark_extra="", # optional
    resources:
        mem_mb=16000,
    threads: CPUS
    wrapper:
        "v1.22.0/bio/gatk/markduplicatesspark"

###############################################################################
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
        "SamTools indexing marked as duplicate merged BAM file"
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
rule unifiedgenotyper:
    # Aim:  Call variants in sequence data. The following parameters comes from the MalariaGEN. There is two difference w/ the 2017's article: --output_mode EMIT_VARIANTS_ONLY &  --genotyping_mode DISCOVERY
    # Use:  java -jar GenomeAnalysisTK.jar \ 
    #       -T UnifiedGenotyper \
    #       -nct {threads.cpus} \ # -nt / --num_threads controls the number of data threads sent to the processor 
    #       -I {sample BAM} \
    #       --alleles {alleles VCF} \ : This option has been removed for the moment.  Alleles against which to genotype (VCF format). Given the sites VCF file is fixed for every sample, and we wish to generalise to future sets of sites/alleles, the VCF file describing sites and alleles should be considered a parameter. This file for A. gambiae (AgamP4) is available at
    #       -R {reference sequence} \
    #       --out {output VCF} \
    #       -A: Annotation \
    #       -XA: eXclude Annotation
    #
    #       Annotations needed as decision tools: see https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
    message:
        "UnifiedGenotyper calling SNVs"
    conda:
        GATK
    input:
        bam = "results/04_Polishing/realigned/{sample}_{aligner}_md_realigned_fixed.bam",
        ref = "resources/genomes/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa",
        index = "results/04_Polishing/realigned/{sample}_{aligner}_md_realigned_fixed.bai",
        intervals = "/Users/loictalignani/snakemake-models/shave/resources/genomes/chrom_list.intervals",
        #alleles = ALLELES
    output:
        vcf="results/05_Variants/{sample}_{aligner}.vcf"
    log:
        "results/11_Reports/unifiedgenotyper/{sample}_{aligner}.log"
    benchmark:
        "benchmarks/unifiedgenotyper/{sample}_{aligner}.tsv"
    threads: CPUS
    shell:
        "gatk3 -T UnifiedGenotyper "                    # Genome Analysis Tool Kit - Broad Institute UnifiedGenotyper
        "-nct {threads} "                               # -nct / --num_cpu_threads_per_data_thread controls the number of CPU threads allocated to each data thread
        "-I {input.bam} "                               # Input indel realigned BAM file
        "-R {input.ref} "                               # Reference sequence in fasta format
        "-L {input.intervals} "                         # interval list (chromosome names)
        "--out {output.vcf} "                           # Output VCF
        "--genotype_likelihoods_model BOTH "            # Genotype likelihoods calculation model to employ -- BOTH is the default option, while INDEL is also available for calling indels and SNP is available for calling SNPs only (SNP|INDEL|BOTH)
        "--genotyping_mode GENOTYPE_GIVEN_ALLELES "                  # Should we output confident genotypes (i.e. including ref calls) or just the variants? (DISCOVERY|GENOTYPE_GIVEN_ALLELES)
        "--heterozygosity 0.01 "                        # Heterozygosity value used to compute prior likelihoods for any locus
        "--indel_heterozygosity 0.001 "                 # Heterozygosity for indel calling
        "--downsampling_type BY_SAMPLE "                # Type of reads downsampling to employ at a given locus. Reads will be selected randomly to be removed from the pile based on the method described here (NONE|ALL_READS| BY_SAMPLE) given locus
        "-dcov 250 "                                    # downsampling coverage
        "--output_mode EMIT_ALL_SITES "                 # Should we output confident genotypes (i.e. including ref calls) or just the variants? (EMIT_VARIANTS_ONLY|EMIT_ALL_CONFIDENT_SITES|EMIT_ALL_SITES)
        "--min_base_quality_score 17 "                  # Minimum base quality required to consider a base for calling
        "-stand_call_conf 30.0 "                        # standard min confidence-threshold for calling. This is the minimum confidence threshold (phred-scaled) at which the program should emit variant sites as called.
        "-contamination 0.05 "                          # Define the fraction of contamination in sequence data (for all samples) to aggressively remove.
        "-A DepthPerAlleleBySample "                    # 
        "-A RMSMappingQuality "                         # MQ: needed as decision tools for hard-filtering. Compares the mapping qualities of the reads supporting the reference allele and the alternate allele. positive value means the mapping qualities of the reads supporting the alternate allele are higher than those supporting the reference allele
        "-A Coverage "                                  # 
        "-A FisherStrand "                              # FS: needed as decision tools for hard-filtering. Strand Bias tells us whether the alternate allele was seen more or less often on the forward or reverse strand than the reference allele. measure strand bias (a type of sequencing bias in which one DNA strand is favored over the other, which can result in incorrect evaluation of the amount of evidence observed for one allele vs. the other
        "-A StrandOddsRatio "                           # SOR: needed as decision tools for hard-filtering. created because FS tends to penalize variants that occur at the ends of exons. Reads at the ends of exons tend to only be covered by reads in one direction
        "-A BaseQualityRankSumTest "                    #
        "-A MappingQualityRankSumTest "                 # MQRankSum: needed as decision tools for hard-filtering
        "-A QualByDepth "                               # QD: needed as decision tools for hard-filtering. Intended to normalize the variant quality in order to avoid inflation caused when there is deep coverage. This is the variant confidence (from the QUAL field) divided by the unfiltered depth of non-hom-ref samples. better to use QD than either QUAL or DP directly.
        "-A ReadPosRankSumTest "                        # ReadPosRankSum: needed as decision tools for hard-filtering. z-approximation from the Rank Sum Test for site position within reads. Compares whether the positions of the reference and alternate alleles are different within the reads.
        "-XA ExcessHet "                                # 
        "-XA InbreedingCoeff "                          #
        "-XA MappingQualityZero "                       #
        "-XA HaplotypeScore "                           #
        "-XA SpanningDeletions "                        #
        "-XA ChromosomeCounts "                         #
        " > {log} 2>&1"

###############################################################################
rule bgzip_vcfs:
    input:
        expand("results/05_Variants/{sample}_{aligner}.vcf", sample=SAMPLE, aligner=ALIGNER),
    output:
        "results/05_Variants/{sample}_{aligner}.vcf.gz",
    params:
        extra="", # optional
    threads: CPUS
    log:
        "results/11_Reports/bgzip/{sample}_{aligner}.vcf.gz.log",
    wrapper:
        "v1.22.0/bio/bgzip"

###############################################################################
rule indexfeaturefile:
    message:
        "Indexing VCFs file for CombineGVCFs"
    conda:
        GATK4
    input:
        vcf ="results/05_Variants/{sample}_{aligner}.vcf.gz",
    output:
        indexvcf = "results/05_Variants/{sample}_{aligner}.vcf.gz.tbi",
    log:
        "results/11_Reports/indexfeaturefile/{sample}_{aligner}.indexvcf.tbi.log",
    shell:
        """
        gatk IndexFeatureFile -I {input.vcf} -O {output.indexvcf} > {log} 2>&1 || true
        """

###############################################################################
rule combinegvcfs:
    message:
        "GATK CombineGVCFs merge VCFs"
    input:
        gvcfs=expand("results/05_Variants/{sample}_{aligner}.vcf.gz", sample=SAMPLE, aligner=ALIGNER),
        ref=REFPATH+REFERENCE,
    output:
        gvcf="results/05_Variants/merged_raw/merged.vcf.gz",
    log:
        "results/11_Reports/combinegvcfs/combinegvcfs.log",
    params:
        extra="",  # optional
        java_opts="",  # optional
    resources:
        mem_mb=16000,
    wrapper:
        "v1.22.0/bio/gatk/combinegvcfs"

###############################################################################
rule variantstotable:
    # Aim:  Call variants in sequence data. The following parameters comes from the MalariaGEN
    # Use:  gatk VariantsToTable\ 
    #       -V input.vcf \
    #       -F CHROM -F POS -F TYPE -GF AD \
    #       -O output.table
    message:
        "VariantsToTable"
    conda:
        GATK4
    input:
        vcf = "results/05_Variants/merged_raw/merged.vcf.gz",
    output:
        table = "results/05_Variants/merged_raw/merged.table",
    log:
        "results/11_Reports/variantstotable/merged_table.log",
    shell:
        """
        gatk VariantsToTable -V {input.vcf} -O {output.table} -F CHROM -F POS -F QUAL -F NS -F DP -F MQ -F FS -F SOR -F MQRankSum -F ReadPosRankSum -GF AD 2> {log}
        """

###############################################################################
rule gatk_filter:
    # Aim: Filter variant calls based on INFO and/or FORMAT annotations.
    # Use: gatk VariantFiltration \
    # -R reference.fasta \
    # -V input.vcf.gz \
    # -O output.vcf.gz
    # --filter-name "my_filter1" \
    # --filter-expression "AB < 0.2" \
    # --filter-name "my_filter2" \
    # --filter-expression "MQ0 > 50"
    message:
        "VariantFiltration Hard-filtering"
    input:
        ref=REFPATH+REFERENCE,  #"resources/genomes/GCA_018104305.1_AalbF3_genomic.fasta",
        vcf="results/05_Variants/merged_raw/merged.vcf.gz",
    output:
        vcf="results/05_Variants/merged_filtered/merged_hardfiltered.vcf.gz",
    params:
        filters={"myfilter": "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"},
        extra="",
        java_opts="",
    resources:
        mem_mb=16000
    log:
        "results/11_Reports/variantfiltration/merged_hardfiltered.log",
    wrapper:
        "v1.22.0/bio/gatk/variantfiltration"

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

###############################################################################
rule multiqc_reports_aggregation:
    # Aim: aggregates bioinformatics analyses results into a single report
    # Use: multiqc [OPTIONS] --output [MULTIQC/] [FASTQC/] [FASTQSCREEN/] [qulimaps/]
    priority: 42
    message:
        "MultiQC reports aggregating"
    conda:
        MULTIQC
    input:
        fastqc = "results/00_Quality_Control/fastqc/",
        trimmed_fastqc = "results/00_Quality_Control/trimmed_fastqc/",
        fastqscreen = "results/00_Quality_Control/fastq-screen/",
        qualimap = "results/00_Quality_Control/qualimap/multisampleBamQcReport.html",
    output:
        multiqc = directory("results/00_Quality_Control/multiqc_old/")
    log:
        "results/11_Reports/multiqc/multiqc_old.log"
    shell:
        "multiqc "                  # Multiqc, searches in given directories for analysis & compiles a HTML report
        "--outdir {output.multiqc} " # -o: Create report in the specified output directory
        "{input.fastqc} "            # Input FastQC files
        "{input.fastqscreen} "       # Input Fastq-Screen
        "{input.qualimap} "          # Input Qualimap
        "--no-ansi "                 # Disable coloured log
        "&> {log}"                   # Log redirection






















