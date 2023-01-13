#######M######I#####V#####E######G######E######C######I#######R######D#########
# Name:                 shave.smk
# Author:               Loïc TALIGNANI
# Affiliation:          IRD_UMR224_MIVEGEC
# Aim:                  Snakefile for SHort-read Alignment pipeline for VEctor 
# Date:                 2022.10.05
# Run:                  snakemake --snakefile workflow/rules/shave.smk --cores X --use-conda
# Latest modification:  2022.11.23
# Done:                 Added RealignerTargetCreator, IndelRealigner, 
#                       UnifiedGenotyper, VariantFiltration, compression and tabix

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
MULTIQC = config["conda"][OS]["multiqc"]            # MultiQC
CUTADAPT = config["conda"][OS]["cutadapt"]          # Cutadapt
SICKLETRIM = config["conda"][OS]["sickle-trim"]     # Sickle-trim
BOWTIE2 = config["conda"][OS]["bowtie2"]            # Bowtie2
BWA = config["conda"][OS]["bwa"]                    # Bwa
SAMTOOLS = config["conda"][OS]["samtools"]          # SamTools
BEDTOOLS = config["conda"][OS]["bedtools"]          # BedTools
BCFTOOLS = config["conda"][OS]["bcftools"]          # BcfTools
GAWK = config["conda"][OS]["gawk"]                  # Gawk
LOFREQ = config["conda"][OS]["lofreq"]              # LoFreq
GATK = config["conda"][OS]["gatk"]                  # GATK 3.8
GATK4 = config["conda"][OS]["gatk4"]                # GATK 4.3.0
PICARD = config["conda"][OS]["picard"]              # Picard 2.24.7
QUALIMAP = config["conda"][OS]["qualimap"]          # Qualimap 2.2.2d

###############################################################################
# PARAMETERS #
LENGTHc = config["cutadapt"]["length"]              # Cutadapt --minimum-length
TRUSEQ = config["cutadapt"]["kits"]["truseq"]       # Cutadapt --adapter Illumina TruSeq
NEXTERA = config["cutadapt"]["kits"]["nextera"]     # Cutadapt --adapter Illumina Nextera
SMALL = config["cutadapt"]["kits"]["small"]         # Cutadapt --adapter Illumina Small

COMMAND = config["sickle-trim"]["command"]          # Sickle-trim command
ENCODING = config["sickle-trim"]["encoding"]        # Sickle-trim --qual-type
QUALITY = config["sickle-trim"]["quality"]          # Sickle-trim --qual-threshold
LENGTH = config["sickle-trim"]["length"]            # Sickle-trim --length-treshold

CONFIG = config["fastq-screen"]["config"]           # Fastq-screen --conf
MAPPER = config["fastq-screen"]["aligner"]          # Fastq-screen --aligner
SUBSET = config["fastq-screen"]["subset"]           # Fastq-screen --subset

ALIGNER = config["aligner"]                         # Aligners ('bwa' or 'bowtie2')
MARKDUP = config["markdup"]                         # Mark Duplicate Program ('picard' or 'samtools')

BWAPATH = config["bwa"]["path"]                     # BWA path to indexes
BT2PATH = config["bowtie2"]["path"]                 # Bowtie2 path to indexes
SENSITIVITY = config["bowtie2"]["sensitivity"]      # Bowtie2 sensitivity preset

REFPATH = config["consensus"]["path"]               # Path to genomes references
REFERENCE = config["consensus"]["reference"]        # Genome reference sequence, in fasta format
ALLELES = config["alleles"]["alleles_target"]          # Alleles against which to genotype (VCF format) 
MINCOV = config["consensus"]["mincov"]              # Minimum coverage, mask lower regions with 'N'
MINAF = config["consensus"]["minaf"]                # Minimum allele frequency allowed
IUPAC = config["consensus"]["iupac"]                # Output variants in the form of IUPAC ambiguity codes

###############################################################################
# FUNCTIONS AND COMMANDS #

# Output Field Separator = <tab>. if the third column is empty, print col 1, col 2 -1 (transform coordinates from 1-based to 0-based for IGV) and col 2.  
# else, print the three column and transform the second one in 0-based coordinates.
# look for "robust quoting for snakemake for more info".
#AWK_CMD_INTERVALS = r"""'BEGIN { OFS = "\t" } { if( $3 == "") { print $1, $2-1, $2 } else { print $1, $2-1, $3}}'"""

############################### O N S T A R T ################################
onstart:
    print("##### Creating profile pipeline #####\n") 
    print("\t Creating jobs output subfolders...\n")
    shell("mkdir -p Cluster_logs/indelrealigner")
    shell("mkdir -p Cluster_logs/realignertargetcreator")
    shell("mkdir -p Cluster_logs/bwa_mapping")

############################# O N S U C C E S S ##############################
onsuccess:
    shutil.rmtree(".snakemake")

################################## A L L #####################################
rule all:
    input:
        multiqc = "results/00_Quality_Control/multiqc/",
        fastqc = "results/00_Quality_Control/fastqc/",
        fastqscreen = "results/00_Quality_Control/fastq-screen/",
        qualimap = "results/00_Quality_Control/qualimap/",
        index_archive = "results/04_Variants/merged_hardfiltered.vcf.gz.tbi",
        archive = "results/04_Variants/merged_hardfiltered.vcf.gz",
        hard_filter = "results/04_Variants/variantfiltration/merged_hardfiltered.vcf",
        vcf = "results/04_Variants/unifiedgenotyper/merged_indels.vcf",
        check = "results/05_Validation/validatesamfile/merged_realign_fixed_ValidateSam.txt",
        stats = "results/05_Validation/realigned/merged_realign_fixed_stats.txt",        
        igv_output = "results/04_Variants/merged_realignertargetcreator.bed",
        index_post_realign = "results/04_Variants/realigned/merged_realign_fixed.bai",
        fixmateinformation = "results/04_Variants/realigned/merged_realign_fixed.bam",
        covstats = "results/03_Coverage/merged_coverage-stats.tsv",
        markduplicate = "results/02_Mapping/merged-mark-dup.bam"

################################ R U L E S ####################################
rule multiqc_reports_aggregation:
    # Aim: aggregates bioinformatics analyses results into a single report
    # Use: multiqc [OPTIONS] --output [MULTIQC/] [FASTQC/] [MULTIQC/]
    priority: 42
    message:
        "MultiQC reports aggregating"
    conda:
        MULTIQC
    input:
        fastqc = "results/00_Quality_Control/fastqc/",
        fastqscreen = "results/00_Quality_Control/fastq-screen/",
        qualimap = "results/00_Quality_Control/qualimap/"
    output:
        multiqc = directory("results/00_Quality_Control/multiqc/")
    log:
        "results/11_Reports/quality/multiqc.log"
    shell:
        "multiqc "                  # Multiqc, searches in given directories for analysis & compiles a HTML report
        "--quiet "                   # -q: Only show log warning
        "--outdir {output.multiqc} " # -o: Create report in the specified output directory
        "{input.fastqc} "            # Input FastQC files
        "{input.fastqscreen} "       # Input Fastq-Screen
        "{input.qualimap} "          # Input Qualimap
        "--no-ansi "                 # Disable coloured log
        "&> {log}"                   # Log redirection

###############################################################################
rule tabix:
    input:
        "results/04_Variants/merged_hardfiltered.vcf.gz",
    output:
        "results/04_Variants/merged_hardfiltered.vcf.gz.tbi",
    log:
        "results/11_Reports/tabix/merged_hardfiltered.log",
    params:
        # pass arguments to tabix (e.g. index a vcf)
        "-p vcf",
    wrapper:
        "v1.21.1/bio/tabix/index"

###############################################################################
rule bgzip:
    input:
        "results/04_Variants/variantfiltration/merged_hardfiltered.vcf",
    output:
        "results/04_Variants/merged_hardfiltered.vcf.gz",
    params:
        extra="", # optional
    threads: CPUS
    log:
        "results/11_Reports/bgzip/merged_hardfiltered.vcf.gz.log",
    wrapper:
        "v1.21.1/bio/bgzip"

###############################################################################
rule hard_filter_calls:
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
    conda:
        GATK4
    input:
        ref="resources/genomes/GCA_018104305.1_AalbF3_genomic.fasta",
        vcf="results/04_Variants/unifiedgenotyper/merged_indels.vcf",
    output:
        vcf="results/04_Variants/variantfiltration/merged_hardfiltered.vcf",
    params:
        filters={"myfilter": "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"},
        extra="",
        java_opts="",
    resources:
        mem_gb=MEM_GB
    log:
        "results/11_Reports/variantfiltration/merged_hardfiltered.log",
    wrapper:
        "0.74.0/bio/gatk/variantfiltration"

###############################################################################
rule unifiedgenotyper:
    # Aim:  Call variants in sequence data. The following parameters comes from the MalariaGEN
    # Use:  java -jar GenomeAnalysisTK.jar \ 
    #       -T UnifiedGenotyper \
    #       -nct {threads.cpus} \ # -nt / --num_threads controls the number of data threads sent to the processor 
    #       -I {sample BAM} \
    #       --alleles {alleles VCF} \ : This option has been removed for the moment.  Alleles against which to genotype (VCF format). Given the sites VCF file is fixed for every sample, and we wish to generalise to future sets of sites/alleles, the VCF file describing sites and alleles should be considered a parameter. This file for A. gambiae (AgamP4) is available at
    #       -R {reference sequence} \
    #       --out {output VCF} \
    message:
        "UnifiedGenotyper calling SNVs"
    conda:
        GATK
    input:
        bam = "results/04_Variants/realigned/merged_realign_fixed.bam",
        ref = "resources/genomes/GCA_018104305.1_AalbF3_genomic.fasta",
        index = "results/04_Variants/realigned/merged_realign_fixed.bai"
        #alleles = ALLELES
    output:
        vcf="results/04_Variants/unifiedgenotyper/merged_indels.vcf"
    log:
        "results/11_Reports/unifiedgenotyper/merged_indels.log"
    benchmark:
        "benchmarks/unifiedgenotyper/merged_indels.tsv"
    threads: CPUS
    shell:
        "gatk3 -T UnifiedGenotyper "                    # Genome Analysis Tool Kit - Broad Institute UnifiedGenotyper
        "-nct {threads} "                               # -nct / --num_cpu_threads_per_data_thread controls the number of CPU threads allocated to each data thread
        "-I {input.bam} "                               # Input indel realigned BAM file
        "-R {input.ref} "                               # Reference sequence in fasta format
        "--out {output.vcf} "                           # Output VCF
        "--genotype_likelihoods_model BOTH "            # Genotype likelihoods calculation model to employ -- BOTH is the default option, while INDEL is also available for calling indels and SNP is available for calling SNPs only (SNP|INDEL|BOTH)
        "--genotyping_mode GENOTYPE_GIVEN_ALLELES "     # Should we output confident genotypes (i.e. including ref calls) or just the variants? (DISCOVERY|GENOTYPE_GIVEN_ALLELES)
        "--heterozygosity 0.015 "                       # Heterozygosity value used to compute prior likelihoods for any locus
        "--heterozygosity_stdev 0.05 "                  # Standard deviation of heterozygosity for SNP and indel calling
        "--indel_heterozygosity 0.001 "                 # Heterozygosity for indel calling
        "--downsampling_type BY_SAMPLE "                # Type of reads downsampling to employ at a given locus. Reads will be selected randomly to be removed from thepile based on the method described here (NONE|ALL_READS| BY_SAMPLE) given locus
        "-dcov 250 "                                    # downsampling coverage
        "--output_mode EMIT_ALL_SITES "                 # Should we output confident genotypes (i.e. including ref calls) or just the variants? (EMIT_VARIANTS_ONLY|EMIT_ALL_CONFIDENT_SITES|EMIT_ALL_SITES)
        "--min_base_quality_score 17 "                  # Minimum base quality required to consider a base for calling
        "-stand_call_conf 0.0 "                         # standard min confidence-threshold for calling
        "-contamination 0.0 "                           # Define the fraction of contamination in sequence data (for all samples) to aggressively remove.
        "-A DepthPerAlleleBySample "                    # 
        "-XA RMSMappingQuality "                        # 
        "-XA Coverage "                                 #
        "-XA ExcessHet "                                #
        "-XA InbreedingCoeff "                          #
        "-XA MappingQualityZero "                       #
        "-XA HaplotypeScore "                           #
        "-XA SpanningDeletions "                        #
        "-XA FisherStrand "                             #
        "-XA StrandOddsRatio "                          #
        "-XA ChromosomeCounts "                         #
        "-XA BaseQualityRankSumTest "                   #
        "-XA MappingQualityRankSumTest "                #
        "-XA QualByDepth "                              #
        "-XA ReadPosRankSumTest "                       #

###############################################################################
rule samtools_flagstat:
    input:
        "results/04_Variants/realigned/merged_realign_fixed.bam",
    output:
        "results/05_Validation/realigned/merged_realign_fixed_bam.flagstat",
    log:
        "results/11_Reports/samtools/flagstat/merged_realign_fixed_bam.log",
    params:
        extra="",  # optional params string
    wrapper:
        "v1.21.1/bio/samtools/flagstat"

###############################################################################
rule samtools_idxstats:
    # Aim: samtools idxstats – reports alignment summary statistics
    #       Retrieve and print stats in the index file corresponding to the input file. Before calling idxstats, the input BAM file should be indexed by samtools index.
    #       If run on a SAM or CRAM file or an unindexed BAM file, this command will still produce the same summary statistics, but does so by reading through the entire file. 
    #       This is far slower than using the BAM indices.
    #       The output is TAB-delimited with each line consisting of reference sequence name, sequence length, # mapped read-segments and # unmapped read-segments. 
    #       It is written to stdout. Note this may count reads multiple times if they are mapped more than once or in multiple fragments.
    input:
        bam="results/04_Variants/realigned/merged_realign_fixed.bam",
        idx="results/04_Variants/realigned/merged_realign_fixed.bai",
    output:
        "results/05_Validation/realigned/merged_realign_fixed_bam.idxstats",
    log:
        "results/11_Reports/samtools/idxstats/merged_realign_fixed_bam.log",
    params:
        extra="",  # optional params string
    wrapper:
        "v1.21.1/bio/samtools/idxstats"

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
        bam = "results/04_Variants/realigned/merged_realign_fixed.bam",
        refpath = "resources/genomes/GCA_018104305.1_AalbF3_genomic.fasta"
    output:
        stats = "results/05_Validation/realigned/merged_realign_fixed_stats.txt"
    log:
        "results/11_Reports/samtools/merged_realign_fixed_stats.log"
    shell:
        """
        samtools stats --threads {resources.cpus} -r {input.refpath} {input.bam} 1> {output.stats} 2> {log}
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
        bam = "results/04_Variants/realigned/merged_realign_fixed.bam",
    output:
        check = "results/05_Validation/validatesamfile/merged_realign_fixed_ValidateSam.txt"
    threads: CPUS
    log:
        "results/11_Reports/validatesamfiles/merged_realign_fixed-mate_sorted_validate_bam.log"
    shell:
        """
        set +e
        picard ValidateSamFile -I {input.bam} -M SUMMARY
        exitcode=$?
        if [ $exitcode -eq 1 ]
        then
            exit 1
        else
            exit 0
        fi
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
        bam = "results/04_Variants/realigned/merged_realign.bam"
    output:
        report = directory("results/00_Quality_Control/qualimap/")
    threads: CPUS
    log:
        "results/11_Reports/qualimap/qualimap.log"
    shell:
        """
        qualimap bamqc -bam {input.bam} -c -nt {threads} -outdir {output.report} -outformat HTML -sd
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
        fixedbam = "results/04_Variants/realigned/merged_realign_fixed.bam"
    output:
        index = "results/04_Variants/realigned/merged_realign_fixed.bai"
    log:
        "results/11_Reports/samtools/merged_realign_fixed_indexed.log"
    threads: 
        CPUS
    shell:
        """
        samtools index -@ {resources.cpus} -b {input.fixedbam} {output.index} &> {log}
        """

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
        realigned = "results/04_Variants/realigned/merged_realign.bam",
    output:
        fixed = "results/04_Variants/realigned/merged_realign_fixed.bam"
    threads: CPUS
    log:
        "results/11_Reports/fixmateinformation/merged_realign_fixed.log"
    shell:
        """
        picard FixMateInformation -I {input.realigned} -O {output.fixed} --ADD_MATE_CIGAR true
        """

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
        bam="results/02_Mapping/merged-mark-dup.bam",
        bai="results/02_Mapping/merged-mark-dup.bai",
        ref="resources/genomes/GCA_018104305.1_AalbF3_genomic.fasta",
        fai="resources/genomes/GCA_018104305.1_AalbF3_genomic.fasta.fai",
        dict="resources/genomes/GCA_018104305.1_AalbF3_genomic.dict",
        target_intervals="results/04_Variants/merged.intervals"
    output:
        bam="results/04_Variants/realigned/merged_realign.bam",
        bai="results/04_Variants/realigned/merged_realign.bai",
    benchmark:
        "benchmarks/indelrealigner/merged.tsv",
    log:
        "results/11_Reports/indelrealigner/merged.log",
    params:
        extra="--defaultBaseQualities 20 --filter_reads_with_N_cigar",  # optional
    threads: CPUS
    resources:
        mem_mb=16000,
    wrapper:
        "v1.19.1/bio/gatk3/indelrealigner"

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
        intervals="results/04_Variants/merged.intervals"
    params:
        cmd = r"""'BEGIN { OFS = "\t" } { if( $3 == "") { print $1, $2-1, $2 } else { print $1, $2-1, $3}}'"""
    output:
        bed = "results/04_Variants/merged_realignertargetcreator.bed"
    log:
        "results/11_Reports/awk/merged_intervals_for_IGV.log"
    threads: 
        CPUS
    shell:
        "awk -F '[:-]' "                # Awk, a program that you can use to select particular records in a file and perform operations upon them
        "{params.cmd} "                 # {AWK_CMD_INTERVALS:q} :q : is asking snakemake to quote the awk command for me. 
        "{input.intervals} "            # Intervals input
        "1> {output.bed} "              # BedGraph output
        "2> {log} "                     # Log redirection

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
        bam="results/02_Mapping/merged-mark-dup.bam",
        bai="results/02_Mapping/merged-mark-dup.bai",
        ref="resources/genomes/GCA_018104305.1_AalbF3_genomic.fasta",
        fai="resources/genomes/GCA_018104305.1_AalbF3_genomic.fasta.fai",
        dict="resources/genomes/GCA_018104305.1_AalbF3_genomic.dict",
    output:
        intervals=temp("results/04_Variants/merged.intervals"),
    benchmark:
        "benchmarks/realignertargetcreator/merged.tsv"
    log:
        "results/11_Reports/realignertargetcreator/merged.log",
    params:
        extra="--defaultBaseQualities 20 --filter_reads_with_N_cigar",  # optional
    resources:
        mem_mb=16000,
    threads: CPUS
    wrapper:
        "v1.19.1/bio/gatk3/realignertargetcreator"

###############################################################################
rule awk_coverage_stats:
    # Aim: computing genomme coverage stats
    # Use: awk {FORMULA} END {{print [RESULTS.tsv] [BEDGRAPH.bed]
    message:
        "Awk compute genome coverage statistics BED"
    conda:
        GAWK
    params:
        mincov = MINCOV
    input:
        genomecov = "results/03_Coverage/merged-genome-cov.bed"
    output:
        covstats = "results/03_Coverage/merged_coverage-stats.tsv"
    log:
        "results/11_Reports/awk/merged_coverage-stats.log"
    threads: 
        CPUS
    shell:
        "awk ' "                                  # Awk, a program that you can use to select particular records in a file and perform operations upon them
        "$4 >= {params.mincov} "                   # Minimum coverage
        "{{supMinCov+=$3-$2}} ; "                  # Genome size >= @ mincov X
        "{{genomeSize+=$3-$2}} ; "                 # Genome size
        "{{totalBases+=($3-$2)*$4}} ; "            # Total bases @ 1 X
        "{{totalBasesSq+=(($3-$2)*$4)**2}} "       # Total bases square @ 1 X
        "END "                                     # End
        "{{print "                                 # Print
        """ "sample_id", "\t", """                 # Sample ID header
        """ "mean_depth", "\t", """                # Mean depth header
        """ "standard_deviation", "\t", """        # Standard deviation header
        """ "cov_percent_@10X" """                 # Coverage percent @ mincov X header
        "ORS "                                     # \n newline
        """ "Merged", "\t", """                     # Sample ID value
        """ int(totalBases/genomeSize), "\t", """  # Mean depth value
        """ int(sqrt((totalBasesSq/genomeSize)-(totalBases/genomeSize)**2)), "\t", """ # Standard deviation value
        """ supMinCov/genomeSize*100 """           # Coverage percent @ mincov X value
        "}}' "                                     #
        "{input.genomecov} "                       # BedGraph coverage input
        "1> {output.covstats} "                    # Mean depth output
        "2> {log}"                                 # Log redirection

###############################################################################
rule bedtools_genome_coverage:
    # Aim: computing genome coverage sequencing
    # Use: bedtools genomecov [OPTIONS] -ibam [MARKDUP.bam] 1> [BEDGRAPH.bed]
    message:
        "BedTools computing genome coverage against reference genome sequence"
    conda:
        BEDTOOLS
    input:
        markdup = "results/02_Mapping/merged-mark-dup.bam",
        index = "results/02_Mapping/merged-mark-dup.bai"
    output:
        genomecov = "results/03_Coverage/merged-genome-cov.bed"
    log:
        "results/11_Reports/bedtools/merged-genome-cov.log"
    shell:
        """
        bedtools genomecov -bga -ibam {input.markdup} 1> {output.genomecov} 2> {log}
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
        markdup = "results/02_Mapping/merged-mark-dup.bam"
    output:
        index = "results/02_Mapping/merged-mark-dup.bai"
    log:
        "results/11_Reports/samtools/merged-mark-dup-index.log"
    shell:
        "samtools index "     # Samtools index, tools for alignments in the SAM format with command to index alignment
        "-@ {resources.cpus} " # --threads: Number of additional threads to use (default: 1)
        "-b "                  # -b: Generate BAI-format index for BAM files (default)
        "{input.markdup} "     # Markdup bam input
        "{output.index} "      # Markdup index output
        "&> {log}"             # Log redirection


###############################################################################
rule mark_duplicates_spark:
    input:
        "results/02_Mapping/merged.bam",
    output:
        bam = "results/02_Mapping/merged-mark-dup.bam",
        metrics="results/02_Mapping/merged-mark-dup_metrics.txt",
    log:
        "results/11_Reports/samtools/merged_picard-mark-dup.log",
    params:
        extra="--remove-sequencing-duplicates",  # optional
        java_opts="",  # optional
        #spark_runner="",  # optional, local by default
        #spark_v1.19.1="",  # optional
        #spark_extra="", # optional
    resources:
        mem_mb=16000,
    threads: CPUS
    wrapper:
        "v1.19.1/bio/gatk/markduplicatesspark"

###############################################################################
rule merge_bams:
    input:
        expand("results/02_Mapping/bams/{sample}_{aligner}_sorted.bam", sample=SAMPLE, aligner=ALIGNER),
    output:
        "results/02_Mapping/merged.bam",
    log:
        "results/11_Reports/logs/picard/mergesamfiles.log",
    params:
        extra="--VALIDATION_STRINGENCY LENIENT",
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=16000,
    wrapper:
        "v1.21.1/bio/picard/mergesamfiles"

###############################################################################
rule samtools_index:
    # Aim: indexing marked as duplicate BAM file for callable_loci rule
    # Use: samtools index -@ [THREADS] -b [MARKDUP.bam] [INDEX.bai] # -b: Generate BAI-format index for BAM files (default)
    message:
        "SamTools indexing realigned fixed BAM file for Picard ValidateSamFile"
    conda:
        SAMTOOLS
    resources:
       cpus = CPUS
    input:
        bam = "results/02_Mapping/{sample}_{aligner}_sorted.bam"
    output:
        index = "results/02_Mapping/{sample}_{aligner}_sorted.bai"
    log:
        "results/11_Reports/samtools/{sample}_{aligner}_sorted_index.log"
    threads: 
        CPUS
    shell:
        """
        samtools index -@ {resources.cpus} -b {input.bam} {output.index} &> {log}
        """

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
rule samtools_view:
    # Aim: Convert or filter SAM/BAM.
    # Use : samtools view -bo aln.bam aln.sam
    message:
        "Samtools view conversion of sample sam in bam format"
    input:
        "results/02_Mapping/{sample}_{aligner}-mapped.sam",
    output:
        bam="results/02_Mapping/bams/{sample}_{aligner}.bam",
    log:
        "results/11_Reports/samtools_view/{sample}_{aligner}.log",
    params:
        extra="",  # optional params string
        region="",  # optional region string
    threads: CPUS
    wrapper:
        "v1.21.1/bio/samtools/view"

###############################################################################
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
        bwapath = BWAPATH,
        reference = REFERENCE,
        extra = r"-R '@RG\tID:{sample}\tSM:{sample}\tCN:SC\tPL:ILLUMINA'" # Manage ReadGroup
    input:
        fwdreads = "results/01_Trimming/sickle/{sample}_sickle-trimmed_R1.fastq.gz",
        revreads = "results/01_Trimming/sickle/{sample}_sickle-trimmed_R2.fastq.gz"
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
        "{params.extra} "                                           #
        "{params.bwapath}{params.reference} "                       # Reference index filename prefix
        "{input.fwdreads} "                                         # Forward input reads
        "{input.revreads} "                                         # Reverse input reads
        "1> {output.mapped} "                                       # SAM output
        "2> {log}"                                                  # Log redirection

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
        fwdreads = "results/01_Trimming/sickle/{sample}_sickle-trimmed_R1.fastq.gz",
        revreads = "results/01_Trimming/sickle/{sample}_sickle-trimmed_R2.fastq.gz"
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
        truseq = TRUSEQ,
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
        "--cores {resources.cpus} "          # -j: Number of CPU cores to use. Use 0 to auto-detect (default: 1)
        "--trim-n "                          # --trim-n: Trim N's on ends of reads
        "--minimum-length {params.length} "  # -m: Discard reads shorter than length
        "--adapter {params.truseq} "         # -a: Sequence of an adapter ligated to the 3' end of the first read
        "-A {params.truseq} "                # -A: 3' adapter to be removed from second read in a pair
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
rule fastqc_quality_control:
    # Aim: reads sequence files and produces a quality control report
    # Use: fastqc [OPTIONS] --output [DIR/] [SAMPLE_1.fastq] ... [SAMPLE_n.fastq]
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
        "mkdir -p {output.fastqc} "     # (*) this directory must exist as the program will not create it
        "2> /dev/null && "              # in silence and then...
        "fastqc "                       # FastQC, a high throughput sequence QC analysis tool
        "--quiet "                      # -q: Supress all progress messages on stdout and only report errors
        "--threads {resources.cpus} "   # -t: Specifies files number which can be processed simultaneously
        "--outdir {output.fastqc} "     # -o: Create all output files in the specified output directory (*)
        "{input.fastq}/*.fastq.gz "     # Input file.fastq
        "&> {log}"                      # Log redirection

###############################################################################