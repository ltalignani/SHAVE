#######M######I#####V#####E######G######E######C######I#######R######D#########
# Name:                 shave.smk
# Author:               Loïc TALIGNANI
# Affiliation:          IRD_UMR224_MIVEGEC
# Aim:                  Snakefile for SHort-read Alignment pipeline for VEctor 
# Date:                 2022.10.05
# Run:                  snakemake --snakefile shave.smk --cores X --use-conda
# Latest modification:  2022.10.07
# Done:                 Added RealignerTargetCreator, IndelRealigner, 
#                       UnifiedGenotyper

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
    shell("mkdir -p Cluster_logs/bowtie2_mapping")

############################# O N S U C C E S S ##############################
onsuccess:
    shutil.rmtree(".snakemake")

################################## A L L #####################################
rule all:
    input:
        multiqc = "results/00_Quality_Control/multiqc/",
        fastqc = "results/00_Quality_Control/fastqc/",
        fastqscreen = "results/00_Quality_Control/fastq-screen/",
        hard_filter = expand("results/04_Variants/variantfiltration/{sample}_{aligner}_{markdup}_{mincov}X_hardfiltered.vcf",
                            sample = SAMPLE, aligner = ALIGNER, markdup = MARKDUP, mincov = MINCOV),
        vcf = expand("results/04_Variants/unifiedgenotyper/{sample}_{aligner}_{markdup}_{mincov}X_indels.vcf",
                            sample = SAMPLE, aligner = ALIGNER, markdup = MARKDUP, mincov = MINCOV),
        callable_loci = expand("results/05_Validation/callableloci/{sample}_{aligner}_{markdup}_{mincov}X_realign_fix-mate_sorted_callable_status.bed",
                            sample = SAMPLE, aligner = ALIGNER, markdup = MARKDUP, mincov = MINCOV),
        stats = expand("results/05_Validation/{sample}_{aligner}_{markdup}_{mincov}X_realign_fix-mate_sorted_stats.txt",
                            sample = SAMPLE, aligner = ALIGNER, markdup = MARKDUP, mincov = MINCOV),        
        igv_output = expand("results/04_Variants/{sample}_{aligner}_{markdup}_{mincov}X_realignertargetcreator.bed",
                            sample = SAMPLE, aligner = ALIGNER, markdup = MARKDUP, mincov = MINCOV),
        index_post_realign = expand("results/05_Validation/{sample}_{aligner}_{markdup}_{mincov}X_realign_fix-mate_sorted.bai",
                            sample = SAMPLE, aligner = ALIGNER, markdup = MARKDUP, mincov = MINCOV),
        realigned_fixed_bam = expand("results/05_Validation/{sample}_{aligner}_{markdup}_{mincov}X_realign_fix-mate_sorted.bam",
                            sample = SAMPLE, aligner = ALIGNER, markdup = MARKDUP, mincov = MINCOV),
        indelqual = expand("results/04_Variants/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual.bam",
                            sample = SAMPLE, aligner = ALIGNER, markdup = MARKDUP, mincov = MINCOV),
        covstats = expand("results/03_Coverage/{sample}_{aligner}_{markdup}_{mincov}X_coverage-stats.tsv",
                            sample = SAMPLE, aligner = ALIGNER, markdup = MARKDUP, mincov = MINCOV),
        markduplicate = expand("results/02_Mapping/{sample}_{aligner}_{markdup}-mark-dup.bam",
                            sample = SAMPLE, aligner = ALIGNER, markdup = MARKDUP)

################################ R U L E S ####################################
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
        "VariantFiltration Hard-filtering for {wildcards.sample} sample ({wildcards.aligner}-{wildcards.mincov})"
    conda:
        GATK4
    input:
        ref="resources/genomes/GCA_018104305.1_AalbF3_genomic.fasta",
        vcf="results/04_Variants/unifiedgenotyper/{sample}_{aligner}_{markdup}_{mincov}X_indels.vcf",
    output:
        vcf="results/04_Variants/variantfiltration/{sample}_{aligner}_{markdup}_{mincov}X_hardfiltered.vcf",
    params:
        filters={"myfilter": "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"},
        extra="",
        java_opts="",
    resources:
        mem_gb=MEM_GB
    log:
        "results/11_Reports/variantfiltration/{sample}_{aligner}_{markdup}_{mincov}X_hardfiltered.log",
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
        "UnifiedGenotyper calling SNVs for {wildcards.sample} sample ({wildcards.aligner}-{wildcards.mincov})"
    conda:
        GATK
    input:
        bam = "results/05_Validation/realigned/{sample}_{aligner}_{markdup}_{mincov}X_realign_fixed-mate_sorted.bam",
        ref = "resources/genomes/GCA_018104305.1_AalbF3_genomic.fasta",
        index = "results/05_Validation/realigned/{sample}_{aligner}_{markdup}_{mincov}X_realign_fixed-mate_sorted.bai"
        #alleles = ALLELES
    output:
        vcf="results/04_Variants/unifiedgenotyper/{sample}_{aligner}_{markdup}_{mincov}X_indels.vcf"
    log:
        "results/11_Reports/unifiedgenotyper/{sample}_{aligner}_{markdup}_{mincov}X.log"
    benchmark:
        "benchmarks/unifiedgenotyper/{sample}_{aligner}_{markdup}_{mincov}X.tsv"
    threads: CPUS
    shell:
        "gatk3 -T UnifiedGenotyper "                    # Genome Analysis Tool Kit - Broad Institute UnifiedGenotyper
        "-nct {threads} "                          # -nct / --num_cpu_threads_per_data_thread controls the number of CPU threads allocated to each data thread
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
rule callable_loci:
    # Aim: Collects statistics on callable, uncallable, poorly mapped, and other parts of the genome.
    # A very common question about a NGS set of reads is what areas of the genome are considered callable. This tool
    # considers the coverage at each locus and emits either a per base state or a summary interval BED file that
    # partitions the genomic intervals into the following callable states:
    # REF_N: The reference base was an N, which is not considered callable the GATK
    # PASS: The base satisfied the min. depth for calling but had less than maxDepth to avoid having EXCESSIVE_COVERAGE
    # NO_COVERAGE: Absolutely no reads were seen at this locus, regardless of the filtering parameters
    # LOW_COVERAGE: There were fewer than min. depth bases at the locus, after applying filters
    # EXCESSIVE_COVERAGE: More than -maxDepth read at the locus, indicating some sort of mapping problem
    # POOR_MAPPING_QUALITY: More than --maxFractionOfReadsWithLowMAPQ at the locus, indicating a poor mapping quality of the reads
    # Important : the bam file need to be indexed, and the index has to be in the same directory.
    # Use: gatk3 -T CallableLoci \
    #     -T CallableLoci \
    #     -R reference.fasta \
    #     -I myreads.bam \
    #     -summary table.txt \
    #     -o callable_status.bed
    message:
        "GATK3 CallableLoci for {wildcards.sample} sample ({wildcards.aligner}"
    conda:
        GATK
    input:
        refpath = "resources/genomes/GCA_018104305.1_AalbF3_genomic.fasta",
        sort = "results/05_Validation/realigned/{sample}_{aligner}_{markdup}_{mincov}X_realign_fixed-mate_sorted.bam",
        index = "results/05_Validation/realigned/{sample}_{aligner}_{markdup}_{mincov}X_realign_fixed-mate_sorted.bai",
    output:
        call = "results/05_Validation/callableloci/{sample}_{aligner}_{markdup}_{mincov}X_realign_fixed-mate_sorted_callable_status.bed",
        summary = "results/05_Validation/callableloci/{sample}_{aligner}_{markdup}_{mincov}X_summary_table.txt"
    threads: 
        CPUS
    log :
        "results/11_reports/callableloci/{sample}_{aligner}_{markdup}_{mincov}X_realign_fixed-mate_sorted_callable_status.log"
    shell:
        "gatk3 -T CallableLoci -R {input.refpath} -I {input.sort} -summary {output.summary} -o {output.call}" #  > {log} 2>&1

###############################################################################
rule validate_sam:
    # Aim: Basic check for bam file validity, as interpreted by the Broad Institute.
    # Use: picard.jar ValidateSamFile \
    #      -I input.bam \
    #      -MODE SUMMARY
    message:
        "Picard ValidateSamFile for {wildcards.sample} sample ({wildcards.aligner})"
    conda:
        PICARD
    input:
        sorted = "results/05_Validation/realigned/{sample}_{aligner}_{markdup}_{mincov}X_realign_fixed-mate_sorted.bam",
        index = "results/05_Validation/realigned/{sample}_{aligner}_{markdup}_{mincov}X_realign_fixed-mate_sorted.bai",
        refpath = "resources/genomes/GCA_018104305.1_AalbF3_genomic.fasta",
    output:
        check = "results/05_Validation/validatesamfile/{sample}_{aligner}_{markdup}_{mincov}X_realign_fixed-mate_sorted_validate_bam.txt"
    threads: CPUS
    log:
        "results/11_Reports/validatesamfiles/{sample}_{aligner}_{markdup}_{mincov}X_realign_fixed-mate_sorted_validate_bam.log"
    shell:
        "scripts/validatesamfile.sh"

###############################################################################
rule samtools_stats:
    # Aim: Collects statistics from BAM files
    # Use: samtools stats -r ref.fa input.bam # -r: Reference sequence (required for GC-depth and mismatches-per-cycle calculation).
    message:
        "SamTools indexing marked as duplicate BAM file {wildcards.sample} sample ({wildcards.aligner})"
    conda:
        SAMTOOLS
    resources:
       cpus = CPUS
    input:
        bam = "results/05_Validation/realigned/{sample}_{aligner}_{markdup}_{mincov}X_realign_fixed-mate_sorted.bam",
        refpath = "resources/genomes/GCA_018104305.1_AalbF3_genomic.fasta"
    output:
        stats = "results/05_Validation/realigned/{sample}_{aligner}_{markdup}_{mincov}X_realign_fixed-mate_sorted_stats.txt"
    log:
        "results/11_Reports/samtools/{sample}_{aligner}_{markdup}_{mincov}X_realign_fixed-mate_sorted_stats.log"
    shell:
        """
        samtools stats --threads {resources.cpus} -r {input.refpath} {input.bam} 1> {output.stats} 2> {log}
        """

###############################################################################
rule samtools_index_post_realign:
    # Aim: indexing marked as duplicate BAM file
    # Use: samtools index -@ [THREADS] -b [MARKDUP.bam] [INDEX.bai] # -b: Generate BAI-format index for BAM files (default)
    message:
        "SamTools indexing realigned fixed-mate sorted BAM file {wildcards.sample} sample ({wildcards.aligner}) for Picard ValidateSamFile"
    conda:
        SAMTOOLS
    resources:
       cpus = CPUS
    input:
        sortedbam = "results/05_Validation/realigned/{sample}_{aligner}_{markdup}_{mincov}X_realign_fixed-mate_sorted.bam"
    output:
        index = "results/05_Validation/realigned/{sample}_{aligner}_{markdup}_{mincov}X_realign_fixed-mate_sorted.bai"
    log:
        "results/11_Reports/samtools/{sample}_{aligner}_{markdup}_{mincov}X_realign_fixed-mate_sorted-index.log"
    threads: 
        CPUS
    shell:
        """
        samtools index -@ {resources.cpus} -b {input.sortedbam} {output.index} &> {log}
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
        "Picard FixMateInformation for {wildcards.sample} sample ({wildcards.aligner})"
    conda:
        PICARD
    input:
        sorted = "results/05_Validation/realigned/{sample}_{aligner}_{markdup}_{mincov}X_realign_fix-mate_sorted.bam",
    output:
        check = "results/05_Validation/realigned/{sample}_{aligner}_{markdup}_{mincov}X_realign_fixed-mate_sorted.bam"
    threads: CPUS
    log:
        "results/11_Reports/fixmateinformation/{sample}_{aligner}_{markdup}_{mincov}X_realign_fixed-mate_sorted_validate_bam.log"
    shell:
        """
        picard FixMateInformation -I {input.sorted} -O {output.check} --ADD_MATE_CIGAR true
        """

###############################################################################
rule samtools_sort_post_realign:
    # Aim: sorting
    # Use: samtools sort -@ [THREADS] -m [MEM] -T [TMPDIR] -O BAM -o [SORTED.bam] [FIXMATE.bam] # -T: Write temporary files to PREFIX.nnnn.bam
    message:
        "SamTools sorting {wildcards.sample} sample reads ({wildcards.aligner})"
    conda:
        SAMTOOLS
    resources:
       cpus = CPUS,
       mem_gb = MEM_GB
    params:
        tmpdir = TMPDIR
    input:
        fixmate = "results/05_Validation/realigned/{sample}_{aligner}_{markdup}_{mincov}X_realign_fix-mate.bam"
    output:
        sorted = "results/05_Validation/realigned/{sample}_{aligner}_{markdup}_{mincov}X_realign_fix-mate_sorted.bam"
    log:
        "results/11_Reports/samtools/{sample}_{aligner}_{markdup}_{mincov}X_realign_fix-mate_sorted.log"
    shell:
        "samtools sort "
        "--threads {resources.cpus} "
        "-m {resources.mem_gb}G "
        "-T {params.tmpdir} "
        "--output-fmt BAM "
        "-o {output.sorted} " 
        "{input.fixmate} "
        "&> {log}"

###############################################################################
rule samtools_fixmate_post_realign:
    # Aim: filling in mate coordinates
    # Use: samtools fixmate -@ [THREADS] -m -O BAM [SORTBYNAMES.bam] [FIXMATE.bam]
    message:
        "SamTools filling in mate coordinates {wildcards.sample} sample reads ({wildcards.aligner})"
    conda:
        SAMTOOLS
    resources:
        cpus = CPUS
    input:
        realignedbynames = "results/05_Validation/realigned/{sample}_{aligner}_{markdup}_{mincov}X_realign_sort-by-names.bam"
    output:
        fixmate = temp("results/05_Validation/realigned/{sample}_{aligner}_{markdup}_{mincov}X_realign_fix-mate.bam")
    log:
        "results/11_Reports/samtools/{sample}_{aligner}_{markdup}_{mincov}X_fixmate.log"
    shell:
        """
        samtools fixmate --threads {resources.cpus} -m --output-fmt BAM {input.realignedbynames} {output.fixmate} &> {log} 
        """

###############################################################################
rule samtools_sortbynames_post_realign:
    # Aim: sorting by names
    # Use: samtools sort -t [THREADS] -n -O BAM -o [SORTBYNAMES.bam] [MAPPED.sam] # -n: Sort by read name (not compatible with samtools index command)
    message:
        "SamTools sorting by names {wildcards.sample} sample reads ({wildcards.aligner})"
    conda:
        SAMTOOLS
    resources:
        cpus = CPUS,
        mem_gb = MEM_GB
    input:
        realigned_bam = "results/05_Validation/realigned/{sample}_{aligner}_{markdup}_{mincov}X_realign.bam",
    output:
        realignedbynames = temp("results/05_Validation/realigned/{sample}_{aligner}_{markdup}_{mincov}X_realign_sort-by-names.bam")
    log:
        "results/11_Reports/samtools/{sample}_{aligner}_{markdup}_{mincov}X_realign_sort-by-names.log"
    shell:
        """
        samtools sort --threads {resources.cpus} -m {resources.mem_gb}G -n --output-fmt BAM -o {output.realignedbynames} {input.realigned_bam} &> {log}
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
        "Awk IGV intervals visualization for {wildcards.sample} sample ({wildcards.aligner}-{wildcards.mincov})"
    input:
        bam="results/04_Variants/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual.bam",
        bai="results/04_Variants/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual.bai",
        ref="resources/genomes/GCA_018104305.1_AalbF3_genomic.fasta",
        fai="resources/genomes/GCA_018104305.1_AalbF3_genomic.fasta.fai",
        dict="resources/genomes/GCA_018104305.1_AalbF3_genomic.dict",
        target_intervals="results/04_Variants/{sample}_{aligner}_{markdup}_{mincov}X.intervals"
    output:
        bam="results/05_Validation/realigned/{sample}_{aligner}_{markdup}_{mincov}X_realign.bam",
        bai="results/05_Validation/realigned/{sample}_{aligner}_{markdup}_{mincov}X_realign.bai",
    benchmark:
        "benchmarks/indelrealigner/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual.tsv",
    log:
        "results/11_Reports/indelrealigner/{sample}_{aligner}_{markdup}_{mincov}X.log",
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
        "Awk IGV intervals visualization for {wildcards.sample} sample ({wildcards.aligner}-{wildcards.mincov})"
    conda:
        GAWK
    input:
        intervals="results/04_Variants/{sample}_{aligner}_{markdup}_{mincov}X.intervals"
    params:
        cmd = r"""'BEGIN { OFS = "\t" } { if( $3 == "") { print $1, $2-1, $2 } else { print $1, $2-1, $3}}'"""
    output:
        bed = "results/04_Variants/{sample}_{aligner}_{markdup}_{mincov}X_realignertargetcreator.bed"
    log:
        "results/11_Reports/awk/{sample}_{aligner}_{markdup}_{mincov}X_min-cov-filt.log"
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
        "RealignerTargetCreator creates a target intervals file for {wildcards.sample} sample ({wildcards.aligner}-{wildcards.mincov})"
    input:
        bam="results/04_Variants/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual.bam",
        bai="results/04_Variants/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual.bai",
        ref="resources/genomes/GCA_018104305.1_AalbF3_genomic.fasta",
        fai="resources/genomes/GCA_018104305.1_AalbF3_genomic.fasta.fai",
        dict="resources/genomes/GCA_018104305.1_AalbF3_genomic.dict",
    output:
        intervals=temp("results/04_Variants/{sample}_{aligner}_{markdup}_{mincov}X.intervals"),
    benchmark:
        "benchmarks/realignertargetcreator/{sample}_{aligner}_{markdup}_{mincov}X.tsv"
    log:
        "results/11_Reports/realignertargetcreator/{sample}_{aligner}_{markdup}_{mincov}X.log",
    params:
        extra="--defaultBaseQualities 20 --filter_reads_with_N_cigar",  # optional
    resources:
        mem_mb=16000,
    threads: CPUS
    wrapper:
        "v1.19.1/bio/gatk3/realignertargetcreator"

###############################################################################
rule samtools_indel_indexing:
    # Aim: indexing indel qualities BAM file
    # Use: samtools index -@ [THREADS] -b [INDELQUAL.bam] [INDEX.bai]
    message:
        "SamTools indexing indel qualities BAM file {wildcards.sample} sample ({wildcards.aligner}-{wildcards.mincov})"
    conda:
        SAMTOOLS
    resources:
       cpus = CPUS
    input:
        indelqual = "results/04_Variants/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual.bam"
    output:
        index = "results/04_Variants/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual.bai"
    log:
        "results/11_Reports/samtools/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual-index.log"
    threads: 
        CPUS
    shell:
        """
        samtools index -@ {resources.cpus} -b {input.indelqual} {output.index} &> {log}
        """

###############################################################################
rule lofreq_indel_qualities:
    # Aim: Indels qualities. Can be used instead of GATK’s BQSR or on non-Illumina data. 
    #      If you have Illumina data and don’t want to use GATK’s BQSR then the easiest thing is to use the --dindel option.
    # Use: lofreq indelqual --dindel -f [MASKEDREF.fasta] -o [INDEL.bam] [MARKDUP.bam]
    # Note: do not realign your BAM file afterwards!
    message:
        "LoFreq insert indels qualities for {wildcards.sample} sample ({wildcards.aligner}-{wildcards.mincov})"
    conda:
        LOFREQ
    input:
        maskedref = "results/04_Variants/{sample}_{aligner}_{markdup}_{mincov}X_masked-ref.fasta",
        markdup = "results/02_Mapping/{sample}_{aligner}_{markdup}-mark-dup.bam"
    output:
        indelqual = "results/04_Variants/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual.bam"
    log:
        "results/11_Reports/lofreq/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual.log"
    threads: 
        CPUS
    shell:
        """
        lofreq indelqual --dindel --ref {input.maskedref} --out {output.indelqual} {input.markdup} &> {log}
        """

###############################################################################
rule bedtools_masking:
    # Aim: masking low coverage regions
    # Use: bedtools maskfasta [OPTIONS] -fi [REFERENCE.fasta] -bed [RANGE.bed] -fo [MASKEDREF.fasta]
    message:
        "BedTools masking low coverage regions for {wildcards.sample} sample ({wildcards.aligner}-{wildcards.mincov})"
    conda:
        BEDTOOLS
    params:
        path = "resources/genomes/GCA_018104305.1_AalbF3_genomic.fasta"
    input:
        lowcovmask = "results/03_Coverage/{sample}_{aligner}_{markdup}_{mincov}X_low-cov-mask.bed"
    output:
        maskedref = "results/04_Variants/{sample}_{aligner}_{markdup}_{mincov}X_masked-ref.fasta"
    log:
        "results/11_Reports/bedtools/{sample}_{aligner}_{markdup}_{mincov}X_masking.log"
    threads: 
        CPUS
    shell:
        """
        bedtools maskfasta -fi {params.path} -bed {input.lowcovmask} -fo {output.maskedref} &> {log}
        """

###############################################################################
rule bedtools_merged_mask:
    # Aim: merging overlaps
    # Use: bedtools merge [OPTIONS] -i [FILTERED.bed] -g [GENOME.fasta]
    message:
        "BedTools merging overlaps for {wildcards.sample} sample ({wildcards.aligner}-{wildcards.mincov})"
    conda:
        BEDTOOLS
    input:
        mincovfilt = "results/03_Coverage/{sample}_{aligner}_{markdup}_{mincov}X_min-cov-filt.bed"
    output:
        lowcovmask = temp("results/03_Coverage/{sample}_{aligner}_{markdup}_{mincov}X_low-cov-mask.bed")
    log:
        "results/11_Reports/bedtools/{sample}_{aligner}_{markdup}_{mincov}X_merging.log"
    threads: 
        CPUS
    shell:
        """
        bedtools merge -i {input.mincovfilt} 1> {output.lowcovmask} 2> {log}
        """

###############################################################################
rule awk_mincovfilt:
    # Aim: minimum coverage filtration
    # Use: awk '$4 < [MINCOV]' [BEDGRAPH.bed] 1> [FILTERED.bed]
    message:
        "Awk minimum coverage filtration for {wildcards.sample} sample ({wildcards.aligner}-{wildcards.mincov})"
    conda:
        GAWK
    params:
        mincov = MINCOV
    input:
        genomecov = "results/03_Coverage/{sample}_{aligner}_{markdup}-genome-cov.bed"
    output:
        mincovfilt = temp("results/03_Coverage/{sample}_{aligner}_{markdup}_{mincov}X_min-cov-filt.bed")
    log:
        "results/11_Reports/awk/{sample}_{aligner}_{markdup}_{mincov}X_min-cov-filt.log"
    threads: 
        CPUS
    shell:
        "awk "                      # Awk, a program that you can use to select particular records in a file and perform operations upon them
        "'$4 < {params.mincov}' "    # Minimum coverage for masking regions in consensus sequence
        "{input.genomecov} "         # BedGraph coverage input
        "1> {output.mincovfilt} "    # Minimum coverage filtered bed output
        "2> {log} "                  # Log redirection

###############################################################################
rule awk_coverage_stats:
    # Aim: computing genomme coverage stats
    # Use: awk {FORMULA} END {{print [RESULTS.tsv] [BEDGRAPH.bed]
    message:
        "Awk compute genome coverage statistics BED {wildcards.sample} sample ({wildcards.aligner}-{wildcards.mincov})"
    conda:
        GAWK
    params:
        mincov = MINCOV
    input:
        genomecov = "results/03_Coverage/{sample}_{aligner}_{markdup}-genome-cov.bed"
    output:
        covstats = "results/03_Coverage/{sample}_{aligner}_{markdup}_{mincov}X_coverage-stats.tsv"
    log:
        "results/11_Reports/awk/{sample}_{aligner}_{markdup}_{mincov}X_coverage-stats.log"
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
        """ "cov_percent_@{wildcards.mincov}X" """ # Coverage percent @ mincov X header
        "ORS "                                     # \n newline
        """ "{wildcards.sample}", "\t", """        # Sample ID value
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
        "BedTools computing genome coverage for {wildcards.sample} sample against reference genome sequence ({wildcards.aligner})"
    conda:
        BEDTOOLS
    input:
        markdup = "results/02_Mapping/{sample}_{aligner}_{markdup}-mark-dup.bam",
        index = "results/02_Mapping/{sample}_{aligner}_{markdup}-mark-dup.bai"
    output:
        genomecov = temp("results/03_Coverage/{sample}_{aligner}_{markdup}-genome-cov.bed")
    log:
        "results/11_Reports/bedtools/{sample}_{aligner}_{markdup}-genome-cov.log"
    shell:
        """
        bedtools genomecov -bga -ibam {input.markdup} 1> {output.genomecov} 2> {log}
        """

###############################################################################
rule samtools_index_markdup:
    # Aim: indexing marked as duplicate BAM file
    # Use: samtools index -@ [THREADS] -b [MARKDUP.bam] [INDEX.bai]
    message:
        "SamTools indexing marked as duplicate BAM file {wildcards.sample} sample ({wildcards.aligner})"
    conda:
        SAMTOOLS
    resources:
       cpus = CPUS
    input:
        markdup = "results/02_Mapping/{sample}_{aligner}_{markdup}-mark-dup.bam"
    output:
        index = "results/02_Mapping/{sample}_{aligner}_{markdup}-mark-dup.bai"
    log:
        "results/11_Reports/samtools/{sample}_{aligner}_{markdup}-mark-dup-index.log"
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
        calmd = "results/02_Mapping/{sample}_{aligner}_sorted_MD.bam",
    output:
        bam = "results/02_Mapping/{sample}_{aligner}_picard-mark-dup.bam",
        metrics="results/02_Mapping/{sample}_{aligner}_picard-markdup_metrics.txt",
    log:
        "results/11_Reports/samtools/{sample}_{aligner}_picard-mark-dup.log",
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
rule samtools_markdup:
    # Aim: marking duplicate alignments
    # Use: samtools markdup -@ [THREADS] -r -s -O BAM [SORTED.bam] [MARKDUP.bam]
    message:
        "SamTools marking duplicate alignments for {wildcards.sample} sample ({wildcards.aligner})"
    conda:
        SAMTOOLS
    resources:
       cpus = CPUS
    input:
        calmd = "results/02_Mapping/{sample}_{aligner}_sorted_MD.bam"
    output:
        markdup = "results/02_Mapping/{sample}_{aligner}_samtools-mark-dup.bam"
    log:
        "results/11_Reports/samtools/{sample}_{aligner}_samtools-mark-dup.log"
    shell:
        "samtools markdup "           # Samtools markdup, tools for alignments in the SAM format with command mark duplicates
        "--threads {resources.cpus} " # -@: Number of additional threads to use (default: 1)
        "-r "                         # -r: Remove duplicate reads
        "-s "                         # -s: Report stats
        "--output-fmt BAM "           # -O: Specify output format: SAM, BAM, CRAM (here, BAM format)
        "{input.calmd} "              # Sorted bam input
        "{output.markdup} "           # Markdup bam output
        "&> {log}"                    # Log redirection

###############################################################################
rule samtools_calmd:
    # Aim: Tool to generate the MD tag for SNP/indel calling w/o lookin at the reference.
    # It does this by carrying information about the reference that the read does not carry, for a particular alignment. 
    # A SNP’s alternate base is carried in the read, but without the MD tag or use of the alignment reference, 
    # it’s impossible to know what the reference base was. Thus, this information is carried in the MD tag.
    # Use: samtools calmd -@ [THREADS] [REF.FA] -b input.bam > output.bam 
    message:
        "SamTools calmd {wildcards.sample} sample reads ({wildcards.aligner})"
    conda:
        SAMTOOLS
    resources:
       cpus = CPUS,
    params:
        refpath = REFPATH,
        reference = REFERENCE
    input:
        sorted = "results/02_Mapping/{sample}_{aligner}_sorted.bam"
    output:
        calmd = temp("results/02_Mapping/{sample}_{aligner}_sorted_MD.bam")
    log:
        "results/11_Reports/samtools/{sample}_{aligner}_sorted.log"
    shell:
       "samtools calmd "                                                 
        "--threads {resources.cpus} "                                   # -@: Number of additional threads to use (default: 1)
        "-b "                                                           # Output compressed BAM
        "{input.sorted} "                                               # bam input (sorted)
        "{params.refpath}{params.reference} "                           # Reference index filename prefix        
        "1> {output.calmd} "                                            # Sorted bam output
        "2> {log}"                                                      # Log redirection

###############################################################################
rule samtools_sorting:
    # Aim: sorting
    # Use: samtools sort -@ [THREADS] -m [MEM] -T [TMPDIR] -O BAM -o [SORTED.bam] [FIXMATE.bam]
    message:
        "SamTools sorting {wildcards.sample} sample reads ({wildcards.aligner})"
    conda:
        SAMTOOLS
    resources:
       cpus = CPUS,
       mem_gb = MEM_GB
    params:
        tmpdir = TMPDIR
    input:
        fixmate = "results/02_Mapping/{sample}_{aligner}_fix-mate.bam"
    output:
        sorted = temp("results/02_Mapping/{sample}_{aligner}_sorted.bam")
    log:
        "results/11_Reports/samtools/{sample}_{aligner}_sorted.log"
    shell:
        "samtools sort "              # Samtools sort, tools for alignments in the SAM format with command to sort alignment file
        "--threads {resources.cpus} "  # -@: Number of additional threads to use (default: 1)
        "-m {resources.mem_gb}G "      # -m: Set maximum memory per thread, suffix K/M/G recognized (default: 768M)
        "-T {params.tmpdir} "          # -T: Write temporary files to PREFIX.nnnn.bam
        "--output-fmt BAM "            # -O: Specify output format: SAM, BAM, CRAM (here, BAM format)
        "-o {output.sorted} "          # Sorted bam output
        "{input.fixmate} "             # Fixmate bam input
        "&> {log}"                     # Log redirection

###############################################################################
rule samtools_fixmate:
    # Aim: This tool fills in mate coordinates, ISIZE and mate related flags fom a name-sorted or a name-collated alignment. 
    # This means that before using samtools-fixmate, samtools-sort with the option of sorting by name must be used.
    # Samtools-fixmate is usually the previous step to using samtools-markdup.
    # Use: samtools fixmate -@ [THREADS] -m -O BAM [SORTBYNAMES.bam] [FIXMATE.bam]
    message:
        "SamTools filling in mate coordinates {wildcards.sample} sample reads ({wildcards.aligner})"
    conda:
        SAMTOOLS
    resources:
        cpus = CPUS
    input:
        sortbynames = "results/02_Mapping/{sample}_{aligner}_sort-by-names.bam"
    output:
        fixmate = temp("results/02_Mapping/{sample}_{aligner}_fix-mate.bam")
    log:
        "results/11_Reports/samtools/{sample}_{aligner}_fixmate.log"
    shell:
        "samtools fixmate "           # Samtools fixmate, tools for alignments in the SAM format with command to fix mate information
        "--threads {resources.cpus} "  # -@: Number of additional threads to use (default: 1)
        "-m "                          # -m: Add mate score tag
        "--output-fmt BAM "            # -O: Specify output format: SAM, BAM, CRAM (here, BAM format)
        "{input.sortbynames} "         # Sortbynames bam input
        "{output.fixmate} "            # Fixmate bam output
        "&> {log}"                     # Log redirection

###############################################################################
rule samtools_sortbynames:
    # Aim: sorting by names
    # Use: samtools sort -t [THREADS] -n -O BAM -o [SORTBYNAMES.bam] [MAPPED.sam]
    message:
        "SamTools sorting by names {wildcards.sample} sample reads ({wildcards.aligner})"
    conda:
        SAMTOOLS
    resources:
        cpus = CPUS,
        mem_gb = MEM_GB
    input:
        mapped = "results/02_Mapping/{sample}_{aligner}-mapped.sam"
    output:
        sortbynames = temp("results/02_Mapping/{sample}_{aligner}_sort-by-names.bam")
    log:
        "results/11_Reports/samtools/{sample}_{aligner}_sort-by-names.log"
    shell:
        "samtools sort "              # Samtools sort, tools for alignments in the SAM format with command to sort alignment file
        "--threads {resources.cpus} "  # -@: Number of additional threads to use (default: 1)
        "-m {resources.mem_gb}G "      # -m: Set maximum memory per thread, suffix K/M/G recognized (default: 768M)
        "-n "                          # -n: Sort by read name (not compatible with samtools index command)
        "--output-fmt BAM "            # -O: Specify output format: SAM, BAM, CRAM (here, BAM format)
        "-o {output.sortbynames} "     # -o: Write final output to FILE rather than standard output
        "{input.mapped} "              # Mapped reads input
        "&> {log}"                     # Log redirection

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
        extra = r"-R '@RG\tID:{sample}\tSM:{sample}'"
    input:
        fwdreads = "results/01_Trimming/sickle/{sample}_sickle-trimmed_R1.fastq.gz",
        revreads = "results/01_Trimming/sickle/{sample}_sickle-trimmed_R2.fastq.gz"
    output:
        mapped = "results/02_Mapping/{sample}_bwa-mapped.sam"
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
        fastqscreen = "results/00_Quality_Control/fastq-screen/"
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
        "--no-ansi "                 # Disable coloured log
        "&> {log}"                   # Log redirection

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