#######M######I#####V#####E######G######E######C######I#######R######D#########
# Name:                 shave.smk
# Author:               Loïc TALIGNANI
# Affiliation:          IRD_UMR224_MIVEGEC
# Aim:                  Snakefile for SHort-read Alignment for VEctor pipeline
# Date:                 2022.10.05
# Run:                  snakemake --snakefile shave.smk --cores X --use-conda
# Latest modification:  2022.10.07
# Done:                 Added HaplotypeCaller, GenotypeGVCFs, VariantFiltration

###############################################################################
# PUBLICATIONS #

###############################################################################
# CONFIGURATION #
configfile: "config/config.yaml"

from snakemake.utils import min_version

min_version("5.18.0")

report: "../report/workflow.rst"

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
PICARD = config["conda"][OS]["picard"]              # Picard 2.18.7

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

BWAPATH = config["bwa"]["path"]                     # BWA path to indexes
BT2PATH = config["bowtie2"]["path"]                 # Bowtie2 path to indexes
SENSITIVITY = config["bowtie2"]["sensitivity"]      # Bowtie2 sensitivity preset

REFPATH = config["consensus"]["path"]               # Path to genomes references
REFERENCE = config["consensus"]["reference"]        # Genome reference sequence, in fasta format
ALLELES = config["alleles"]["known_sites"]          # Known allele sites, in VCF format 
MINCOV = config["consensus"]["mincov"]              # Minimum coverage, mask lower regions with 'N'
MINAF = config["consensus"]["minaf"]                # Minimum allele frequency allowed
IUPAC = config["consensus"]["iupac"]                # Output variants in the form of IUPAC ambiguity codes

###############################################################################
# FUNCTIONS #

###############################################################################
rule all:
    input:
        multiqc = "results/00_Quality_Control/multiqc/",
        index_archive = expand("results/04_Variants/{sample}_{aligner}_{mincov}X_variant-filt.bgz.tbi",
                            sample = SAMPLE, aligner = ALIGNER, mincov = MINCOV),
        archive = expand("results/04_Variants/{sample}_{aligner}_{mincov}X_variant-filt.vcf.bgz",   
                            sample = SAMPLE, aligner = ALIGNER, mincov = MINCOV),
        consensus = expand("results/06_Consensus/{sample}_{aligner}_{mincov}X_consensus.fasta",
                            sample = SAMPLE, aligner = ALIGNER, mincov = MINCOV),

###############################################################################
rule sed_rename_headers:
    # Aim: rename all fasta header with sample name
    # Use: sed 's/[OLD]/[NEW]/' [IN] > [OUT]
    message:
        "Sed rename header for {wildcards.sample} sample consensus fasta ({wildcards.aligner}-{wildcards.mincov})"
    input:
        constmp = "results/06_Consensus/{sample}_{aligner}_{mincov}X_consensus.fasta.tmp"
    output:
        consensus = "results/06_Consensus/{sample}_{aligner}_{mincov}X_consensus.fasta"
    log:
        "results/11_Reports/sed/{sample}_{aligner}_{mincov}X_fasta-header.log"
    shell:
        "sed " # Sed, a Stream EDitor used to perform basic text transformations on an input stream
        "'s/^>.*$/>{wildcards.sample}_{wildcards.aligner}_{wildcards.mincov}/' "
        "{input.constmp} "       # Input file
        "1> {output.consensus} " # Output file
        "2> {log}"               # Log redirection

###############################################################################
rule bcftools_consensus:
    # Aim: create consensus
    # Use: bcftools consensus -f [REFERENCE] [VARIANTS.vcf.gz] -o [CONSENSUS.fasta]
    message:
        "BcfTools consensus for {wildcards.sample} sample ({wildcards.aligner}-{wildcards.mincov})"
    conda:
        BCFTOOLS
    params:
        iupac = IUPAC
    input:
        maskedref = "results/04_Variants/{sample}_{aligner}_{mincov}X_masked-ref.fasta",
        archive = "results/04_Variants/{sample}_{aligner}_{mincov}X_variant-filt.vcf.bgz",
        index = "results/04_Variants/{sample}_{aligner}_{mincov}X_variant-filt.bgz.tbi"
    output:
        constmp = temp("results/06_Consensus/{sample}_{aligner}_{mincov}X_consensus.fasta.tmp")
    log:
        "results/11_Reports/bcftools/{sample}_{aligner}_{mincov}X_consensus.log"
    shell:
        "bcftools "                       # Bcftools, tools for variant calling and manipulating VCFs and BCFs
        "consensus "                      # Create consensus sequence by applying VCF variants to a reference fasta file
        "--fasta-ref {input.maskedref} "  # -f: reference sequence in fasta format
        "{params.iupac} "                 # -I, --iupac-codes: output variants in the form of IUPAC ambiguity codes
        "{input.archive} "                # SNVs and Indels filtered VCF archive file
        "--output {output.constmp} "      # -o: write output to a file (default: standard output)
        "2> {log}"                        # Log redirection

###############################################################################
rule tabix_tabarch_indexing:
    # Aim: tab archive indexing
    # Use: tabix [OPTIONS] [TAB.bgz]
    message:
        "Tabix tab archive indexing for {wildcards.sample} sample ({wildcards.aligner}-{wildcards.mincov})"
    conda:
        SAMTOOLS
    input:
        archive = "results/04_Variants/{sample}_{aligner}_{mincov}X_variant-filt.vcf.bgz"
    output:
        index = "results/04_Variants/{sample}_{aligner}_{mincov}X_variant-filt.bgz.tbi"
    log:
        "results/11_Reports/tabix/{sample}_{aligner}_{mincov}X_variant-archive-index.log"
    shell:
        "tabix "             # Tabix, indexes a TAB-delimited genome position file in.tab.bgz and creates an index file
        "{input.archive} "   # The input data file must be position sorted and compressed by bgzip
        "-f "                # overwrite index if already existing   
        "1> {output.index} " # Tabix output TBI index formats
        "2> {log}"           # Log redirection

###############################################################################
rule bgzip_variant_archive:
    # Aim: Variant block compressing
    # Use: bgzip [OPTIONS] -c -@ [THREADS] [INDEL.vcf] 1> [COMPRESS.vcf.bgz]
    message:
        "Bgzip variant block compressing for {wildcards.sample} sample ({wildcards.aligner}-{wildcards.mincov})"
    conda:
        SAMTOOLS
    resources:
        cpus = CPUS
    input:
        variantfilt = "results/04_Variants/variantfiltration/{sample}_{aligner}_{mincov}X_hardfiltered.vcf.gz"                  #results/04_Variants/lofreq/{sample}_{aligner}_{mincov}X_variant-filt.vcf"
    output:
        archive = "results/04_Variants/{sample}_{aligner}_{mincov}X_variant-filt.vcf.bgz"
    log:
        "results/11_Reports/bgzip/{sample}_{aligner}_{mincov}X_variant-archive.log"
    shell:
        "bgzip "                      # Bgzip, block compression/decompression utility
        "--stdout "                   # -c: Write to standard output, keep original files unchanged
        "--threads {resources.cpus} " # -@: Number of threads to use (default: 1)
        "{input.variantfilt} "        # VCF input file, gzip suuported, no streaming supported
        "1> {output.archive} "        # VCF output file, gzip supported (default: standard output)
        "2> {log}"                    # Log redirection

###############################################################################
rule hard_filter_calls:
    # Aim: Perform joint genotyping on one or more samples pre-called with HaplotypeCaller.
    # In any case, the input samples must possess genotype likelihoods produced by HaplotypeCaller with `-ERC GVCF` or `-ERC BP_RESOLUTION`.
    # Use: gatk --java-options "-Xmx4g" GenotypeGVCFs \
    # -R Homo_sapiens_assembly38.fasta \
    # -V input.g.vcf.gz \
    # -O output.vcf.gz
    message:
        "GenotypeGVCFs calling genotypes for {wildcards.sample} sample ({wildcards.aligner}-{wildcards.mincov})"
    input:
        ref=REFPATH,
        vcf="results/04_Variants/unifiedgenotyper/{sample}_{aligner}_{mincov}X_indels.vcf",
    output:
        vcf="results/04_Variants/variantfiltration/{sample}_{aligner}_{mincov}X_hardfiltered.vcf.gz",
    params:
        filters={"myfilter": "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"},
        extra="",
        java_opts="",
    resources:
        mem_gb=MEM_GB
    log:
        "results/11_Reports/variantfiltration/{sample}_{aligner}_{mincov}X_hardfiltered.log",
    wrapper:
        "0.74.0/bio/gatk/variantfiltration"

###############################################################################
rule unifiedgenotyper:
    # Aim:  Call variants in sequence data. The following parameters comes from the MalariaGEN
    # Use:  java -jar GenomeAnalysisTK.jar \ 
    #       -T UnifiedGenotyper \
    #       -I {sample BAM} \
    #       --alleles {alleles VCF} \
    #       -R {reference sequence} \
    #       --out {output VCF} \
    message:
        "UnifiedGenotyper calling SNVs for {wildcards.sample} sample ({wildcards.aligner}-{wildcards.mincov})"
    #conda:
    #    GATK
    input:
        bam = "results/05_Validation/realigned/{sample}_{aligner}_{mincov}X_realign_fix-mate.bam",
        ref = REFPATH,
        alleles = ALLELES
    output:
        vcf="results/04_Variants/unifiedgenotyper/{sample}_{aligner}_{mincov}X_indels.vcf"
    log:
        "results/11_Reports/unifiedgenotyper/{sample}_{aligner}_{mincov}X.log"
    shell:
        "gatk3 "                                        # Genome Analysis Tool Kit - Broad Institute
        "-T UnifiedGenotyper "                          # UnifiedGenotyper
        "-I {input.bam} "                               # Input indel realigned BAM file
        "--alleles {input.alleles} "                    # Alleles against which to genotype (VCF format). Given the sites VCF file is fixed for every sample, and we wish to generalise to future sets of sites/alleles, the VCF file describing sites and alleles should be considered a parameter. This file for A. gambiae (AgamP4) is available at
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