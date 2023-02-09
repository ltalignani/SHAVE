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
MULTIQC = config["conda"][OS]["multiqc"]            # MultiQC 1.14
SAMTOOLS = config["conda"][OS]["samtools"]          # SamTools
GAWK = config["conda"][OS]["gawk"]                  # Gawk
GATK = config["conda"][OS]["gatk"]                  # GATK 3.6
GATK4 = config["conda"][OS]["gatk4"]                # GATK 4.3.0
PICARD = config["conda"][OS]["picard"]              # Picard 2.27.5
QUALIMAP = config["conda"][OS]["qualimap"]          # Qualimap 2.2.0

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
        indexvcf = expand("results/05_Variants/{sample}_{aligner}.vcf.gz.tbi", sample=SAMPLE, aligner=ALIGNER),
        table = "results/05_Variants/merged_raw/merged.table",
        combinegvcfs = "results/05_Variants/merged_raw/merged.vcf.gz",
        bgzip_vcfs = expand("results/05_Variants/{sample}_{aligner}.vcf.gz", sample=SAMPLE, aligner=ALIGNER),
        vcf = expand("results/05_Variants/{sample}_{aligner}.vcf", sample=SAMPLE, aligner=ALIGNER),

################################ R U L E S ####################################
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
