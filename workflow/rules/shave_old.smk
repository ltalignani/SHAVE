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
MINCOV = config["consensus"]["mincov"]              # Minimum coverage, mask lower regions with 'N'
MINAF = config["consensus"]["minaf"]                # Minimum allele frequency allowed
IUPAC = config["consensus"]["iupac"]                # Output variants in the form of IUPAC ambiguity codes

#FILTERS = config["filtering"["hard"]]["snvs"]

###############################################################################
# FUNCTIONS #

###############################################################################
rule all:
    input:
        multiqc = "results/00_Quality_Control/multiqc/",
        covstats = expand("results/03_Coverage/{sample}_{aligner}_{mincov}X_coverage-stats.tsv",
                          sample = SAMPLE, aligner = ALIGNER, mincov = MINCOV),
        index_archive = expand("results/04_Variants/{sample}_{aligner}_{mincov}X_variant-filt.gz.tbi",
                           sample = SAMPLE, aligner = ALIGNER, mincov = MINCOV),        
        archive = expand("results/04_Variants/{sample}_{aligner}_{mincov}X_variant-filt.vcf.gz",
                           sample = SAMPLE, aligner = ALIGNER, mincov = MINCOV),       
        index = expand("results/04_Variants/{sample}_{aligner}_{mincov}X_indel-qual.bai",
                           sample = SAMPLE, aligner = ALIGNER, mincov = MINCOV),        
        stats = expand("results/05_Validation/{sample}_{aligner}_mark-dup.txt",
                           sample = SAMPLE, aligner = ALIGNER, mincov = MINCOV),
        check = expand("results/05_Validation/validatesamfile/{sample}_{aligner}_check_mark_dup_bam.txt",
                           sample = SAMPLE, aligner = ALIGNER, mincov = MINCOV)

###############################################################################
rule tabix_tabarch_indexing:
    # Aim: tab archive indexing
    # Use: tabix [OPTIONS] [TAB.bgz]
    message:
        "Tabix tab archive indexing for {wildcards.sample} sample ({wildcards.aligner}-{wildcards.mincov})"
    conda:
        SAMTOOLS
    input:
        archive = "results/04_Variants/{sample}_{aligner}_{mincov}X_variant-filt.vcf.gz"
    output:
        index = "results/04_Variants/{sample}_{aligner}_{mincov}X_variant-filt.gz.tbi"
    log:
        "results/11_Reports/tabix/{sample}_{aligner}_{mincov}X_variant-archive-index.log"
    shell:
        "tabix "             # Tabix, indexes a TAB-delimited genome position file in.tab.bgz and creates an index file
        "{input.archive} "   # The input data file must be position sorted and compressed by bgzip
        "-f "                # overwrite index if already existing   
        "1> {output.index} " # Tabix output TBI index formats
        "2> {log}"           # Log redirection

###############################################################################
rule bcftools_variant_filt_archive:
    # Aim: Variant block compressing
    # Use: bgzip [OPTIONS] -c -@ [THREADS] [INDEL.vcf] 1> [COMPRESS.vcf.bgz]
    message:
        "Bgzip variant block compressing for {wildcards.sample} sample ({wildcards.aligner}-{wildcards.mincov})"
    conda:
        BCFTOOLS
    resources:
        cpus = CPUS
    input:
        variantfilt = "results/04_Variants/variantfiltration/{sample}_{aligner}_{mincov}X_hardfiltered.vcf"                  #results/04_Variants/lofreq/{sample}_{aligner}_{mincov}X_variant-filt.vcf"
    output:
        archive = "results/04_Variants/{sample}_{aligner}_{mincov}X_variant-filt.vcf.gz"
    log:
        "results/11_Reports/bcftools/{sample}_{aligner}_{mincov}X_variant-archive.log"
    shell:
        "bcftools "                         # bcftools,  a set of utilities that manipulate variant calls in the Variant Call Format (VCF).
        "view "                             # view : subset, filter and convert VCF and BCF files
        "--threads {resources.cpus} "       # -@: Number of threads to use (default: 1)
        "{input.variantfilt} "              # VCF input file,
        "-Oz -o {output.archive} "          # -O[z|b]: output-type -o: VCF output file,
        "&> {log}"  

###############################################################################
rule hard_filter_calls:
    # Aim: Perform joint genotyping on one or more samples pre-called with HaplotypeCaller.
    # In any case, the input samples must possess genotype likelihoods produced by HaplotypeCaller with `-ERC GVCF` or `-ERC BP_RESOLUTION`.
    # Use: gatk --java-options "-Xmx4g" GenotypeGVCFs \
    # -R Homo_sapiens_assembly38.fasta \
    # -V input.g.vcf.gz \
    # -O output.vcf.gz
    message:
        "Hard-filtering for {wildcards.sample} sample ({wildcards.aligner}-{wildcards.mincov})"
    input:
        ref="resources/genomes/GCA_018104305.1_AalbF3_genomic.fasta",
        vcf="results/04_Variants/unifiedgenotyper/{sample}_{aligner}_{mincov}X_indels.vcf",
    output:
        vcf="results/04_Variants/variantfiltration/{sample}_{aligner}_{mincov}X_hardfiltered.vcf",
    params:
        filters={"myfilter": "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"},
        extra="",
        java_opts="",
    resources:
        mem_gb=MEM_GB
    threads:
        CPUS
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
    conda:
        GATK
    input:
        bam="results/05_Validation/realigned/{sample}_{aligner}_{mincov}X_realigned_fix-mate.bam",
        ref=REFERENCE,
        alleles=""
    output:
        vcf="results/04_Variants/unifiedgenotyper/{sample}_{aligner}_{mincov}X_indels.vcf"
    log:
        "results/11_Reports/unifiedgenotyper/{sample}_{aligner}_{mincov}X.log"
    params:
        extra=""  # optional
    threads: CPUS
    resources:
        mem_gb = MEM_GB
    shell:
        "gatk3 "                                        # Genome Analysis Tool Kit - Broad Institute
        "-T UnifiedGenotyper "                          # UnifiedGenotyper
        "-I {input.bam} "                               # Input indel realigned BAM file
        "--alleles {input.alleles} "                    # Alleles against which to genotype (VCF format). Given the sites VCF file is fixed for every sample, and we wish to generalise to future sets of sites/alleles, the VCF file describing sites and alleles should be considered a parameter. This file for A. gambiae (AgamP4) is available at
        "-R {input.ref} "                               # Reference sequence in fasta format
        "--out {output VCF} "                           # Output VCF
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
rule samtools_fixmate_after_realignment:
    # Aim: filling in mate coordinates
    # Use: samtools fixmate -@ [THREADS] -m -O BAM [SORTBYNAMES.bam] [FIXMATE.bam]
    message:
        "SamTools filling in mate coordinates {wildcards.sample} sample reads ({wildcards.aligner})"
    conda:
        SAMTOOLS
    resources:
        cpus = CPUS
    input:
        bam = "results/05_Validation/realigned/{sample}_{aligner}_{mincov}X_realigned.bam"
    output:
        fixmate = temp("results/05_Validation/{sample}_{aligner}_{mincov}X_realigned_fix-mate.bam")
    log:
        "results/11_Reports/samtools/{sample}_{aligner}_{mincov}X_fixmate.log"
    shell:
        "samtools fixmate "           # Samtools fixmate, tools for alignments in the SAM format with command to fix mate information
        "--threads {resources.cpus} "  # -@: Number of additional threads to use (default: 1)
        "-m "                          # -m: Add mate score tag
        "--output-fmt BAM "            # -O: Specify output format: SAM, BAM, CRAM (here, BAM format)
        "{input.bam} "         # Sortbynames bam input
        "{output.fixmate} "            # Fixmate bam output
        "&> {log}"                     # Log redirection

###############################################################################
rule indelrealigner:
    # Aim:  takes a coordinate-sorted and indexed BAM and a target intervals file generated by RealignerTargetCreator. 
    #       IndelRealigner then performs local realignment on reads coincident with the target intervals using consenses 
    #       from indels present in the original alignment.
    # Use:  java -Xmx8G -Djava.io.tmpdir=/tmp -jar GenomeAnalysisTK.jar \
    #       -T IndelRealigner \
    #       -R human_g1k_v37_decoy.fasta \
    #       -targetIntervals 7156_realignertargetcreator.intervals \
    #       -known INDEL_chr10_1Mb_b37_1000G_phase3_v4_20130502.vcf \ 
    #       -I 7156_snippet.bam \
    #       -o 7156_snippet_indelrealigner.bam
    #      
    message:
        "Awk IGV intervals visualization for {wildcards.sample} sample ({wildcards.aligner}-{wildcards.mincov})"
    input:
        bam="results/04_Variants/{sample}_{aligner}_{mincov}X_indel-qual.bam",
        bai="results/04_Variants/{sample}_{aligner}_{mincov}X_indel-qual.bai",
        ref=REFERENCE,
        known="",
        known_idx="",
        target_intervals="results/04_Variants/{sample}_{aligner}_{mincov}X.intervals"
    output:
        bam="results/05_Validation/realigned/{sample}_{aligner}_{mincov}X_realigned.bam",
        bai="results/05_Validation/realigned/{sample}_{aligner}_{mincov}X_realigned.bai",
        java_temp=temp(directory("results/04_Variants/indelrealigner/{sample}_{aligner}_{mincov}X")),
    log:
        "results/11_Reports/indelrealigner/{sample}_{aligner}_{mincov}X.log"
    params:
        extra=""  # optional
    threads: CPUS
    resources:
        mem_gb = MEM_GB
    wrapper:
        "bio/gatk/indelrealigner"

###############################################################################
rule awk_intervals:
    # Aim: View intervals on IGV
    # Use: awk -F '[:-]' 'BEGIN { OFS = "\t" } { if( $3 == "") { print $1, $2-1, $2 } else { print $1, $2-1, $3}}' \ 
    #      7156_realignertargetcreator.intervals > 7156_realignertargetcreator.bed
    message:
        "Awk IGV intervals visualization for {wildcards.sample} sample ({wildcards.aligner}-{wildcards.mincov})"
    conda:
        GAWK
    input:
        intervals="results/04_Variants/{sample}_{aligner}_{mincov}X.intervals"
    output:
        bed = "results/03_Coverage/{sample}_{aligner}_{mincov}X_realignertargetcreator.bed"
    log:
        "results/11_Reports/awk/{sample}_{aligner}_{mincov}X_min-cov-filt.log"
    shell:
        "awk "                          # Awk, a program that you can use to select particular records in a file and perform operations upon them
        "-F '[:-]' 'BEGIN { OFS = '\t' } { if( $3 == "") { print $1, $2-1, $2 } else { print $1, $2-1, $3}}' "
        "{input.intervals} "            # Intervals input
        "1> {output.bed} "              # BedGraph output
        "2> {log} "                     # Log redirection

###############################################################################
rule realignertargetcreator:
    # Aim: Local realignment around indels. Takes a coordinate-sorted and indexed BAM and a VCF of known indels and creates a target intervals file.
    # Use: java -jar GenomeAnalysisTK.jar \
    #           -T RealignerTargetCreator \
    #           -R human_g1k_v37_decoy.fasta \
    #           -L 10:96000000-97000000 \
    #           -known INDEL_chr10_1Mb_b37_1000G_phase3_v4_20130502.vcf \
    #           -I 7156_snippet.bam \
    #           -o 7156_realignertargetcreator.intervals
    message:
        "RealignerTargetCreator creates a target intervals file for {wildcards.sample} sample ({wildcards.aligner}-{wildcards.mincov})"
    input:
        bam="results/04_Variants/{sample}_{aligner}_{mincov}X_indel-qual.bam",
        ref=REFERENCE,
        known="",
    output:
        intervals="results/04_Variants/{sample}_{aligner}_{mincov}X.intervals",
        java_temp=temp(directory("results/04_Variants/indelrealigner/{sample}_{aligner}_{mincov}X")),
    log:
        "results/11_Reports/realignertargetcreator/{sample}_{aligner}_{mincov}X.log",
    params:
        extra="", # optional
    resources:
        mem_gb= MEM_GB,
    threads: CPUS
    wrapper:
        "bio/gatk3/realignertargetcreator"

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
        indelqual = "results/04_Variants/{sample}_{aligner}_{mincov}X_indel-qual.bam"
    output:
        index = "results/04_Variants/{sample}_{aligner}_{mincov}X_indel-qual.bai"
    log:
        "results/11_Reports/samtools/{sample}_{aligner}_{mincov}X_indel-qual-index.log"
    shell:
        "samtools index "      # Samtools index, tools for alignments in the SAM format with command to index alignment
        "-@ {resources.cpus} " # Number of additional threads to use (default: 0)
        "-b "                  # -b: Generate BAI-format index for BAM files (default)
        "{input.indelqual} "   # Sorted bam input
        "{output.index} "      # Markdup bam output
        "&> {log}"             # Log redirection

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
        maskedref = "results/04_Variants/{sample}_{aligner}_{mincov}X_masked-ref.fasta",
        markdup = "results/02_Mapping/{sample}_{aligner}_mark-dup.bam"
    output:
        indelqual = "results/04_Variants/{sample}_{aligner}_{mincov}X_indel-qual.bam"
    log:
        "results/11_Reports/lofreq/{sample}_{aligner}_{mincov}X_indel-qual.log"
    shell:
        "lofreq "                   # LoFreq, fast and sensitive inference of SNVs and Indels
        "indelqual "                # Insert indel qualities into BAM file (required for indel predictions)
        "--dindel "                 # Add Dindel's indel qualities Illumina specifics (need --ref and clashes with -u)
        "--ref {input.maskedref} "  # -f: Reference (masked) sequence used for mapping (only required for --dindel)
        "--out {output.indelqual} " # -o: Indel BAM file output (default: standard output)
        "{input.markdup} "          # Markdup BAM input
        "&> {log}"                  # Log redirection

###############################################################################
rule bedtools_masking:
    # Aim: masking low coverage regions
    # Use: bedtools maskfasta [OPTIONS] -fi [REFERENCE.fasta] -bed [RANGE.bed] -fo [MASKEDREF.fasta]
    message:
        "BedTools masking low coverage regions for {wildcards.sample} sample ({wildcards.aligner}-{wildcards.mincov})"
    conda:
        BEDTOOLS
    params:
        path = REFPATH,
        reference = REFERENCE
    input:
        lowcovmask = "results/03_Coverage/{sample}_{aligner}_{mincov}X_low-cov-mask.bed"
    output:
        maskedref = "results/04_Variants/{sample}_{aligner}_{mincov}X_masked-ref.fasta"
    log:
        "results/11_Reports/bedtools/{sample}_{aligner}_{mincov}X_masking.log"
    shell:
        "bedtools maskfasta "                        # Bedtools maskfasta, mask a fasta file based on feature coordinates
        "-fi {params.path}{params.reference}.fasta " # Input FASTA file
        "-bed {input.lowcovmask} "                   # BED/GFF/VCF file of ranges to mask in -fi
        "-fo {output.maskedref} "                    # Output masked FASTA file
        "&> {log}"                                   # Log redirection

###############################################################################
rule bedtools_merged_mask:
    # Aim: merging overlaps
    # Use: bedtools merge [OPTIONS] -i [FILTERED.bed] -g [GENOME.fasta]
    message:
        "BedTools merging overlaps for {wildcards.sample} sample ({wildcards.aligner}-{wildcards.mincov})"
    conda:
        BEDTOOLS
    input:
        mincovfilt = "results/03_Coverage/{sample}_{aligner}_{mincov}X_min-cov-filt.bed"
    output:
        lowcovmask = temp("results/03_Coverage/{sample}_{aligner}_{mincov}X_low-cov-mask.bed")
    log:
        "results/11_Reports/bedtools/{sample}_{aligner}_{mincov}X_merging.log"
    shell:
        "bedtools merge "        # Bedtools merge, merges overlapping BED/GFF/VCF entries into a single interval
        "-i {input.mincovfilt} "  # -i: BED/GFF/VCF input to merge
        "1> {output.lowcovmask} " # merged output
        "2> {log}"                # Log redirection

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
        genomecov = "results/03_Coverage/{sample}_{aligner}_genome-cov.bed"
    output:
        mincovfilt = temp("results/03_Coverage/{sample}_{aligner}_{mincov}X_min-cov-filt.bed")
    log:
        "results/11_Reports/awk/{sample}_{aligner}_{mincov}X_min-cov-filt.log"
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
        genomecov = "results/03_Coverage/{sample}_{aligner}_genome-cov.bed"
    output:
        covstats = "results/03_Coverage/{sample}_{aligner}_{mincov}X_coverage-stats.tsv"
    log:
        "results/11_Reports/awk/{sample}_{aligner}_{mincov}X_coverage-stats.log"
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
        markdup = "results/02_Mapping/{sample}_{aligner}_mark-dup.bam",
        index = "results/02_Mapping/{sample}_{aligner}_mark-dup.bai"
    output:
        genomecov = temp("results/03_Coverage/{sample}_{aligner}_genome-cov.bed")
    log:
        "results/11_Reports/bedtools/{sample}_{aligner}_genome-cov.log"
    shell:
        "bedtools genomecov "    # Bedtools genomecov, compute the coverage of a feature file among a genome
        "-bga "                   # Report depth in BedGraph format, regions with zero coverage are also reported
        "-ibam {input.markdup} "  # The input file is in BAM format, must be sorted by position
        "1> {output.genomecov} "  # BedGraph output
        "2> {log} "               # Log redirection

###############################################################################
rule samtools_stats:
    # Aim: Collects statistics from BAM files
    # Use: samtools stats -r ref.fa input.bam
    message:
        "SamTools indexing marked as duplicate BAM file {wildcards.sample} sample ({wildcards.aligner})"
    conda:
        SAMTOOLS
    resources:
       cpus = CPUS
    input:
        bam = "results/02_Mapping/{sample}_{aligner}_mark-dup.bam",
        refpath = "resources/genomes/GCA_018104305.1_AalbF3_genomic.fasta"
    output:
        stats = "results/05_Validation/{sample}_{aligner}_mark-dup.txt"
    log:
        "results/11_Reports/samtools/{sample}_{aligner}_mark-dup.log"
    shell:
        "samtools stats "               # Samtools stats, collects statistics from BAM files. The output can be visualized using plot-bamstats.
        "--threads {resources.cpus} "   # -@: Number of additional threads to use (default: 1)
        "-r {input.refpath}"          # -r: Reference sequence (required for GC-depth and mismatches-per-cycle calculation).
        "{input.bam} "                  # mark-dup bam input
        "1> {output.stats}"             # stats output
        "2> {log}"                      # Log redirection

###############################################################################
rule validate_sam:
    # Aim: Basic check for bam file validity, as interpreted by the Broad Institute
    # Use: picard.jar ValidateSamFile \
    #      -I=input.bam \
    #      MODE=SUMMARY
    message:
        "SamTools indexing marked as duplicate BAM file {wildcards.sample} sample ({wildcards.aligner})"
    conda:
        PICARD
    input:
        markdup = "results/02_Mapping/{sample}_{aligner}_mark-dup.bam",
        refpath = "resources/genomes/GCA_018104305.1_AalbF3_genomic.fasta",
    output:
        check = "results/05_Validation/validatesamfile/{sample}_{aligner}_check_mark_dup_bam.txt"
    log:
       "results/11_Reports/validatesamfiles/{sample}_{aligner}_mark-dup.log"
    shell:
        "picard ValidateSamFile " # Validates a SAM/BAM/CRAM file.<p>This tool reports on the validity of a SAM/BAM/CRAM file relative to the SAM format specification
        "-I {input.markdup} "       # Input SAM/BAM/CRAM required
        "-R {input.refpath} "
        "-MODE SUMMARY "            # This mode outputs a summary table listing the numbers of all 'errors' and 'warnings'
        "--QUIET true "
        "&> {log}"

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
        markdup = "results/02_Mapping/{sample}_{aligner}_mark-dup.bam"
    output:
        index = "results/02_Mapping/{sample}_{aligner}_mark-dup.bai"
    log:
        "results/11_Reports/samtools/{sample}_{aligner}_mark-dup-index.log"
    shell:
        "samtools index "     # Samtools index, tools for alignments in the SAM format with command to index alignment
        "-@ {resources.cpus} " # --threads: Number of additional threads to use (default: 1)
        "-b "                  # -b: Generate BAI-format index for BAM files (default)
        "{input.markdup} "     # Markdup bam input
        "{output.index} "      # Markdup index output
        "&> {log}"             # Log redirection

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
        markdup = "results/02_Mapping/{sample}_{aligner}_mark-dup.bam"
    log:
        "results/11_Reports/samtools/{sample}_{aligner}_mark-dup.log"
    shell:
        "samtools markdup "          # Samtools markdup, tools for alignments in the SAM format with command mark duplicates
        "--threads {resources.cpus} " # -@: Number of additional threads to use (default: 1)
        "-r "                         # -r: Remove duplicate reads
        "-s "                         # -s: Report stats
        "--output-fmt BAM "           # -O: Specify output format: SAM, BAM, CRAM (here, BAM format)
        "{input.calmd} "             # Sorted bam input
        "{output.markdup} "           # Markdup bam output
        "&> {log}"                    # Log redirection

###############################################################################
rule samtools_calmd:
    # Aim: Generate the MD tag
    # Use: samtools calmd -@ [THREADS] [REF.FA] -b input.bam > output.bam 
    message:
        "SamTools calmd {wildcards.sample} sample reads ({wildcards.aligner})"
    conda:
        SAMTOOLS
    resources:
       cpus = CPUS,
    params:
        bwapath = BWAPATH,
        reference = REFERENCE
    input:
        sorted = "results/02_Mapping/{sample}_{aligner}_sorted.bam"
    output:
        calmd = "results/02_Mapping/{sample}_{aligner}_sorted_MD.bam"
    log:
        "results/11_Reports/samtools/{sample}_{aligner}_sorted.log"
    threads:
        CPUS
    shell:
        "samtools calmd "                                               # Samtools calmd, tools to generate the MD tag for SNP/indel calling w/o lookin at the reference
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
    # Aim: filling in mate coordinates
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
        mapped = temp("results/02_Mapping/{sample}_bwa-mapped.sam")
    log:
        "results/11_Reports/bwa/{sample}.log"
    shell:
        "bwa mem "                              # BWA-MEM algorithm, performs local alignment.
        "-M "                                   # Mark shorter split hits as secondary (for Picard compatibility). 
        "-T 0 "                                 # Don’t output alignment with score lower than INT. This option only affects output.
        "-R readgroup_string "                  # Complete read group header line. ’\t’ can be used in STR and will be converted to a TAB in the output SAM. The read group ID will be attached to every read in the output. An example is ’@RG\tID:foo\tSM:bar’
        "-t {resources.cpus} "                  # -t: Number of threads (default: 12)
        "-v 1 "                                 # -v: Verbosity level: 1=error, 2=warning, 3=message, 4+=debugging
        "{params.extra} "                       # 
        "{params.bwapath}{params.reference} "   # Reference index filename prefix
        "{input.fwdreads} "                     # Forward input reads
        "{input.revreads} "                     # Reverse input reads
        "1> {output.mapped} "                   # SAM output
        "2> {log}"                              # Log redirection

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
       "sickle "                # Sickle, a windowed adaptive trimming tool for FASTQ files using quality
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
        #report("results/00_Quality/multiqc/multiqc.html", caption="../report/multiqc.rst", category="Quality control"),
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
