###############################################################################
rule tabix:
    input:
        "results/05_Variants/merged.vcf.gz",
    output:
        "results/05_Variants/merged.vcf.gz.tbi",
    log:
        "results/11_Reports/tabix/merged_tabix.log",
    params:
        # pass arguments to tabix (e.g. index a vcf)
        "-p vcf",
    wrapper:
        "v1.21.2/bio/tabix/index"

###############################################################################
rule bgzip:
    input:
        "results/05_Variants/merged.vcf",
    output:
        "results/05_Variants/merged.vcf.gz",
    params:
        extra="", # optional
    threads: CPUS
    log:
        "results/11_Reports/bgzip/merged_hardfiltered.vcf.gz.log",
    wrapper:
        "v1.21.2/bio/bgzip"

###############################################################################
rule bcftools_merge:
    # Aim: Merge multiple VCF/BCF files from non-overlapping sample sets to create one multi-sample file. 
    #      For example, when merging file A.vcf.gz containing samples S1, S2 and S3 and file B.vcf.gz containing samples S3 and S4, 
    #      the output file will contain five samples named S1, S2, S3, 2:S3 and S4.
    message:
        "bcftools merging all VCFs"
    input:
        calls=expand("results/05_Variants/{sample}_{aligner}.vcf", sample=SAMPLE, aligner=ALIGNER),
    output:
        temp("results/05_Variants/merged.vcf"),
    log:
        "results/11_Reports/bcftools/mergevcfs.log",
    params:
        uncompressed_bcf=False,
        extra="",  # optional parameters for bcftools concat (except -o)
    threads: CPUS
    resources:
        mem_mb=16000,
    wrapper:
        "v1.21.2/bio/bcftools/merge"###############################################################################
