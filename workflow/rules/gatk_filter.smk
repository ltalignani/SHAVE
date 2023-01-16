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
        ref="resources/genomes/GCA_018104305.1_AalbF3_genomic.fasta",
        vcf="results/05_Variants/merged.vcf.gz",
    output:
        vcf="results/05_Variants/merged_hardfiltered.vcf.gz",
    params:
        filters={"myfilter": "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"},
        extra="",
        java_opts="",
    resources:
        mem_mb=16000
    log:
        "results/11_Reports/variantfiltration/merged_hardfiltered.log",
    wrapper:
        "v1.21.2/bio/gatk/variantfiltration"