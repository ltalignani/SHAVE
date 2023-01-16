###############################################################################
rule bedtools_genome_coverage:
    # Aim: computing genome coverage sequencing
    # Use: bedtools genomecov [OPTIONS] -ibam [MARKDUP.bam] 1> [BEDGRAPH.bed]
    message:
        "BedTools computing genome coverage against reference genome sequence"
    conda:
        BEDTOOLS
    input:
        markdup = "results/02_Mapping/{sample}_{aligner}_sorted-mark-dup-fx.bam",
        index = "results/02_Mapping/{sample}_{aligner}_sorted-mark-dup-fx.bai"
    output:
        genomecov = "results/03_Coverage/{sample}_{aligner}_sorted-mark-dup-genome-cov.bed"
    log:
        "results/11_Reports/bedtools/{sample}_{aligner}_sorted-mark-dup-genome-cov.log"
    shell:
        """
        bedtools genomecov -bga -ibam {input.markdup} 1> {output.genomecov} 2> {log}
        """
