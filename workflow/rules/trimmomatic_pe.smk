###############################################################################
rule trimmomatic_pe:
    input:
        r1="resources/reads/{sample}_R1.fastq.gz",
        r2="resources/reads/{sample}_R2.fastq.gz"
    output:
        r1="results/01_Trimming/trimmomatic/{sample}_trimmomatic-trimmed_R1.fastq.gz",
        r2="results/01_Trimming/trimmomatic/{sample}_trimmomatic-trimmed_R2.fastq.gz",
        # reads where trimming entirely removed the mate
        r1_unpaired="results/01_Trimming/trimmomatic/{sample}.1.unpaired.fastq.gz",
        r2_unpaired="results/01_Trimming/trimmomatic/{sample}.2.unpaired.fastq.gz"
    log:
        "results/11_Reports/trimmomatic/{sample}.log"
    params:
        # list of trimmers (see manual)
        trimmer=["TRAILING:3"],
        # optional parameters
        extra="",
        compression_level="-9"
    threads:
        CPUS
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=16000
    wrapper:
        "v1.21.2/bio/trimmomatic/pe"


rule trimmomatic:
    input:
        reads = get_fastq,
        adapters = config['trimmomatic']["adapters"]
    output:
        forward_reads   = WORKING_DIR + "trimmed/{samples}_forward.fastq.gz",
        reverse_reads   = WORKING_DIR + "trimmed/{samples}_reverse.fastq.gz",
        forwardUnpaired = temp(WORKING_DIR + "trimmed/{samples}_forward_unpaired.fastq.gz"),
        reverseUnpaired = temp(WORKING_DIR + "trimmed/{samples}_reverse_unpaired.fastq.gz")
    log:
        RESULT_DIR + "logs/trimmomatic/{samples}.log"
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
    threads: 10
    conda:
        "../envs/trimmomatic.yaml"
    message:
        "Trimming reads for {wildcards.samples}"
    shell:
        """
        trimmomatic PE \
        -threads 10 \
        -phred33 \
         {input.reads} \
         {output.forward_reads} {output.forwardUnpaired} {output.reverse_reads} {output.reverseUnpaired} \
         ILLUMINACLIP:{input.adapters}:{params.seedMisMatches}:{params.palindromeClipTreshold}:{params.simpleClipThreshhold} \
         LEADING:3 \
         TRAILING:3 \
         SLIDINGWINDOW:4:15 \
         MINLEN:40 &>{log}
        """