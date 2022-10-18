Variants where called based on the `GATK best practices workflow`_:
Reads were mapped onto {{ snakemake.config["consensus"]["species"] }} build {{ snakemake.config["consensus"]["reference"] }}  with `BWA mem`_, and both optical and PCR duplicates were removed with Samtools_, 
followed by indel qualities estimation with LoFreq_, and BAM indel indexing with Samtools_. 
The GATK_ HaplotypeCaller was used to call variants per-sample, including summarized evidence for non-variant sites (GVCF_ approach).
Then, genotyping was done in a joint way over GVCF_ files of all samples, using GATK_ GenotypeGVCFs.
{% if snakemake.config["filtering"]["vqsr"] %}
Genotyped variants were filtered with the GATK_ VariantRecalibrator approach.
{% else %}
Genotyped variants were filtered using hard thresholds.
For SNVs, the criterion ``{{ snakemake.config["filtering"]["hard"]["snvs"] }}`` was used, for Indels the criterion ``{{ snakemake.config["filtering"]["hard"]["indels"] }}`` was used.
{% endif %}
Finally, SnpEff_ was used to predict and report variant effects.
In addition, quality control was performed with FastQC_, Samtools_, and aggregated into an interactive report via MultiQC_.

.. _GATK: https://software.broadinstitute.org/gatk/
.. _BWA mem: http://bio-bwa.sourceforge.net/
.. _GATK best practices workflow: https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145
.. _LoFreq: http://csb5.github.io/lofreq/
.. _GVCF: https://gatkforums.broadinstitute.org/gatk/discussion/4017/what-is-a-gvcf-and-how-is-it-different-from-a-regular-vcf
.. _SnpEff: http://snpeff.sourceforge.net
.. _MultiQC: http://multiqc.info/
.. _Samtools: http://samtools.sourceforge.net/
.. _FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
