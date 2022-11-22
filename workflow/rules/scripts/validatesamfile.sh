#!/usr/bin/env bash

echo "Picard ValidateSamFile on sample ${snakemake_wildcards[sample]} 2> "${snakemake_log[0]}"

picard ValidateSamFile -I {input.sorted} -R {input.refpath} -O {output.check} -MODE VERBOSE