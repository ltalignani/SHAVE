#!/bin/bash

###### Required for local computing ######
# For Mac with Apple Silicon processors, create a new empty osx-64 specific environment :
# conda deactivate base
# CONDA_SUBDIR=osx-64 conda create -n shave
# conda activate shave
# conda config --env --set subdir osx-64

##### Colors ######
red="\033[1;31m"   # red
green="\033[1;32m" # green
ylo="\033[1;33m"   # yellow
blue="\033[1;34m"  # blue
nc="\033[0m"       # no color

###### About ######
echo ""
echo -e "${green}------------------------------------------------------------------------${nc}"
echo -e "${green}#####${nc} ${red}ABOUT${nc} ${green}#####${nc}"
echo -e "${green}-----------------${nc}"
echo ""
echo -e "${blue}Name${nc} __________________ Start_shave.sh"
echo -e "${blue}Author${nc} ________________ Loïc Talignani"
echo -e "${blue}Affiliation${nc} ___________ UMR_MIVEGEC"
echo -e "${blue}Aim${nc} ___________________ Bash script for ${red}SH${nc}ort-read ${red}A${nc}lignment pipeline for ${red}VE${nc}ctor v.1"
echo -e "${blue}Date${nc} __________________ 2022.10.05"
echo -e "${blue}Run${nc} ___________________ bash Start_shave1.sh"
echo -e "${blue}Latest Modification${nc} ___ "


###### Hardware ######
echo ""
echo -e "${green}------------------------------------------------------------------------${nc}"
echo -e "${green}#####${nc} ${red}OPERATING SYSTEM${nc} ${green}#####${nc}"
echo -e "${green}--------------------${nc}"
echo ""

# Operating System
case "$OSTYPE" in
  linux*)   os="Linux" ;;
  bsd*)     os="BSD" ;;
  darwin*)  os="OSX" ;;
  solaris*) os="Solaris" ;;
  msys*)    os="Windows" ;;
  cygwin*)  os="Windows" ;;
  *)        os="Unknown (${OSTYPE})" ;;
esac

echo -e "${blue}Operating system${nc} _______ ${red}${os}${nc}" # Print operating system


###### Hardware ######
echo ""
echo -e "${green}------------------------------------------------------------------------${nc}"
echo -e "${green}###########${nc} ${red}HARDWARE CHARACTERISTICS${nc} ${green}###########${nc}"
echo -e "${green}------------------------------------------------------------------------${nc}"
echo ""

if [[ ${os} == "OSX" ]]
then
    model_name=$(sysctl -n machdep.cpu.brand_string) # Get chip model name
    physical_cpu=$(sysctl -n hw.physicalcpu)         # Get physical cpu
    logical_cpu=$(sysctl -n hw.logicalcpu)           # Get logical cpu
    mem_size=$(sysctl -n hw.memsize)                 # Get memory size (bit)
    ram_size=$(expr ${mem_size} \/ $((1024**3)) )    # / 1024**3 = Gb
elif [[ ${os} == "Linux" ]]
then
    model_name=$(lscpu | grep -o -E "Model name: +.+" | sed -r "s/Model name: +//")                           # Get chip model name
    physical_cpu=$(lscpu | grep -o -E "^CPU\(s\): +[0-9]+" | sed -r "s/CPU\(s\): +//")                        # Get physical cpu
    threads_cpu=$(lscpu | grep -o -E "^Thread\(s\) per core: +[0-9]+" | sed -r "s/Thread\(s\) per core: +//") # Get thread(s) per core
    logical_cpu=$(expr ${physical_cpu} \* ${threads_cpu})                                                     # Calcul logical cpu
    mem_size=$(grep -o -E "MemTotal: +[0-9]+" /proc/meminfo | sed -r "s/MemTotal: +//")                       # Get memory size (Kb)
    ram_size=$(expr ${mem_size} \/ $((1024**2)) )                                                             # / 1024**2 = Gb
else
    "Please, use 'OSX' or 'Linux' operating system"
    exit 1
fi

echo -e "${ylo}Brand(R)${nc} | ${ylo}Type(R)${nc} | ${ylo}Model${nc} | ${ylo}@ Speed GHz${nc}" # Print header chip model name
echo -e "${blue}Chip Model Name${nc} _______ ${model_name}"                     # Print chip model name
echo -e "${blue}Physical CPUs${nc} _________ ${red}${physical_cpu}${nc} cores" # Print physical cpu
echo -e "${blue}Logical CPUs${nc} __________ ${red}${logical_cpu}${nc} threads" # Print logical cpu
echo -e "${blue}System Memory${nc} _________ ${red}${ram_size}${nc} Gb of RAM"  # Print RAM size

###### Settings ######
echo ""
echo -e "${green}------------------------------------------------------------------------${nc}"
echo -e "${green}###################${nc} ${red}SETTINGS${nc} ${green}###################${nc}"
echo -e "${green}------------------------------------------------------------------------${nc}"
echo ""

workdir=$(cd "$(dirname "${BASH_SOURCE[0]}" )" && pwd)                                          # Get working directory
fastq=$(expr $(ls -l ${workdir}/resources/reads/*.fastq.gz | wc -l))                            # Count fastq.gz files
samples=$(expr ${fastq} \/ 2)                                                                   # / 2 = samples (paired-end)
max_threads=$(grep -o -E "cpus: [0-9]+" ${workdir}/config/config.yaml | sed "s/cpus: //")       # Get user config for max threads
max_memory=$(grep -o -E "mem_gb: [0-9]+" ${workdir}/config/config.yaml | sed "s/mem_gb: //")    # Get user config for max memory
reference=$(grep -o -E "reference: '.+'" ${workdir}/config/config.yaml | sed "s/reference: //") # Get user config genome reference
aligner=$(grep -o -E "^aligner: '[a-z]+'" ${workdir}/config/config.yaml | sed "s/aligner: //")  # Get user config aligner
min_cov=$(grep -o -E "mincov: [0-9]+" ${workdir}/config/config.yaml | sed "s/mincov: //")       # Get user config minimum coverage
min_af=$(grep -o -E "minaf: [0-1]\.[0-9]+" ${workdir}/config/config.yaml | sed "s/minaf: //")   # Get user config minimum allele frequency
time_stamp_start=$(date +"%Y-%m-%d %H:%M")                                                      # Get analyzes starting time
SECONDS=0                                                                                       # Initialize SECONDS counter

echo -e "${blue}Working Directory${nc} _____ ${workdir}/"                                                     # Print working directory
echo -e "${blue}Samples Processed${nc} _____ ${red}${samples}${nc} samples (${ylo}${fastq}${nc} fastq files)" # Print samples number
echo -e "${blue}Maximum Threads${nc} _______ ${red}${max_threads}${nc} of ${ylo}${logical_cpu}${nc} threads available" # Print max threads
echo -e "${blue}Maximum Memory${nc} ________ ${red}${max_memory}${nc} of ${ylo}${ram_size}${nc} Gb available" # Print max memory
echo -e "${blue}Genome Reference${nc} ______ ${red}${reference}${nc}"                                         # Print user config genome reference
echo -e "${blue}Aligner${nc} _______________ ${ylo}${aligner}${nc}"                                           # Print user config aligner
echo -e "${blue}Min. Coverage${nc} _________ ${red}${min_cov}${nc}x"                                          # Print user config minimum coverage
echo -e "${blue}Min. Allele Frequency${nc} _ ${red}${min_af}${nc}"                                            # Print user config minimum Al.Freq.
echo -e "${blue}Start Time${nc} ____________ ${time_stamp_start}"                                             # Print analyzes starting time


###### Installations ######
echo ""
echo -e "${green}------------------------------------------------------------------------${nc}"
echo -e "${green}#########${nc} ${red}APPS INSTALLATIONS WITH CONDA${nc} ${green}########${nc}"
echo -e "${green}------------------------------------------------------------------------${nc}"
echo ""

# Snakemake
#snake_ver="7.8.2"
snake_ver="6.12.3"
if ls ~/miniconda3/bin/snakemake 2> /dev/null
then
    echo ""
#   source ${HOME}/miniconda3/etc/profile.d/conda.sh
else

    conda install -c conda-forge -c bioconda snakemake==${snake_ver} --yes
fi

# Rename
if ls ~/miniconda3/bin/rename 2> /dev/null
then
    echo ""
else
    conda install -c bioconda rename --yes
fi

# gatk3
gatk_ver="3.8"
if ls ~/miniconda3/bin/gatk3 2> /dev/null
then
    echo ""
#   source ${HOME}/miniconda3/etc/profile.d/conda.sh
else

    conda install -c conda-forge -c bioconda gatk==${gatk_ver} --yes
fi

# Picard tools
picard_ver="2.27.4"
if ls ~/miniconda3/bin/picard 2> /dev/null
then
    echo ""
else
    conda install -c bioconda picard==${picard_ver} --yes
fi

###### Rename samples ######
echo ""
echo -e "${green}------------------------------------------------------------------------${nc}"
echo -e "${green}##############${nc} ${red}RENAME FASTQ FILES${nc} ${green}##############${nc}"
echo -e "${green}------------------------------------------------------------------------${nc}"
echo ""

# Rename fastq files to remove "_001" Illumina pattern.
## De/comment (#) if you want keep Illumina barcode-ID and/or Illumina line-ID
#rename "s/_S\d+_/_/" ${workdir}/resources/reads/*.fastq.gz                 # Remove barcode-ID like {_S001_}
rename "s/_L\d+_/_/" ${workdir}/resources/reads/*.fastq.gz                  # Remove line-ID ID like {_L001_}
rename "s/_001.fastq.gz/.fastq.gz/" ${workdir}/resources/reads/*.fastq.gz   # Remove end-name ID like {_001}.fastq.gz
#rename 's/(\w_).*_(R[1-2]).*(.fastq.gz)/$1$2$3/' *.fastq.gz                 # Keep only expr. in ( )

###### Call snakemake pipeline ######
echo ""
echo -e "${green}------------------------------------------------------------------------${nc}"
echo -e "${green}###########${nc} ${red}SNAKEMAKE PIPELINE START${nc} ${green}###########${nc}"
echo -e "${green}------------------------------------------------------------------------${nc}"
echo ""

echo -e "${blue}Unlocking working directory:${nc}"
echo ""
# Specify working directory (relative paths in the snakefile will use this as their origin).
# The workflow definition in form of a snakefile.
# Set or overwrite values in the workflow config object.
# Re-run all jobs the output of which is recognized as incomplete.
# Remove a lock on the working directory.
snakemake \
    --directory ${workdir}/ \
    --snakefile ${workdir}/workflow/rules/shave.smk \
    --config os=${os} \
    --rerun-incomplete \
    --unlock

echo ""
echo -e "${blue}Conda environments list:${nc}"
echo ""
# Specify working directory (relative paths in the snakefile will use this as their origin).
# The workflow definition in form of a snakefile.
# Use at most N CPU cores/jobs in parallel. If N is omitted or ‘all’, the limit is set to the number of available CPU cores.
# Set or overwrite values in the workflow config object.
# Re-run all jobs the output of which is recognized as incomplete.
# List all conda environments and their location on disk.
snakemake \
    --directory ${workdir}/ \
    --snakefile ${workdir}/workflow/rules/shave.smk \
    --cores ${max_threads} \
    --config os=${os} \
    --rerun-incomplete \
    --list-conda-envs

echo ""
echo -e "${blue}Conda environments update:${nc}"
echo ""
# Specify working directory (relative paths in the snakefile will use this as their origin).
# The workflow definition in form of a snakefile.
# Use at most N CPU cores/jobs in parallel. If N is omitted or ‘all’, the limit is set to the number of available CPU cores.
# Re-run all jobs the output of which is recognized as incomplete.
# Set or overwrite values in the workflow config object.
# Cleanup unused conda environments.
snakemake \
    --directory ${workdir}/ \
    --snakefile ${workdir}/workflow/rules/shave.smk \
    --cores ${max_threads} \
    --config os=${os} \
    --rerun-incomplete \
    --conda-cleanup-envs

echo ""
echo -e "${blue}Conda environments setup:${nc}"
echo ""
# Specify working directory (relative paths in the snakefile will use this as their origin).
# The workflow definition in form of a snakefile.
# Use at most N CPU cores/jobs in parallel. If N is omitted or ‘all’, the limit is set to the number of available CPU cores.
# Set or overwrite values in the workflow config object.
# Re-run all jobs the output of which is recognized as incomplete.
# If defined in the rule, run job in a conda environment.
# If specified, only creates the job-specific conda environments then exits. The –use-conda flag must also be set.
# If mamba package manager is not available, or if you still prefer to use conda, you can enforce that with this setting (default: 'mamba').
snakemake \
    --directory ${workdir}/ \
    --snakefile ${workdir}/workflow/rules/shave.smk \
    --cores ${max_threads} \
    --config os=${os} \
    --rerun-incomplete \
    --use-conda \
    --conda-create-envs-only \
    --conda-frontend conda # Default "mamba", recommended because much faster, but : "Library not loaded: @rpath/libarchive.13.dylib"

echo ""
echo -e "${blue}Dry run:${nc}"
echo ""
# Specify working directory (relative paths in the snakefile will use this as their origin).
# The workflow definition in form of a snakefile.
# Use at most N CPU cores/jobs in parallel. If N is omitted or ‘all’, the limit is set to the number of available CPU cores.
# Set or overwrite values in the workflow config object.
# Re-run all jobs the output of which is recognized as incomplete.
# If defined in the rule, run job in a conda environment.
# Tell the scheduler to assign creation of given targets (and all their dependencies) highest priority.
# Do not execute anything, and display what would be done. If very large workflow, use –dry-run –quiet to just print a summary of the DAG of jobs.
# Do not output any progress or rule information.
snakemake \
    --directory ${workdir}/ \
    --snakefile ${workdir}/workflow/rules/shave.smk \
    --cores ${max_threads}\
    --config os=${os} \
    --rerun-incomplete \
    --use-conda \
    --conda-frontend conda \
    --prioritize multiqc_reports_aggregation \
    --dry-run \
    --quiet

echo ""
echo -e "${blue}Let's run!${nc}"
echo ""
# Specify working directory (relative paths in the snakefile will use this as their origin).
# The workflow definition in form of a snakefile.
# Use at most N CPU cores/jobs in parallel. If N is omitted or ‘all’, the limit is set to the number of available CPU cores.
# Define a global maximum number of threads available to any rule. Rules requesting more threads will have their values reduced to the maximum.
# Set or overwrite values in the workflow config object.
# Re-run all jobs the output of which is recognized as incomplete.
# Go on with independent jobs if a job fails.
# If defined in the rule, run job in a conda environment.
# Tell the scheduler to assign creation of given targets (and all their dependencies) highest priority.
# Print out the shell commands that will be executed.
snakemake \
    --directory ${workdir}/ \
    --snakefile ${workdir}/workflow/rules/shave.smk \
    --cores ${max_threads} \
    --max-threads ${max_threads} \
    --config os=${os} \
    --rerun-incomplete \
    --keep-going \
    --use-conda \
    --conda-frontend conda \
    --prioritize multiqc_reports_aggregation \
    --printshellcmds


###### Create usefull graphs, summary and logs ######
echo ""
echo -e "${green}------------------------------------------------------------------------${nc}"
echo -e "${green}###########${nc} ${red}SNAKEMAKE PIPELINE LOGS${nc} ${green}############${nc}"
echo -e "${green}------------------------------------------------------------------------${nc}"
echo ""

mkdir ${workdir}/results/10_Graphs/ 2> /dev/null

graph_list="dag rulegraph filegraph"
extention_list="pdf png"

for graph in ${graph_list} ; do
    for extention in ${extention_list} ; do
	snakemake \
	    --directory ${workdir}/ \
            --snakefile ${workdir}/workflow/rules/shave.smk \
            --${graph} | \
	    dot -T${extention} > \
		${workdir}/results/10_Graphs/${graph}.${extention} ;
    done ;
done

snakemake \
    --directory ${workdir} \
    --snakefile ${workdir}/workflow/rules/shave.smk \
    --summary > ${workdir}/results/11_Reports/files_summary.txt

cp ${workdir}/config/config.yaml ${workdir}/results/11_Reports/config.yaml

echo "                        Brand(R) | Type(R) | Model | @ Speed GHz" >> ${workdir}/results/11_Reports/settings.log                       # Log header chip model name
echo "Chip Model Name _______ ${model_name}" >> ${workdir}/results/11_Reports/settings.log                                                  # Log chip model name
echo "Physical CPUs _________ ${physical_cpu} cores" >> ${workdir}/results/11_Reports/settings.log                                         # Log physical cpu
echo "Logical CPUs __________ ${logical_cpu} threads" >> ${workdir}/results/11_Reports/settings.log                                         # Log logical cpu
echo "System Memory _________ ${ram_size} Gb of RAM" >> ${workdir}/results/11_Reports/settings.log                                          # Log RAM size

echo "Working Directory _____ ${workdir}/" >> ${workdir}/results/11_Reports/settings.log                                                    # Log working directory
echo "Samples Processed _____ ${samples} samples (${fastq} fastq files)" >> ${workdir}/results/11_Reports/settings.log                      # Log samples number
echo "Maximum Threads _______ ${max_threads} of ${logical_cpu} threads available" >> ${workdir}/results/11_Reports/settings.log             # Log max threads
echo "Maximum Memory ________ ${max_memory} of ${ram_size} Gb available" >> ${workdir}/results/11_Reports/settings.log                      # Log max memory
echo "Genome Reference ______ ${reference}" >> ${workdir}/results/11_Reports/settings.log                                                   # Log user config genome reference
echo "Aligner _______________ ${aligner}" >> ${workdir}/results/11_Reports/settings.log                                                     # Log user config aligner
echo "Min. Coverage _________ ${min_cov}" >> ${workdir}/results/11_Reports/settings.log                                                     # Log user config minimum coverage
echo "Min. Allele Frequency _ ${min_af}" >> ${workdir}/results/11_Reports/settings.log                                                      # Log user config snvs cov min
echo "Start Time ____________ ${time_stamp_start}" >> ${workdir}/results/11_Reports/settings.log                                            # Log analyzes starting time

###### copy multiqc_report.html to results/ dir root ######
echo ""
echo -e "${blue}------------------------------------------------------------------------${nc}"
echo -e "${blue}#############${nc} ${red}COPY QC REPORT TO ROOT${nc} ${blue}############${nc}"
echo -e "${blue}------------------------------------------------------------------------${nc}"
echo ""

# and copy multiqc_report.html to results/ dir root
cp ${workdir}/results/00_Quality_Control/multiqc/multiqc_report.html ${workdir}/results/All_readsQC_reports.html

###### Concatenate all coverage stats ######
echo ""
echo -e "${blue}------------------------------------------------------------------------${nc}"
echo -e "${blue}###########${nc} ${red}CONCATENATE COVERAGE STATS${nc} ${blue}##########${nc}"
echo -e "${blue}------------------------------------------------------------------------${nc}"
echo ""

cat ${workdir}/results/03_Coverage/*coverage-stats.tsv > ${workdir}/results/All_genome_coverages.tsv

awk "NR==1 || NR%2==0" ${workdir}/results/All_genome_coverages.tsv > ${workdir}/results/GENCOV.tmp \
    && mv ${workdir}/results/GENCOV.tmp ${workdir}/results/All_genome_coverages.tsv


###### End managment ######
echo ""
echo -e "${blue}------------------------------------------------------------------------${nc}"
echo -e "${blue}##################${nc} ${red}SCRIPT END${nc} ${blue}###################${nc}"
echo -e "${blue}------------------------------------------------------------------------${nc}"
echo ""

find ${workdir}/results/ -type f -empty -delete                 # Remove empty file (like empty log)
find ${workdir}/results/ -type d -empty -delete                 # Remove empty directory

time_stamp_end=$(date +"%Y-%m-%d %H:%M")                        # Get date / hour ending analyzes
elapsed_time=${SECONDS}                                         # Get SECONDS counter
minutes=$((${elapsed_time}/60))                                 # / 60 = minutes
seconds=$((${elapsed_time}%60))                                    # % 60 = seconds

echo -e "${red}End Time${nc} ______________ ${time_stamp_end}"                                                             # Print analyzes ending time
echo -e "${red}Processing Time${nc} _______ ${ylo}${minutes}${nc} minutes and ${ylo}${seconds}${nc} seconds elapsed"       # Print total time elapsed

echo ""
echo -e "${blue}------------------------------------------------------------------------${nc}"
echo ""
