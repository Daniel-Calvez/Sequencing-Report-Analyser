#! /bin/bash

Help()
{
   # Display Help
   echo "Add description of the script functions here."
   echo "Syntax: pipeline_report_script.sh [-r REFERENCE|-s SOURCE|-o OUTPUT|-c CPUNUMBER|-h]"
   echo "options:"
   echo "-r    The fasta reference file"
   echo "-s    The source file"
   echo "-o    The output folder"
   echo "-c    The cpu number to use (1 by default)"
   echo "-m    If you want to fass a folder containing multiple fastq files"
   echo "-h    The helper"
}
multipleFiles=0
### Docker and conda need to be installed ###
while getopts "hr:s:o:c:m" option; do
   case $option in
        h) # display Help
            Help
            exit;;
        r) # store reference
            reference=$OPTARG;;
        s) # store source
            source=$OPTARG;;
        o) # store output directory path
            Output_Dir=$OPTARG;;
        c) # store number of CPU used
            CPU_Number=$OPTARG;;
        m)
            multipleFiles=1;;
        \?)
            echo "Error: Invalid option"
            exit;;
   esac
done
 
if [ -z "$reference" ];
then echo "Error: reference undefined"
exit
fi

if [ -z "$source" ];
then echo "Error: source undefined"
exit
fi

if [ -z "$Output_Dir" ];
then echo "Error: output directory undefined"
exit
fi

if [ -z "$CPU_Number" ];
then echo "Error: CPU number undefined"
exit
fi

mkdir -p "$Output_Dir"
# Get the folder from reference
refFolder=$(dirname "$reference")

refFile=$(basename "$reference")

tempPath="$PWD"
cd "$Output_Dir" || return
Output_Dir="$PWD"
cd "$tempPath" || return



tempPath="$PWD"
cd "$refFolder" || return
refFolder="$PWD"
cd "$tempPath" || return

# Base case: if there is only one patient
if [ $multipleFiles == 0 ];
then 
    mkdir -p "$Output_Dir/python_parse_fastq_results"
    mkdir -p "$Output_Dir/results/deepvariant_results/intermediate_results_dir"
    mkdir -p "$Output_Dir/results/Bam_Sam"
    mkdir -p "$Output_Dir/results/Bam_Sam"
    mkdir -p "$Output_Dir/python_parse_sam_results"
    mkdir -p "$Output_Dir/results/Consensus_Sequence/"
    touch "$Output_Dir/results/noCoverageZones.bed"
    touch "$Output_Dir/python_vcf_result.vcf"

    Result_Dir="${Output_Dir}/results"

# unzip fastq
gzip -k -d "$source"
# script generation sam 
minimap2 -a "${reference}" "${source}" > "${Result_Dir}/Bam_Sam/alignment.sam"
# script generation bam
samtools view -b "${Result_Dir}/Bam_Sam/alignment.sam" -o "${Result_Dir}/Bam_Sam/bamAlignment.bam"
# script sort bam
samtools sort "${Result_Dir}/Bam_Sam/bamAlignment.bam" -o "${Result_Dir}/Bam_Sam/bamAlignment_sorted.bam"
# script generation bai
samtools index "${Result_Dir}/Bam_Sam/bamAlignment_sorted.bam" -o "${Result_Dir}/Bam_Sam/bamAlignment_sorted.bai"
# script generation fai
samtools faidx "$reference"

# script docker: 
source=$2 docker pull google/deepvariant


## script docker deepvariant --------------------------------------------------  
sudo docker run \
 -v "${refFolder}":"/input" \
 -v "${Output_Dir}":"/output" \
 google/deepvariant /opt/deepvariant/bin/run_deepvariant \
 --model_type ONT_R104 \
 --ref=/input/"$refFile" \
 --reads=/output/results/Bam_Sam/bamAlignment_sorted.bam \
 --output_vcf=/output/results/deepvariant_results/results.output.vcf.gz \
 --output_gvcf=/output/results/deepvariant_results/results.output.g.vcf.gz \
 --intermediate_results_dir=/output/results/deepvariant_results/intermediate_results_dir
## fin du script -------------------------------------------

gzip -k -d "$Output_Dir/results/deepvariant_results/results.output.vcf.gz"

# script remove lock
sudo chmod -R 777 "${Output_Dir}/results/deepvariant_results"

# script for bed file generation and getting zones with 0 coverage
bedtools genomecov -bga -ibam "$Output_Dir/results/Bam_Sam/bamAlignment_sorted.bam" | grep -w 0$ > "$Output_Dir/results/noCoverageZones.bed"

#Script for regular bed file generation
bedtools genomecov -bga -ibam "$Output_Dir/results/Bam_Sam/bamAlignment_sorted.bam" > "$Output_Dir/results/Coverage.bed"

# script consensus (replaces no coverage zones by N)
bcftools consensus -m "${Output_Dir}/results/noCoverageZones.bed" -f "$reference" "${Output_Dir}/results/deepvariant_results/results.output.vcf.gz" > "${Output_Dir}/results/Consensus_Sequence/consensus.fa"

# Script analyse python sam
python "$PWD/python/parse_sam.py" -i "$Output_Dir/results/Bam_Sam/alignment.sam" -o "$Output_Dir/python_parse_sam_results" -c "$CPU_Number"

# Script analyse python vcf
python "$PWD/python/filter_vcf.py" -r -b "$Output_Dir/results/Bam_Sam/bamAlignment_sorted.bam" -v "$Output_Dir/results/deepvariant_results/results.output.vcf" -f "$reference" -o "$Output_Dir/python_vcf_result.vcf" -c "$CPU_Number"

# Script analyse python fastq
python "$PWD/python/parse_fastq.py" -i "$source" -o "$Output_Dir/python_parse_fastq_results" -c "$CPU_Number"

# Generate report
Rscript "${PWD}/R/Render_report.R" -s "${Output_Dir}" -n Rapport -p "${PWD}/R/Rapport_sequencage.Rmd" -w "${PWD}/R/Rapport_sequencage_html.Rmd" -r "$reference" 

#If there are multiple patients
else 
nbFiles=1
for entry in "$source"/*
do
filename=$(basename -- "$entry")
extension="${filename##*.}"
filename="${filename%.*}"
echo "$extension"
if [ "$extension" != "fastq" ];
then continue
fi
    tempOutputDir="$Output_Dir/$filename"
    mkdir -p "$tempOutputDir"
    mkdir -p "$tempOutputDir/python_parse_fastq_results"
    mkdir -p "$tempOutputDir/results/deepvariant_results/intermediate_results_dir"
    mkdir -p "$tempOutputDir/results/Bam_Sam"
    mkdir -p "$tempOutputDir/results/Bam_Sam"
    mkdir -p "$tempOutputDir/python_parse_sam_results"
    mkdir -p "$tempOutputDir/results/Consensus_Sequence"
    touch "$tempOutputDir/results/noCoverageZones.bed"
    touch "$tempOutputDir/python_vcf_result.vcf"
    
    tempSource=$entry
    
    Result_Dir="${tempOutputDir}/results"

# unzip fastq
gzip -k -d "$tempSource"
# script generation sam 
minimap2 -a "${reference}" "${tempSource}" > "${Result_Dir}/Bam_Sam/alignment.sam"
# script generation bam
samtools view -b "${Result_Dir}/Bam_Sam/alignment.sam" -o "${Result_Dir}/Bam_Sam/bamAlignment.bam"
# script sort bam
samtools sort "${Result_Dir}/Bam_Sam/bamAlignment.bam" -o "${Result_Dir}/Bam_Sam/bamAlignment_sorted.bam"
# script generation bai
samtools index "${Result_Dir}/Bam_Sam/bamAlignment_sorted.bam" -o "${Result_Dir}/Bam_Sam/bamAlignment_sorted.bai"
# script generation fai
samtools faidx "$reference"

# script docker: 
docker pull google/deepvariant

## script docker deepvariant --------------------------------------------------  
sudo docker run \
 -v "${refFolder}":"/input" \
 -v "${tempOutputDir}":"/output" \
 google/deepvariant /opt/deepvariant/bin/run_deepvariant \
 --model_type ONT_R104 \
 --ref=/input/"$refFile" \
 --reads=/output/results/Bam_Sam/bamAlignment_sorted.bam \
 --output_vcf=/output/results/deepvariant_results/results.output.vcf.gz \
 --output_gvcf=/output/results/deepvariant_results/results.output.g.vcf.gz \
 --intermediate_results_dir=/output/results/deepvariant_results/intermediate_results_dir
## fin du script -------------------------------------------

gzip -k -d "$tempOutputDir/results/deepvariant_results/results.output.vcf.gz"

# script remove lock
sudo chmod -R 777 "${tempOutputDir}/results/deepvariant_results"

# script for bed file generation and getting zones with 0 coverage
bedtools genomecov -bga -ibam "$tempOutputDir/results/Bam_Sam/bamAlignment_sorted.bam" | grep -w 0$ > "$tempOutputDir/results/noCoverageZones.bed"

#Script for regular bed file generation
bedtools genomecov -bga -ibam "$tempOutputDir/results/Bam_Sam/bamAlignment_sorted.bam" > "$tempOutputDir/results/Coverage.bed"

# script consensus (replaces no coverage zones by N)
bcftools consensus -m "${tempOutputDir}/results/noCoverageZones.bed" -f "$reference" "${tempOutputDir}/results/deepvariant_results/results.output.vcf.gz" > "${tempOutputDir}/results/Consensus_Sequence/consensus.fa"

# Script analyse python sam
python "$PWD/python/parse_sam.py" -i "$tempOutputDir/results/Bam_Sam/alignment.sam" -o "$tempOutputDir/python_parse_sam_results" -c "$CPU_Number"

# Script analyse python vcf
python "$PWD/python/filter_vcf.py" -r -b "$tempOutputDir/results/Bam_Sam/bamAlignment_sorted.bam" -v "$tempOutputDir/results/deepvariant_results/results.output.vcf" -f "$reference" -o "$tempOutputDir/python_vcf_result.vcf" -c "$CPU_Number"

# Script analyse python fastq
python "$PWD/python/parse_fastq.py" -i "$tempSource" -o "$tempOutputDir/python_parse_fastq_results" -c "$CPU_Number"

# Generate report for each single patient
Rscript "${PWD}/R/Render_report.R" -s "${tempOutputDir}" -n Rapport -p "${PWD}/R/Rapport_sequencage.Rmd" -w "${PWD}/R/Rapport_sequencage_html.Rmd" -r "$reference" 

nbFiles=$((nbFiles+1))
done

# Generate a report for all patients
Rscript "${PWD}/R/Render_report.R" -s "${Output_Dir}" -n Rapport_Global -p "${PWD}/R/Rapport_sequencage_multi.Rmd" -w "${PWD}/R/Rapport_sequencage_html_multi.Rmd" -r "$reference"

fi
