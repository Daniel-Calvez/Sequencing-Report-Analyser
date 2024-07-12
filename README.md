# Sequencing-Report-Analyser
A tool to generate reports from a fastq file and a reference
---

# Setup

## Prerequisite

To run the script you will need conda and docker

<a href="https://conda.io/projects/conda/en/latest/user-guide/install/index.html">Link to install Conda</a>

<a href="https://docs.docker.com/engine/install/">Link to install Docker</a>

## Conda setup

To use this tool you can manually install python, R and all the necessary packages but it's easier to use Conda

Move into the folder with the project in it: ```cd Sequencing-Report-Analyser```

Create the conda environment with: ```conda env create --file env_pipeline_report.yml```

Activate the conda environment with: ```conda activate env_pipeline_report```

---

# Using the tool

First, you need to move to the folder containing the project: ```cd Sequencing-Report-Analyser```

Then you need to launch the script: <br>
```bash pipeline_report_script.sh -r [referenceFilePath] -s [fastqFile] -o [OutputFolder] -c [CpuNumber] [-m]```

[referenceFilePath]: The path of the fasta reference file<br>
[fastqFile]: The path of the fastqFile from the sequencer<br>
[OutputFolder]: The output folder that will contain the created files along with the<br>
[CpuNumber]: The number of cpu used for the scripts (optional, default is 1)<br>
[-m]: Option to enable multi-patient analysis, if you use it, pass a folder containing all fastq files to analyse instead of a single file <br>
The pdf and html reports are in the output folder
