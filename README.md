# manuscript

## Table of Contents

- [Installation](#installation)
- [Pipeline](#pipeline)
- [Usage](#usage)
- [Configuration](#configuration)
- [Examples](#examples)
- [License](#license)

## Installation

### To install, clone this repository:
```bash
git clone https://github.com/alankang25/manuscript.git
cd manuscript
```

### This pipeline uses a conda environment to manage dependencies. To install with Conda:
```bash
# 1. Create environment (first time only)
conda env create -f environment.yml

# 2. Activate new environment
conda activate manuscript-pipeline
```

## Pipeline
### Below is an outline of the pipeline which takes in an ENCODE formatted TSV file, ATAC-seq BigWig and peak files pre- and post- treatment to list the importance of various chromatin features in explaining accessibility changes across treatment groups. This output is organized in a .csv file named output.csv. In this example, data from BRG1/BRM Associated Factors Complex (BAF) inhibited GM12878 cells will be used.

![Pipeline diagram showing data flow and model steps](docs/pipeline_diagram.svg)

## Usage
### FRiP and Peak Count Filtering

```bash
# For this demo, a tsv with the following filters from ENCODE was downloaded:
# biosample=GM12878, Output Type=IDR thresholded peaks, File format=bed, Assay Title=TF ChIP-seq, Status=Released

    python FRiP_filter.py -i "INPUT_ENCODE_TSV" -t "CHIP-SEQ_TYPE"
```

### To reset the ENCODE BED download directory
```bash
    python Reset_bed.py
```

### Random Forest and Linear Regression Feature Analysis
```bash
    python Feature_analysis.py -c "PATH_TO_CONTROL_FILES" -t "PATH_TO_TREATMENT_FILES"
```

## Configuration
### If needed, there are some configurations available for the pipeline.
```bash
# While filtering for the highest FRiP value in a target, you can also set a minimum peak number cutoff. For example: 
    python FRiP_filter.py -i "INPUT_ENCODE_TSV" -t "CHIP-SEQ_TYPE" -m 1000 
#will remove all BED files with less than 1000 peaks total. This way, you can filter out files that have low signal.


#
