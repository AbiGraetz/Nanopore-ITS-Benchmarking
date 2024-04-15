# Quality control for nanopore ITS sequences

# Requirements
- Nextflow (v23.10.0 or later) 
  - [Install directly](https://www.nextflow.io/docs/latest/install.html)
  - or via conda:
    - `conda create -n nextflow-23 bioconda::nextflow=23.10.1`
- A container environment: either Docker or Singularity.
  - [Install Docker](https://docs.docker.com/engine/install/)

# Simple quickstart

1. Ensure you have nextflow running correctly:
    ```shell
    nextflow -version
    ```
2. Run the pipeline
    ```shell
    nextflow run main.nf --input 'data/*.fq.gz' --outdir 'output/'
    ```
    The `--input` parameter indicates which gzipped fastq files to analyze. Wilcards (`'*'`) allows for multiple files to be selected.
    The `--outdir` parameter indicates where the output files should be placed
3. 