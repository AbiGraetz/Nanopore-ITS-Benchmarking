# Quality control for nanopore ITS sequences

## Report
[Quality report and read statistics for all samples combined](docs/report.md)

# Running the pipeline
## Requirements
- Nextflow (v23.10.0 or later) 
  - [Install directly](https://www.nextflow.io/docs/latest/install.html)
  - or via conda:
    - `conda create -n nextflow-23 bioconda::nextflow=23.10.1`
- A container environment: either Docker or Singularity.
  - [Install Docker](https://docs.docker.com/engine/install/)

When running with docker or singularity, nextflow will download dependencies automatically.

## Simple quickstart (docker)

1. Ensure you have nextflow running correctly:
    ```shell
    nextflow -version
    ```
2. Run the pipeline
    ```shell
    nextflow run main.nf --input 'data/' --outdir 'output/'
    ```
    The `--input` parameter indicates where the fastq files to analyze are located. 
    The pipeline is will take all files matching `data/Fastq/*.fq.gz`

    The `--outdir` parameter indicates where the output files should be placed

    The pipeline with run with docker by default
3. To resume a pipeline run that has stopped (or failed), use the `-resume` flag
    ```shell
    nextflow run main.nf --input 'data/' --outdir 'output/' -resume
    ```
   
## Configuring the pipeline
### Run with singularity
Use the `singularity` profile using the `-p` flag.
```shell
nextflow run main.nf --input 'data/' --outdir 'output/' -p singularity
```

### Multithreading
The ability to use multiple cpus for compute intensive processes can be customised using the following configuration files:
- [conf/envs/cluster.config](conf/envs/cluster.config) - when running on my laptop (2 cpus for intensive tasks)
- [conf/envs/cluster.config](conf/envs/cluster.config) - when running on my cluster (8 cpus for intensive tasks)

To run the pipeline with one of these files use the `-c` flag:
```shell
nextflow run main.nf --input 'data/' --outdir 'output/' -c conf/envs/local.config
```

### Primer sequences
The current primers used in this pipeline are **NS5** forward and **LR6** reverse primers ([See UNITE](https://unite.ut.ee/primers.php)).
These have been customized by editing python files required by [Porechop](https://github.com/rrwick/Porechop).

The [modules/local/porechop-custom-adapters/adapters.py](modules/local/porechop-custom-adapters/adapters.py) file is used to customise which primers are trimmed.