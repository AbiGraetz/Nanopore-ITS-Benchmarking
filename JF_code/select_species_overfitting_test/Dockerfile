FROM condaforge/mambaforge:latest

RUN mamba create -y -n its-var -c conda-forge -c bioconda \
        fastp seqkit vsearch mafft emboss \
        r-base r-ape \
    && mamba clean -afy

SHELL ["conda", "run", "-n", "its-var", "/bin/bash", "-c"]

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "its-var"]
CMD ["bash"]
