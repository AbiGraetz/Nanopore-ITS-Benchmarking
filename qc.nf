process plotQualityProfile {
  tag "$meta.id"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
          'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.30.0--r43hf17093f_0' :
          'biocontainers/bioconductor-dada2:1.30.0--r43hf17093f_0' }"
  
  publishDir "${params.outdir}/QC/${meta.stage_dir}/${meta.id}/qualityProfile", mode: 'copy'

  input:
      tuple val(meta), path(fq_files)
  output:
      tuple val(meta), path("*.png"), emit: plots

  script:
  def r_files = "$fq_files".tokenize().collect{"\'$it\'"}.join(',')
  """
  #!/usr/bin/env Rscript
  library(dada2)
  library(ggplot2)

  plot <- plotQualityProfile(c($r_files), aggregate = TRUE)
  ggsave('${meta.id}.png')
  """
}

process nanoplot {
  tag "$meta.id"

  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://depot.galaxyproject.org/singularity/nanoplot:1.41.6--pyhdfd78af_0' :
      'biocontainers/nanoplot:1.41.6--pyhdfd78af_0' }"
  
  publishDir "${params.outdir}/QC/${meta.stage_dir}/${meta.id}/nanoplot", mode: 'copy'

  input:
      tuple val(meta), path(ontfiles, name: "*.fastq.gz")

  output:
      tuple val(meta), path("*.html")                , emit: html
      tuple val(meta), path("*.png") , optional: true, emit: png
      tuple val(meta), path("*.txt")                 , emit: txt
      tuple val(meta), path("*.log")                 , emit: log

  script:
  def args = task.ext.args ?: ''
  """
  NanoPlot \\
      $args \\
      -t $task.cpus \\
      --fastq $ontfiles
  """
}