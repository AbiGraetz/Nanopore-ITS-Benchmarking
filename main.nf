#!/usr/bin/env nextflow
include {
  nanoplot;
  nanoplot as nanoplot_agg;
  plotQualityProfile;
  plotQualityProfile as plotQualityProfile_agg;
} from './qc.nf'

params.input = "$baseDir/data/"
params.outdir = 'output'

def collectWithId(id, ch) {
  ch.collect { meta, read -> read }
    .map { [ [id: id], it ] }
}

def removeFileEndings(file, String extension, String... rest) {
  for (ext in [extension, *rest]) {
    if (file.endsWith(ext)) {
      return file.replaceAll("$ext\$", "")
    }
  }
  return file
}

workflow {
  ch_raw_reads = Channel.fromPath("$params.input/Fastq/*.fq.gz", checkIfExists: true)

  ch_reads = ch_raw_reads
    .map{ fastq -> 
      [ [id: removeFileEndings(fastq.name, ".fastq.gz", ".fq.gz")], fastq]
    }
  
  plotQuality('raw_reads', ch_reads)
}

workflow plotQuality {
  take:
    stage_name
    reads // [ Map{id, stage_dir}, FastqFile (1+) ]
    
  main:
      rs = reads.map{meta, reads -> [meta + [stage_dir: stage_name], reads]} 
      
      // sample level quality plots
      nanoplot(rs)
      plotQualityProfile(rs)

      combined_reads = collectWithId('all', reads)
        .map{meta, reads -> [meta + [stage_dir: stage_name], reads]}

      // aggregated quality plots
      plotQualityProfile_agg(combined_reads)
      nanoplot_agg(combined_reads)
}