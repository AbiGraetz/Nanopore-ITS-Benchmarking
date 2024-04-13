#!/usr/bin/env nextflow
params.input = '$baseDir/data/'
params.outdir = 'output'

include {
  plotQuality;
  plotQuality as plotQuality_pc;
} from './qc/qc_wf.nf'
include { PORECHOP_PORECHOP } from './modules/nf-core/porechop/porechop'

def removeFileEndings(file, String extension, String... rest) {
  for (ext in [extension, *rest]) {
    if (file.endsWith(ext)) {
      return file.replaceAll("$ext\$", "")
    }
  }
  return file
}

workflow {
  ch_raw_reads = Channel.fromPath('$params.input/Fastq/*.fq.gz', checkIfExists: true)

  ch_reads = ch_raw_reads
    .map{ fastq -> 
      [ [id: removeFileEndings(fastq.name, '.fastq.gz', '.fq.gz')], fastq]
    }
  
  plotQuality('raw_reads', ch_reads)

  porechop_out = PORECHOP_PORECHOP(ch_reads)
  plotQuality_pc('porechop', porechop_out.reads)

  // Find primers

  // filter by quality, length etc
}