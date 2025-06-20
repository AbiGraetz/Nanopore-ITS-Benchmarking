#!/usr/bin/env nextflow
params.input = '$baseDir/data/'
params.outdir = 'output'

include {
  plotQuality;
  plotQuality as plotQuality_pc;
  plotQuality as plotQuality2;
  plotQuality as plotQuality3;
  plotQuality as plotQuality4;
  plotQuality as plotQuality_pc_custom;
} from './qc/qc_wf.nf'
include { PORECHOP_PORECHOP as porechop_adapters } from './modules/nf-core/porechop/porechop'
include { porechop_custom_adapters as porechop_primers } from './modules/local/porechop-custom-adapters'
include { 
  chopper as chopper_15;
  chopper as chopper_17;
  chopper as chopper_20;
} from './qc/chopper.nf'

def removeFileEndings(file, String extension, String... rest) {
  for (ext in [extension, *rest]) {
    if (file.endsWith(ext)) {
      return file.replaceAll("$ext\$", "")
    }
  }
  return file
}

workflow {
  ch_raw_reads = Channel.fromPath("$params.input", checkIfExists: true)

  ch_reads = ch_raw_reads
    .map{ fastq -> 
      [ [id: removeFileEndings(fastq.name, '.fastq.gz', '.fq.gz')], fastq]
    }
  
  plotQuality('raw_reads', ch_reads)


  porechop_out = porechop_adapters(ch_reads).reads // trim ONT adapters and barcodes
    | porechop_primers // trim primers (customized)

  plotQuality_pc('porechop', porechop_out.reads)

  q15 = chopper_15(
    params.readQuality,
    porechop_out.reads)
  
  plotQuality2('Qmin15', q15.reads)

  q17 = chopper_17(
    [minQualityPhred: 17, minLength: params.readQuality.minLength, maxLength: params.readQuality.maxLength],
    porechop_out.reads)
  
  plotQuality3('Qmin17', q17.reads)

  q20 = chopper_20(
    [minQualityPhred: 20, minLength: params.readQuality.minLength, maxLength: params.readQuality.maxLength],
    porechop_out.reads)
  
  plotQuality4('Qmin20', q20.reads)
}