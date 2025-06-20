include {
  nanoplot;
  nanoplot as nanoplot_agg;
  plotQualityProfile;
  plotQualityProfile as plotQualityProfile_agg;
} from './qc.nf'

def collectWithId(id, ch) {
  ch.collect { meta, read -> read }
    .map { [ [id: id], it ] }
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