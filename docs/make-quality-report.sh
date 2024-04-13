#!/bin/bash
qc_dir=${1:-"../output/QC/raw_reads"}

if [ ! -d "$qc_dir" ]; then
    echo "directory $qc_dir not found" 1>&2
    exit -1
fi

write_plots_for_sample () {
    local sample_dir=$1
    local sample=$2

    echo "| Avg. Read Quality vs Length | Base Quality by Position |"
    echo "| -- | -- |"
    echo "<img src=\"$sample_dir/nanoplot/LengthvsQualityScatterPlot_dot.png\"  width=\"400\" /> | <img src=\"$sample_dir/qualityProfile/$sample.png\" width=\"400\"  /> |"
    echo ""
}

write_stats_for_sample() {
    local sample_dir=$1

    echo "| NanoStats | |"
    echo "| -- | -- |"

    local stats=$(grep -E 'number_of_reads|number_of_bases|median_read_length|mean_read_length|read_length_stdev|mean_qual|median_qual|Reads >' \
        "$sample_dir/nanoplot/NanoStats.txt"
    )
    echo "$stats" | while read stat; do
        echo $(echo "$stat" | sed 's/\t/ | /')
    done
    echo ""
}

samples=$(ls -1 "$qc_dir" | sort)
echo "$samples" | while read sample; do
    sample_dir="$qc_dir/$sample"
    echo "# $sample"

    write_plots_for_sample "$sample_dir" "$sample"
    write_stats_for_sample "$sample_dir"
done