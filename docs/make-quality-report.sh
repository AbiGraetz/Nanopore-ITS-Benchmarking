#!/bin/bash
qc_dir=${1:-"../output/QC/raw_reads"}

if [ ! -d "$qc_dir" ]; then
    echo "directory $qc_dir not found" 1>&2
    exit -1
fi

write_plots_for_sample () {
    local sample_name=$1
    local img_dir="$qc_dir/$sample_name"

    echo "# $sample_name"
    echo "| Avg. Read Quality vs Length | Base Quality by Position |"
    echo "| -- | -- |"
    echo "<img src=\"$img_dir/nanoplot/LengthvsQualityScatterPlot_dot.png\"  width=\"400\" /> | <img src=\"$img_dir/qualityProfile/$sample_name.png\" width=\"400\"  /> |"
}

samples=$(ls -1 "$qc_dir" | sort)
echo "$samples" | while read sample; do
    write_plots_for_sample "$sample"
done