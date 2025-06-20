#!/usr/bin/env bash
qc_dir=${1:-"../output/QC"}

main() {
    if [ ! -d "$qc_dir" ]; then
        echo "directory $qc_dir not found" 1>&2
        exit 1
    fi

    declare -A stages=( ["raw_reads"]="Raw reads" \
        ["porechop"]="Porechop (adapters, barcodes and primers)" \
        ["Qmin15"]="Filtering min Q15" ["Qmin17"]="Filtering min Q17"\
        ["Qmin20"]="Filtering min Q20")
    stage_order=(raw_reads porechop Qmin15 Qmin17 Qmin20)
    
    for stage in "${stage_order[@]}"; do
        mkdir -p ./sample
        write_stats_table_for_all_samples "$stage" > "./sample/$stage.md"

        echo "# ${stages[$stage]}"
        write_plots_for_sample "$qc_dir/$stage" "all"
        write_stats_for_sample "$qc_dir/$stage/all"
        echo "[Sample details](./sample/$stage.md)"
        echo ""
    done
}

write_plots_for_sample () {
    local stage_dir=$1
    local sample=$2
    local sample_dir="$stage_dir/$sample"

    echo "| Avg. Read Quality vs Length | Base Quality by Position |"
    echo "| -- | -- |"
    echo "| <img src=\"$sample_dir/nanoplot/LengthvsQualityScatterPlot_dot.png\"  width=\"500\" /> | <img src=\"$sample_dir/qualityProfile/$sample.png\" width=\"500\"  /> |"
    echo ""
}

write_stats_for_sample() {
    local sample_dir=$1

    echo "| NanoStats | |"
    echo "| -- | -- |"

    declare -A stats_dict="($(awk -vRS="\n" -vFS="\t" -vORS=" " 'NR> 1 {print "[\""$1"\"]=\""$2"\"" }' "$sample_dir/nanoplot/NanoStats.txt"))"

    metrics=(number_of_reads number_of_bases median_read_length mean_read_length \
        read_length_stdev mean_qual median_qual "Reads >Q5:" "Reads >Q7:" "Reads >Q10:" "Reads >Q12:" "Reads >Q15:")
    for metric in "${metrics[@]}"
    do
        echo "| $metric | ${stats_dict[$metric]} |"
    done
    echo ""
}

write_stats_table_for_all_samples() {
    local stage_name=${1:-raw_reads}
    local samples=$(ls -1 "$qc_dir/$stage_name" | sort)

    echo "# $stage_name (sample details)"

    metrics=(number_of_reads number_of_bases median_read_length mean_read_length \
            read_length_stdev mean_qual median_qual "Reads >Q5:" "Reads >Q7:" "Reads >Q10:" "Reads >Q12:" "Reads >Q15:")

    echo "| Sample |$(IFS='|' ; echo "${metrics[*]}")|"
    local len=${#metrics[@]}
    echo "|$(for ((i=0;i<=len;i++));do printf " -- |"; done)"


    echo "$samples" | while read -r sample; do
        sample_dir="$qc_dir/$stage_name/$sample"
        declare -A stats_dict="($(awk -vRS="\n" -vFS="\t" -vORS=" " 'NR> 1 {print "[\""$1"\"]=\""$2"\"" }' "$sample_dir/nanoplot/NanoStats.txt"))"
        echo "| $sample | $(for metric in "${metrics[@]}";do printf '%s | ' "${stats_dict[$metric]}"; done)"
    done
}

write_plots_and_stats_for_all_samples() {
    local stage_name=${1:-raw_reads}
    local samples=$(ls -1 "$qc_dir/$stage_name" | sort)
    echo "# $stage_name (sample details)"

    echo "$samples" | while read -r sample; do
        echo "## $sample"

        write_plots_for_sample "$qc_dir/$stage_name/" "$sample"
        write_stats_for_sample "$qc_dir/$stage_name/$sample"
    done
}

main > report.md
