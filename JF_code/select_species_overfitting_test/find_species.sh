#!/usr/bin/env bash

set -euo pipefail

IMG="its-var-img"        # Docker image
MAX_JOBS=20              


echo $MAX_JOBS

DATA_DIR="$(realpath "${1:-./fastas}")"
OUT_DIR="$(realpath "${2:-./results}")"
mkdir -p "$OUT_DIR"

echo -e "sample\tmean_p_distance" > "$OUT_DIR/p_distance_mean.tsv"

shopt -s nullglob
job_count=0

process_one() {
  local fa="$1"
  local base="$(basename "$fa")"
  # local prefix="${base%%.*}"
  local prefix="${base%.fasta}"
  echo "processing $prefix"

  echo "▶  [$$] processing $base ..."   

  docker run --rm \
    -v "$DATA_DIR":/data:ro \
    -v "$OUT_DIR":/out \
    "$IMG" bash -c "
      set -e
      conda run -n its-var mafft --thread 1 --auto /data/$base > /out/${prefix}_align.fasta
      conda run -n its-var distmat -sequence /out/${prefix}_align.fasta -nucmethod 2 -outfile /out/${prefix}_distmat.txt
  "

}

for fa in "$DATA_DIR"/*.fasta "$DATA_DIR"/*.fa "$DATA_DIR"/*.fna "$DATA_DIR"/*.uniq.fasta; do
  process_one "$fa" &
  # echo "processing $fa" >> name.txt &
  (( ++job_count ))

  if (( job_count % MAX_JOBS == 0 )); then
    echo "test"
    wait -n               
  fi
done

wait                       
echo -e "\n 全部完成！結果位於：$OUT_DIR"

