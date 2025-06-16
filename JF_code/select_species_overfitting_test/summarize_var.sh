#!/usr/bin/env bash

INPUT_FOLDER=$1
OUTFILE="variance_summary.tsv"

echo -e "File\tVariance" > "$OUTFILE"

for f in ${INPUT_FOLDER}/*.txt; do
  var=$(awk '
    BEGIN { row=0 }
    /^[^#]/ && NF>1 {
        row++
        for(i=2; i<=NF; i++){
        if($i ~ /^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/){
            d = $i         
            sum  += d
            sum2 += d*d
            cnt++
        }
        }
    }
    END {
        if(cnt>0){
        mean = sum/cnt
        var  = sum2/cnt - mean*mean
        printf "pairs=%d\tmean=%.4f\tvar=%.4f\n", cnt, mean, var
        } else {
        print "NA"
        }
    }
    ' "$f")

  echo -e "${f}\t${var}" >> "$OUTFILE"
done

echo "saved in $OUTFILE"
