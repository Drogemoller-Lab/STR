#/bin/bash

awk 'BEGIN{OFS="\t"} $1=="chr5" && $3=="gene" && $4 <= 73782133 && $5 >= 73782133' gencode.v47.basic.annotation.gtf
