#!/bin/bash
sed 's/ /\t/g; s/;/\t/g; s/"/\t/g' genes.gtf >genes_temp.gtf
awk 'BEGIN{FS="\t"; OFS="\t"}{if($3 == "gene") print $26 ,$1, $4, $5}' genes_temp.gtf >gene_order.txt
