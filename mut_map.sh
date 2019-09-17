#!/bin/bash

## Mapping mutation from maf file to sc bam data to get which cell has which ref/var allele
## Usage: bash run_mut_map.sh <scRNA.bam> <maf> <annotation> <output.dir>


BAM=$1;shift
MAF=$1;shift
ANNOTATION=$1;shift
OUTPUT_DIR=$1;shift
SAMPLE=`echo $BAM | cut -f8 -d "/"`
SCRIPT_DIR=/diskmnt/Projects/Users/simonmo/tools/lijun_mapping_pipeline
#step1
cd ${OUTPUT_DIR}
perl ${SCRIPT_DIR}/10Xmapping.pl --mapq 20 --bam ${BAM} --maf ${MAF} --out ${SAMPLE}_mapping.txt
#step2
perl ${SCRIPT_DIR}/parse_scrna_bc.pl ${SAMPLE}_mapping.txt ${SAMPLE}_mapping_parse.txt
#step3
cat ${ANNOTATION} | tail -n+2 | awk '{print "'"${SAMPLE}"'" "\t" $1 "\t" $2}' > ${SAMPLE}_annotations.txt #You need to close the ' string, insert the shell ${SAMPLE} (inside " in case there are special characters), then reopen the ' string.
mv ${SAMPLE}_mapping_parse.txt ${SAMPLE}.txt
perl ${SCRIPT_DIR}/add_celltype_bc.pl ${SAMPLE}_annotations.txt ${SAMPLE}.txt ${SAMPLE}_map_ct_bc.txt
#step4
cat ${ANNOTATION} | tail -n+2 | cut -f2 | sort -u | while read CELL_TYPE
do
	perl ${SCRIPT_DIR}/refvar_heatmap_cell.pl ${CELL_TYPE} ${SAMPLE}_map_ct_bc.txt ${MAF} ${SAMPLE}_mapping_heatmap_${CELL_TYPE}.txt
done

