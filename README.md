# 10Xmapping

Song Cao

scao@wustl.edu

mapping reference allele and variant allele to 10X barcodes from an input of scRNA-seq bam and the unique 10X barcode corresponds to each individual cell

usage: 

step 1: extracting reads contains the reference allele and variant allele of a somatic variant

perl 10Xmapping.pl --mapq --bam --maf --out

mapq:  the mapping quality; Default 0

bam: the scRNA-seq bam path

maf: the maf file containing a list of somatic variants

out: the output file

step 2: get the corresponding table between ref allele, var allele and 10X barcode 

perl  parse_scrna_bc.pl $fin $fout

fin: the input file from step 1's output

fout: the output file

step 3: output number of reads supporting reference and variant alleles and VAF

perl rc.pl $fin $fout

fin: the input file from step 1's output

fout: the output file
