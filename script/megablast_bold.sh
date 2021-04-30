blastn -task megablast -query $1 -db ~/projects/PhyloRefDB/BOLD/bold.fasta -num_alignments 100 -outfmt 5 -out ${1//.fa*}_bold.xml -num_threads 3
gzip ${1//.fa*}_bold.xml
