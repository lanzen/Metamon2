#!/bin/bash
# v 44 2019-03-17 -Jon modified, added shebang, Anders testing the effect of not removing same length in front and other

fprimer=GCTTGWCTCAAAGATTAAGCC #cutadapt also cuts any subsequent/preceeding bases (so will remove random 12)
fprimlength=21
rprimer_rc=TCYAAGGAAGGCAGCAGG #Reverse complement of CCTGCTGCCTTCCTTRGA! (most bases removed before!)

#Make directory where the merged read files will be placed (fastq and fasta)
mkdir merged_reads

#Iterate over all subdirectories, assume they have the sample name, enter and
#uncompress fastq files
for d in $*; do #The $* means all *arguments*, $0, $1 etc signifies arguments given when executing script.
    echo $d
    cd $d
    pwd
    i=${d//\/} #Is d with slashes removed
    echo $i

 
    # Merge the reads (non-staggering since adaptors removed) max 5 differences

    echo $d >> ../readprep.log

    vsearch --fastq_mergepairs *_R1*.fastq.gz --reverse *_R2*.fastq.gz  --fastq_allowmergestagger --fastqout ${i}_m.fastq 2>> ../readprep.log

    cutadapt -g $fprimer -m 100 -O $fprimlength --discard-untrimmed -o ${i}_trimmed1.fastq ${i}_m.fastq -u 12 >> ../readprep.log 
    cutadapt -a $rprimer_rc -m 100 --discard-untrimmed -o ${i}_trimmed.fastq ${i}_trimmed1.fastq >>../readprep.log

    rm  ${i}_trimmed1.fastq

    #Jon commented out try no length trimming 
    #vsearch --fastq_filter ${i}_trimmed.fastq --fastaout ${i}_merged_QF.fasta --fastq_maxee 1 --fastq_trunclen $crop_length 2>>../readprep.log    

    vsearch --fastq_filter ${i}_trimmed.fastq --fastaout ${i}_merged_QF.fasta --fastq_maxee 1 --fastq_minlen 330 --fastq_maxlen 450  2>>../readprep.log
    
    #Fix read names to include "barcode" label
    fixreads.py ${i}_merged_QF.fasta $i > ../merged_reads/${i}_reads_fixed.fasta

    rm ${i}_merged_QF.fasta

    cd ..
done

echo "Preparing FASTQC reports"

cat */*_m.fastq > all_merged.fastq
fastqc all_merged.fastq
rm all_merged.fastq

cat */*_trimmed.fastq > all_trimmed.fastq
fastqc all_trimmed.fastq
rm all_trimmed.fastq

rm */*.fastq
rm *fastqc.zip

for h in *.html; do gnome-open $h; done
