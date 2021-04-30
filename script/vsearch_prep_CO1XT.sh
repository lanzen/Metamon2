# v1.1 2020-07-08
fprimer=GGWACWRGWTGRACWNTNTAYCCYCC

rprimer=TANACYTCNGGRTGNCCRAARAAYCA

fprimlength=26

#Reverse primer - reverse complement!
rprimer_rc=TGRTTYTTYGGNCAYCCNGARGTNTA

mkdir merged_reads

for d in $*; do
    cd $d
    i=${d//\/}
    
#    gunzip -kf *.fastq.gz

    echo "Sample: $d"
    echo "Sample: $d" >> ../readprep.log
    echo `pwd`

    # Merge the reads (non-staggering since adaptors removed) max 40 difference\s (defualt 5 - but low qual V3 and >40 does not improve anyway)

    vsearch --fastq_mergepairs *_R1*.fastq.gz --reverse *_R2*.fastq.gz --fastq_maxdiffs 40 --fastq_allowmergestagger --fastqout ${i}_m.fastq 2>> ../readprep.log


    # Allow up to 2 mismatches and one missing base in primers. Max length 450, min 300

    cutadapt -g $fprimer -O 26 --max-n=0 --discard-untrimmed -o ${i}_trimmed1.fastq ${i}_m.fastq >> ../readprep.log

    cutadapt -a $rprimer_rc -O 26 --discard-untrimmed -o ${i}_trimmed.fastq ${i}_trimmed1.fastq >>../readprep.log

    # Min overlap of 215 corresponds to 333 nt and max of 257 to 291 for 300 nt reads, but since some reads are shorter 274 occurs too using flash
    vsearch --fastq_filter ${i}_trimmed.fastq --fastaout ${i}_merged_QF.fasta --fastq_maxee 1 --fastq_minlen 274 --fastq_maxlen 333  2>>../readprep.log

       
    #Fix read names to include "barcode" label
    python ~/script/drive5/fixreads.py ${i}_merged_QF.fasta $i > ../merged_reads/${i}_reads_fixed.fasta

    rm *.fasta

    cd ..

done
ls

echo "Preparing FASTQC reports"

cat */*_m.fastq > all_merged.fastq
fastqc all_merged.fastq
rm all_merged.fastq */*_m.fastq

cat */*_trimmed.fastq > all_trimmed.fastq
fastqc all_trimmed.fastq
rm all_trimmed.fastq */*_trimmed.fastq

rm */*.fastq
rm *fastqc.zip

