filename=$1

# remove adaptors
if [ ! -f "$filename"_R1_raw_trimmo_paired_truseq3-2_2_10_5_1.fastq ]; then
echo "...removing adaptors" > run_"$filename"_unstranded.log
java -jar /home/Program_NGS_sl-pw-srv01/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 6 \
-trimlog trimmolog1.txt \
"$filename"_R1_raw.fastq \
"$filename"_R2_raw.fastq \
"$filename"_R1_raw_trimmo_paired_truseq3-2_2_10_5_1.fastq \
"$filename"_R1_raw_trimmo_unpaired_truseq3-2_2_10_5_1.fastq \
"$filename"_R2_raw_trimmo_paired_truseq3-2_2_10_5_1.fastq \
"$filename"_R2_raw_trimmo_unpaired_truseq3-2_2_10_5_1.fastq \
ILLUMINACLIP:/home/Program_NGS_sl-pw-srv01/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa:2:10:5:1 2> trimmolog2.txt
else
echo "...adaptors already removed" > run_"$filename"_unstranded.log
fi 



#  map with tophat
echo "...mapping with tophat" >> run_"$filename"_unstranded.log
mkdir tophat_results
tophat --phred64-quals \
--max-multihits 1 \
--num-threads 6 \
--library-type fr-unstranded \
--GTF /home/Reference_genomes/Arabidopsis_thaliana_19genomes_Ler0/consolidated_annotation.Ler_0.gtf \
--transcriptome-index /home/Reference_genomes/Arabidopsis_thaliana_19genomes_Ler0/consolidated_annotation.Ler_0 \
--output-dir tophat_results/"$filename"_trimmo_19genomes_ler_nomixed_unstranded/ \
--no-mixed \
/home/Reference_genomes/Arabidopsis_thaliana_19genomes_Ler0/ler_0.v7 \
"$filename"_R1_raw_trimmo_paired_truseq3-2_2_10_5_1.fastq "$filename"_R2_raw_trimmo_paired_truseq3-2_2_10_5_1.fastq 2> tophatlog1.txt

# rename and sort
mv tophat_results/"$filename"_trimmo_19genomes_ler_nomixed_unstranded/accepted_hits.bam \
tophat_results/"$filename"_trimmo_19genomes_ler_nomixed_unstranded/"$filename"_trimmo_paired_2_10_5_1_tophat_19genomes_ler_nomixed_unstranded.bam

echo "...sorting mapped reads" >> run_"$filename"_unstranded.log
samtools sort tophat_results/"$filename"_trimmo_19genomes_ler_nomixed_unstranded/"$filename"_trimmo_paired_2_10_5_1_tophat_19genomes_ler_nomixed_unstranded.bam \
tophat_results/"$filename"_trimmo_19genomes_ler_nomixed_unstranded/"$filename"_trimmo_paired_2_10_5_1_tophat_19genomes_ler_nomixed_unstranded_sorted 2>> run_"$filename"_unstranded.log


# mark and remove duplicates 
echo "...removing duplicates" >> run_"$filename"_unstranded.log
java -Xmx4g -jar /home/Program_NGS_sl-pw-srv01/picard-tools-1.103/MarkDuplicates.jar \
INPUT=tophat_results/"$filename"_trimmo_19genomes_ler_nomixed_unstranded/"$filename"_trimmo_paired_2_10_5_1_tophat_19genomes_ler_nomixed_unstranded_sorted.bam \
OUTPUT=tophat_results/"$filename"_trimmo_19genomes_ler_nomixed_unstranded/"$filename"_trimmo_paired_2_10_5_1_tophat_19genomes_ler_nomixed_unstranded_sorted_rmdup_picard.bam \
METRICS_FILE=dup.txt VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true 2> markdup_stderr.txt


# indexing
samtools index tophat_results/"$filename"_trimmo_19genomes_ler_nomixed_unstranded/"$filename"_trimmo_paired_2_10_5_1_tophat_19genomes_ler_nomixed_unstranded_sorted_rmdup_picard.bam

# word count raw data
raw_line=$(wc -l < "$filename"_R1_raw.fastq) 
raw_count=$(echo "$raw_line/4" | bc -l)
echo "number of raw reads: $raw_count" >> run_"$filename"_unstranded.log

# get flagstat
echo "number of clean reads mapped:" >> run_"$filename"_unstranded.log
samtools flagstat tophat_results/"$filename"_trimmo_19genomes_ler_nomixed_unstranded/"$filename"_trimmo_paired_2_10_5_1_tophat_19genomes_ler_nomixed_unstranded_sorted_rmdup_picard.bam >> run_"$filename"_unstranded.log

# estimate genome average and normalise
echo "...normalising reads" >> run_"$filename"_unstranded.log
genomeCoverageBed -split -bg -ibam tophat_results/"$filename"_trimmo_19genomes_ler_nomixed_unstranded/"$filename"_trimmo_paired_2_10_5_1_tophat_19genomes_ler_nomixed_unstranded_sorted_rmdup_picard.bam \
-g /home/Reference_genomes/Arabidopsis_thaliana_19genomes_Ler0/chrom_size_19genomes_Ler0.txt \
> tophat_results/"$filename"_trimmo_19genomes_ler_nomixed_unstranded/"$filename"_trimmo_paired_2_10_5_1_tophat_19genomes_ler_nomixed_unstranded_sorted_rmdup_picard.bedgraph


sum=$(samtools depth tophat_results/"$filename"_trimmo_19genomes_ler_nomixed_unstranded/"$filename"_trimmo_paired_2_10_5_1_tophat_19genomes_ler_nomixed_unstranded_sorted_rmdup_picard.bam | awk '{sum+=$3;cnt++}END{printf "%.0f", sum}')
sum_norm=$(echo "$sum/117706417" | bc -l)
echo "genome normalised coverage: $sum_norm" >> run_"$filename"_unstranded.log

export MYVAR=$sum_norm
perl -e 'print $ENV{MYVAR}."\n"'

# normalise read counts by genome-wide coverage
perl -ne 'chomp($_); @a=split(/\t/,$_);print $a[0]."\t".$a[1]."\t".$a[2]."\t".$a[3]/$ENV{MYVAR}."\n";' \
tophat_results/"$filename"_trimmo_19genomes_ler_nomixed_unstranded/"$filename"_trimmo_paired_2_10_5_1_tophat_19genomes_ler_nomixed_unstranded_sorted_rmdup_picard.bedgraph \
> tophat_results/"$filename"_trimmo_19genomes_ler_nomixed_unstranded/"$filename"_trimmo_paired_2_10_5_1_tophat_19genomes_ler_nomixed_unstranded_sorted_rmdup_picard_genomenorm.bedgraph

# convert bedgraph to bigwig
bedGraphToBigWig tophat_results/"$filename"_trimmo_19genomes_ler_nomixed_unstranded/"$filename"_trimmo_paired_2_10_5_1_tophat_19genomes_ler_nomixed_unstranded_sorted_rmdup_picard_genomenorm.bedgraph \
/home/Reference_genomes/Arabidopsis_thaliana_19genomes_Ler0/chrom_size_19genomes_Ler0.txt \
tophat_results/"$filename"_trimmo_19genomes_ler_nomixed_unstranded/"$filename"_trimmo_paired_2_10_5_1_tophat_19genomes_ler_nomixed_unstranded_sorted_rmdup_picard_genomenorm.bw

echo "...calling FPKMs" >> run_"$filename"_unstranded.log

cufflinks --output-dir cufflinks_results/cufflinks_"$filename"_trimmo_19genomes_ler_nomixed_unstranded/ \
--GTF /home/Reference_genomes/Arabidopsis_thaliana_19genomes_Ler0/consolidated_annotation.Ler_0_mod2.gff3 \
--num-threads 8 \
--frag-bias-correct /home/Reference_genomes/Arabidopsis_thaliana_19genomes_Ler0/ler_0.v7.fas \
--multi-read-correct \
--library-type fr-unstranded \
tophat_results/"$filename"_trimmo_19genomes_ler_nomixed_unstranded/"$filename"_trimmo_paired_2_10_5_1_tophat_19genomes_ler_nomixed_unstranded_sorted_rmdup_picard.bam

# obtain raw reads using HTseq-count
echo "...calling raw reads" >> run_"$filename"_unstranded.log
samtools sort -n tophat_results/"$filename"_trimmo_19genomes_ler_nomixed_unstranded/"$filename"_trimmo_paired_2_10_5_1_tophat_19genomes_ler_nomixed_unstranded_sorted_rmdup_picard.bam \
tophat_results/"$filename"_trimmo_19genomes_ler_nomixed_unstranded/"$filename"_trimmo_paired_2_10_5_1_tophat_19genomes_ler_nomixed_unstranded_sorted_rmdup_picard_name_sorted 2>> run_"$filename"_unstranded.log

htseq-count -r name -s no -f bam -t exon -i gene_id \
tophat_results/"$filename"_trimmo_19genomes_ler_nomixed_unstranded/"$filename"_trimmo_paired_2_10_5_1_tophat_19genomes_ler_nomixed_unstranded_sorted_rmdup_picard_name_sorted.bam \
/home/Reference_genomes/Arabidopsis_thaliana_19genomes_Ler0/consolidated_annotation.Ler_0_mod2.gtf \
> tophat_results/"$filename"_trimmo_19genomes_ler_nomixed_unstranded/"$filename"_trimmo_paired_2_10_5_1_tophat_19genomes_ler_nomixed_unstranded_sorted_rmdup_picard_name_sorted_htseq_count.ct
 
# combine FPKM/TPM/RAW
echo "...combining FPKM/TPM/Raw" >> run_"$filename"_unstranded.log
cp cufflinks_results/cufflinks_"$filename"_trimmo_19genomes_ler_nomixed_unstranded/genes.fpkm_tracking .
cut -f1,4,7,10 genes.fpkm_tracking > genes.fpkm_tracking_fpkm
cp tophat_results/"$filename"_trimmo_19genomes_ler_nomixed_unstranded/"$filename"_trimmo_paired_2_10_5_1_tophat_19genomes_ler_nomixed_unstranded_sorted_rmdup_picard_name_sorted_htseq_count.ct genes.raw
 
Rscript -e 'genes.fpkm <- read.table("genes.fpkm_tracking_fpkm", sep="\t", quote="", header=T, colClasses = c(rep("character", 3), "numeric"))' \
-e 'genes.fpkm$TPM <-genes.fpkm$FPKM/sum(genes.fpkm$FPKM)*10^6' \
-e 'genes.raw <- read.table("genes.raw", sep="\t", quote="", header=F, colClasses = c("character", "numeric"))' \
-e 'names(genes.raw) <- c("tracking_id", "Raw")' \
-e 'genes.fpkm.raw <- merge(genes.fpkm, genes.raw, by = "tracking_id")' \
-e 'write.table(genes.fpkm.raw, "genes_fpkm_tpm_raw.txt", sep="\t", quote=F, row.names=F)' 

mv genes_fpkm_tpm_raw.txt "$filename"_trimmo_paired_2_10_5_1_tophat_19genomes_ler_nomixed_unstranded_sorted_rmdup_picard_combined_read.txt
rm genes.fpkm_tracking
rm genes.fpkm_tracking_fpkm
rm genes.raw




