#######################
#minimap2
#! /bin/bash
fq=$1
sampleid=$2
workdir=$3
fqdir=$4
ref=$5

clean_fq=$workdir/clean.$sampleid.fq.gz
sam=$workdir/${sampleid}.sam
bam=$workdir/${sampleid}.bam
clean_fastq_plots=$workdir/${sampleid}.clean-fastq-plots
bam_sort=$workdir/${sampleid}_minimap2_sorted.bam
bam_sort_tag=$workdir/${sampleid}_minimap2_sorted_tag.bam

echo "$(date) 1. Start QC: $sampleid"
gunzip -c $fq |/biosoft/NanoFilt -q 7  -l 1000 |gzip > $clean_fq
echo "$(date) 1. Finish QC: $sampleid"

echo "$(date) 2. Start NanoPlot: $sampleid"
/biosoft/NanoPlot  --fastq $clean_fq \
   -o $clean_fastq_plots \
   --maxlength 40000 \
   -t 20 \
   --plots hex dot
echo "$(date) 2. Finish NanoPlot: $sampleid"

echo "$(date) 3. Start mapping: $sampleid"
/biosoft/minimap2/minimap2 --MD -ax map-ont -L -t 20 $ref $clean_fq > $sam
echo "$(date) 3. Finish mapping: $sampleid"

echo "$(date) 3. Start bam: $sampleid"
/biosoft/samtools-1.9/samtools view -@ 20 -bS $sam -o $bam
echo "$(date) 3. Finish bam: $sampleid"

echo "$(date) 3. Start Delete sam: $sampleid"
rm $sam
echo "$(date) 3. Finish Delete sam: $sampleid"

echo "$(date) 3. Start sort: $sampleid"
/biosoft/samtools-1.9/samtools sort  -@ 20 $bam -o $bam_sort
echo "$(date) 3. Finish sort: $sampleid"

echo "$(date) 4. Start add tag: $sampleid"
(/biosoft/samtools-1.9/samtools view -H "$bam_sort"; /biosoft/samtools-1.9/samtools view "$bam_sort"|awk 'BEGIN{OFS="\t"}{$1="tumor"$1; print $0}')|/biosoft/samtools-1.9/samtools view -S -b - -o "$bam_sort_tag"
echo "$(date) 4. Finish add tag: $sampleid"

echo "$(date) 5. Start Delete bam: $sampleid"
rm $bam
echo "$(date) 5. Finish Delete bam: $sampleid"

echo "$(date) 6. Start Delete bam_sort: $sampleid"
rm $bam_sort
echo "$(date) 6. Finish Delete bam_sort: $sampleid"

