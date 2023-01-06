#!/bin/bash

source activate QC

#move all samples into one file called raw


cd ../../storehouse/gbaoteng/bombus/

mkdir raw

mkdir trimmed

mkdir IndexReads

mkdir AlignedReads

find -iname '*.fastq.gz' -exec mv '{}' raw/ \;

cat list.txt | parallel -j 4 "fastp --thread 2 -c -r -j -h -p -i raw/{}_R1_001.fastq.gz -I raw/{}_R2_001.fastq.gz -o trimmed/{}_R1_001_trimmed.fastq.gz -O trimmed/{}_R2_001_trimmed.fastq.gz"

echo "QC and trimming process done"

#echo "Going ahead to check for rRNA contamination"
#mkdir rrnaContam

#bwa index bombus_rRNA.fasta

#cat list.txt | parallel -j 1 "bwa mem -t 8 -P reference/GCF_000188095.3_BIMP_2.2_genomic.fna trimmed/{}_R1_001_trimmed.fastq.gz trimmed/{}_R2_001_trimmed.fastq.gz | samtools view - @ 4 -bS -o rnaContam/{}_rrna.bam && samtools flagstat -@ 8 rrnaConta/{}_rrna.bam > rrnaContam/{}_rrna.out"

#mkdir microContam
#cat list.txt | parallel -j -2 "kranken2 --db minikraken_8GB_20200312 --threads 4 --use-names --memory-mapping --output microContam/{}.out --report miroContam/{}.report --paired starAligned/{}Unmapped.out.mate1 starAligned/{}Unmapped.out.mate2"


STAR --runThreadN 16 --runMode genomeGenerate --genomeDir starIndex --genomeFastaFiles reference/GCF_000188095.3_BIMP_2.2_genomic.fna --sjdbGTFfile reference/GCF_000188095.3_BIMP_2.2_genomic.gff

cat list.txt | parallel -j 2 "gunzip trimmed/{}_R*_001.fastq.gz && STAR --runThreadN 32 --genomeDir starIndex --readFilesIn trimmed/{}_R1_001_trimmed.fastq trimmed/{}_R2_001_trimmed.fastq --outFilterIntronMotifs RemoveNoncanonical --outFileNamePrefix AlignedReads/{} --outSAMtype BAM SortedByCoordinate && pigz trimmed/{}_*R"

#bwa mem -t 8 reference/GCF_000188095.3_BIMP_2.2_genomic.fna


#cat list.txt | parallel -j 1 "parallel gunzip ::: trimmed/{}_R*.gz && STAR --runThreadN 64 --genomeLoad LoadAndKeep --genomeDir IndexReads/ --readFilesIn trimmed/{}_R1 trimmed/{}_R2 --outFilterIntronMotifs RemoveNoncanonical --outFileNamePrefix AlignedReads/{} --outSAMtype BAM SortedByCoordinate && pigz trimmed/{}_*R"

mkdir CountFiles

featureCounts -T 16 -a reference/GCF_000188095.3_BIMP_2.2_genomic.gff -g "ID" -o CountFiles/CountFiles.txt AlignedReads/*.bam

cat list.txt | parallel -j 2 "stringtie -p 16 -G reference/GCF_000188095.3_BIMP_2.2_genomic.gff -o gft/{}.gff -l {} AlignedReads/{}Aligned.sortedByCoord.out.bam"
stringtie --merge -p 16 -G ../reference/GCF_000188095.3_BIMP_2.2_genomic.gff -o Stringtie_merged.gff *.gff 
cat list.txt | parallel -j 2 "stringtie -e -B -p 8 -G gft/Stringtie_merged.gff -o CountFiles/{}.gff AlignedReads/{}Aligned.sortedByCoord.out.bam"

$ scp gbaoteng@134.129.113.53/cd../../storehouse/gbaoteng/bombus/Aligned.sortedByCoord.out.bam* /path/to/directory_on_laptop


