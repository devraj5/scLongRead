screen -S AD_splice_validation

# working directory ---> /data/project/rizzardilab/drbhatta/long_read/AD_EDA/AD_splice_validation

srun --mem-per-cpu=50G --cpus-per-task=16 --time=8:00:00 --partition=largemem --job-name=AD_EDA --pty bash


module load Anaconda3/2023.07-2
micromamba activate splice_validation

## I am using sample 1238 (Ctrl) to extract the splice junction information. There were 2 replicates each with 2 runs.
# Downloading SRA files for sample 1238 rep1 (GSM6619552) Ran on Illumina NovaSeq 6000
prefetch SRR21834823 SRR21834824

# Downloading SRA files for sample 1238 rep2 (GSM6619553) Ran on Illumina NovaSeq 6000
prefetch SRR21834821 SRR21834822

# Converting SRA to FASTQ format
mkdir -p raw_data/
fasterq-dump SRR21834823 --split-files --threads 16 -O raw_data/
fasterq-dump SRR21834824 --split-files --threads 16 -O raw_data/
fasterq-dump SRR21834821 --split-files --threads 16 -O raw_data/
fasterq-dump SRR21834822 --split-files --threads 16 -O raw_data/

##I got just *_4.fastq files - did --include-technical --split-files and got 4 files for each replicate *_1.fastq, *_2.fastq, *_3.fastq, *_4.fastq
##later from 10x page found that *_4.fastq files  are the actual insert reads (read 2 - 90 bp) and *_3.fastq files are the UMI+CB (read 1 - 28 bp)

## so SRR21834823_4.fastq is read 2  and SRR21834823_3.fastq is read 1 for rep1
## SRR21834824_4.fastq is read 2  and SRR21834824_3.fastq is read 1 for rep1
## SRR21834821_4.fastq is read 2  and SRR21834821_3.fastq is read 1 for rep2
## SRR21834822_4.fastq is read 2  and SRR21834822_3.fastq is read 1 for rep2

## Reference genome

fasta: "/data/project/rizzardilab/drbhatta/long_read/references/GRCh38.primary_assembly.genome.fa"
gtf: "/data/project/rizzardilab/drbhatta/long_read/references/gencode.v32.annotation.gtf"

## STAR index   
STAR --runMode genomeGenerate \
 --genomeDir /data/project/rizzardilab/drbhatta/long_read/references/GRCh38_index \
 --genomeFastaFiles /data/project/rizzardilab/drbhatta/long_read/references/GRCh38.primary_assembly.genome.fa \
 --sjdbGTFfile /data/project/rizzardilab/drbhatta/long_read/references/gencode.v32.annotation.gtf \
 --sjdbOverhang 89 \
 --runThreadN 16


## Alignment
mkdir -p starsolo_output/rep1_run1 starsolo_output/rep1_run2 \
         starsolo_output/rep2_run1 starsolo_output/rep2_run2

# Aligning replicate 1 run 1
STAR --genomeDir /data/project/rizzardilab/drbhatta/long_read/references/GRCh38_index \
     --readFilesIn raw_data/SRR21834823_4.fastq raw_data/SRR21834823_3.fastq \
     --soloType CB_UMI_Simple \
     --soloCBwhitelist /data/project/rizzardilab/drbhatta/737K-arc-v1.txt \
     --soloBarcodeReadLength 28 \
     --soloCBstart 1 --soloCBlen 16 \
     --soloUMIstart 17 --soloUMIlen 12 \
     --soloCellFilter EmptyDrops_CR \
     --soloCellReadStats Standard \
     --soloFeatures Gene GeneFull SJ \
     --soloMultiMappers EM \
     --outFileNamePrefix starsolo_output/rep1_run1/1238_rep1_run1_ \
     --outSAMtype BAM SortedByCoordinate \
     --runThreadN 16

# Aligning replicate 1 run 2
STAR --genomeDir /data/project/rizzardilab/drbhatta/long_read/references/GRCh38_index \
     --readFilesIn raw_data/SRR21834824_4.fastq raw_data/SRR21834824_3.fastq \
     --soloType CB_UMI_Simple \
     --soloCBwhitelist /data/project/rizzardilab/drbhatta/737K-arc-v1.txt \
     --soloBarcodeReadLength 28 \
     --soloCBstart 1 --soloCBlen 16 \
     --soloUMIstart 17 --soloUMIlen 12 \
     --soloCellFilter EmptyDrops_CR \
     --soloCellReadStats Standard \
     --soloFeatures Gene GeneFull SJ \
     --soloMultiMappers EM \
     --outFileNamePrefix starsolo_output/rep1_run2/1238_rep1_run2_ \
     --outSAMtype BAM SortedByCoordinate \
     --runThreadN 16


# Aligning replicate 2 run 1
STAR --genomeDir /data/project/rizzardilab/drbhatta/long_read/references/GRCh38_index \
     --readFilesIn raw_data/SRR21834821_4.fastq raw_data/SRR21834821_3.fastq \
     --soloType CB_UMI_Simple \
     --soloCBwhitelist /data/project/rizzardilab/drbhatta/737K-arc-v1.txt \
     --soloBarcodeReadLength 28 \
     --soloCBstart 1 --soloCBlen 16 \
     --soloUMIstart 17 --soloUMIlen 12 \
     --soloCellFilter EmptyDrops_CR \
     --soloCellReadStats Standard \
     --soloFeatures Gene GeneFull SJ \
     --soloMultiMappers EM \
     --outFileNamePrefix starsolo_output/rep2_run1/1238_rep2_run1_ \
     --outSAMtype BAM SortedByCoordinate \
     --runThreadN 16


# Aligning replicate 2 run 2
STAR --genomeDir /data/project/rizzardilab/drbhatta/long_read/references/GRCh38_index \
     --readFilesIn raw_data/SRR21834822_4.fastq raw_data/SRR21834822_3.fastq \
     --soloType CB_UMI_Simple \
     --soloCBwhitelist /data/project/rizzardilab/drbhatta/737K-arc-v1.txt \
     --soloBarcodeReadLength 28 \
     --soloCBstart 1 --soloCBlen 16 \
     --soloUMIstart 17 --soloUMIlen 12 \
     --soloCellFilter EmptyDrops_CR \
     --soloCellReadStats Standard \
     --soloFeatures Gene GeneFull SJ \
     --soloMultiMappers EM \
     --outFileNamePrefix starsolo_output/rep2_run2/1238_rep2_run2_ \
     --outSAMtype BAM SortedByCoordinate \
     --runThreadN 16


# Indexing the BAM files
samtools index starsolo_output/rep1_run1/1238_rep1_run1_Aligned.sortedByCoord.out.bam
samtools index starsolo_output/rep1_run2/1238_rep1_run2_Aligned.sortedByCoord.out.bam
samtools index starsolo_output/rep2_run1/1238_rep2_run1_Aligned.sortedByCoord.out.bam
samtools index starsolo_output/rep2_run2/1238_rep2_run2_Aligned.sortedByCoord.out.bam

# Creating short read junction files list (comma-separated)
SJ_FILES="starsolo_output/rep1_run1/1238_rep1_run1_SJ.out.tab,starsolo_output/rep1_run2/1238_rep1_run2_SJ.out.tab,starsolo_output/rep2_run1/1238_rep2_run1_SJ.out.tab,starsolo_output/rep2_run2/1238_rep2_run2_SJ.out.tab"

# Creating SR_BAM files list (a "file of file names" or fofn)
echo "starsolo_output/rep1_run1/1238_rep1_run1_Aligned.sortedByCoord.out.bam" >> sqanti3_optimized/sr_bam.fofn
echo "starsolo_output/rep1_run2/1238_rep1_run2_Aligned.sortedByCoord.out.bam" >> sqanti3_optimized/sr_bam.fofn
echo "starsolo_output/rep2_run1/1238_rep2_run1_Aligned.sortedByCoord.out.bam" >> sqanti3_optimized/sr_bam.fofn
echo "starsolo_output/rep2_run2/1238_rep2_run2_Aligned.sortedByCoord.out.bam" >> sqanti3_optimized/sr_bam.fofn


## Adapted from https://github.com/SziKayLeung/LOGen/wiki
### I started with the GTF transcript model file from my pipeline rather than the BAM file because:
### I think the transcript model GTF represents assembled transcripts that have already gone through correction, alignment, and transcript calling in the pipeline.
### SQANTI highly recommends to collapse the long read data before running the tool.

mkdir -p LOGen

## This is the long read data
## If i give the gtf file for transcript models, i dont have to give the fasta file. But if i give the fasta file, i need to convert the gtf file to fasta file.
# Creating FASTA from long read transcript models 
gffread /data/project/rizzardilab/drbhatta/long_read/scnanoseq_sqanti/work/9f/d4b3e6872ff1bc54c545571cc8cf5b/OUT/OUT.transcript_models.gtf \
  -g /data/project/rizzardilab/drbhatta/long_read/references/GRCh38.primary_assembly.genome.fa \
  -w LOGen/OUT.long_read_transcripts.fasta

# Aligning FASTA with pbmm2 
pbmm2 align --preset ISOSEQ \
  --sort \
  /data/project/rizzardilab/drbhatta/long_read/references/GRCh38.primary_assembly.genome.fa \
  LOGen/OUT.long_read_transcripts.fasta \
  LOGen/merged_mapped.bam \
  --log-level TRACE \
  --log-file LOGen/merged.log

# Collapsing isoforms with Iso-Seq3
isoseq3 collapse LOGen/merged_mapped.bam \
  LOGen/merged_collapse.gff 



micromamba deactivate

conda activate SQANTI3.env

mkdir -p sqanti3_optimized

# Running SQANTI3 on collapsed long read data with all the short read data
python /data/project/rizzardilab/drbhatta/long_read/SQANTI3/SQANTI3-5.2.2/sqanti3_qc.py \
  LOGen/merged_collapse.gff \
  /data/project/rizzardilab/drbhatta/long_read/references/gencode.v32.annotation.gtf \
  /data/project/rizzardilab/drbhatta/long_read/references/GRCh38.primary_assembly.genome.fa \
  --coverage $SJ_FILES \
  --polyA_motif_list /data/project/rizzardilab/drbhatta/long_read/SQANTI3/mouse_and_human.polyA_motif.txt \
  --CAGE_peak /data/project/rizzardilab/drbhatta/long_read/SQANTI3/SQANTI3-5.2.2/data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed \
  --polyA_peak /data/project/rizzardilab/drbhatta/long_read/SQANTI3/SQANTI3-5.2.2/data/atlas.clusters.2.0.GRCh38.96.bed \
  --SR_bam sqanti3_optimized/sr_bam.fofn \
  --dir sqanti3_optimized \
  --output optimized_sqanti3_results \
  --force_id_ignore \
  --cpus 16 \
  --report both

########################################################

# SQANTI3 using on uncollapsed long read data using one of the short read data
mkdir -p sqanti3_output/rep1_run1 sqanti3_output/rep1_run2 sqanti3_output/rep2_run1 sqanti3_output/rep2_run2

# SQANTI3 on the long read data using one of the  STARsolo SJ.out.tab files
python /data/project/rizzardilab/drbhatta/long_read/SQANTI3/SQANTI3-5.2.2/sqanti3_qc.py \
  --coverage starsolo_output/rep1_run1/1238_rep1_run1_SJ.out.tab \
  --polyA_motif_list /data/project/rizzardilab/drbhatta/long_read/SQANTI3/mouse_and_human.polyA_motif.txt \
  --dir sqanti3_output/rep1_run1 \
  --output 1238_rep1_run1_sc_sqanti3_results \
  /data/project/rizzardilab/drbhatta/long_read/scnanoseq_sqanti/work/9f/d4b3e6872ff1bc54c545571cc8cf5b/OUT/OUT.transcript_models.gtf \
  /data/project/rizzardilab/drbhatta/long_read/references/gencode.v32.annotation.gtf \
  /data/project/rizzardilab/drbhatta/long_read/references/GRCh38.primary_assembly.genome.fa

########################################################

# SQANTI3 on the uncollapsed long read data without short read data

mkdir -p sqanti3_output/rep1_run1/just_gtf

python /data/project/rizzardilab/drbhatta/long_read/SQANTI3/SQANTI3-5.2.2/sqanti3_qc.py \
  --dir sqanti3_output/rep1_run1/just_gtf \
  --output 1238_rep1_run1_sc_sqanti3_results_just_gtf \
  /data/project/rizzardilab/drbhatta/long_read/scnanoseq_sqanti/work/9f/d4b3e6872ff1bc54c545571cc8cf5b/OUT/OUT.transcript_models.gtf \
  /data/project/rizzardilab/drbhatta/long_read/references/gencode.v32.annotation.gtf \
  /data/project/rizzardilab/drbhatta/long_read/references/GRCh38.primary_assembly.genome.fa


###########################################################

# SQANTI3 on the collapsed long read data with without short read data
python /data/project/rizzardilab/drbhatta/long_read/SQANTI3/SQANTI3-5.2.2/sqanti3_qc.py \
  LOGen/merged_collapse.gff \
  /data/project/rizzardilab/drbhatta/long_read/references/gencode.v32.annotation.gtf \
  /data/project/rizzardilab/drbhatta/long_read/references/GRCh38.primary_assembly.genome.fa \
  --dir LOGen/sqanti3_output/single_cell \
  --output merged_sqanti3_results_single_cell


## using gff file from wf-single-cell 

## I did not collapse the long read data because i assume the gff file is already collapsed from the wf-single-cell pipeline
## converting gff to gtf
gffread /data/project/rizzardilab/drbhatta/long_read/wf-single_AD_SC_LR/1238_Ctrl_1/1238_Ctrl_1/1238_Ctrl_1.transcriptome.gff -T -o 1238_Ctrl_1.transcriptome.gtf

screen -S AD_splice_validation

srun --mem-per-cpu=50G --cpus-per-task=16 --time=8:00:00 --partition=largemem --job-name=AD_EDA --pty bash


module load Anaconda3/2023.07-2
conda activate SQANTI3.env

mkdir -p sqanti3_optimized/sqanti_wf_single_cell

# I was having problem with the incompatibility of the gtf file with SQANTI, TSS ratio calculation not possible due to some transcript having no strand information
# looked at similar issues in github and found that latest version has solve the issue so i updated my SQANTI to latest version
#### Navigating to SQANTI3 directory
cd /data/project/rizzardilab/drbhatta/long_read/SQANTI3/

# Backing up the current version
mv SQANTI3-5.2.2 SQANTI3-5.2.2-backup

# Cloning the latest version from GitHub
git clone https://github.com/ConesaLab/SQANTI3.git SQANTI3-latest

# check the version of SQANTI3
python /data/project/rizzardilab/drbhatta/long_read/SQANTI3/SQANTI3-latest/sqanti3_qc.py --version

SJ_FILES="starsolo_output/rep1_run1/1238_rep1_run1_SJ.out.tab,starsolo_output/rep1_run2/1238_rep1_run2_SJ.out.tab,starsolo_output/rep2_run1/1238_rep2_run1_SJ.out.tab,starsolo_output/rep2_run2/1238_rep2_run2_SJ.out.tab"

 # Running SQANTI3 using gtf and fasta file from wf-single-cell
python /data/project/rizzardilab/drbhatta/long_read/SQANTI3/SQANTI3-latest/sqanti3_qc.py \
   1238_Ctrl_1.transcriptome.gtf \
  /data/project/rizzardilab/drbhatta/long_read/references/refdata-gex-GRCh38-2024-A/genes/genes.gtf \
  /data/project/rizzardilab/drbhatta/long_read/references/refdata-gex-GRCh38-2024-A/fasta/genome.fa \
  --coverage $SJ_FILES \
  --polyA_motif_list /data/project/rizzardilab/drbhatta/long_read/SQANTI3/SQANTI3-latest/data/polyA_motifs/mouse_and_human.polyA_motif.txt \
  --CAGE_peak /data/project/rizzardilab/drbhatta/long_read/SQANTI3/SQANTI3-latest/data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed \
  --polyA_peak /data/project/rizzardilab/drbhatta/long_read/SQANTI3/SQANTI3-5.2.2-backup/data/atlas.clusters.2.0.GRCh38.96.bed \
  --SR_bam sqanti3_optimized/sr_bam.fofn \
  --dir sqanti3_optimized/sqanti_wf_single_cell \
  --output sqanti_wf_single_cell_results \
  --force_id_ignore \
  --cpus 16  

  # Running SQANTI3 using gtf and fasta file from gencode
python /data/project/rizzardilab/drbhatta/long_read/SQANTI3/SQANTI3-latest/sqanti3_qc.py \
   1238_Ctrl_1.transcriptome.gtf \
  /data/project/rizzardilab/drbhatta/long_read/references/gencode.v32.annotation.gtf \
  /data/project/rizzardilab/drbhatta/long_read/references/GRCh38.primary_assembly.genome.fa \
  --coverage $SJ_FILES \
  --polyA_motif_list /data/project/rizzardilab/drbhatta/long_read/SQANTI3/SQANTI3-latest/data/polyA_motifs/mouse_and_human.polyA_motif.txt \
  --CAGE_peak /data/project/rizzardilab/drbhatta/long_read/SQANTI3/SQANTI3-latest/data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed \
  --polyA_peak /data/project/rizzardilab/drbhatta/long_read/SQANTI3/SQANTI3-5.2.2-backup/data/atlas.clusters.2.0.GRCh38.96.bed \
  --SR_bam sqanti3_optimized/sr_bam.fofn \
  --dir sqanti3_optimized/sqanti_wf_single_cell \
  --output sqanti_wf_single_cell_results \
  --force_id_ignore \
  --cpus 16 

## The result was not much different. its just that the gtf file needs to be collapsed first before i run sqanti.
