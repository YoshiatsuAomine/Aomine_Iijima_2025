#WSL Ubuntu 24.04 LTS
#Anaconda3 Python 3.12.2

conda install wget

conda install -c bioconda fastp

conda install -c milaboratories mixcr

#-------

namelist=("CD4_Br_S6" "CD4_SK_S8" "CD4_SP_S1" "CD4_Stem_S7" "CD4_TG_S10" "CD8_Br_S1" "CD8_SK_S4" "CD8_SP_S5" "CD8_Stem_S3" "CD8_TG_S2" "Na_SP_CD4_S32" "Na_SP_CD8_S33" "Na_TG_CD4_S30" "Na_TG_CD8_S31")
namelist=("Tet_TG+Br_CD8_S1_L001")
namelist=("CD4_Br_TG_S6_L001")

cd ./projectfile

mkdir ./01fastq

mkdir ./02fastp

cd ./02fastp

mkdir report

cd ../01fastq

for ID in "${namelist[@]}"; do
          fastp -i ./${ID}_R1_001.fastq.gz \
                -I ./${ID}_R2_001.fastq.gz \
                -o ../02fastp/${ID}_1.fq.gz \
                -O ../02fastp/${ID}_2.fq.gz \
                --detect_adapter_for_pe \
                --length_required 30 \
                -j ../02fastp/report/${ID}.json \
                -h ../02fastp/report/${ID}.html
          done

#MiXCR
mixcr activate-license
E-GQXAPOLWVNRPQHUSJFPMDPZLWXEVINQUJOEMPDOIGIIVSGCO

mkdir ./03mixcr

cd ./02fastp

for ID in "${namelist[@]}"; do
    mixcr analyze takara-mouse-rna-tcr-smarter\
    ${ID}_1.fq.gz ${ID}_2.fq.gz ../03mixcr_SMARTer/${ID}_output
    done
#or
for ID in "${namelist[@]}"; do
    mixcr analyze generic-amplicon\
    --species mmu\
    --rna\
    --rigid-left-alignment-boundary\
    --floating-right-alignment-boundary C\
    --assemble-clonotypes-by CDR3\
    ${ID}_1.fq.gz ${ID}_2.fq.gz ../03mixcr/${ID}_output
    done


#root/projectfile/rawdata/ some fastq.gz of data
#get reference data
cd ./projectfile
mkdir ref
#https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.27/
#download GCF_000001635.27_GRCm39_genomic.gtf and GCF_000001635.27_GRCm39_genomic.fna


#Trimmomatic or Trim_galore
#Trimmomatic
wget https://github.com/usadellab/Trimmomatic/files/5854859/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip

mkdir 02_trimmed

cd ./01_rawdata

#threads 18 on hikidaLAB CPU
files=(`ls -1 ./`)
    for ID in "${files[@]}"
        do
        trimmomatic SE -threads 18 -phred33 \
        -trimlog ../02_trimmed/${ID}_log.txt \
        -summary ../02_trimmed/${ID}_sum.txt \
        ${ID} \
        ../02_trimmed/trim_${ID} \
        ILLUMINACLIP:../Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:50
        done

#FastQC
cd ../projectfile

mkdir 03_FastQC

cd ./02_trimmed

files=(`ls -1 *.gz ./`)
for ID in "${files[@]}"
    do
    fastqc -t 10 --nogroup -o ../03_FastQC -f fastq ${ID}
    done

# trim_galore
mkdir 02_trimmed

cd ./01_rawdata

files=(`ls -1 ./`)
for ID in "${files[@]}"
    do
    trim_galore --fastqc --cores 4 --output_dir ../02_trimmed  ${ID}
    done


#STAR
cd ../

cd ./00_ref

gunzip GCF_000001635.27_GRCm39_genomic.fna.gz
zcat GCF_000001635.27_GRCm39_genomic.gtf.gz | grep -v _PAR_Y > GCF_000001635.27_GRCm39_genomic.gtf

STAR --runMode genomeGenerate -- genomeDir star_rsem -- runThreadN 18 \
--sjdbOverhang 99 --genomeFastaFiles GCF_000001635.27_GRCm39_genomic.fna \
--sjdbGTFfile GCF_000001635.27_GRCm39_genomic.gtf

tail -n 2 ./star_rsem/Log.out

cd ../

mkdir 04_trimmed_log_sum

cd ./02_trimmed

files=(`ls -1 *.gz ./`)
for ID in "${files[@]}"
    do
    mkdir ../03_mapping/${ID}
    STAR --genomeDir ../00_ref/star_rsem --runThreadN 18 \
    --outFileNamePrefix ../03_mapping/${ID}/ --quantMode TranscriptomeSAM \
    --outSAMtype BAM SortedByCoordinate \
    --readFilesIn ${ID} --readFilesCommand gunzip -c
    done

cd ../05_mapping

files=(`ls -1 *.gz ./`)
for ID in "${files[@]}"
    do
    cat ${ID}/Log.final.out
    done

#RSEM
cd ./RSEM-1.3.3

rsem-prepare-reference --gtf ../projectfile/ref/GCF_000001635.27_GRCm39_genomic.gtf -p 10 \
../projectfile/ref/GCF_000001635.27_GRCm39_genomic.fna ../projectfile/ref/star_rsem/mg

cd ../projectfile/ref/star_rsem

ls mg.*
#return  [mg.chrlist  mg.grp  mg.idx.fa  mg.n2g.idx.fa  mg.seq  mg.ti  mg.transcripts.fa]

cd ./projectfile

mkdir rsem

cd ./mapping

files=(`ls -1 *.gz ./`)
for ID in "${files[@]}"
    do
    rsem-calculate-expression --alignments -p 16 \
    --forward-prob=0 --append-names --estimate-rspd \
    ${ID}/Aligned.toTranscriptome.out.bam ../ref/star_rsem/mg ../rsem/${ID}
    done

cd ./mapping

files=(`ls -1 *.gz ./`)
for ID in "${files[@]}"
    do
    ls ../rsem/${ID}.*
    done

cd ../rsem

files=(`ls -1 *.genes.results`)

rsem-generate-data-matrix \
${files[@]} \
>240501iijimasensei_STAR_genes_count_matrix.tsv

cd ../rsem

files=(`ls -1 *.isoforms.results`)
rsem-generate-data-matrix \
${files[@]} \
>240501iijimasensei_STAR_isoforms_count_matrix.tsv
