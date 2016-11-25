#!/bin/bash

#SBATCH -C new
#SBATCH --ntasks=25
##SBATCH --cpus-per-task=8
#SBATCH -J hets
#SBATCH -t 14-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=berebolledo@gmail.com

source /galaxy/home/boris/.bash_profile
subject=${1}


###############################################################################

reference="/galaxy/home/boris/boris1/mtprojectreference/mtprojectref"

f1=${subject}_1.fastq.gz
f2=${subject}_2.fastq.gz
mkdir -p mapped

bwa7 mem \
    -t 8 \
    -M   \
    -R "@RG\tID:${subject}\tLB:amelia\tPL:ILLUMINA\tSM:${subject}"\
    ${reference}      \
    ${f1}             \
    ${f2} 2>/dev/null|\
    samtools view -Sb - > mapped/${subject}.mapped.bam 2>/dev/null

###############################################################################

cd mapped
file=${subject}.mapped.bam
name=${subject}


#regex="[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9-]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*."
regex="null"
ref="/home/boris/boris1/mtprojectreference/rCRS.fasta"
picardTools="/galaxy/home/boris/boris4/software/bin/picardTools"

###############################################################################

echo "----determine read length started..."

readLength=`samtools view ${file} |head -100|awk '{sum+=length($10)} END{print sum/NR}'`

echo "----determine read length finished..."
###############################################################################

echo "----sort by coordinate started..."
java -jar ${picardTools}/SortSam.jar \
    ${picard}                        \
    I=${file}                        \
    O=srt.${file}                    \
    SO=coordinate
echo "----sort by coordinate finished..."
###############################################################################

echo "----markDup started..."
java -jar ${picardTools}/MarkDuplicates.jar \
    ${picard}                               \
    I=srt.${file}                           \
    O=marked.srt.${file}                    \
    METRICS_FILE=srt.${file}.metrics        \
    READ_NAME_REGEX=${regex}
echo "----markDup finished..."
###############################################################################

echo "----Selection of proper started..."
bamtools filter                   \
    -region chrM                  \
    -isPaired true                \
    -isProperPair true            \
    -isMateMapped true            \
    -isPrimaryAlignment true      \
    -in marked.srt.${file}        \
    -out proper.marked.srt.${file}
echo "----Selection of proper finished..."
###############################################################################

echo "----sort by name started..."
java -jar ${picardTools}/SortSam.jar  \
    ${picard}                         \
    I=proper.marked.srt.${file}       \
    O=query.proper.marked.srt.${file} \
    SO=queryname
echo "----sort by name finished..."
###############################################################################

echo "----dechim and rlen started..."
rm_chim_in_pair.py query.proper.marked.srt.${file} ${readLength}
echo "----dechim and rlen finished..."
###############################################################################

echo "----realignment started..."
samtools view -b dechim.rlen.query.proper.marked.srt.${file}|               \
bamleftalign -f ${ref}|                                                     \
samtools view -b - > realigned.dechim.rlen.query.proper.marked.srt.${file}
echo "----realignment finished..."
###############################################################################

echo "----sort by coordinate started..."
java -jar ${picardTools}/SortSam.jar                            \
    ${picard}                                                   \
    I=realigned.dechim.rlen.query.proper.marked.srt.${file}     \
    O=srt.realigned.dechim.rlen.query.proper.marked.srt.${file} \
    SO=coordinate
echo "----sort by coordinate finished..."
###############################################################################

echo "----index bam started..."
samtools index srt.realigned.dechim.rlen.query.proper.marked.srt.${file}
echo "----index bam finished..."
###############################################################################

echo "----major sequence started..."
#get_major_from_bam.py srt.realigned.dechim.rlen.query.proper.marked.srt.${file}
name2="srt.realigned.dechim.rlen.query.proper.marked.srt.${file}"
angsd              \
    -i ${name2}    \
    -doHetPlas 2   \
    -out ${name2}  \
    -nThreads 4    \
    -minQ 30       \
    -minMapQ 20    \
    -r chrM        \
    -nLines 10000  \
    -doCounts 1    \
    -dumpCounts 3  \
    -doQsDist 1    \
    -howOften 10

mitoMajorFromThorGL.py ${name2}.hetGL
echo "----major sequence finished..."
###############################################################################

echo "----recalculate NM started..."
samtools fillmd               \
    -b ${name2}               \
    ${name2}.hetGL.major.fa   \
    2>/dev/null > md.${name2}
echo "----recalculate NM finished..."
###############################################################################

echo "----sort by name 2 started..."
java -jar ${picardTools}/SortSam.jar \
    ${picard}                        \
    I=md.${name2}                    \
    O=query.md.${name2}              \
    SO=queryname
echo "----sort by name 2 finished..."
###############################################################################

echo "----NM selection started..."
nm-ratio.select.py query.md.${name2}
echo "----NM selection finished..."
###############################################################################

mv nm-ratio.query.md.${name2} ${name}.nvcReady.bam
###############################################################################

echo "----sort by coordinate started..."
java -jar ${picardTools}/SortSam.jar \
    ${picard}                        \
    I=${name}.nvcReady.bam           \
    O=srt.${name}.nvcReady.bam       \
    SO=coordinate
echo "----sort by coordinate finished..."
###############################################################################

echo "----cleaning started..."
mkdir -p nvcReady
mv srt.${name}.nvcReady.bam nvcReady
mkdir -p fasta
mv *.fa fasta/
rm -fr *srt.${file}.* *.nvcReady.*
echo "----cleaning finished..."
###############################################################################

echo "----NVC started..."
cd nvcReady
samtools index srt.${name}.nvcReady.bam

naive_variant_caller.py             \
    -b srt.${name}.nvcReady.bam     \
    -i srt.${name}.nvcReady.bam.bai \
    -o ${name}.vcf                  \
    -r ${ref}                       \
    -s -q 30 -m 20                  \
    -t uint32                       \
    --region chrM &>/dev/null
echo "----NVC finished..."
###############################################################################

echo "----Allele counts started..."
allele-counts.py        \
    -i ${name}.vcf      \
    -o ${name}.counts   \
    -n -s
echo "----Allele counts finished..."
