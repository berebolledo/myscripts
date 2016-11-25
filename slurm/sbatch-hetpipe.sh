#!/bin/bash

#SBATCH -C new
##SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH -J pipeline
#SBATCH -t 14-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=berebolledo@gmail.com
#SBATCH -e slurm-%j.err


# -----------
# Definitions
# -----------

source /galaxy/home/boris/.bash_profile

reference="/nfs/brubeck.bx.psu.edu/scratch1/boris/mtprojectreference/mtprojectref.fa"
rCRS="/nfs/brubeck.bx.psu.edu/scratch1/boris/mtprojectreference/rCRS.fasta"
picardTools="/nfs/brubeck.bx.psu.edu/scratch4/boris/software/bin/picardTools"
picard="TMP_DIR=tmp VALIDATION_STRINGENCY=SILENT VERBOSITY=ERROR QUIET=true"

regex="null"
#regex="[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9-]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*."

subject=${1} # SRA id
f1=${subject}_1.fastq.gz
f2=${subject}_2.fastq.gz

readLength=$(zcat ${f1} | awk 'NR%2==0'| \
           awk 'NR%2==1'| head -1000 |   \
           awk '{sum+=length($1)} END{print sum/NR}')

mkdir -p mapped

# ---------
# Alignment
# ---------

bwa mem  \
    -t 4 \
    -M   \
    -R "@RG\tID:${subject}\tLB:amelia\tPL:ILLUMINA\tSM:${subject}"\
    ${reference}  \
    ${f1}         \
    ${f2}         \
    | samtools view -Sb - | bamtools filter -region chrM -out mapped/${subject}.chrM.mapped.bam




# --------------------
# Sort by coordinate I
# --------------------

cd mapped
current=${subject}.chrM.mapped.bam

java -jar ${picardTools}/SortSam.jar \
    ${picard}                        \
    I=${current}                     \
    O=srt.${current}                 \
    SO=coordinate

if [ $? -eq 0 ] ; then
    rm -f ${current}
    tmp_name=${current}
    current=srt.${tmp_name}
fi

# ---------------
# Mark duplicates
# ---------------

java -jar ${picardTools}/MarkDuplicates.jar \
    ${picard}                               \
    I=${current}                            \
    O=marked.${current}                     \
    METRICS_FILE=${current}.metrics         \
    READ_NAME_REGEX=${regex}

if [ $? -eq 0 ] ; then
    rm -f ${current} ${current}.metrics
    tmp_name=${current}
    current=marked.${tmp_name}
fi

# -----------------------------
# Select properly aligned reads
# -----------------------------

bamtools filter                   \
    -region chrM                  \
    -isPaired true                \
    -isProperPair true            \
    -isMateMapped true            \
    -isPrimaryAlignment true      \
    -in ${current}                \
    -out proper.${current}

if [ $? -eq 0 ] ; then
    rm -f ${current}
    tmp_name=${current}
    current=proper.${tmp_name}
fi


# --------------
# Sort by name I
# --------------

java -jar ${picardTools}/SortSam.jar  \
    ${picard}                         \
    I=${current}                      \
    O=query.${current}                \
    SO=queryname

if [ $? -eq 0 ] ; then
    rm -f ${current}
    tmp_name=${current}
    current=query.${tmp_name}
fi


# -----------------------------------------------
# Remove chimeric reads and select by read length
# -----------------------------------------------

rm_chim_in_pair.py ${current} ${readLength}

if [ $? -eq 0 ] ; then
    rm -f ${current}
    tmp_name=${current}
    current=dechim.rlen.${tmp_name}
fi

# -----------
# Realignment
# -----------

samtools view -b ${current} | bamleftalign -f ${rCRS} |  \
samtools view -o realigned.${current} -b - 


if [ $? -eq 0 ] ; then
    rm -f ${current}
    tmp_name=${current}
    current=realigned.${tmp_name}
fi

# ---------------------
# Sort by coordinate II
# ---------------------

java -jar ${picardTools}/SortSam.jar  \
    ${picard}                         \
    I=${current}                      \
    O=srt.${current}                  \
    SO=coordinate

if [ $? -eq 0 ] ; then
    rm -f ${current}
    tmp_name=${current}
    current=srt.${tmp_name}
fi


# --------------
# Index BAM file
# --------------

samtools index ${current}

# -------------------------------
# Calculate major allele sequence
# -------------------------------

angsd               \
    -i ${current}   \
    -doHetPlas 2    \
    -out ${current} \
    -nThreads 4     \
    -minQ 30        \
    -minMapQ 20     \
    -r chrM         \
    -nLines 10000   \
    -doCounts 1     \
    -dumpCounts 3   \
    -doQsDist 1     \
    -howOften 10

mitoMajorFromThorGL.py ${current}.hetGL
mkdir -p fasta
mv ${current}.hetGL.major.fa fasta/
rm -f ${current}.*

# ---------------------------------------------------------------
# Recalculate NM tag (SAM). Number of mismatches to the reference
# ---------------------------------------------------------------


samtools fillmd                      \
    -b ${current}                    \
    fasta/${current}.hetGL.major.fa  \
    2>/dev/null > md.${current}


if [ $? -eq 0 ] ; then
    rm -f ${current}
    tmp_name=${current}
    current=md.${tmp_name}
fi


# ---------------
# Sort by name II
# ---------------

java -jar ${picardTools}/SortSam.jar \
    ${picard}                        \
    I=${current}                     \
    O=query.${current}               \
    SO=queryname

if [ $? -eq 0 ] ; then
    rm -f ${current}
    tmp_name=${current}
    current=query.${tmp_name}
fi

#----------------------------------------------
# Select reads on the (N)number of (M)ismatches
# ---------------------------------------------

nm-ratio.select.py ${current}

if [ $? -eq 0 ] ; then
    rm -f ${current}
    tmp_name=${current}
    current=nm-ratio.${tmp_name}
fi

#-------
# Rename
# ------

mv ${current} ${subject}.nvcReady.bam
current=${subject}.nvcReady.bam 


# ----------------------
# Sort by coordinate III
# ----------------------


java -jar ${picardTools}/SortSam.jar \
    ${picard}                        \
    I=${current}                     \
    O=srt.${current}                 \
    SO=coordinate

if [ $? -eq 0 ] ; then
    rm -f ${current}
    tmp_name=${current}
    current=srt.${tmp_name}
fi

# ------------------------
# Allele count calculation
# ------------------------ 

mkdir -p nvcReady
mv ${current} nvcReady
cd nvcReady

samtools index ${current}

naive_variant_caller.py  \
    -b ${current}        \
    -i ${current}.bai    \
    -o ${subject}.vcf    \
    -r ${rCRS}           \
    -s -q 30 -m 20       \
    -t uint32            \
    --region chrM


allele-counts.py         \
    -i ${subject}.vcf    \
    -o ${subject}.counts \
    -n -s

