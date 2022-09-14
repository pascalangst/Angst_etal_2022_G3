#!/bin/bash

# to make the multiple sequence alignment
# in my pipeline used in the folder "fasta_references" and all "masked" ortholg-alignments

# create directories for the output 
mkdir -p ../fasta_MSA/{FASTA,FASTA_PEP}

# loop to create the multiple sequence alignment
# input are the realigned and masked reference orthologs and the in-group sequences

for fasta in *.aligned_masked.fas

do
newname=${fasta//.aligned_masked.fas/.fasta}

# add the sequences of the "ingroup" samples to the alignments.
# Make sure that they are in capital letters
mafft --add ../fasta_ingroup/$newname --keeplength --preservecase  $fasta > ../fasta_MSA/$newname

done


cd ../fasta_MSA


# final input preparation: mask all the positions of the in-group samples, already masked in the Oc4 sequences
for fasta in *.fasta

do
ref=${fasta//.fasta/_Oc}
gene=${fasta//.fasta/""}

# The positions in the reference sequence with a N should also be N in the ingroup samples:
# define positions with N in the references and take only the one of Oc
python ../originals/generate_masked_ranges.py $fasta  > ranges.txt

grep $ref < ranges.txt > rangesOc.bed

# the actual masking for each sample
samples=$(bcftools query -l VCF.vcf)

cp $fasta $gene.allmasked.fasta

for sample in $samples
do
sed "s/_Oc/|$sample/g" rangesOc.bed > ranges$sample.bed

bedtools maskfasta -fi $gene.allmasked.fasta -bed ranges$sample.bed -fo tmp.allmasked.fasta

mv tmp.allmasked.fasta $gene.allmasked.fasta

done

# define sequence names for popgenome and reduce them to the population
genename=${gene//_/""}

sed -i 's/_//g' $gene.allmasked.fasta
sed -i "s/$genename|/""/g" $gene.allmasked.fasta
sed -i "s/$genename/""/g" $gene.allmasked.fasta

# move the masked multiple sequence alignments to the input folder for popgenome
mv $gene.allmasked.fasta FASTA/

# mask gaps from in all sequences (for translation)

# define positions with -
python ../originals/generate_masked_ranges_gap.py FASTA/$gene.allmasked.fasta > rangesCombined.txt

# generate one range that contains the positions of all samples
samples=$(cat <(echo Oc) <(echo Ocl) <(bcftools query -l VCF.vcf))

for sample in $samples
do
sed -i "s/$sample/combined/g" rangesCombined.txt
done
bedtools sort -i rangesCombined.txt > rangesCombined.bed
bedtools merge -i rangesCombined.bed > rangesCombinedFinal.bed

cp FASTA/$gene.allmasked.fasta $gene.allmaskedgap.fasta

for sample in $samples
do
# we are using the same range for all samples, it must be renamed though
sed "s/combined/$sample/g" rangesCombinedFinal.bed > ranges$sample.bed
# masking with gap character
bedtools maskfasta -fi $gene.allmaskedgap.fasta -bed ranges$sample.bed -fo tmp.allmaskedgap.fasta -mc -

mv tmp.allmaskedgap.fasta $gene.allmaskedgap.fasta
done

# define sequence names and reduce them to the population
sed -i 's/_//g' $gene.allmaskedgap.fasta
sed -i "s/$genename|/""/g" $gene.allmaskedgap.fasta
sed -i "s/$genename/""/g" $gene.allmaskedgap.fasta

mv $gene.allmaskedgap.fasta FASTA_PEP/
done

