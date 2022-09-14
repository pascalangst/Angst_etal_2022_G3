# get CDS as multiple sequence alignments

samples=$(bcftools query -l VCF.vcf)

mkdir originals
mkdir tmp
mkdir fasta_ingroup
mkdir -p fasta_references/{fasta_MdSP,fasta_Ocl}

# ingroup #
###########

# prepare the reference genome for gatk
# samples is a vector of samplenames in the VCF

for sample in $samples
do

	# new vcf containing only one sample
	GenomeAnalysisTK -T SelectVariants -R Reference.fasta -o $sample.vcf -V VCF.vcf -sn $sample

	# new "reference" for the one sample 
	GenomeAnalysisTK -T FastaAlternateReferenceMaker -R Reference.fasta -o $sample.fasta -V $sample.vcf --use_IUPAC_sample $sample

    # Replace blanks by dashes.
	sed -E -e "s/([A,T,C,G])   ([A,T,C,G])/\1---\2/g;s/([A,T,C,G])  ([A,T,C,G])/\1--\2/g;s/([A,T,C,G]) ([A,T,C,G]) ([A,T,C,G])/\1-\2-\3/g;s/([A,T,C,G]) ([A,T,C,G])/\1-\2/g;s/^ ([A,T,C,G])/-\1/g;s/([A,T,C,G]) $/\1-/g" $sample.fasta > $sample.gap2-.fasta

	# cut the gene sequences out of the genome-fasta according to the gff for the sample
	perl originals/gff2fasta.pl $sample.gap2-.fasta GFF.gff3 $sample

	# prepare new sample fasta (only containing coding sequences)
	sed -i -E "s/(>.*)/\1|$sample/g" $sample.cds.fasta
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $sample.cds.fasta > $sample.cds.oneline.fasta
    sed -i 1d $sample.cds.oneline.fasta

	# create a fasta file for each of the genes
	awk -F "|" '/^>/ {close(F) ; F = $1".fasta"} {print >> F}' $sample.cds.oneline.fasta

	# move all unneeded files to the tmp folder
	mv $sample* tmp/

done

for f in *FUN*; do mv "$f" "$(echo "fasta_ingroup/$f" | sed -e 's/>//')";done

for sample in FUN_00*; do awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $sample > $sample.oneline.fasta; done
paste FUN_00*oneline.fasta > all.fasta
sed -E "s/    //g" all.fasta > all2.fasta
sed -E 's/(>[^>]*)>.*/\1/g' all2.fasta > all.fasta
sed -i 1d all.fasta


##########################################################################################
# Reference plus outgroup #
###########################

# Prepare the reference for cutting the orthologs

# cut the gene sequences out of the reference according to the gff
perl originals/gff2fasta.pl Reference.fasta GFF.gff3 Reference

# prepare reference fasta (only containing coding sequences)
sed -i -E "s/(>.*)/\1|Oc/g" Reference.cds.fasta
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < Reference.cds.fasta  > Reference.cds.oneline.fasta
# cut orthologs
for i in `awk '{if ($3 == 1) print $5}' orthologs.csv | sed -e 's/-T1//'`;do awk 'BEGIN{RS=">"} /'$i'/{print ">" $0}' Reference.cds.oneline.fasta;done > Reference.orthologs.fasta

# create a fasta file for each of the genes
awk -F "|" '/^>/ {close(F) ; F = $1".fasta"} {print >> F}' Reference.orthologs.fasta

# move the gene-fastas to the fasta_references folder
for f in *FUN*; do cp "$f" "$(echo "fasta_references/fasta_Oc/$f" | sed -e 's/>//')";done

# prepare the outgroup fasta (only containing coding sequences)
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < Outgroup.fasta > cds_Ocl.oneline.fna

# cut orthologs
for i in `awk '{if ($3 == 1) print $4}' orthologs.csv`;do awk 'BEGIN{RS=">"} /'$i'/{print ">" $0}' cds_Ocl.oneline.fna;done > Ocl.orthologs.fasta
# translate gene IDs from Rozella to Md
paste <(awk '{if ($3 == 1) print $4}' orthologs.csv) <(awk '{if ($3 == 1) print $5}' orthologs.csv | sed -e 's/-T1//') | while read a b; do sed -i "s/$a/$b/" Ocl.orthologs.fasta; done

# add the genes to the fastas of the reference genes and move them to the common folder
awk -F "|" '/^>/ {close(F) ; F = $1".fasta"} {print >> F}' Ocl.orthologs.fasta
for f in *FUN*; do mv "$f" "$(echo "fasta_references/$f" | sed -e 's/>//;s/.fasta//')";done

# move the gene-fastas to the fasta_references folder and the rest to the tmp folder
awk -F "|" '/^>/ {close(F) ; F = $1".fasta"} {print >> F}' Ocl.orthologs.fasta
for f in *FUN*; do mv "$f" "$(echo "fasta_references/fasta_Ocl/$f" | sed -e 's/>//')";done
mv Reference* tmp/
mv *Ocl* tmp/


##########################################################################################
# multiple sequence alignment #
###############################

cd fasta_references

ls FUN_* | parallel -j 6 "Rscript ../originals/EstimateKaKs2_modified.Rscript {}"
ls *aligned* | parallel -j 6 "Rscript ../originals/mask_alignment.Rscript {}"
ls *aligned_masked* | parallel -j 6 "Rscript ../originals/mask_paml.Rscript {}"
bash ../originals/multiple_sequence_alignment.sh 

cd ../fasta_MSA/FASTA

for sample in FUN_00*; do awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $sample > $sample.oneline.fasta; done
paste FUN_00*oneline.fasta > all.fasta
sed -E "s/    //g" all.fasta > all2.fasta
sed -E 's/(>[^>]*)>.*/\1/g' all2.fasta > all.fasta
sed -i 1d all.fasta
