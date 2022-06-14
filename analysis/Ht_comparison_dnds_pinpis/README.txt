# First bash code, followed by R code down below.

# dn/ds

join -1 2 -2 1 <(join <(sort groups.txt) <(sort PAML_DnDs_Estimates_aligned_masked_FIOER33-pacbio_BEOM2.tsv) | sort -k 2) <(sort PAML_DnDs_Estimates_aligned_masked_Oc4_Ocl.tsv) > joined_tables.txt

awk 'NF==14 {print $0}' joined_tables.txt > joined_tables_cleaned.txt 

join -a1 <(sort joined_tables_cleaned.txt) <(sort single_copy_busco_sequences_Oc.txt) > joined_tables_Oc_busco_marked.txt
join -a1 -1 2 -2 1 <(sort -k 2 joined_tables_cleaned.txt) <(sort single_copy_busco_sequences_ht.txt) > joined_tables_Ht_busco_marked.txt

grep BUSCO joined_tables_Ht_busco_marked.txt >  joined_tables_Ht_busco.txt
grep -v BUSCO joined_tables_Ht_busco_marked.txt >  joined_tables_Ht_non-busco.txt

grep BUSCO joined_tables_Oc_busco_marked.txt >  joined_tables_Oc_busco.txt
grep -v BUSCO joined_tables_Oc_busco_marked.txt >  joined_tables_Oc_non-busco.txt

# dn/ds Ht is in column 6 
# dn/ds Oc is in column 12

#Ht busco: 0.2176
#Ht non-busco: 0.2606
#
#Oc busco: 0.2958
#Oc non-busco: 0.3149

# pin/pis

join -1 2 -2 1 <(join <(sort groups.txt) <(sort selectionStats_Ht_ingroup.txt) | sort -k 2) <(sort selectionStats_Oc_ingroup.txt) > joined_tables_pi.txt

awk 'NF==16 {print $0}' joined_tables_pi.txt > joined_tables_pi_cleaned.txt 

join -a1 <(sort joined_tables_pi_cleaned.txt ) <(sort single_copy_busco_sequences_Oc.txt) > joined_tables_pi_Oc_busco_marked.txt
join -a1 -1 2 -2 1 <(sort -k 2 joined_tables_pi_cleaned.txt ) <(sort single_copy_busco_sequences_ht.txt) > joined_tables_pi_Ht_busco_marked.txt

grep BUSCO joined_tables_pi_Ht_busco_marked.txt >  joined_tables_pi_Ht_busco.txt
grep -v BUSCO joined_tables_pi_Ht_busco_marked.txt >  joined_tables_pi_Ht_non-busco.txt

grep BUSCO joined_tables_pi_Oc_busco_marked.txt >  joined_tables_pi_Oc_busco.txt
grep -v BUSCO joined_tables_pi_Oc_busco_marked.txt >  joined_tables_pi_Oc_non-busco.txt

# pin/pis Ht is in column (6/7)/(8/9)
# pin/pis Oc is in column (13/14)/(15/16) 

#Ht busco: 
#Ht non-busco: 
#
#Oc busco: 
#Oc non-busco: 




join -a1 <(sort selectionStats_Oc_ingroup.txt) <(sort single_copy_busco_sequences_Oc.txt) > pi_Oc_busco_marked.txt
join -a1 <(sort selectionStats_Ht_ingroup.txt) <(sort single_copy_busco_sequences_Ht.txt) > pi_Ht_busco_marked.txt

grep BUSCO  pi_Ht_busco_marked.txt >   pi_Ht_busco.txt
grep -v BUSCO  pi_Ht_busco_marked.txt >   pi_Ht_non-busco.txt

grep BUSCO  pi_Oc_busco_marked.txt >   pi_Oc_busco.txt
grep -v BUSCO  pi_Oc_busco_marked.txt >   pi_Oc_non-busco.txt















############## R code

joined_tables_pi_Ht_busco <- read.table("/Volumes/USB PASCAL/dnds_files/joined_tables_pi_Ht_busco.txt", quote="\"", comment.char="")
median(na.rm = T, (joined_tables_pi_Ht_busco$V6/joined_tables_pi_Ht_busco$V7)/(joined_tables_pi_Ht_busco$V8/joined_tables_pi_Ht_busco$V9))

joined_tables_pi_Ht_nonbusco <- read.table("/Volumes/USB PASCAL/dnds_files/joined_tables_pi_Ht_non-busco.txt", quote="\"", comment.char="")
median(na.rm = T, (joined_tables_pi_Ht_nonbusco$V6/joined_tables_pi_Ht_nonbusco$V7)/(joined_tables_pi_Ht_nonbusco$V8/joined_tables_pi_Ht_nonbusco$V9))

joined_tables_pi_Oc_busco <- read.table("/Volumes/USB PASCAL/dnds_files/joined_tables_pi_Oc_busco.txt", quote="\"", comment.char="")
median(na.rm = T, ((joined_tables_pi_Oc_busco$V13/joined_tables_pi_Oc_busco$V14)/(joined_tables_pi_Oc_busco$V15/joined_tables_pi_Oc_busco$V16))[is.finite((joined_tables_pi_Oc_busco$V13/joined_tables_pi_Oc_busco$V14)/(joined_tables_pi_Oc_busco$V15/joined_tables_pi_Oc_busco$V16))])

joined_tables_pi_Oc_nonbusco <- read.table("/Volumes/USB PASCAL/dnds_files/joined_tables_pi_Oc_non-busco.txt", quote="\"", comment.char="")
median(na.rm = T, ((joined_tables_pi_Oc_nonbusco$V13/joined_tables_pi_Oc_nonbusco$V14)/(joined_tables_pi_Oc_nonbusco$V15/joined_tables_pi_Oc_nonbusco$V16))[is.finite((joined_tables_pi_Oc_nonbusco$V13/joined_tables_pi_Oc_nonbusco$V14)/(joined_tables_pi_Oc_nonbusco$V15/joined_tables_pi_Oc_nonbusco$V16))])


(sum(joined_tables_pi_Ht_busco$V6)/sum(joined_tables_pi_Ht_busco$V7))/(sum(joined_tables_pi_Ht_busco$V8)/sum(joined_tables_pi_Ht_busco$V9))
(sum(joined_tables_pi_Ht_nonbusco$V6)/sum(joined_tables_pi_Ht_nonbusco$V7))/(sum(joined_tables_pi_Ht_nonbusco$V8)/sum(joined_tables_pi_Ht_nonbusco$V9))

(sum(joined_tables_pi_Oc_busco$V6)/sum(joined_tables_pi_Oc_busco$V7))/(sum(joined_tables_pi_Oc_busco$V8)/sum(joined_tables_pi_Oc_busco$V9))
(sum(joined_tables_pi_Oc_nonbusco$V6)/sum(joined_tables_pi_Oc_nonbusco$V7))/(sum(joined_tables_pi_Oc_nonbusco$V8)/sum(joined_tables_pi_Oc_nonbusco$V9))




joined_tables_pi_Ht_busco$Ht_ratio<-(joined_tables_pi_Ht_busco$V6/joined_tables_pi_Ht_busco$V7)/(joined_tables_pi_Ht_busco$V8/joined_tables_pi_Ht_busco$V9)
joined_tables_pi_Ht_busco$Oc_ratio<-(joined_tables_pi_Ht_busco$V13/joined_tables_pi_Ht_busco$V14)/(joined_tables_pi_Ht_busco$V15/joined_tables_pi_Ht_busco$V16)

#joined_tables_pi_Ht_busco <- joined_tables_pi_Ht_busco[!(joined_tables_pi_Ht_busco$Ht_ratio==Inf | joined_tables_pi_Ht_busco$Oc_ratio==Inf | is.na(joined_tables_pi_Ht_busco$Ht_ratio) | is.na(joined_tables_pi_Ht_busco$Oc_ratio)),]
joined_tables_pi_Ht_busco[ is.na(joined_tables_pi_Ht_busco$Ht_ratio),]$Ht_ratio <- 0 
joined_tables_pi_Ht_busco[is.na(joined_tables_pi_Ht_busco$Oc_ratio),]$Oc_ratio <- 0

joined_tables_pi_Oc_nonbusco$Ht_ratio<-(joined_tables_pi_Oc_nonbusco$V6/joined_tables_pi_Oc_nonbusco$V7)/(joined_tables_pi_Oc_nonbusco$V8/joined_tables_pi_Oc_nonbusco$V9)
joined_tables_pi_Oc_nonbusco$Oc_ratio<-(joined_tables_pi_Oc_nonbusco$V13/joined_tables_pi_Oc_nonbusco$V14)/(joined_tables_pi_Oc_nonbusco$V15/joined_tables_pi_Oc_nonbusco$V16)

#joined_tables_pi_Oc_nonbusco <- joined_tables_pi_Oc_nonbusco[!(joined_tables_pi_Oc_nonbusco$Ht_ratio==Inf | joined_tables_pi_Oc_nonbusco$Oc_ratio==Inf | is.na(joined_tables_pi_Oc_nonbusco$Ht_ratio) | is.na(joined_tables_pi_Oc_nonbusco$Oc_ratio)),]
joined_tables_pi_Oc_nonbusco[ is.na(joined_tables_pi_Oc_nonbusco$Ht_ratio),]$Ht_ratio <- 0 
joined_tables_pi_Oc_nonbusco[ is.na(joined_tables_pi_Oc_nonbusco$Oc_ratio),]$Oc_ratio <- 0



median(joined_tables_pi_Ht_busco$Ht_ratio)
median(joined_tables_pi_Ht_busco$Oc_ratio)
median(joined_tables_pi_Oc_nonbusco$Ht_ratio)
median(joined_tables_pi_Oc_nonbusco$Oc_ratio)












pi_Ht_busco <- read.table("/Volumes/USB PASCAL/dnds_files/pi_Ht_busco.txt", quote="\"", comment.char="")

pi_Ht_nonbusco <- read.table("/Volumes/USB PASCAL/dnds_files/pi_Ht_non-busco.txt", quote="\"", comment.char="", skip = 1)

pi_Oc_busco <- read.table("/Volumes/USB PASCAL/dnds_files/pi_Oc_busco.txt", quote="\"", comment.char="")

pi_Oc_nonbusco <- read.table("/Volumes/USB PASCAL/dnds_files/pi_Oc_non-busco.txt", quote="\"", comment.char="", skip = 1)


(sum(pi_Ht_busco$V5)/sum(pi_Ht_busco$V6))/(sum(pi_Ht_busco$V7)/sum(pi_Ht_busco$V8))
(sum(pi_Ht_nonbusco$V5)/sum(pi_Ht_nonbusco$V6))/(sum(pi_Ht_nonbusco$V7)/sum(pi_Ht_nonbusco$V8))

(sum(pi_Oc_busco$V5)/sum(pi_Oc_busco$V6))/(sum(pi_Oc_busco$V7)/sum(pi_Oc_busco$V8))
(sum(pi_Oc_nonbusco$V5)/sum(pi_Oc_nonbusco$V6))/(sum(pi_Oc_nonbusco$V7)/sum(pi_Oc_nonbusco$V8))




pi_Ht_busco$Ht_ratio<-(pi_Ht_busco$V5/pi_Ht_busco$V6)/(pi_Ht_busco$V7/pi_Ht_busco$V8)
pi_Ht_nonbusco$Ht_ratio<-(pi_Ht_nonbusco$V5/pi_Ht_nonbusco$V6)/(pi_Ht_nonbusco$V7/pi_Ht_nonbusco$V8)

#pi_Ht_busco <- pi_Ht_busco[!(pi_Ht_busco$Ht_ratio==Inf | is.na(pi_Ht_busco$Ht_ratio)),]
#pi_Ht_nonbusco <- pi_Ht_nonbusco[!(pi_Ht_nonbusco$Ht_ratio==Inf | is.na(pi_Ht_nonbusco$Ht_ratio)),]
pi_Ht_busco[pi_Ht_busco$Ht_ratio==Inf |is.na(pi_Ht_busco$Ht_ratio),]$Ht_ratio <- 0
pi_Ht_nonbusco[pi_Ht_nonbusco$Ht_ratio==Inf |is.na(pi_Ht_nonbusco$Ht_ratio),]$Ht_ratio <- 0

pi_Oc_busco$Oc_ratio<-(pi_Oc_busco$V5/pi_Oc_busco$V6)/(pi_Oc_busco$V7/pi_Oc_busco$V8)
pi_Oc_nonbusco$Oc_ratio<-(pi_Oc_nonbusco$V5/pi_Oc_nonbusco$V6)/(pi_Oc_nonbusco$V7/pi_Oc_nonbusco$V8)

#pi_Oc_busco <- pi_Oc_busco[!(pi_Oc_busco$Oc_ratio==Inf | is.na(pi_Oc_busco$Oc_ratio)),]
#pi_Oc_nonbusco <- pi_Oc_nonbusco[!(pi_Oc_nonbusco$Oc_ratio==Inf | is.na(pi_Oc_nonbusco$Oc_ratio)),]
pi_Oc_busco[pi_Oc_busco$Oc_ratio==Inf |is.na(pi_Oc_busco$Oc_ratio),]$Oc_ratio <- 0
pi_Oc_nonbusco[pi_Oc_nonbusco$Oc_ratio==Inf |is.na(pi_Oc_nonbusco$Oc_ratio),]$Oc_ratio <- 0



median(pi_Ht_busco$Ht_ratio, na.rm = T)
median(pi_Ht_nonbusco$Ht_ratio, na.rm = T)
median(pi_Oc_busco$Oc_ratio, na.rm = T)
median(pi_Oc_nonbusco$Oc_ratio, na.rm = T)

kruskal.frame<-as.data.frame(rbind(cbind(pi_Ht_busco$Ht_ratio, "HtB"),
      cbind(pi_Ht_nonbusco$Ht_ratio, "HtnB"),
      cbind(pi_Oc_busco$Oc_ratio, "OcB"),
      cbind(pi_Oc_nonbusco$Oc_ratio, "OcnB")
))
colnames(kruskal.frame) <- c("ratio", "group")
kruskal.test(kruskal.frame$ratio~kruskal.frame$group)


(sum(pi_Ht_busco$V5)/sum(pi_Ht_busco$V6))/(sum(pi_Ht_busco$V7)/sum(pi_Ht_busco$V8))
(sum(pi_Ht_nonbusco$V5)/sum(pi_Ht_nonbusco$V6))/(sum(pi_Ht_nonbusco$V7)/sum(pi_Ht_nonbusco$V8))

(sum(pi_Oc_busco$V5)/sum(pi_Oc_busco$V6))/(sum(pi_Oc_busco$V7)/sum(pi_Oc_busco$V8))
(sum(pi_Oc_nonbusco$V5)/sum(pi_Oc_nonbusco$V6))/(sum(pi_Oc_nonbusco$V7)/sum(pi_Oc_nonbusco$V8))



