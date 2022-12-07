####2022-09-02####
####2022-10-23 Zhiyuan Cai####

####download SRAdata####
prefetch SRRnum

####SRA-to-fastq####
fastq-dump --split-3 --readids SRRnum 


####kneaddata####
kneaddata -t 16 --input SRRnum_1.fastq --input SRRnum_2.fastq --output kneaddata_demo_output --trimmomatic trimmomatic-0.39-2/ --trimmomatic-options 'ILLUMINACLIP:trimmomatic-0.39-2/adapters/TruSeq3-PE.fa:2:40:15 SLIDINGWINDOW:4:20 MINLEN:50' --reference-db database/kneaddata/ --bowtie2-options '--very-sensitive --dovetail' --remove-intermediate-output

####metaphlan2####
metaphlan2.py --input_type fastq <(cat SRRnum_kneaddata_paired_?.fastq) >  SRRnum.tsv 

####megahit####
megahit -1 SRRnum_1.fastq -2 kSRRnum_2.fastq -o megahitSRRnum --continue

####PlasFlow####
PlasFlow.py --input SRRnum.fasta --output SRRnum.tsv --threshold 0.7 


####filter contig(>1kb)####
cat plasmidsSRRnum.fasta | awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } '| paste - - |awk 'length($7) >=1000 {print $1 "\n" $7}' > plasmid_1kbSRRnum.fasta ;done &

####PRODIGAL v2.6.3 [February, 2016]####
prodigal  -p meta -a protein_seqSRRnum.fasta -m -d nucleotide_seqSRRnum.fasta  -o genesSRRnum.gff -f gff -s poteintialSRRnum.stat -i plasmid_1kbSRRnum.fasta 


####cd-hit-est####
cd-hit-est -i nucleotide_seqSRRnum.fasta -o plasmid_nrgeneSRRnum.fasta -aS 0.9 -c 0.95 -G 0 -M 0 -T 9 -g 1 


####salmon-1.5.2####
####index####
salmon  index -t plasmid_1kbSRRnum.fasta  -p 32 -k 31 -i plasmid_indexSRRnum 
salmon  index -t nucleotide_seqSRRnum.fasta  -p 32 -k 31 -i plasmid_gene_indexSRRnum 

####quant####
salmon quant -i indexSRRnum -l A -1 SRRnum_kneaddata_paired_1.fastq -2 SRRnum_kneaddata_paired_2.fastq -o plasmid_salmonSRRnum --meta 
salmon quant -i indexSRRnum -l A -1 SRRnum_kneaddata_paired_1.fastq -2 SRRnum_kneaddata_paired_2.fastq -o plasmid_gene_salmonSRRnum --meta


####blast####
makeblastdb -in plasmiddb/plasmids.fna  -dbtype nucl -parse_seqids  -out plasmiddb

blastn -num_threads 16 -query plasmid_1kbSRRnum.fasta  -db plasmiddb -out blast_plasmid_1kbSRRnum.tsv -outfmt '6 qseqid sseqid pident length qcovs qcovhsp qcovus mismatch gapopen qstart qend sstart send evalue bitscore'


####eggNOG-annotation####
emapper.py  -i protein_seqSRRnum.fasta -o NOGSRRnum --cpu 8 --no_file_comments --override

####CAZY-annotation####
run_dbcan.py protein_seqSRRnum.fasta  protein --db_dir CAZyDB/ --out_dir cazyout_SRRnum


####MMUPHin####
fit_adjust_batch <- adjust_batch(feature_abd = crc_ddcount,
                                 batch = "batch",
                                 covariates = "group",
                                 data = csif)
CRC_abd_adj <- fit_adjust_batch$feature_abd_adj

fit_lm_meta <- lm_meta(feature_abd = CRC_abd_adj,
                       batch = "batch",
                       exposure = "group",
                       data = csif,
                       control = list(verbose = FALSE))

####Boruta####
seed(2022)
boruta.train <- Boruta(as.factor(group) ~ .,data = CRC_ddcount_adj, doTrace = 2)

####randomForest####
prei <- randomForest(group~.,data= fold_train,ntree=500,proximity=TRUE,importance=TRUE)



