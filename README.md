# Master-thesis

In this master thesis, the variant calling step was adapted in order to detect and collect INDEL information. The renewed script is called 'snp_indel_calling' . This script makes use of three other scripts that were adapted: 

- Filtering_indels collects all the INDELs called by SAMtools

- filterSAMfile_INDEL collects  INDELs from dbSNP that are present in the analyzed SAMfile

- snpIndexBuilder_new generates  an index for every INDEL and SNP

In addition, the old assembly script is replaced with a new version called 'Assembly_new'. This script enables the insertion of SNPs and INDELs in the  transcript sequences. These sequences are also translated into protein sequences and stored in an SQLite table. 
