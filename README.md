# Master-thesis

In this master thesis, the PROTEOFORMER pipeline was adapted in order to detect, collect  and include INDEL information in the generated custom protein databases. Therefore, the variant calling step and the assembly step were changed. 


The renewed variant calling script is called 'snp_indel_calling'. This script makes use of three other scripts that were adapted: 

- Filtering_indels collects all the INDELs called by SAMtools

- filterSAMfile_INDEL collects INDELs from dbSNP that are present in the analyzed SAMfile

- snpIndexBuilder_new generates an index for every INDEL and SNP

In addition, the old assembly script is replaced with a new version called 'Assembly_new'. This script enables the insertion of SNPs and INDELs in the  transcript sequences. These sequences are also translated into protein sequences and stored in an SQLite table. 


More information of the PROTEOFORMER pipeline is available at http://www.biobix.be/proteoformer/
