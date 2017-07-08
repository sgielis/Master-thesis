# Master-thesis

In this master thesis, the PROTEOFORMER pipeline was adapted in order to detect, collect  and include INDEL information in the generated custom protein databases. Therefore, the variant calling step and the assembly step were changed. 


The renewed variant calling script is called 'snp_indel_calling'. This script makes use of three other scripts that were adapted: 

- Filtering_indels collects all the INDELs called by SAMtools. This script needs to be in the same directory as the snp_indel_calling script. 

- filterSAMfile_INDEL collects INDELs from dbSNP that are present in the analyzed SAMfile. This file replaces the old filterSAMfile script and needs to be in the folder mentioned in the --tooldir argument. Also the original filterSAMfile needs to be in this folder. 

- snpIndexBuilder_new generates an index for every INDEL and SNP. This file replaces the old snpINDEXBuilder script and needs to be in the folder mentioned in the --tooldir argument.

- The original splitVCFaltRecords.pl also needs to be in the folder mentioned in the -- tooldir argument.


In addition, the old assembly script is replaced with a new version called 'Assembly_new'. This script enables the insertion of SNPs and INDELs in the  transcript sequences. These sequences are also translated into protein sequences and stored in an SQLite table. (It is possible that the name of this SQLite db needs to be adapted after running the script in order to be in agreement with consecutive scripts.)


More information of the PROTEOFORMER pipeline is available at http://www.biobix.be/proteoformer/


# Adaptations:
The new Assembly script is adapted in order to include INDELs from COSMIC. This script is called Assembly_cosmic_indels.py.
