#!/usr/bin/python
__author__ = "Sofie Gielis"

'''
Arguments:
    -w   | --workdir						 Path to the working directory. (optional argument, default = current directory)
    -sq  | --sqlitedb						 The SQLite database holding all RIBO-pipeline results. (mandatory argument)
    -in  | --indel 					         (optional argument, default = 'NO', others can be 'samtools', 'samtools_dbSNP')
    -snp | --snp     						 The SNP calling algorithm applied. (optional argument, default 'NO', others can be 'samtools', 'samtools_dbSNP')
    -t   | --tmp_folder						 Folder where temporary files are stored. (optional argument, default = workdir/tmp)
    -i   | --tis_id 						 List of the TIS IDs you want to analyze. (mandatory argument)
    -o   | --sqlite_out                      The SQLite database holding the output. (optional argument, default = sqlitedb = same as input SQLite DB)
'''


import sqlite3
import os
import re
import sys
import json, ast
import timeit
import math
import multiprocessing
import io

import shutil
import collections
from collections import OrderedDict
from string import maketrans
import csv
import traceback
import getopt
from multiprocessing import Pool
def main():

    print(" ")
    print("-----------------------------------")
    print(" ASSEMBLY OF TRANSLATION PRODUCTS")
    print("-----------------------------------")
    print("	")

    # Start time:
    startTime = timeit.default_timer()

    # Catch command line with getopt:
    try:
        myopts, args = getopt.getopt(sys.argv[1:],"w:sq:snp:in:t:i:o",["workdir=","sqlitedb=","snp=","indel=",'tmp_folder=','tis_id=','sqlite_out='])
    except getopt.GetoptError as err:
        print err
        sys.exit()

    # Catch arguments:
    # o == option in command line
    # a == argument passed to o
    for o, a in myopts:
        if o in ('-w','--workdir'):
            workdir=a
        elif o in ('-sq','--sqlitedb'):
            sqlitedb=a
        elif o in ('-snp','--snp'):
            snp=a
        elif o in ('-in','--indel'):
            indel=a
        elif o in ('-t','--tmp_folder'):
            tmp_folder=a
        elif o in ('-i','--tis_id'):
            tis_id=a
        elif o in ('-o','--sqlite_out'):
            sqlite_out=a

    try:
        workdir
    except:
        workdir=''
    try:
        sqlitedb
    except:
        sqlitedb=''
    try:
        snp
    except:
        snp=''
    try:
        indel
    except:
        indel=''
    try:
        tmp_folder
    except:
        tmp_folder=''
    try:
        tis_id
    except:
        tis_id=''
    try:
        sqlite_out
    except:
        sqlite_out=''

    # Check for correct arguments
    if(workdir == ''):
        workdir = os.getcwd()
    if(workdir !=''):
        os.chdir(workdir)
    if(sqlitedb==''):
        print "Error: do not forget to fill in the SQLite DB argument!"
        sys.exit()
    if(tis_id==''):
        print "Error: do not forget to fill in the TIS IDs argument!"
        sys.exit()
    tis_id_list = []
    if(tis_id=='all'):
        tis_id_list = get_all_tis_ids(sqlitedb)
    elif(tis_id.isdigit()):
        tis_id_list.append(tis_id)
    elif(bool(re.search(',',tis_id))):
        tis_id_list = re.split(',', tis_id)
    else:
        print "Error: TIS IDs argument should be a number, a comma separated enumeration of numbers or 'all'!"
        sys.exit()
    if(indel==''):
        indel='NO'
    if (snp==''):
        snp='NO'
    if (tmp_folder==''):
        if (not os.path.isdir(workdir + "/tmp")):
            os.system("mkdir " + workdir + "/tmp")
        tmp_folder=workdir+'/tmp'
    if (sqlite_out==''):
        sqlite_out=sqlitedb

    # Print used arguments and parameters:
    print("Used arguments and parameters:")
    print(" Working directory:                         " + workdir)
    print(" SQLite DB:                                 " + sqlitedb)
    print(" TMP folder:                                " + tmp_folder)
    print(" SNP					                       " + snp)
    print(" INDELs:    		                           " + indel)
    print(" TIS IDs that need to be analyzed:          " + str(tis_id_list))
    print(" SQLite out DB:							   " + sqlite_out)
    print("")


    # Get the input variables:
    ensembldb, igenomes_root, species, ens_version, nr_of_cores = get_arguments(sqlitedb)
    print('Getting the arguments from the arguments table:')
    print(' ensembldb:					' + ensembldb)
    print(' igenomes_root:				' + igenomes_root)
    print(' species:					' + species)
    print(' ens_version:				' + str(ens_version))
    print(' nr_of_cores:				' + str(nr_of_cores))
    print(' ')

    # Conversion of species terminology:
    speciesLatin = "Mus_musculus" if species == "mouse" else \
        "Homo_sapiens" if species == "human" else \
        "Arabidopsis_thaliana" if species == "arabidopsis" else \
        "Drosophila_melanogaster" if species == "fruitfly" else ""
    speciesShort = "mmu" if species == "mouse" else \
        "hsa" if species == "human" else \
        "ath" if species == "arabidopsis" else \
        "dme" if species == "fruitfly" else ""
    # Determine the assembly:
    assembly = "GRCm38" if species == "mouse" and ens_version >= 70 else \
        "NCBIM37" if species == "mouse" and ens_version < 70 else \
        "GRCh38" if species == "human" and ens_version > 75 else \
        "GRCh37" if species == "human" and ens_version <= 75 else \
        "TAIR10" if species == "arabidopsis" else \
        "BDGP5" if species == "fruitfly" and ens_version < 79 else \
        "BDGP6" if species == "fruitfly" and ens_version >= 79 else ""

    ###########################################
    ## Assembly of the translation products: ##
    ###########################################

    # Get the chromosomes with their sizes:
    print('Get the chromosomes with their sizes.')
    chrs={}
    chromosomeSizesFile = igenomes_root+"/"+speciesLatin+"/Ensembl/"+assembly+"/Annotation/Genes/ChromInfo.txt"
    if os.path.isfile(chromosomeSizesFile):
        chrs = get_chrs_sizes(chromosomeSizesFile, species)
    else:
        chromosomeSizesFile = igenomes_root+"/"+speciesLatin+"/Ensembl/"+assembly+"/Sequence/WholeGenomeFasta/GenomeSize.xml"
        if os.path.isfile(chromosomeSizesFile):
            chrs = get_chrsXml_sizes(chromosomeSizesFile, species)
        else:
            print "ERROR: chromosome sizes file could not be found."
            sys.exit()
    print(chrs)
    print(" ")

    # Create binary chromosomes if they don't exist:
    print("Creating binary chromosomes if they don't exist.")
    BIN_chrom_dir=workdir+'/chromosomes_BIN'
    if os.path.isdir(BIN_chrom_dir):
        print('File already exists')
    else:
        BIN_chromosomes=create_BIN_chromosomes(nr_of_cores,BIN_chrom_dir,chrs,igenomes_root,speciesLatin,assembly)
    print(" ")

    # Loop over all selected TIS IDs:
    for id in tis_id_list:
        analysis_ID=id
        print("Processing TIS ID "+analysis_ID+".")
        # Construct the translation products for the selected TIS ID.
        construct_translation_product(species,assembly,chrs,ensembldb,workdir,sqlitedb,snp,indel,analysis_ID,nr_of_cores)
        # Add SNP and INDEL info to the TIS_overview table.
        print(" Update TIS_overview.")
        update_TIS_overview(analysis_ID,snp,sqlitedb,indel)
    print(" ")

    # Calculate the run time:
    stopTime = timeit.default_timer()
    runtime = math.ceil(stopTime - startTime)
    m, s = divmod(runtime, 60)
    h, m = divmod(m, 60)
    runtime = "%d:%02d:%02d" % (h, m, s)
    print "\n   ------ Program COMPLETED ------"
    print "   The program run time: " + runtime + "\n\n"

    return


####################
## The functions: ##
####################

### get_arguments ###
def get_arguments(path_to_sqlitedb):
    '''
    :param path_to_sqlitedb: the location of the SQLite database containing the table arguments
    :return: the name of the Ensembl database, the location of the iGenomes root folder, the species that is studied, the Ensembl version and the number of cores used
    '''
    try:
        con=sqlite3.connect(path_to_sqlitedb)
    except:
        print "Could not connect to the SQLite database."
        sys.exit()

    #Init
    ens_db=''
    igenomes_root=''
    species=''
    ensembl_version=0
    nr_of_cores=0

    with con:
        cur = con.cursor()

        if cur.execute("SELECT value FROM arguments WHERE variable='ens_db';"):
            ens_db = str(cur.fetchone()[0])
        else:
            print "ERROR: could not fetch the Ensembl DB from the arguments table in the SQLite DB."
            sys.exit()

        if cur.execute("SELECT value FROM arguments WHERE variable='igenomes_root';"):
            igenomes_root = str(cur.fetchone()[0])
        else:
            print "ERROR: could not fetch the iGenomes root from the arguments table in the SQLite DB."
            sys.exit()

        if cur.execute("SELECT value FROM arguments WHERE variable='species';"):
            species = str(cur.fetchone()[0])
        else:
            print "ERROR: could not fetch the species from the arguments table in the SQLite DB."
            sys.exit()

        if cur.execute("SELECT value FROM arguments WHERE variable='ensembl_version';"):
            ensembl_version = int(cur.fetchone()[0])
        else:
            print "ERROR: could not fetch the Ensembl version from the arguments table in the SQLite DB."
            sys.exit()

        if cur.execute("SELECT value FROM arguments WHERE variable='nr_of_cores';"):
            nr_of_cores = int(cur.fetchone()[0])
        else:
            print "ERROR: could not fetch the number of cores from the arguments table in the SQLite DB."
            sys.exit()

    return ens_db, igenomes_root, species, ensembl_version, nr_of_cores

### get_all_tis_ids ###
def get_all_tis_ids(path_to_sqlitedb):
    '''
    :param path_to_sqlitedb: the location of the SQLite database containing the table TIS_overview
    :return: a list of all the TIS IDs from table TIS_overview
    '''
    try:
        con=sqlite3.connect(path_to_sqlitedb)
    except:
        print "Error: could not connect to the SQLite database."
        sys.exit()
    with con:
        cur=con.cursor()
        if cur.execute("SELECT ID FROM TIS_overview"):
            tis_ids_list = [ item[0] for item in cur.fetchall()]
    return tis_ids_list

### get_chrs_sizes ###
def get_chrs_sizes(chrSizeFile, species):
    '''
    :param chrSizeFile: file which contains the sizes of every chromosome
    :param species: the species that is studied
    :return: dictionary containing the chromosomes as keys and the length of the chromosomes as values
    '''

    #Init
    chrs = {}

    #open file
    try:
        FR = open(chrSizeFile, 'r')
    except:
        print "ERROR with opening chromosome sizes file."

    #parse
    for line in FR:
        parts = re.split('\W+',line)
        if species=="fruitfly" and parts[0]=="M":
            chrs["dmel_mitochondrion_genome"] = parts[1]
        else:
            chrs[parts[0]] = parts[1]

    return chrs

### get_chrsXml_sizes ###
def get_chrsXml_sizes(chrSizeFile, species):
    '''
    :param chrSizeFile: file which contains the sizes of every chromosome
    :param species: the species that is studied
    :return: dictionary containing the chromosomes as keys and the length of the chromosomes as values
    '''
    #Init
    chrs = {}

    #open file
    try:
        FR = open(chrSizeFile, 'r')
    except:
        print "ERROR with opening chromosome sizes xml file."

    #parse
    for line in FR:
        pattern = re.compile('contigName=\"(\w{1,2})\" totalBases\W\"(\d+)\"')
        m = pattern.search(line)
        if m:
            if species=="fruitfly" and m.group(1)=="M":
                chrs["dmel_mitochondrion_genome"] = m.group(2)
            else:
                chrs[m.group(1)] = m.group(2)

    return chrs

### create_BIN_chromosomes ###
def create_BIN_chromosomes(nr_of_cores,BIN_chrom_dir,chrs,igenomes_root,speciesLatin,assembly):
    os.mkdir(BIN_chrom_dir)
    #from multiprocessing import Pool
    pool = multiprocessing.Pool(processes=nr_of_cores)
    for key in chrs:
       CHR = igenomes_root + "/" + speciesLatin + "/Ensembl/" + assembly + "/Sequence/Chromosomes/" + key + ".fa"
       file = open(CHR, "r")
       with io.FileIO("CHR_BIN" + key, "w") as newfile:
          for line in file:
             if line.startswith('>'):
                continue
             else:
                line = line.replace('\n', '')
                newfile.write(line)
       shutil.move('CHR_BIN' + key, BIN_chrom_dir)
    pool.close()
    pool.join()

### get_chrs_sequence ###
def get_chrs_sequence(path_to_ensembldb, assembly, chrs, species):
    '''
    :param path_to_ensembldb: the location of the Ensembl database
    :param assembly: the assembly that is used
    :param chrs: a dictionary with the chromosomes as keys and the corresponding length as values
    :param species: the species that is studied
    :return: a dictionary with the chromosomes as keys and the corresponding seq_region_IDs as values
    '''

    #Init:
    chrs_seqRegionID= {}

    try:
        con = sqlite3.connect(path_to_ensembldb)
    except:
        print "Error: could not connect to Ensembl database."
        sys.exit()

    with con:
        cur = con.cursor()
        if cur.execute("SELECT coord_system_id FROM coord_system WHERE name = 'chromosome' AND version = '"+assembly+"';"):
            coordSystemID = cur.fetchone()[0]
        else:
            print "ERROR: could not fetch the coord_system_id from the Ensembl database."
            sys.exit()

    # Get chrs with their seq_region_id:
    for key in chrs:
        seqRegionID=""
        if (species =="fruitfly"):
            if (key =="M"):
                chrom = "dmel_mitochondrion_genome"
            else:
                chrom = key
        else:
            chrom = key
        cur = con.cursor()
        if cur.execute("SELECT seq_region_id FROM seq_region where coord_system_id = '"+str(coordSystemID)+"' and name = '"+chrom+"';"):
            seqRegionID =[ item[0] for item in cur.fetchall()]
        else:
            print "ERROR: could not fetch the seq_region_id from the Ensembl database."
            sys.exit()
        chrs_seqRegionID[chrom] = seqRegionID

    return chrs_seqRegionID

### get_transcripts_per_chromosome ###
def get_transcripts_per_chromosome(path_to_ensembldb,chromosome,analysis_id):
    '''
    :param path_to_ensembldb: the location of the Ensembl database
    :param chromosome: the number of the chromosome that needs to be studied
    :param analysis_id: a number representing the TIS ID that is being analyzed
    :return: a dictionary with the transcript_start as key and the info of the transcript as value
    '''

    transcript_starts = {}
    try:
        con = sqlite3.connect(path_to_ensembldb)
    except:
        print "Error: could not connect to Ensembl database."
        sys.exit()

    with con:
        cur = con.cursor()
        if cur.execute("SELECT transcript_id||'_'||start as transcript_start,transcript_id,biotype,chr,strand,start,dist_to_transcript_start,dist_to_aTIS,annotation,aTIS_call,start_codon,peak_shift,count,Rltm_min_Rchx FROM TIS_"+str(analysis_id)+" WHERE chr = '"+str(chromosome)+"';"):
            transcript_starts1 = cur.fetchall()
            for row in transcript_starts1:
                transcript_start=""
                #print(row)
                transcript_start = row[0]
                #print (transcript_start)
                transcript_starts[transcript_start]=row

        else:
            print "ERROR: could not fetch the information from the Ensembl database."
            sys.exit()

    # The dictionary transcrip_starts contains the desired values, but these values are not named.
    # Therefore a new dictionary is made containing a name for every value.
    transcript_per_chromosome = {}
    for key in transcript_starts:
        transcript_info = {}
        if transcript_starts.has_key(key):
            # print transcript_starts[key]
            transcript_start = transcript_starts[key][0]
            transcript_id = transcript_starts[key][1]
            biotype = transcript_starts[key][2]
            chr = transcript_starts[key][3]
            strand = transcript_starts[key][4]
            start = transcript_starts[key][5]
            dist_to_transcript_start = transcript_starts[key][6]
            dist_to_aTIS = transcript_starts[key][7]
            annotation = transcript_starts[key][8]
            aTIS_call = transcript_starts[key][9]
            start_codon = transcript_starts[key][10]
            peak_shift = transcript_starts[key][11]
            count = transcript_starts[key][12]
            Rltm_min_Rchx = transcript_starts[key][13]
            transcript_info['transcript_start'] = transcript_start
            transcript_info['transcript_id'] = transcript_id
            transcript_info['biotype'] = biotype
            transcript_info['chr'] = chr
            transcript_info['strand'] = strand
            transcript_info['start'] = start
            transcript_info['dist_to_transcript_start'] = dist_to_transcript_start
            transcript_info['dist_to_aTIS'] = dist_to_aTIS
            transcript_info['annotation'] = annotation
            transcript_info['aTIS_call'] = aTIS_call
            transcript_info['start_codon'] = start_codon
            transcript_info['peak_shift'] = peak_shift
            transcript_info['count'] = count
            transcript_info['Rltm_min_Rchx'] = Rltm_min_Rchx
        transcript_per_chromosome[transcript_start] = transcript_info
    transcript_per_chromosome = ast.literal_eval(json.dumps(transcript_per_chromosome))
    return transcript_per_chromosome

### get_reverse_complement ###
def get_reverse_complement(sequence):
    '''
    :param sequence: a mRNA sequence
    :return: the reverse complement of the mRNA sequence
    '''
    reverse=sequence[::-1]
    complement=maketrans('AGTCagtc','TCAGtcag')
    reverse_complement=reverse.translate(complement)
    return reverse_complement

### get_SNPS_per_chromosome ###
def get_SNPs_per_chromosome(path_to_sqlitedb,chromosome,snp):
    '''
    :param path_to_sqlitedb: the location of the SQLite database holding the table with SNP information
    :param chromosome: the number of the chromosome that needs to be studied
    :param snp: the SNP algorithm that is used: 'NO', 'samtools' or 'samtools_dbSNP'
    :return: a dictionary containing the following info for every SNP on this chromosome: position, reference SNP, alternative SNP, AF and the chromosome number itself
    '''

    all_SNPs_for_this_chromosome = {}

    # Catch all the SNPs for this chromosome:
    try:
        con = sqlite3.connect(path_to_sqlitedb)
    except:
        print "Error: could not connect to SQLite database."
        sys.exit()

    with con:
        cur = con.cursor()
        if snp=='samtools':
            if cur.execute("Select id,chr,pos,ref,alt,case when new='m' then 0.5 else af end as af from snp_"+snp+"  where chr = '" + str(chromosome) + "' and new<>'m';"):
                snp_info = cur.fetchall()
                for row in snp_info:
                    id_snp = ""
                    #print(row)
                    id_snp = row[0]
                    #print (id_snp)
                    all_SNPs_for_this_chromosome[id_snp] = ast.literal_eval(json.dumps(row))
                    #print(all_SNPs_for_this_chromosome)
            else:
                print "ERROR: could not fetch the information from the SQLite database."
                sys.exit()
        elif snp=="samtools_dbSNP":
            name = re.split('_',snp)[0]
            if cur.execute("Select id,chr,pos,ref,alt,case when new='m' then 0.5 else af end as af from snp_"+name+"  where chr = '" + str(chromosome) + "';"):
                snp_info = cur.fetchall()
                for row in snp_info:
                    id_snp = ""
                    #print(row)
                    id_snp = row[0]
                    #print (id_snp)
                    all_SNPs_for_this_chromosome[id_snp] = ast.literal_eval(json.dumps(row))
            else:
                print "ERROR: could not fetch the information from the SQLite database."
                sys.exit()

    # Put all the information in a dictionary called all_snp_info_for_this_chromosome.
    all_snp_info_for_this_chromosome={}
    position_list = []
    for key in all_SNPs_for_this_chromosome:
        #print(key)
        info_snp={}
        if all_SNPs_for_this_chromosome.has_key(key):
            #print(all_SNPs_for_this_chromosome[key]) #prints value for every key
            pos = all_SNPs_for_this_chromosome[key][2]
            #print(pos)
            position_list.append(pos)
            #print(position_list.count(pos))
            if position_list.count(pos)==2:
                #print(pos)
                alt1 = all_SNPs_for_this_chromosome[key][4]
                alt2=all_snp_info_for_this_chromosome[pos]['alt'][0]
                #print(alt1,alt2)
                all_snp_info_for_this_chromosome[pos]['alt']=[alt1,alt2]
                #print(all_snp_info_for_this_chromosome)

            elif position_list.count(pos)==3:
                #print(pos)
                alt1 = all_SNPs_for_this_chromosome[key][4]
                alt2 = all_snp_info_for_this_chromosome[pos]['alt'][0]
                alt3=all_snp_info_for_this_chromosome[pos]['alt'][1]
                #print(alt1, alt2,alt3)
                all_snp_info_for_this_chromosome[pos]['alt'] = [alt1, alt2, alt3]

            else:
                chr = all_SNPs_for_this_chromosome[key][1]
                pos = all_SNPs_for_this_chromosome[key][2]
                ref = all_SNPs_for_this_chromosome[key][3]
                alt = all_SNPs_for_this_chromosome[key][4]
                af = all_SNPs_for_this_chromosome[key][5]
                info_snp['chr'] = chr
                info_snp['pos'] = pos
                info_snp['ref'] = ref
                info_snp['alt'] = alt
                info_snp['af'] = af
                all_snp_info_for_this_chromosome[pos]=info_snp

    return all_snp_info_for_this_chromosome


### get_INDELS_per_chromosome
def get_INDELS_per_chromosome(path_to_sqlitedb,chromosome,indel):
    '''
    :param path_to_sqlitedb: path to the SQLite database holding the table with INDEL information
    :param chromosome: the number of the chromosome that needs to be studied
    :param indel: the algorithm that is used: 'NO', 'samtools' or 'samtools_dbSNP'
    :return: A dictionary with the position of the INDEL as key and the following information as value: chromosome number, position, reference INDEL, alternative INDEL, allelic frequency
    '''

    all_INDELS_for_this_chromosome = {}

    # Catch all the INDELs for this chromosome:
    try:
        con = sqlite3.connect(path_to_sqlitedb)
    except:
        print "Error: could not connect to SQLite database."
        sys.exit()

    with con:
        cur = con.cursor()
        if indel=='samtools':
            if cur.execute("Select id,chr,pos,ref,alt,case when new='m' then 0.5 else af end as af from indel_samtools  where chr = '" + str(chromosome) + "' and new<>'m';"):
                INDEL_info = cur.fetchall()
                for row in INDEL_info:
                    id_INDEL = ""
                    #print(row) #This gives you the info for every position
                    id_INDEL = row[0]
                    #print (id_INDEL) #This gives you id
                    all_INDELS_for_this_chromosome[id_INDEL] = ast.literal_eval(json.dumps(row))
            else:
                print "ERROR: could not fetch the information from the SQLite database."
                sys.exit()
        elif indel=='samtools_dbSNP':
            name = re.split('_', indel)[0]
            if cur.execute("Select id,chr,pos,ref,alt,case when new='m' then 0.5 else af end as af from indel_samtools where chr = '" + str(chromosome) + "';"):
                INDEL_info = cur.fetchall()
                for row in INDEL_info:
                    id_INDEL = ""
                    # print(row)#This gives you the info for every position
                    id_INDEL = row[0]
                    # print (id_INDEL) #This gives you the id
                    all_INDELS_for_this_chromosome[id_INDEL]=ast.literal_eval(json.dumps(row))
            else:
                print "ERROR: could not fetch the information from the SQLite database."
                sys.exit()

    # Put all the information in a dictionary called all_INDEL_info_for_this_chromosome.
    all_INDEL_info_for_this_chromosome = {}
    position_list=[]
    for key in all_INDELS_for_this_chromosome:
        info_indel = {}
        if all_INDELS_for_this_chromosome.has_key(key):
            #print(all_INDELS_for_this_chromosome[key]) #prints value for every key
            pos = all_INDELS_for_this_chromosome[key][2]
            position_list.append(pos)
            #print('aantal',position_list.count(pos))

            if position_list.count(pos) == 2:
                #print(pos)
                alt1 = all_INDELS_for_this_chromosome[key][4]
                alt2 = all_INDEL_info_for_this_chromosome[pos]['alt'][0]
                all_INDEL_info_for_this_chromosome[pos]['alt'] = [alt1, alt2]

            elif position_list.count(pos) == 3:
                #print(pos)
                alt1 = all_INDELS_for_this_chromosome[key][4]
                alt2 = all_INDEL_info_for_this_chromosome[pos]['alt'][0]
                alt3 = all_INDEL_info_for_this_chromosome[pos]['alt'][1]
                all_INDEL_info_for_this_chromosome[pos]['alt'] = [alt1, alt2, alt3]

            else:
                #print(pos)
                chr = all_INDELS_for_this_chromosome[key][1]
                pos = all_INDELS_for_this_chromosome[key][2]
                ref = all_INDELS_for_this_chromosome[key][3]
                alt = all_INDELS_for_this_chromosome[key][4]
                af = all_INDELS_for_this_chromosome[key][5]
                info_indel['chr'] = chr
                info_indel['pos'] = pos
                info_indel['ref'] = ref
                info_indel['alt'] = [alt]
                info_indel['af'] = af
                all_INDEL_info_for_this_chromosome[pos] = info_indel
                #print(all_INDEL_info_for_this_chromosome)

    return all_INDEL_info_for_this_chromosome

### get_exons_for_transcript ###
def get_exons_for_transcript(path_to_ensembldb,transcript_id):
    '''
    :param path_to_ensembldb: path to the Ensembl database
    :param transcript_id: the Ensembl transcript ID
    :return: a dictionary with rank as key and another dictionary as value containing info about the exons
    '''
    exons_for_this_transcript = {}

    # Catch the exons for this transcript:
    try:
        con = sqlite3.connect(path_to_ensembldb)
    except:
        print "Error: could not connect to Ensembl database."
        sys.exit()

    with con:
        cur = con.cursor()
        if cur.execute("select transcript.transcript_id, transcript.stable_id, exon_transcript.rank, seq_region.name chr,exon.seq_region_start,exon.seq_region_end,exon.seq_region_strand,exon.phase,exon.end_phase,exon.stable_id exon_stable_id from transcript inner join exon_transcript on transcript.transcript_id = exon_transcript.transcript_id inner join exon on exon.exon_id = exon_transcript.exon_id inner join seq_region on seq_region.seq_region_id = transcript.seq_region_id where transcript.transcript_id = '" + str(transcript_id) + "';"):
            exons = cur.fetchall()
            for row in exons:
                rank=""
                #print(row)#Gives you all the info for the exon
                rank = row[2]
                #print (rank)
                exons_for_this_transcript[rank]=row

        else:
            print "ERROR: could not fetch the information from the Ensembl database."
            sys.exit()
    exons_for_this_transcript = ast.literal_eval(json.dumps(exons_for_this_transcript))

    # Put all the information in a dictionary all_exons_for_this_transcript.
    all_exons_for_this_transcript={}
    for key in exons_for_this_transcript:
        exon_info ={}
        if exons_for_this_transcript.has_key(key):
            transcript_id = exons_for_this_transcript[key][0]
            tr_stable_id = exons_for_this_transcript[key][1]
            rank = exons_for_this_transcript[key][2]
            chr = exons_for_this_transcript[key][3]
            seq_region_start = exons_for_this_transcript[key][4]
            seq_region_end = exons_for_this_transcript[key][5]
            seq_region_strand = exons_for_this_transcript[key][6]
            phase = exons_for_this_transcript[key][7]
            end_phase = exons_for_this_transcript[key][8]
            exon_stable_id = exons_for_this_transcript[key][9]
            exon_info['transcript_id'] = transcript_id
            exon_info['transcript_stable_id'] = tr_stable_id
            exon_info['rank'] = rank
            exon_info['chr'] = chr
            exon_info['seq_region_start'] = seq_region_start
            exon_info['seq_region_end'] = seq_region_end
            exon_info['seq_region_strand'] = seq_region_strand
            exon_info['phase'] = phase
            exon_info['end_phase'] = end_phase
            exon_info['exon_stable_id'] = exon_stable_id
        all_exons_for_this_transcript[rank] = exon_info

    return all_exons_for_this_transcript

### fetch_SNPs_per_exon ###
def fetch_SNPs_per_exon(exon_start,exon_end,all_SNPs_for_this_chromosome,CDS_tmp_length,strand):
    '''
    :param exon_start: the start position of the exon (= seq_region_start from the Ensembl db)
    :param exon_end: the end position of the exon (= seq_region_end from the Ensembl db)
    :param all_SNPs_for_this_chromosome: a dictionary containing the position of every SNP as key and the info (alt, ref, pos, af, chr) of every SNP as value for the examined chromosome
    :param CDS_tmp_length: the length of the CDS up until now
    :param strand: 1 (= sense) or -1 (=antisense)
    :return: a dictionary containing all the SNPs for this exon
    '''
    all_SNPs_for_this_exon = {}
    for key in all_SNPs_for_this_chromosome:
        fetch_SNPs_for_this_exon = {}
        if (exon_start <= int(key) and int(key) <= exon_end):
            af = all_SNPs_for_this_chromosome[key]['af']
            pos = all_SNPs_for_this_chromosome[key]['pos']
            if strand ==1:
                sequence_position=pos-exon_start+1+CDS_tmp_length
                ref = all_SNPs_for_this_chromosome[key]['ref']
                alt = all_SNPs_for_this_chromosome[key]['alt']
                pos_exon=pos-exon_start+1
            if strand==-1:
                alt=[]
                for bases in all_SNPs_for_this_chromosome[key]['alt']:
                    reversed_base=get_reverse_complement(bases)
                    alt.append(reversed_base)
                ref = get_reverse_complement(all_SNPs_for_this_chromosome[key]['ref'])
                sequence_position = exon_end - (pos + len(ref) - 1) + 1 + CDS_tmp_length
                pos_exon=exon_end - (pos + len(ref) - 1) + 1
            fetch_SNPs_for_this_exon['ref'] = ref
            fetch_SNPs_for_this_exon['alt'] = alt
            fetch_SNPs_for_this_exon['af'] = af
            fetch_SNPs_for_this_exon['pos'] = pos
            fetch_SNPs_for_this_exon['pos_exon']=pos_exon
            fetch_SNPs_for_this_exon['type']='SNP'
            fetch_SNPs_for_this_exon['sequence_position'] = sequence_position
            all_SNPs_for_this_exon[key] = fetch_SNPs_for_this_exon
    return all_SNPs_for_this_exon



### fetch_INDELS_per_exon ###
def fetch_INDELS_per_exon(exon_start,exon_end,all_INDELS_for_this_chromosome,CDS_tmp_length,strand):
    '''
    :param exon_start: the start position of the exon (= seq_region_start from the Ensembl db)
    :param exon_end: the end position of the exon (= seq_region_end from the Ensembl db)
    :param all_INDELS_for_this_chromosome: a dictionary containing the position of every INDEL as key and the info (alt, ref, pos, af, chr) of every INDEL as value for the examined chromosome
    :param CDS_tmp_length: the length of the CDS up until now
    :param strand: 1 (=sense) or -1 (=antisense)
    :return: a dictionary containing all the INDELS for this exon
    '''
    all_INDELS_for_this_exon = {}
    for key in all_INDELS_for_this_chromosome:
        fetch_INDELS_for_this_exon={}
        if (exon_start <= int(key) and int(key) <= exon_end):
            af = all_INDELS_for_this_chromosome[key]['af']
            pos = all_INDELS_for_this_chromosome[key]['pos']
            if strand ==1:
                sequence_position=pos-exon_start+1+CDS_tmp_length
                ref = all_INDELS_for_this_chromosome[key]['ref']
                alt = all_INDELS_for_this_chromosome[key]['alt']
                pos_exon=pos-exon_start+1

            if strand==-1:
                alt=[]
                for bases in all_INDELS_for_this_chromosome[key]['alt']:
                    reversed_base=get_reverse_complement(bases)
                    alt.append(reversed_base)
                ref = get_reverse_complement(all_INDELS_for_this_chromosome[key]['ref'])
                sequence_position = exon_end - (pos + len(ref) - 1) + 1 + CDS_tmp_length
                pos_exon=exon_end - (pos + len(ref) - 1) + 1
            fetch_INDELS_for_this_exon['ref'] = ref
            fetch_INDELS_for_this_exon['alt'] = alt
            fetch_INDELS_for_this_exon['af'] = af
            fetch_INDELS_for_this_exon['pos'] = pos
            fetch_INDELS_for_this_exon['pos_exon']=pos_exon
            fetch_INDELS_for_this_exon['type']='INDEL'
            fetch_INDELS_for_this_exon['sequence_position'] = sequence_position
            all_INDELS_for_this_exon[key] = fetch_INDELS_for_this_exon
    return all_INDELS_for_this_exon

### merge_INDELS_and_SNPS_per_exon ###
def merge_INDELS_and_SNPs_per_exon(all_SNPs_for_this_exon, all_INDELS_for_this_exon):
    '''
    :param all_SNPs_for_this_exon: dictionary containing all SNPs for this exon
    :param all_INDELS_for_this_exon: dictionary containing all INDELS for this exon
    :return: dictionary containing all SNPs and INDELS
    '''
    INDELS_and_SNPs_per_exon = {}
    INDELS_and_SNPs_per_exon.update(all_SNPs_for_this_exon)
    INDELS_and_SNPs_per_exon.update(all_INDELS_for_this_exon)

    return INDELS_and_SNPs_per_exon

### get_sequence ###
def get_sequence(chromosome,exon_start,exon_end):
    '''
    :param chromosome: the number of the chromosome
    :param exon_start: the start position of the exon (= seq_region_start from the Ensembl db)
    :param exon_end: the end position of the exoxn (= seq_region_end from the Ensembl db)
    :return: the exon sequence
    '''
    length = exon_end - exon_start + 1
    offset = exon_start - 1
    with open('chromosomes_BIN/CHR_BIN' + str(chromosome), 'rb') as inputfile:
        inputfile.seek(offset)
        sequence = inputfile.read(length)
        inputfile.close()
    return sequence

### translate_sequence ###
def translate_sequence(mrna_sequence):
    '''
    :param mrna_sequence: the mRNA sequence you want to translate
    :return: the protein sequence = the translated mRNA-sequence
    '''

    translated_sequence = ''
    position_stop_in_sequence=''

    # The codon table represents the corresponding amino acid for every codon.
    CodonTable={'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
    'AAN':'X','ATN':'X','ACN':'T','AGN':'X','ANN':'X','ANA':'X','ANG':'X','ANC':'X',
    'ANT':'X','TAN':'X','TTN':'X','TCN':'S','TGN':'X','TNN':'X','TNA':'X','TNG':'X',
    'TNC':'X','TNT':'X','CAN':'X','CTN':'L','CCN':'P','CGN':'R','CNN':'X','CNA':'X',
    'CNG':'X','CNC':'X','CNT':'X','GAN':'X','GTN':'V','GCN':'A','GGN':'G','GNN':'X',
    'GNA':'X','GNG':'X','GNC':'X','GNT':'X','NAN':'X','NAA':'X','NAG':'X','NAC':'X',
    'NAT':'X','NTN':'X','NTA':'X','NTG':'X','NTC':'X','NTT':'X','NGN':'X','NGA':'X',
    'NGG':'X','NGC':'X','NGT':'X','NCN':'X','NCA':'X','NCG':'X','NCC':'X','NCT':'X',
    'NNN':'X','NNA':'X','NNG':'X','NNC':'X','NNT':'X'}

    # After the translation, we only want to retain the protein sequences comming from mRNA-sequences having a stop codon.
    # Protein sequences comming from mRNA-sequences not containin a stop codon are eliminated. Therefore this function will
    # give an extra parameter together with the protein sequence: the complete_protein_sequence parameter. Complete_protein_sequence='NO'
    # corresponds to sequences lacking the stop codon; complete_protein_sequence ='YES' corresponds to sequences having a stop codon.
    complete_protein_sequence='NO'

    for position_in_sequence in range(0,len(mrna_sequence)-2,3):
        codon=mrna_sequence[position_in_sequence:position_in_sequence+3]

        if re.match(r'TAA|TAG|TGA',codon):
            complete_protein_sequence='YES'
            position_stop_in_sequence = position_in_sequence # This is the start position of the stop codon.
            translated_sequence+=CodonTable[codon]
            break #Stop the translation when a stop codon is encountered.
        else:
            translated_sequence+=CodonTable[codon]

    return(translated_sequence,complete_protein_sequence,position_stop_in_sequence)

### store in SQLite database ###
def store_in_db(chrs, workdir, path_to_sqlitedb, snp, analysis_ID, indel):
    '''
    :param chrs: dictionary containing the chromosomes as keys and the length of the chromosomes as values
    :param workdir: path to the working directory
    :param path_to_sqlitedb: path to the SQLite database where you want to store the information
    :param snp: the SNP algorithm that is used: 'NO', 'samtools' or 'samtools_dbSNP'
    :param analysis_ID: a number representing the TIS ID that is being analysed
    :param indel:
    :return:
    '''
    try:
        # Make DB connection:
        try:
            con = sqlite3.connect(path_to_sqlitedb)
        except:
            print "Could not connect to "+path_to_sqlitedb
            sys.exit()
        with con:
            cur = con.cursor()
            if snp=="NO":
                if indel =='samtools':
                    tableName = "TIS_" + str(analysis_ID) + "_indel" + indel + "_transcripts"
                if indel=='dbSNP'
                    tableName = "TIS_" + str(analysis_ID) + "_indelsamtools_indel" + indel + "_transcripts"

                else:
                    tableName = "TIS_" + str(analysis_ID) + "_" + "transcripts"
            if snp=='samtools':
                if indel=="samtools":
                    tableName = "TIS_"+str(analysis_ID)+"_snpsamtools_indelsamtools_transcripts"
                if indel=="dbSNP"
                    tableName ="TIS_"+str(analysis_ID)+"_snpsamtools_indelsamtools_indeldbSNP_transcripts"
                else:
                    tableName ="TIS_"+str(analysis_ID)+"_snpsamtools_transcripts"

            if snp=='dbSNP':
                if indel=="samtools":
                    tableName = "TIS_"+str(analysis_ID)+"_snpsamtools_snpdbSNP_indelsamtools_transcripts"
                if indel=="dbSNP"
                    tableName ="TIS_"+str(analysis_ID)+"_snpsamtools_snpdbSNP_indelsamtools_indeldbSNP_transcripts"
                else:
                    tableName ="TIS_"+str(analysis_ID)+"_snpsamtools_snpdbSNP_transcripts"


            # Remove possible existing table:
            dropQuery = "DROP TABLE IF EXISTS "+tableName
            cur.execute(dropQuery)
            # Create new table:
            createQuery = "CREATE TABLE IF NOT EXISTS '"+tableName+"' ("\
                    "'tr_stable_id' varchar(128) NOT NULL default '',"\
                    "'chr' char(50) NOT NULL default '',"\
                    "'strand' int(2) NOT NULL default '',"\
                    "'start' int(10) NOT NULL default '',"\
                    "'start_codon' varchar(128) NOT NULL default '',"\
                    "'stop' int(10) NOT NULL default '',"\
                    "'starts_list' varchar(512) NOT NULL default '',"\
                    "'ends_list' varchar(512) NOT NULL default '',"\
                    "'dist_to_transcript_start' int(10) NOT NULL default '',"\
                    "'dist_to_aTIS' int(10) NOT NULL default 'NA' ,"\
                    "'annotation' varchar(128) NOT NULL default 'NA',"\
                    "'aTIS_call' varchar(128) NOT NULL default 'NA',"\
                    "'peak_shift' int(2) NOT NULL default '',"\
                    "'count'float default NULL,"\
                    "'Rltm_min_Rchx' decimal(11,8) NOT NULL default '0',"\
                    "'coverage' decimal (11,8) NOT NULL default '0',"\
                    "'FPKM' decimal (11,8) NOT NULL default '0',"\
                    "'SNP' varchar(256) NOT NULL default '0',"\
                    "'INDEL' varchar(256) NOT NULL default '0',"\
                    "'tr_seq' TEXT NOT NULL default '',"\
                    "'aa_seq' TEXT NOT NULL default '')"
            cur.execute(createQuery)
            #Store info from chromosomal csv files in SQLite DB:
            for chr in chrs:
                chr=str(chr)
                try:
                    os.system("sqlite3 -separator , " + path_to_sqlitedb + " \".import " + workdir + "/tmp/" + chr + "_tmp.csv " + tableName + "\"")
                except:
                    print ("CSV to SQLite failed for chromosome "+str(chr))
            # Remove tmp files
            for chr in chrs:
                chr=str(chr)
                os.system("rm -rf "+workdir+"/tmp/"+chr+"_tmp.csv")
    except:
        traceback.print_exc()
    return

### Construct translation product ###
def construct_translation_product(species,assembly,chrs,ensembldb,workdir,sqlitedb,snp,indel,analysis_ID,nr_of_cores):
    '''
    :param species:
    :param assembly:
    :param chrs:
    :param ensembldb:
    :param workdir:
    :param sqlitedb:
    :param snp:
    :param indel:
    :param analysis_ID:
    :return:
    '''
    # Get seq_region_IDs for every chromosome
    chrs_seqRegionID = get_chrs_sequence(ensembldb, assembly, chrs,species)

    #Open each chromosome in separated core:
    pool = multiprocessing.Pool(processes=nr_of_cores)
    [pool.apply_async(assembly_per_chr, args=(species,assembly,chrs,ensembldb,workdir,sqlitedb,snp,indel,analysis_ID,chr)) for chr in chrs]

    print (" Waiting for all processes to finish...")
    pool.close()
    pool.join()

    # Store in SQLite database.
    print(' Store the results in the SQLite database.')
    store_in_db(chrs,workdir,sqlitedb,snp,analysis_ID,indel)
    return

def assembly_per_chr(species,assembly,chrs,ensembldb,workdir,sqlitedb,snp,indel,analysis_ID,chr):
    #for chr in chrs:
    print (chr)
    sys.stdout.flush()

    # Output  db-file per process
    chr_tmp_file=open(workdir+"/tmp/"+str(chr)+"_tmp.csv","wb")
    csvWriter = csv.writer(chr_tmp_file)

    # Get SNPs and INDELs for this chromosome
    print(' Get SNPs and INDELs for this chromosome')
    if ( snp != 'NO' ):
        all_SNPs_for_this_chromosome = get_SNPs_per_chromosome(sqlitedb,chr,snp) #command line: snp = samtools
        #print(all_SNPs_for_this_chromosome)
    else:
        print(' ---> There are no SNPs incorporated')
        all_SNPs_for_this_chromosome={}

    if ( indel != 'NO' ):
        all_INDELS_for_this_chromosome= get_INDELS_per_chromosome(sqlitedb,chr,indel)
        #print(all_INDELS_for_this_chromosome)
    else:
        print(' ---> There are no INDELS incorporated')
        all_INDELS_for_this_chromosome={}

    # Get all the transcripts for this chromosome
    print(' Get all the transcripts for this chromosome')
    transcript_per_chromosome=get_transcripts_per_chromosome(sqlitedb,chr,analysis_ID)
    #print(transcript_per_chromosome)
    print(' ')

    #Iterate over the transcripts
    for key in transcript_per_chromosome:

        #Init
        all_INDELS_and_SNPs_for_this_sequence={}
        sort_all_INDELS_and_SNPs_for_this_sequence={}
        TIS = 0
        ref_seq = ''

        # Get the start site of this transcript.
        # This will be used to determine the first exon of the transcript sequence.
        start_site = transcript_per_chromosome[key]['start']

        #Test dit eerst voor 1 transcript, later mag dit weg
        #if key == '11017976_2009580':
        transcript_start=key

        # Get all the exons of the transcript.
        transcript_ID = transcript_per_chromosome[key]['transcript_id']
        all_exons = get_exons_for_transcript(ensembldb, transcript_ID)
        # Determine the strand: 1 = sense; -1 = antisense.
        strand=all_exons[1]['seq_region_strand']

        # Order the exons ascending with their rank.
        exons_for_transcript = {}
        for key in all_exons:
            # print key
            value = all_exons[key]
            # print(value)
            exons_for_transcript[key] = value
        exons_for_transcript = OrderedDict(sorted(exons_for_transcript.items()))
        #print('exons for transcript',exons_for_transcript)

        ## Iterate over the exons.
        # While iterating over the exons, information about the ORF structure will be stored in a dictionary called exon_in_orf.
        # Together with extra transcript information, this dictionary will be stored in a bigger dictionary called tmp_orf_structure.
        CDS_tmp_length = 0
        exon_in_orf={}
        tmp_orf_structure = {}
        length_exon_structure=0

        for rank in exons_for_transcript:
            #Init
            all_SNPs_for_this_exon={}
            all_INDELS_for_this_exon={}
            INDELS_and_SNPs_per_exon={}
            sort_INDELS_and_SNPs_for_this_exon = {}

            # Get the begin and end positons of the exon.
            # This will be used to determine the first exon of the transcript sequence and to catch the exon sequence.
            exon_start = exons_for_transcript[rank]['seq_region_start']
            exon_end = exons_for_transcript[rank]['seq_region_end']

            # Get information about the exons. Later on, this information will be stored in the dictionary tmp_orf_structure.
            tr_stable_id=exons_for_transcript[rank]['transcript_stable_id']
            exon_stable_id=exons_for_transcript[rank]['exon_stable_id']
            end_phase=exons_for_transcript[rank]['end_phase']
            start_phase=exons_for_transcript[rank]['phase']

            if TIS==1:
                # Get the sequence of this exon and append it to the sequences of previous exons.
                sequence = get_sequence(chr, exon_start, exon_end)
                if strand==1:
                    exon_sequence=sequence
                if strand==-1:
                    exon_sequence=get_reverse_complement(sequence)
                ref_seq+=exon_sequence

                # Get all the SNPs and INDELS for this exon.
                all_SNPs_for_this_exon = fetch_SNPs_per_exon(exon_start, exon_end, all_SNPs_for_this_chromosome,CDS_tmp_length,strand)
                all_INDELS_for_this_exon = fetch_INDELS_per_exon(exon_start, exon_end, all_INDELS_for_this_chromosome,CDS_tmp_length,strand)

                # Put all the variants into one dictionary.
                INDELS_and_SNPs_per_exon = merge_INDELS_and_SNPs_per_exon(all_SNPs_for_this_exon, all_INDELS_for_this_exon)

                CDS_tmp_length = exon_end + 1 + CDS_tmp_length - exon_start

                ## Collect information about the ORF structure.
                orf_info={}
                orf_info['exon_start'] = exon_start
                orf_info['exon_end'] = exon_end
                orf_info['end_phase'] = end_phase
                orf_info['start_phase'] = start_phase
                orf_info['tr_stable_id'] = tr_stable_id
                orf_info['e_stable_id'] = exon_stable_id
                orf_info['rank']=rank
                # Calculate the length of the exon.
                length_exon = exon_end - exon_start + 1
                orf_info['length_exon']=length_exon
                # Store the acquired information in the dictionary exon_in_orf.
                exon_in_orf[rank] = orf_info
                length_exon_structure += length_exon

            elif (exon_start <= start_site and start_site <= exon_end):
                # Change the value of TIS from 0 to 1. In this way, only the exons following the start exon will be considered.
                TIS = 1

                if strand==1:
                    # Get the sequence of this exon.
                    sequence = get_sequence(chr, start_site, exon_end)
                    exon_sequence=sequence
                    # Get all the SNPs and INDELS for this exon.
                    all_SNPs_for_this_exon = fetch_SNPs_per_exon(start_site, exon_end, all_SNPs_for_this_chromosome,CDS_tmp_length, strand)
                    all_INDELS_for_this_exon = fetch_INDELS_per_exon(start_site, exon_end,all_INDELS_for_this_chromosome, CDS_tmp_length,strand)
                    # Put all the variants into one dictionary.
                    INDELS_and_SNPs_per_exon = merge_INDELS_and_SNPs_per_exon(all_SNPs_for_this_exon, all_INDELS_for_this_exon)
                    # Calculate the CDS length up until now.
                    CDS_tmp_length = exon_end - start_site + 1

                    ## Collect information about the ORF structure.
                    orf_info = {}
                    orf_info['exon_start'] = exon_start
                    orf_info['exon_end'] = exon_end
                    orf_info['end_phase'] = end_phase
                    orf_info['start_phase'] = start_phase
                    orf_info['tr_stable_id'] = tr_stable_id
                    orf_info['e_stable_id'] = exon_stable_id
                    orf_info['rank']=rank
                    # Calculate the length of the exon.
                    length_exon = exon_end - exon_start + 1
                    orf_info['length_exon'] = length_exon
                    # Determine the start codon.
                    StartCodon=sequence[0:3]
                    # Store the acquired information in the dictionary exon_in_orf.
                    exon_in_orf[rank]=orf_info
                    #Add transcript information to the dictionary tmp_orf_structure.
                    tmp_orf_structure['chr']=chr
                    tmp_orf_structure['strand']=strand
                    tmp_orf_structure['start_codon']=StartCodon
                    tmp_orf_structure['transcript_id']=transcript_ID

                if strand==-1:
                    # Get the sequence of this exon.
                    sequence = get_sequence(chr, exon_start, start_site)
                    exon_sequence=get_reverse_complement(sequence)
                    # Get all the SNPs and INDELS for this exon.
                    all_SNPs_for_this_exon = fetch_SNPs_per_exon(exon_start, start_site, all_SNPs_for_this_chromosome,CDS_tmp_length, strand)
                    all_INDELS_for_this_exon = fetch_INDELS_per_exon(exon_start, start_site,all_INDELS_for_this_chromosome, CDS_tmp_length,strand)
                    # Put all the variants into one dictionary.
                    INDELS_and_SNPs_per_exon = merge_INDELS_and_SNPs_per_exon(all_SNPs_for_this_exon, all_INDELS_for_this_exon)
                    # Calculate the CDS length up until now.
                    CDS_tmp_length = start_site - exon_start + 1

                    ## Collect information about the ORF structure.
                    orf_info={}
                    orf_info['exon_start'] = exon_start
                    orf_info['exon_end'] = exon_end
                    orf_info['end_phase'] = end_phase
                    orf_info['start_phase'] = start_phase
                    orf_info['tr_stable_id'] = tr_stable_id
                    orf_info['e_stable_id'] = exon_stable_id
                    orf_info['rank']=rank
                    # Calculate the length of the exon.
                    length_exon = exon_end - exon_start+ 1
                    orf_info['length_exon'] = length_exon
                    # Determine the start codon.
                    StartCodon = exon_sequence[0:3]
                    # Store the acquired information in the dictionary exon_in_orf.
                    exon_in_orf[rank]=orf_info
                    # Add transcript information to the dictionary tmp_orf_structure.
                    tmp_orf_structure['chr'] = chr
                    tmp_orf_structure['strand'] = strand
                    tmp_orf_structure['start_codon'] = StartCodon
                    tmp_orf_structure['transcript_id'] = transcript_ID

                # Put the obtained sequence in the variable ref_seq (=reference sequence).
                ref_seq+=exon_sequence
                # Calculate the length of the exon structure.
                length_exon_structure+=length_exon

            else:
                pass #You do nothing. In this way the exons before the TIS will not be included in the following steps.

            # Collect the information in the dictionary tmp_orf_structure.
            tmp_orf_structure['orf_structure'] = exon_in_orf
            tmp_orf_structure['length_exon_structure']=length_exon_structure

            #Collect all the INDELs and SNPs of all the exons for this sequence in ascending order.
            all_INDELS_and_SNPs_for_this_sequence.update(INDELS_and_SNPs_per_exon)
            for genomeposition in all_INDELS_and_SNPs_for_this_sequence:
                #print(genomeposition)
                variant_info = all_INDELS_and_SNPs_for_this_sequence[genomeposition]
                sequence_position = variant_info['sequence_position']
                sort_all_INDELS_and_SNPs_for_this_sequence[sequence_position] = variant_info
            sort_all_INDELS_and_SNPs_for_this_sequence = OrderedDict(sorted(sort_all_INDELS_and_SNPs_for_this_sequence.items()))

        #print('sort all INDELS and SNPs for this sequence',sort_all_INDELS_and_SNPs_for_this_sequence)


        ## Collect all the information in one dictionary: all_possible_sequences_for_this_transcript.
        # Init
        old_offset=0
        reference_sequence={}
        start_sequence={}
        variants={}
        all_possible_sequences_for_this_transcript = {}
        # Store in the dictionary.
        start_sequence['sequence']=ref_seq
        start_sequence['offset']=old_offset
        start_sequence['variant info']=variants
        start_sequence['orf_structure'] = tmp_orf_structure
        reference_sequence['posibility1']=start_sequence
        all_possible_sequences_for_this_transcript.update(reference_sequence)
        #print(all_possible_sequences_for_this_transcript)

        # Loop over all the INDELs and SNPs and build them into the sequence.
        for position_in_sequence in sort_all_INDELS_and_SNPs_for_this_sequence:
            #print(position_in_sequence)
            variant = sort_all_INDELS_and_SNPs_for_this_sequence[position_in_sequence]
            type = variant['type']
            allelic_frequency = variant['af']
            alternative=variant['alt']
            reference=variant['ref']
            genomic_position=variant['pos']
            #print (alternative)
            #print(variant)

            all_sequences={}
            i=0
            original_position_in_sequence=position_in_sequence


            if type == 'SNP':
                # Allele frequence equals 1: only the alternative need to be retained.
                if allelic_frequency >= 0.90:
                    # Iterate over all the sequences present in the dictionary all_possible_sequences_for_this_transcript.
                    for sequences in list(all_possible_sequences_for_this_transcript):
                        for base in alternative:
                            # Use a counter to determine the number of possible sequences.
                            i += 1

                            # Get the existing sequence and its offset from the dictionary all_possible_sequences_for_this_transcript.
                            seq = all_possible_sequences_for_this_transcript[sequences]
                            old_sequence = seq['sequence']
                            old_offset = seq['offset']

                            # The positions of the variants in the dictionary sort_all_INDELS_and_SNPs_for_this_sequence are given at offset 0.
                            # Due to possibility of a change in the offset, a new position has to be calculated:
                            position_in_sequence=original_position_in_sequence+old_offset

                            # A new sequence including the alternative SNP is build up.
                            new_sequence = old_sequence[:position_in_sequence - 1] + base + old_sequence[position_in_sequence:]

                            ## Make a dictionary containing the information about all the variants.
                            # Collect variant information of variants that are already build in.
                            previous_variants = seq['variant info']  # This gives the information of the variants that is already build in this sequence.
                            # Collect variant information of the new variant that has to build in.
                            this_variant = {}
                            this_variant['position in sequence'] = position_in_sequence
                            this_variant['alternative'] = base
                            this_variant['reference'] = reference
                            this_variant['allelic frequency'] = allelic_frequency
                            this_variant['genomic position'] = genomic_position
                            this_variant['type'] = 'SNP'
                            this_variant['reference/alternative'] = 'alternative'
                            # The new variant and the previous variants are put together.
                            all_variants = {}
                            all_variants[position_in_sequence] = this_variant
                            all_variants.update(previous_variants)
                            all_variants = OrderedDict(sorted(all_variants.items()))
                            # The variant information, the new sequence and the final offset are collected into one dictionary.
                            new_seq = {}
                            possibility = {}
                            new_seq['variant info'] = all_variants
                            new_seq['sequence'] = new_sequence
                            new_seq['offset'] = old_offset
                            new_seq['orf_structure']=tmp_orf_structure
                            possibility['possibility' + str(i)] = new_seq

                            # Collect all the possible sequences generated within this loop.
                            all_sequences.update(possibility)

                    all_possible_sequences_for_this_transcript = all_sequences

                #Allele frequence equals 0.5: both the reference and the alternative need to be retained.
                elif (0.40 <= allelic_frequency <= 0.60):
                    #Iterate over all the sequences present in the dictionary all_possible_sequences_for_this_transcript.
                    for sequences in list(all_possible_sequences_for_this_transcript):
                        # Use a counter to determine the number of possible sequences.
                        i+=1

                        # Get the existing sequence and its offset from the dictionary all_possible_sequences_for_this_transcript.
                        seq = all_possible_sequences_for_this_transcript[sequences]
                        old_sequence = seq['sequence']
                        old_offset = seq['offset']

                        # Put the reference sequence and the current offset in the dictionary old_seq.
                        old_seq = {}
                        old_seq['sequence'] = old_sequence
                        old_seq['offset'] = old_offset

                        # The positions of the variants in the dictionary sort_all_INDELS_and_SNPs_for_this_sequence are given at offset 0.
                        # Due to possibility of a change in the offset, a new position has to be calculated:
                        position_in_sequence = original_position_in_sequence + old_offset

                        ## Make a dictionary containing the information about all the variants.
                        # Collect variant information of variants that are already build in.
                        previous_variants = seq['variant info']
                        # Collect variant information of the new variant that has to be build in.
                        this_variant = {}
                        this_variant['position in sequence'] = position_in_sequence
                        this_variant['alternative'] = alternative
                        this_variant['reference'] = reference
                        this_variant['allelic frequency'] = allelic_frequency
                        this_variant['genomic position'] = genomic_position
                        this_variant['type'] = 'SNP'
                        this_variant['reference/alternative']='reference'
                        # The new variant and the previous variants are put together.
                        all_variants={}
                        all_variants[position_in_sequence] = this_variant
                        all_variants.update(previous_variants)
                        all_variants = OrderedDict(sorted(all_variants.items()))
                        # The reference sequence, the current offset and the complete variant information are collected in the dictionary possibility.
                        old_seq['variant info']=all_variants
                        old_seq['orf_structure']=tmp_orf_structure
                        possibility = {}
                        possibility['possibility' + str(i)] = old_seq

                        # Collect all the possible sequences generated within this loop.
                        all_sequences.update(possibility)

                        ### Add the alternative sequences.
                        #iterate over the alternative SNPs.
                        for base in alternative:
                            # Use a counter to determine the number of possible sequences.
                            i+=1

                            # Get the existing sequence and its offset from the dictionary all_possible_sequences_for_this_transcript.
                            new_sequence = old_sequence[:position_in_sequence - 1] + base + old_sequence[position_in_sequence:]

                            ## Make a dictionary containing the information about all the variants.
                            # Collect variant information of variants that are already build in.
                            previous_variants = seq['variant info']
                            # Collect variant information of the new variant that has to build in.
                            this_variant = {}
                            this_variant['position in sequence'] = position_in_sequence
                            this_variant['alternative'] = base
                            this_variant['reference'] = reference
                            this_variant['allelic frequency'] = allelic_frequency
                            this_variant['genomic position'] = genomic_position
                            this_variant['type'] = 'SNP'
                            this_variant['reference/alternative'] = 'alternative'
                            # The new variant and the previous variants are put together.
                            all_variants = {}
                            all_variants[position_in_sequence] = this_variant
                            all_variants.update(previous_variants)
                            all_variants = OrderedDict(sorted(all_variants.items()))
                            # The new sequence, the current offset and the complete variant information are collected in the dictionary possibility.
                            new_seq = {}
                            possibility = {}
                            new_seq['sequence']=new_sequence
                            new_seq['offset'] = old_offset
                            new_seq['variant info']=all_variants
                            new_seq['orf_structure']=tmp_orf_structure
                            possibility['possibility'+str(i)]=new_seq

                            # Collect all the possible sequences generated within this loop.
                            all_sequences.update(possibility)

                    all_possible_sequences_for_this_transcript = all_sequences

            else:
                # Allele frequence equals 1: only the alternative need to be retained.
                if allelic_frequency >= 0.90:
                    # Iterate over all the sequences present in the dictionary all_possible_sequences_for_this_transcript.
                    for sequences in list(all_possible_sequences_for_this_transcript):
                        # Iterate over all the possible alternative INDELs.
                        for bases in alternative:
                            # Use a counter to determine the number of possible sequences.
                            i += 1

                            # Get the existing sequence and its offset from the dictionary all_possible_sequences_for_this_transcript.
                            seq = all_possible_sequences_for_this_transcript[sequences] #gives all the information about a sequence in the dictionary all_possible_sequences_for_this_transcript: variant info, sequence, offset.
                            old_sequence = seq['sequence']#gives the current sequence
                            old_offset = seq['offset']#gives the current offset

                            # The positions of the variationary sort_all_INDELS_and_SNPs_for_this_sequence are given at offset 0.
                            # Due to possibility of a change in the offnts in the dicset, a new position has to be calculated:
                            position_in_sequence=original_position_in_sequence+old_offset

                            # A new sequence including the alternative INDEL is build up.
                            new_sequence = old_sequence[:position_in_sequence - 1] + bases + old_sequence[position_in_sequence-1+len(reference):]  # python lists start at index 0

                            ## Make a dictionary containing the information about all the variants.
                            # Collect variant information of variants that are already build in.
                            previous_variants = seq['variant info']  # This gives the information of the variants that is already build in this sequence.
                            # Collect variant information of the new variant that has to build in.
                            this_variant = {}
                            this_variant['position in sequence'] = position_in_sequence
                            this_variant['alternative'] = bases
                            this_variant['reference'] = reference
                            this_variant['allelic frequency'] = allelic_frequency
                            this_variant['genomic position'] = genomic_position
                            this_variant['type'] = 'INDEL'
                            this_variant['reference/alternative'] = 'alternative'
                            # The new variant and the previous variants are put together.
                            all_variants = {}
                            all_variants[position_in_sequence] = this_variant
                            all_variants.update(previous_variants)
                            all_variants = OrderedDict(sorted(all_variants.items()))
                            # The variant information, the new sequence and the final offset are collected into one dictionary.
                            new_seq = {}
                            possibility = {}
                            new_seq['variant info'] = all_variants
                            new_seq['sequence'] = new_sequence
                            new_seq['offset'] = old_offset + len(bases) - len(reference)
                            new_seq['orf_structure']=tmp_orf_structure
                            possibility['possibility' + str(i)] = new_seq

                            # Collect all the possible sequences generated within this loop.
                            all_sequences.update(possibility)
                    all_possible_sequences_for_this_transcript = all_sequences

                # Allele frequence equals 0.5: both the reference and the alternative variant need to be retained.
                elif (0.40 <= allelic_frequency <= 0.60):
                    # Iterate over all the sequences present in the dictionary all_possible_sequences_for_this_transcript.
                    for sequences in list(all_possible_sequences_for_this_transcript):
                        ## Add reference sequence
                        ## Put the sequence and the offset of the reference sequence in a dictionary old_seq.
                        # Use a counter to determine the number of possible sequences.
                        i+=1

                        # Get the existing sequence and its offset from the dictionary all_possible_sequences_for_this_transcript.
                        seq = all_possible_sequences_for_this_transcript[sequences]
                        old_sequence = seq['sequence']
                        old_offset = seq['offset']
                        # The positions of the variants in the dictionary sort_all_INDELS_and_SNPs_for_this_sequence are given at offset 0.
                        # Due to possibility of a change in the offset, a new position has to be calculated:
                        position_in_sequence = original_position_in_sequence + old_offset

                        # Put the reference sequence and the current offset in the dictionary old_seq.
                        old_seq = {}
                        old_seq['sequence'] = old_sequence
                        old_seq['offset'] = old_offset

                        ## Make a dictionary containing the information about all the variants.
                        # Collect variant information of variants that are already build in.
                        previous_variants = seq['variant info']
                        # Collect variant information of the new variant that has to be build in.
                        this_variant = {}
                        this_variant['position in sequence'] = position_in_sequence
                        this_variant['alternative'] = alternative
                        this_variant['reference'] = reference
                        this_variant['allelic frequency'] = allelic_frequency
                        this_variant['genomic position'] = genomic_position
                        this_variant['type'] = 'INDEL'
                        this_variant['reference/alternative'] = 'reference'
                        # The new variant and the previous variants are put together.
                        all_variants = {}
                        all_variants[position_in_sequence] = this_variant
                        all_variants.update(previous_variants)
                        all_variants = OrderedDict(sorted(all_variants.items()))

                        # The reference sequence, the current offset and the complete variant information are collected in the dictionary possibility.
                        old_seq['variant info'] = all_variants
                        old_seq['orf_structure']=tmp_orf_structure
                        possibility = {}
                        possibility['possibility' + str(i)] = old_seq

                        # Collect all the possible sequences generated within this loop.
                        all_sequences.update(possibility)

                        ### Add the alternative sequences.
                        # Iterate over all the posible alternative INDELs.
                        for bases in alternative:
                            # Use a counter to determine the number of possible sequences.
                            i+=1

                            # Get the existing sequence and its offset from the dictionary all_possible_sequences_for_this_transcript.
                            seq = all_possible_sequences_for_this_transcript[sequences]
                            old_sequence = seq['sequence']
                            old_offset = seq['offset']

                            # A new sequence including the alternative INDEL is build up.
                            new_sequence=old_sequence[:position_in_sequence-1]+bases+old_sequence[position_in_sequence-1+len(reference):]

                            ## Make a dictionary containing the information about all the variants.
                            # Collect variant information of variants that are already build in.
                            previous_variants = seq['variant info']
                            # Collect variant information of the new variant that has to build in.
                            this_variant = {}
                            this_variant['position in sequence'] = position_in_sequence
                            this_variant['alternative'] = bases
                            this_variant['reference'] = reference
                            this_variant['allelic frequency'] = allelic_frequency
                            this_variant['genomic position'] = genomic_position
                            this_variant['type'] = 'INDEL'
                            this_variant['reference/alternative'] = 'alternative'

                            # The new variant and the previous variants are put together.
                            all_variants = {}
                            all_variants[position_in_sequence] = this_variant
                            all_variants.update(previous_variants)
                            all_variants = OrderedDict(sorted(all_variants.items()))
                            # The new sequence, the current offset and the complete variant information are collected in the dictionary possibility.
                            new_seq = {}
                            possibility = {}
                            new_seq['sequence'] = new_sequence
                            new_seq['offset'] = old_offset + len(bases) - len(reference)
                            new_seq['variant info'] = all_variants
                            new_seq['orf_structure']=tmp_orf_structure
                            possibility['possibility' + str(i)] = new_seq

                            # Collect all the possible sequences generated within this loop.
                            all_sequences.update(possibility)

                    all_possible_sequences_for_this_transcript = all_sequences

        #print('All possible sequences for this transcript:',all_possible_sequences_for_this_transcript)


        # In the following loop, proteins containing a stop codon and showing non-synonymous SNPs are retained and translated.
        # All the non-synonymous proteins will be stored in a dictionary called: non_synonymous_proteins:
        non_synonymous_proteins = {}
        # We use a counter to count the number of final possible protein sequences.
        counter=0

        # Order the exons ascending.
        exons_in_orf = {}
        for key in exon_in_orf:
            #print key
            value = exon_in_orf[key]
            #print(value)
            exons_in_orf[key] = value
        exons_in_orf = OrderedDict(sorted(exons_in_orf.items()))
        #print('exons in orf',exons_in_orf)

        # Translation of the mRNA-sequences into protein sequences:
        for possible_sequence in all_possible_sequences_for_this_transcript:
            #print('Possible sequence:',possible_sequence)
            #print(all_possible_sequences_for_this_transcript[possible_sequence])
            keep_protein='NO'

            # Control if first codon is (near-)cognate and replace near-cognate start with the cognate methionine.
            mrna_sequence = all_possible_sequences_for_this_transcript[possible_sequence]['sequence']
            start_codon = mrna_sequence[0:3]
            near_cognate = re.compile(r"[ACTG]TG|A[ACTG]G|AT[ACTG]")
            if re.match(near_cognate, start_codon):
                mrna_sequence = "ATG" + mrna_sequence[3:]
                keep_protein='YES'

            # Translation of the mRNA-sequences into protein sequences:
            protein_sequence, complete_protein_sequence,position_stop_in_sequence = translate_sequence(mrna_sequence)
            #print('complete protein',complete_protein_sequence)

            # Only complete proteins (=derived from a mRNA sequence containing a stop codon) will be retained.
            if complete_protein_sequence=='YES':
                # Calculate the length of the mRNA sequence from start codon up to and including the stop codon.
                length_mrna_seq = len(protein_sequence) * 3

                # Determine how long the translated sequence spans over the exon structure.
                orf_structure={}
                length_exon_structure=all_possible_sequences_for_this_transcript[possible_sequence]['orf_structure']['length_exon_structure']

                orf_structure['length_mrna_sequence']=length_mrna_seq
                orf_structure['length_exon_structure']=length_exon_structure
                # Catch the part of the mRNA sequence from start to stop codon.
                mrna_sequence=mrna_sequence[0:position_stop_in_sequence+3]


                # We will use a variable stop_codon_in_INDEL to determine if there are sequences whith a stop codon lying in an INDEL.
                stop_codon_in_INDEL = 'NO'

                ### Make underscore seperated sequences of start and stop coordinates.
                # These sequences also contain the start site and the genomic position of the stop codon.
                # The start site is already given, but the stop site needs to be determined.
                # There are two possibilities for the position of a stop codon: inside or outside an INDEL.
                # The genomic position of a stop codon inside an INDEL will be determined in the next loop.
                # The genomic position of a stop codon not lying in an INDEL will be determined later on.
                if strand==1:
                    # Init
                    pos_shift_stop_codon = 0
                    stop_codon_is_determined='NO'
                    mismatch =0

                    # Loop over all the INDELs and SNPs.
                    SNP_INDELS = all_possible_sequences_for_this_transcript[possible_sequence]['variant info']
                    for snp_indels in SNP_INDELS:
                        #print(snp_indels)
                        variant_type = SNP_INDELS[snp_indels]['type']
                        if variant_type == 'INDEL':
                            reference_alternative = SNP_INDELS[snp_indels]['reference/alternative']
                            # Only for the alternative sequence, you need to calculate the genomic stop position by using this loop.
                            # For the reference sequence, you don't need to take into account an offset.
                            if reference_alternative == 'alternative':
                                # Only take into account the INDELs that are situated before the stop codon.
                                if snp_indels<=length_mrna_seq:
                                    INDEL_position=SNP_INDELS[snp_indels]['position in sequence']
                                    reference_INDEL = SNP_INDELS[snp_indels]['reference']
                                    alternative_INDEL = SNP_INDELS[snp_indels]['alternative']
                                    genomic_position_INDEL=SNP_INDELS[snp_indels]['genomic position']

                                    # Determine if the position of the last base of the stop codon is lying inside the INDEL.

                                    if INDEL_position <= length_mrna_seq <= INDEL_position + len(alternative_INDEL)-1:
                                        # Add an extra letter 'X' to the shortest indel. This is needed for INDELS where the alternative is composed by the reference plus
                                        # an extra sequence. (For example ATT and ATTCC). This is important for the next step in which the position is determined. (Comparing
                                        # ATT and ATTCC will give an error due to the impossibility to compare strings with different lengths. In contrast, comparing ATTX and ATTC will give the first difference in
                                        # bases and will then stop the loop which avoids the length problem.

                                        length_shortest_indel = min(len(reference_INDEL), len(alternative_INDEL))
                                        if len(reference_INDEL) == length_shortest_indel:
                                            reference_INDEL += 'X'
                                        if len(alternative_INDEL) == length_shortest_indel:
                                            alternative_INDEL += 'X'

                                        stop_codon_in_INDEL='YES'

                                        # Determine the first position in the reference and alternative INDEL where these two differ in base.
                                        for i in range(len(reference_INDEL)):
                                            if reference_INDEL[i] != alternative_INDEL[i]:
                                                mismatch=i
                                                break  # You cannot compare sequences of different lengths. Therefore you have to stop afther finding the first difference.
                                            # This is no problem because you only need this first difference.


                                        if length_mrna_seq<INDEL_position+mismatch:
                                            # If the stop codon is lying before the mismatch, the genomic position will be determined later on.
                                            stop_codon_is_determined='NO'
                                        else:
                                            stop_codon_is_determined='YES'
                                            position_stop_codon=genomic_position_INDEL+mismatch-1
                                    else:
                                        # If an INDEL is incorporated before a stop codon, an offset will be calculated.
                                        # This offset will be used later on to determine the reference genomic position of the stop codon.
                                        pos_shift_stop_codon += (len(reference_INDEL) - len(alternative_INDEL))


                        if variant_type == 'SNP':
                            pos_shift_stop_codon += 0
                    #print('Position shift', pos_shift_stop_codon)

                    ## Collect all the exon start and stop coordinates and put them in a list.
                    all_end_positions = []
                    all_start_positions = []
                    # Loop over the exons incorporated in the mRNA-transcripts.
                    # The start site is lying in the first exon and represents the end of the sequence.In order
                    # to incorporate this start site in the all_end_positions list and not the end of the first exon, a parameter first_exon is used and set to YES.
                    first_exon = 'YES'
                    intron = 0
                    length_including_this_exon = 0
                    only_one_exon='YES'


                    for exons in exons_in_orf:
                        end = exons_in_orf[exons]['exon_end']
                        start = exons_in_orf[exons]['exon_start']

                        if first_exon == 'YES':
                            # Append the start site to the list with start positions.
                            all_start_positions.append(str(start_site))
                            # Set first_exon to 'NO'. In this way, you only take this loop once.
                            first_exon = 'NO'
                            # Calculate the length of the transcript up until now.
                            length_including_this_exon += (end - start_site+1)
                            # If this length is smaller than the length of the mRNA sequence, the mRNA sequence overspans at least one more exon.
                            if length_including_this_exon<length_mrna_seq:
                                # Append the exon end position to the list with end positions.
                                all_end_positions.append(str(end))
                                # Start with the calculation of the intron length.
                                intron += -end
                                only_one_exon='NO' # By setting this variable to 'NO', it is clear that the mRNA sequence overspans at least one more exon.

                        else:
                            if only_one_exon=='NO':
                                # Append the exon start position to the list of start positions.
                                all_start_positions.append(str(start))

                                # The last exon contains the stop codon. You want to incorporate the position of this stopcodon
                                # in the all_end_positions list and not the end position of this exon.
                                # Therefore you add the end positions to the list as long as these positions are smaller than
                                # the actual stop site. This will be the case for all exons except the last one. Here, the position of the
                                # stopcodon will be included. Therefore, calculate the length of the transcrip up unil now.
                                # If this length is smaller than the length of the mRNA sequence, the mRNA sequence overspans at least one more exon.
                                length_including_this_exon += (end - start+1)
                                #print('Length including this exon', length_including_this_exon)
                                if length_including_this_exon < length_mrna_seq:
                                    # Append the exon end position to the list with end positions.
                                    all_end_positions.append(str(end))
                                    # Continue with the calculation of the intron length.
                                    intron += (start - end) - 1
                                else:
                                    # Complete the calculation of the intron length.
                                    intron += start - 1
                                    break
                            else:
                                #print('Last exon')
                                break

                    # Calculate the genomic position of the stop codon.
                    if stop_codon_is_determined=='NO':
                        genomic_position_stop_codon = start_site +length_mrna_seq - 1 + intron + pos_shift_stop_codon
                    else:
                        genomic_position_stop_codon=position_stop_codon
                    #print('Genomic position stop codon',genomic_position_stop_codon,pos_shift_stop_codon)

                    # Append this genomic position to the list of all exon end positions.
                    all_end_positions.append(str(genomic_position_stop_codon))

                    # Join all the ends in an underscore separated sequence called stops_list.
                    stops_list = ""
                    separator = "_"
                    stops_list = separator.join(all_end_positions)
                    # Join all the starts in an underscore separated sequence called starts_list.
                    starts_list = ""
                    starts_list = separator.join(all_start_positions)
                    #print('Stop list',stops_list)
                    #print('Start list',starts_list)

                if strand==-1:
                    # Init
                    only_one_exon='YES'
                    pos_shift_stop_codon = 0
                    stop_codon_is_determined = 'NO'
                    mismatch = 0

                    # Loop over all the INDELs and SNPs.
                    SNP_INDELS = all_possible_sequences_for_this_transcript[possible_sequence]['variant info']
                    for snp_indels in SNP_INDELS:
                        #print(snp_indels)
                        variant_type = SNP_INDELS[snp_indels]['type']
                        if variant_type == 'INDEL':
                            reference_alternative = SNP_INDELS[snp_indels]['reference/alternative']
                            # Only for the alternative sequence, you need to calculate the genomic stop position by using this loop.
                            # For the reference sequence, you don't need to take into account an offset.
                            if reference_alternative == 'alternative':
                                # Only take into account the INDELs that are situated before the stop codon.
                                if snp_indels <= length_mrna_seq:
                                    INDEL_position = SNP_INDELS[snp_indels]['position in sequence']
                                    reference_INDEL = SNP_INDELS[snp_indels]['reference']
                                    alternative_INDEL = SNP_INDELS[snp_indels]['alternative']
                                    # Calculate the genomic position of the last base of the INDEL. This represents the genomic position of the first base of this INDEL in antisense direction.
                                    genomic_position_INDEL_antisense = SNP_INDELS[snp_indels]['genomic position']+len(reference_INDEL)-1


                                    # Determine if the stop position is lying inside the INDEL.
                                    if INDEL_position <= length_mrna_seq <= INDEL_position + len(alternative_INDEL) - 1:
                                        # Add an extra letter 'X' to the shortest indel. This is needed for INDELS where the alternative is composed by the reference plus
                                        # an extra sequence. (For example ATT and ATTCC). This is important for the next step in which the position is determined. (Comparing
                                        # ATT and ATTCC will give an error due to the impossibility to compare strings with different lengths. In contrast, comparing ATTX and ATTC will give the first difference in
                                        # bases and will then stop the loop which avoids the length problem.
                                        length_shortest_indel = min(len(reference_INDEL), len(alternative_INDEL))
                                        if len(reference_INDEL) == length_shortest_indel:
                                            reference_INDEL += 'X'
                                        if len(alternative_INDEL) == length_shortest_indel:
                                            alternative_INDEL += 'X'


                                        stop_codon_in_INDEL='YES'

                                        # Determine the first position in the reference and alternative INDEL where these two differ in base.
                                        for i in range(len(reference_INDEL)):
                                            if reference_INDEL[i] != alternative_INDEL[i]:
                                                mismatch = i
                                                break  # You cannot compare sequences of different lengths. Therefore you have to stop afther finding the first difference.
                                            # This is no problem because you only need this first difference.


                                        if length_mrna_seq < INDEL_position + mismatch:
                                            # If the stop codon is lying before the mismatch, the genomic position will be determined later on.
                                            stop_codon_is_determined = 'NO'
                                        else:
                                            stop_codon_is_determined = 'YES'
                                            position_stop_codon = genomic_position_INDEL_antisense - mismatch + 1
                                    else:
                                        # If an INDEL is incorporated before a stop codon, an offset will be calculated.
                                        # This offset will be used later on to determine the reference genomic position of the stop codon.
                                        pos_shift_stop_codon += (len(reference_INDEL) - len(alternative_INDEL))

                        if variant_type == 'SNP':
                            pos_shift_stop_codon += 0

                    ## Collect all the start and stop coordinates and put them in a list.
                    # Collect all the exon start and stop coordinates and put them in a list.
                    all_end_positions=[]
                    all_start_positions=[]
                    # Loop over the exons incorporated in the mRNA-transcripts.
                    # The start site is lying in the first exon and represents the end of the sequence.In order
                    # to incorporate this start site in the all_end_positions list and not the end of the first exon, a parameter first_exon is used and set to YES.
                    first_exon = 'YES'
                    intron=0

                    length_including_this_exon=0
                    for exons in exons_in_orf:
                        end = exons_in_orf[exons]['exon_end']
                        start=exons_in_orf[exons]['exon_start']

                        if first_exon=='YES':
                            # The start site is appended to the list of end positions because the sequence is antisense.
                            all_end_positions.append(str(start_site))
                            # Set first_exon to 'NO'. In this way, you only take this loop once.
                            first_exon = 'NO'
                            # Calculate the length of the transcript up until now.
                            length_including_this_exon+=(start_site-start+1)
                            # If this length is smaller than the length of the mRNA sequence, the mRNA sequence overspans at least one more exon.
                            if length_including_this_exon<length_mrna_seq:
                                # Append the exon start to the list with start positions.
                                all_start_positions.append(str(start))
                                # Start the calculation of the intron length.
                                intron += start
                                only_one_exon='NO' # By setting this variable to 'NO', it is clear that the mRNA sequence overspans at least one more exon.

                        else:
                            if only_one_exon=='NO':
                                # Append the exon end to the list of end positions.
                                all_end_positions.append(str(end))

                                # The last exon contains the stop codon. You want to incorporate the position of this stopcodon
                                # in the all_start_positions list and not the start position of this exon.
                                # Therefore you add the start positions to the list as long as these positions are bigger than
                                # the actual start_site. This will be the case for all exons except the last one. Here, the position of the
                                # stopcodon will be included. Therefore, calculate the length of the transcrip up unil now.
                                # If this length is smaller than the length of the mRNA sequence, the mRNA sequence overspans at least one more exon.
                                length_including_this_exon += (end - start+1)
                                #print(' Length including this exon',length_including_this_exon)
                                if length_including_this_exon<length_mrna_seq:
                                    # Append the exon start to the list of start sites.
                                    all_start_positions.append(str(start))
                                    # Continue the calculation of the intron length.
                                    intron+=(start-end)-1
                                else:
                                    # Complete the calculation of the intron length.
                                    intron+=-end-1
                                    break
                            else:
                                    #print('Last exon')
                                    break

                    # Calculate the genomic position of the stop codon.
                    if stop_codon_is_determined=='NO':
                        genomic_position_stop_codon = start_site - length_mrna_seq + 1 - intron-pos_shift_stop_codon
                    else:
                        genomic_position_stop_codon=position_stop_codon
                    # Append this genomic position to the list with all start positions.
                    all_start_positions.append(str(genomic_position_stop_codon))

                    # Join all the ends in an underscore separated sequence called stops_list.
                    stops_list=""
                    separator="_"
                    stops_list=separator.join(reversed(all_end_positions))
                    # Join all the starts in an underscore separated sequence called starts_list.
                    starts_list=""
                    starts_list=separator.join(reversed(all_start_positions))


                ### Select the non-synonymous SNPs.
                # You only want to retain non-redundant protein sequences. If the reference and the alternative protein
                # sequences are the same, you keep only the reference sequence, except when the allelic frequency is one. Then you keep
                # the alternative sequence. In all cases, only the non-synonymous SNPs will be included in the description line of proteins.

                # Initialize the lists, dictionaries and variables you need later.
                final_protein = {}
                non_synonymous_SNPs = []
                all_INDELS=[]
                all_non_synonymous_SNPs=""
                all_INDELS_in_final_protein=""
                #keep_protein = 'YES' #Not needed anymore. Needs to be reviewed.

                # Collect the information about SNPs and INDELs.
                SNP_INDELS = all_possible_sequences_for_this_transcript[possible_sequence]['variant info']

                # If there are no SNPs or INDELs, there is only one sequence possible: the reference sequence.
                # In order to retain this sequence, we give the variable 'keep_protein' the value YES.
                if SNP_INDELS=={}:
                    #keep_protein='YES'
                    pass

                else:
                    # If there are SNPs or INDELs present, we will loop over every codon and see if the SNP or INDEL
                    # lies inside this codon. When this is the case, the codon in the reference sequence will be compared
                    # to the corresponding codon in the alternative sequence. Only when these codons are leading to two
                    # different amino acids, both the reference and the alternative sequence are retained.If the the reference
                    # and alternative sequence are equal, we will only retain the reference.
                    for i in range(0,len(mrna_sequence)-2,3):
                        # Loop over al the codons in the sequence. This can be the reference or an alternative sequence.
                        alternative_codon = mrna_sequence[i:i + 3]# This codon is named alternative in order to be able to compare it with the reference codon.

                        # Loop over al the SNPs and INDELs.
                        for snp_indels in SNP_INDELS:
                            variant_type=SNP_INDELS[snp_indels]['type']
                            if variant_type=='INDEL':
                                reference_alternative = SNP_INDELS[snp_indels]['reference/alternative']
                                if reference_alternative=='reference':
                                    AF = SNP_INDELS[snp_indels]['allelic frequency']
                                    # If the allelic frequency is 1, the alternative sequence is the only biological relevant sequence.
                                    # Therefore you want to retain this alternative sequence and forget the reference sequence.
                                    #if AF>0.95:
                                        #keep_protein='NO'
                                        #pass
                                    # If the allelic frequency is 0.5, the reference sequence is also biological relevant.
                                    # Therefore you keep this reference sequence.
                                    #else:
                                        #keep_protein='YES'
                                # You keep all sequences(reference and alternatives) when af=0.5, but you only want to include  the INDEL information
                                # in the description line of the proteins if the alternative INDEL is present.
                                if reference_alternative=='alternative':
                                    # Sequences including INDELs are always non-synonymous due to the event of frame shiftings.
                                    # Therefore you want to retain al sequences with INDEls.

                                    # Add the INDEL to the description line of the proteins.

                                    position_of_INDEL = SNP_INDELS[snp_indels]['position in sequence']
                                    if i + 1 <= position_of_INDEL <= i + 3:
                                        #keep_protein = 'YES'
                                        # The INDEL description is composed as follows: position of INDEL in mRNA sequence_reference INDEL_alternative INDEL_allelic frequency.
                                        reference_INDEL = SNP_INDELS[snp_indels]['reference']
                                        alternative_INDEL = SNP_INDELS[snp_indels]['alternative']
                                        AF = SNP_INDELS[snp_indels]['allelic frequency']
                                        description_for_this_INDEL = str(position_of_INDEL)+ '_'+reference_INDEL+'_'+alternative_INDEL+'_'+str(AF)
                                        all_INDELS.append(description_for_this_INDEL)

                            # SNPs can be synonymous. If the SNPs are not leading to a change in amino acids, only the reference
                            # sequences are kept.
                            if variant_type=='SNP':
                                reference_alternative=SNP_INDELS[snp_indels]['reference/alternative']
                                if reference_alternative=='reference':
                                    AF = SNP_INDELS[snp_indels]['allelic frequency']
                                    # If the allelic frequency is 1, the alternative sequence is the only biological relevant sequence.
                                    # Therefore you want to retain this alternative sequence and forget the reference sequence.
                                    #if AF>0.95:
                                        #keep_protein='NO'
                                        #pass
                                    # If the allelic frequency is 0.5, the reference sequence is also biological relevant.
                                    # Therefore you keep this reference sequence.
                                    #else:
                                        #keep_protein='YES'
                                if reference_alternative=='alternative':
                                    # Loop over the SNPs. If a SNP is lying in a codon, the amino acids comming from the alternative
                                    # and the reference codon are compared. Only if these amino acids are different, the SNP
                                    # will be included in the description line of the proteins.
                                    position_of_SNP = SNP_INDELS[snp_indels]['position in sequence']
                                    if i+1<= position_of_SNP<= i+3:
                                        reference_SNP = SNP_INDELS[snp_indels]['reference']
                                        alternative_SNP=SNP_INDELS[snp_indels]['alternative']
                                        AF=SNP_INDELS[snp_indels]['allelic frequency']
                                        #before_SNP=mrna_sequence[i:position_of_SNP-1]
                                        #after_SNP=mrna_sequence[position_of_SNP:i+3]
                                        reference_codon=mrna_sequence[i:position_of_SNP-1]+reference_SNP+mrna_sequence[position_of_SNP:i+3]

                                        #print('reference codon',reference_codon,'alternative codon',alternative_codon)
                                        reference_amino_acid=translate_sequence(reference_codon)[0]
                                        alternative_amino_acid=translate_sequence(alternative_codon)[0]
                                        #print(reference_amino_acid,alternative_amino_acid)

                                        if reference_amino_acid!=alternative_amino_acid:
                                            #keep_protein='YES'
                                            description_for_this_SNP=str(position_of_SNP)+'_'+reference_SNP+'_'+alternative_SNP+'_'+str(AF)
                                            non_synonymous_SNPs.append(description_for_this_SNP)
                                        # If the allelic frequency is 1, you automatically retain the alternative sequence and forget the reference.
                                        #if AF>=0.95:
                                            #keep_protein='YES'

                ### Put the non-redundant sequences in a csv file.
                if keep_protein=='YES':
                    ### Collect all the information you want to store in a SQLite database.
                    tr_stable_id = exons_for_transcript[exons]['transcript_stable_id']
                    dist_to_transcript_start = transcript_per_chromosome[transcript_start]['dist_to_transcript_start']
                    dist_to_aTIS = transcript_per_chromosome[transcript_start]['dist_to_aTIS']
                    annotation = transcript_per_chromosome[transcript_start]['annotation']
                    aTIS_call = transcript_per_chromosome[transcript_start]['aTIS_call']
                    start_codon = transcript_per_chromosome[transcript_start]['start_codon']
                    peak_shift = transcript_per_chromosome[transcript_start]['peak_shift']
                    count = transcript_per_chromosome[transcript_start]['count']
                    Rltm_min_Rchx = transcript_per_chromosome[transcript_start]['Rltm_min_Rchx']

                    # Collect al the non-synonymous SNPs.
                    separator=":"
                    all_non_synonymous_SNPs=separator.join(non_synonymous_SNPs)

                    # INDELs are always non-synonymous, so you want to collect all the INDELS that are incorporated.
                    separator=":"
                    all_INDELS_in_final_protein=separator.join(all_INDELS)




                    # Write the information in a csv file. In a later step, the contents of all generated csv files will be collected in a SQLite database.
                    csvWriter.writerow((tr_stable_id,chr,strand,start_site,start_codon,genomic_position_stop_codon,starts_list,stops_list,dist_to_transcript_start,dist_to_aTIS,annotation,aTIS_call,peak_shift,count,Rltm_min_Rchx,"","",all_non_synonymous_SNPs,all_INDELS_in_final_protein,start_codon+mrna_sequence[3:],protein_sequence))

                    ## Collect all the information in one dictionary:  non_synonymous_proteins.
                    counter += 1
                    possible_protein_sequence = {}
                    variant_information = all_possible_sequences_for_this_transcript[possible_sequence]['variant info']
                    final_offset = all_possible_sequences_for_this_transcript[possible_sequence]['offset']
                    final_protein['variant info'] = variant_information
                    final_protein['offset'] = final_offset
                    final_protein['protein sequence'] = protein_sequence
                    final_protein['orf_structure']=orf_structure
                    possible_protein_sequence['possibility' + str(counter)] = final_protein
                    non_synonymous_proteins.update(possible_protein_sequence)

                    if stop_codon_in_INDEL=='YES':
                        print('Following sequence has the stop codon coordinate lying in the INDEL:')
                        #print(' mRNA sequence:')
                        #print(mrna_sequence)
                        #print(' Protein sequence:')
                        #print(protein_sequence)
                        print(' Genomic position stop codon:')
                        print(genomic_position_stop_codon)
                        print(' INDEL information:')
                        print(all_INDELS_in_final_protein)
                        print(' SNP information:')
                        print(all_non_synonymous_SNPs)
                        print(' Transcript information:')
                        print('transcript_start:',transcript_start,'tr_stable_id',tr_stable_id)
                        print('')
                        print('')

    # Close the tmp file.
    chr_tmp_file.close()

    print(" Finished translating chromosome "+str(chr))
    return


### Update TIS_overview ###
def update_TIS_overview(analysis_id,snp, path_to_sqlitedb,indel):
    '''
    :param analysis_id: a number representing the TIS ID that is being analysed
    :param snp: the SNP algorithm that is used: 'NO', 'samtools' or 'samtools_dbSNP'
    :param path_to_sqlitedb: path to the SQLite database
    :param include_INDELs: 'YES' if you want to include INDELs, 'NO' if you don't
    :return:
    '''
    try:
        #Make DB connection
        try:
            con = sqlite3.connect(path_to_sqlitedb)
        except:
            print "Could not connect to "+path_to_sqlitedb
            sys.exit()
        with con:
            cur = con.cursor()
            addColumn="ALTER table TIS_overview add indel varchar(20) default''"
            try:
                cur.execute(addColumn)
            except:
                pass
            update1="UPDATE TIS_overview set SNP = '"+snp+"' where ID = "+analysis_id;
            update2="UPDATE TIS_overview set indel = '"+indel+"' where ID = "+analysis_id;
            cur.execute(update1)
            cur.execute(update2)
    except:
        traceback.print_exc()
    return


################################# MAG WEG: #############################################



#print('Getting the analysis IDs that need to be processed')
#Get the analysis IDss that need to be processed
# tis_id="1,2"
# tis_id_list = []
# if(tis_id=='all'):
#     tis_id_list = get_all_tis_ids('SQLite/results.db')
# elif(tis_id.isdigit()):
#     tis_id_list.append(tis_id)
# elif(bool(re.search(',',tis_id))):
#     tis_id_list = re.split(',', tis_id)
# else:
#     print "Error: TIS IDs argument should be a number, a comma separated list of numbers or 'all'!"
#     sys.exit()
# #print(tis_id_list)


#######Set Main##################
if __name__ == "__main__":		#
    try:						#
        main()					#
    except Exception, e:		#
        traceback.print_exc()	#
#################################
