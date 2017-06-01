#!/usr/bin/python

#Haal alle indels uit het bestand snp.chr19.vcf en stop dit in een niev vcf bestand
#Je gebruikt 2 ifs omdat in de header ook indel staat en die wil je niet in je uiteindelijke vcf file

import io
import re

file = open("all_variants.txt", "r")
with io.FileIO("indel_variants.txt", "w") as newfile:
    for line in file:
        if "INDEL" in line:
            if "DP4" in line:
                newfile.write(line)
#Je krijgt nu een apart bestand met enkel de indels
#Hiervan willen we een tabel maken met volgende kolommen:
# id               chr         pos         ref         alt         type (indel) dp format
# remove unnecessary variant information, but keep the read depth and allelic frequency

for line in open("indel_variants.txt",'r'):
  columns = line.strip().split('\t')
  #print columns[0], columns[1], columns[3], columns[4], columns[7]
  column7 = columns[7].strip().split(';')
  #print column7

  #We willen de readt depht weten
  m=re.search("DP=([\d]+)",line)
  if m:
      found = m.groups(1)[0]
      #print("dp",found)
  n = re.search("AF1=([0-9.]+)", line)

  #We willen de allelic frequency weten
  if n:
      found2 = n.groups(1)[0]
      #print("af1", found2)
  #print columns[0], columns[1], columns[3], columns[4], found, found2
  #De resultaten van dit chromosoom staan dus nu in indelsnp.txt
  s = ";"
  seq = (columns[0], columns[1], columns[3], columns[4], found, found2)
  print s.join(seq)

