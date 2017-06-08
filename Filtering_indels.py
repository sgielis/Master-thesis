#!/usr/bin/python

import io
import re

file = open("all_variants.txt", "r")
with io.FileIO("indel_variants.txt", "w") as newfile:
    for line in file:
        if "INDEL" in line:
            if "DP4" in line:
                newfile.write(line)

# remove unnecessary variant information, but keep the read depth and allelic frequency

for line in open("indel_variants.txt",'r'):
  columns = line.strip().split('\t')
  #print columns[0], columns[1], columns[3], columns[4], columns[7]
  column7 = columns[7].strip().split(';')
  #print column7

  
  m=re.search("DP=([\d]+)",line)
  if m:
      found = m.groups(1)[0]
      #print("dp",found)
  n = re.search("AF1=([0-9.]+)", line)

  
  if n:
      found2 = n.groups(1)[0]
      #print("af1", found2)
  
  s = ";"
  seq = (columns[0], columns[1], columns[3], columns[4], found, found2)
  print s.join(seq)

