#!/usr/bin/perl -w
#
# Extract mismatches (SNPs and INDELs), aka alternative bases, and their corresponding reference bases from mapped reads in a SAM file.
# Mismatches are found based on the CIGAR string and MD tag (for more info: http://samtools.sourceforge.net/SAMv1.pdf):
#   - does the MD tag contain an A, T, C or G?
#       yes: there are one or more mismatches
#       no: no mismatches, go to the next read
#   - find the position of the mismatch in the read based on the numbers of matching bases that precede the mismatch
#     using the number in the MD tag right before the A, T, C or G
#   - check the CIGAR string for clipped bases (indicated with an S) and calculate the genomic position of the mismatch
# Input = SAM file
# Output = mismatches.txt file with the following columns:
#   chromosome name | position | reference base | alternative base(s)
#
# Remark: the output will contain duplicates, so don't forget to remove those (eg sort mismatches.txt | uniq > unique_mismatches.txt).
#

use strict;


my $startTime = time;

my $usage = "\nUsage: $0 <infile.sam>\n\n";
my $samFile = shift or die $usage;

open(SAM, "<$samFile") or die "Couldn't open the sam file $samFile\n".$!."\n";
open(SNP, ">>mismatches_indel.txt") or die "Couldn't open mismatch file\n".$!."\n";

print "Reading sam file...\n";

while (<SAM>){
    
    chomp;
    my $samFileLine = $_;
    
    if ($samFileLine !~ /^@/){
        
        my @lineArray = split('\t', $samFileLine);
        # the columns of interest are:
        # column 3: chromosome
        # column 4: position
        # column 6: CIGAR string
        # column 10: read sequence
        # column 19: MD tag
        my $chr = $lineArray[2];
        my $pos = $lineArray[3];
        my $cigar = $lineArray[5];
        my $read = $lineArray[9];
        my $md = $lineArray[-1]; # the number of columns depends on whether duplicate reads were removed with picard or not, but the MD tag will be the last column either way
        $md =~ s/MD:Z://;
        
        my $validCigar = 1;
        if ($md =~ /[ATCG]/ or $cigar =~ /[ID]/){ # is there a SNP in the MD tag? 
            
            # assign a genomic reference position to each base of the mapped read
            #my @genomicPositions;
            my $readPos = $pos;
            my $cigarCopy = $cigar;
            my $counter = 0;
            my $offset = 0;
            my $number_of_D=0;
            while ($cigarCopy !~ /^$/){
                if ($cigarCopy =~ /^([0-9]+[MNIDS])/){
                    my $cigar_part = $1;
                    #print SNP "$cigar_part\n";
                    if ($cigar_part =~ /(\d+)M/){
                        for (my $t = 0; $t < $1; $t++){
                            #push(@genomicPositions, $readPos);
                            $readPos++;
                            $counter++;
                        }
                    # this part will have to deal with INDELS:
                    } elsif ($cigar_part =~ /(\d+)I/){
                        my $length_indel=$1;
                        my $genomic_pos=$readPos-1;
                        
                         
                        my $indelpos=$counter-1 + $offset;
                        my $alt_indel= substr $read, $indelpos,$length_indel+1;
                        my $ref_indel= substr $read, $indelpos,1;
                        my $begin= substr $read, 1, 2;
                       # my $readPos= $readPos+$length_indel;
                        $offset=$offset+$length_indel;
                                               

                        print SNP "$chr;$genomic_pos;$ref_indel;$alt_indel\n";
                       
                    } elsif ($cigar_part =~ /(\d+)D/){
                        my $length_indel=$1;
                        #print SNP "$length_indel\n";
                       
                       $number_of_D+=1;
                       #print SNP "$md\n";
                       # my @dig;
                       # @dig=split /\^/,$md;
                       # print SNP "@dig\n";
                       my $str = $md;
                       my @fields = split /\^/, $str;
                       #print SNP "@fields\n";
                       
                       my $md_part = $fields[$number_of_D];
                       my $ref=substr $md_part,0,$length_indel;
                       #print SNP "$readPos\n";
                       my $alternative_indel=substr $read,$counter-1+$offset,1;
                       my $ref_indel= (substr $read, $counter-1+$offset,1).$ref;
                       my $genomic_pos=$readPos-1;

                       print SNP "$chr;$genomic_pos;$ref_indel;$alternative_indel\n";
                       $readPos += $1;
                      
                      
                    } elsif ($cigar_part =~ /(\d+)S/){
                        for (my $t = 0; $t < $1; $t++){
                            $readPos++;
                            $counter++;
                        }
                    } elsif ($cigar_part =~ /(\d+)N/){
                        $readPos += $1;
                        #$counter+=$1;
                    }
                    $cigarCopy =~ s/$cigar_part//;
                } else {
                    $validCigar = 0;
                    last;
                }
            }
                  
        }
        
    }
    
}

close SAM;
close SNP;

my $runTime = time - $startTime;
printf("runtime: %02d:%02d:%02d", $runTime/3600, ($runTime % 3600)/60, ($runTime % 3600) % 60);
print "\n";