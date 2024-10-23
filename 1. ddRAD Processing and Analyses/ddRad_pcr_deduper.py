#!/usr/bin/env python2

import itertools
import gzip
import re
import sys
import os
from optparse import OptionParser


""" Parsing command line options """
parser = OptionParser(prog="PCRdups_remover", usage="%prog [options", version="%prog 2.0.0")
parser.add_option('-d', '--dbr', action="store", type="string", dest="dbr_motif", default="NNNNNNNNNNNNN",
                  help="DBR motif, e.g.: NNNNNNNNNNNNN, or 13 x N, note this will ignore one putative degenerate base at 5' end\n")
parser.add_option('-g', '--gDNAlen', action="store", type="int", dest="gDNA_len", default=30,
                  help="number of base pairs to compare after dbr+r2cutsite, def. 30\n")
parser.add_option('-1', '--r1', action="store", type="string", dest="f1",
                  help="name of the R1 file\n")
parser.add_option('-2', '--r2', action="store", type="string", dest="f2",
                  help="name of the R2 file\n")
parser.add_option('-o', '--outfolder', action="store", type="string", dest="outfolder", default="output",
                  help="name of the subfolder to write results, def. output\n")
(options, args) = parser.parse_args()


""" Translation of ambiguities to regular expressions for DBR comparison """
ambiguities = {"A":"A",
               "T":"T",
               "C":"C",
               "G":"G",
               "R":"[AG]",
               "Y":"[CT]",
               "K":"[GT]",
               "M":"[AC]",
               "S":"[CG]",
               "W":"[AT]",
               "B":"[CGT]",
               "D":"[AGT]",
               "H":"[ACT]",
               "V":"[ACG]",
               "N":"[ACGT]"}


""" Checking for mandatory options """
if not all([options.f1, options.f2]):
    print "\n\tMust include --r1 (-1) and --r2 (-2) options\n"
    sys.exit()


""" Generate regular expression from DBR motif """
dbr_motif = options.dbr_motif
dbr_len = len(options.dbr_motif)
dbr_regex = ""
dbr_motif = list(dbr_motif)
for nuc in dbr_motif:
    dbr_regex += ambiguities[nuc]


""" Pass other options to variables """
gDNA_len = options.gDNA_len
f1_basename = options.f1.split(".fq")[0]
f2_basename = options.f2.split(".fq")[0]


""" Load fastq files 4 lines at a time """
if ".gz" in options.f1: # if files are compressed
    f1 = gzip.open(options.f1, "rb")
    f2 = gzip.open(options.f2, "rb")
else: # if files are not compressed
    f1 = open(options.f1, "r")
    f2 = open(options.f2, "r")
k1 = itertools.izip(*[iter(f1)]*4)
k2 = itertools.izip(*[iter(f2)]*4)


""" Create output folder if not already created """
outfolder = "./" + options.outfolder + "/"
if not os.path.exists(outfolder):
    os.makedirs(outfolder)


""" Start processing of files """
print "\n\n\t==================================="
print     "\tRemover of PCR-duplicates for ddRAD"
print     "\t===================================\n"

print "\n\tProcessing "+ options.f1 + " and " + options.f2 + "..."

pcr_dups = {} # Dictionary to contain duplicated portions of reads
pcr_dups_min1 = {} 
pcr_dups_pls1 = {}
dbr = re.compile(dbr_regex)
seq_len = dbr_len + gDNA_len

r1_good = []
r2_good = []
r1_bad = []
r2_bad = []
r1_dups = []
r2_dups = []

while 1:
    try: r2 = list(k2.next())
    except StopIteration: break
    r1 = list(k1.next())
    seq = r2[1][1:dbr_len+1] #best match is bases 2-14
    seqout = r2[1][dbr_len-2:]#Ideally this should fully exclude the DBR region, but might include the cut site
    qualout = r2[3][dbr_len-2:]
    if dbr.match(seq): #""" Checking that read is at least 14bp long """
        seq_uniq = r2[1][1:seq_len+1]
        seq_uniq_plus1 = r2[1][2:seq_len+2]
        seq_uniq_minus1 = r2[1][0:seq_len]
        #seq_uniq_plus2 = r2[1][3:seq_len+3] # if decide to shift to longer search window
        if (seq_uniq in pcr_dups): 
        	pcr_dups[seq_uniq] += 1
        	pcr_dups_pls1[seq_uniq_plus1] = 1
        	pcr_dups_min1[seq_uniq_minus1] = 1
        	r2[1] = seqout
        	r2[3] = qualout
        	r1_dups.append("".join(r1))
        	r2_dups.append("".join(r2))
        elif (seq_uniq in pcr_dups_pls1):
        	pcr_dups[seq_uniq] = 1
        	pcr_dups_pls1[seq_uniq_plus1] = 1
        	pcr_dups_min1[seq_uniq_minus1] = 1
        	r2[1] = seqout
        	r2[3] = qualout
        	r1_dups.append("".join(r1))
        	r2_dups.append("".join(r2))
        elif (seq_uniq in pcr_dups_min1):
        	pcr_dups[seq_uniq] = 1
        	pcr_dups_pls1[seq_uniq_plus1] = 1
        	pcr_dups_min1[seq_uniq_minus1] = 1
        	r2[1] = seqout
        	r2[3] = qualout
        	r1_dups.append("".join(r1))
        	r2_dups.append("".join(r2))
        else:
        	pcr_dups[seq_uniq] = 1
        	pcr_dups_pls1[seq_uniq_plus1] = 1
        	pcr_dups_min1[seq_uniq_minus1] = 1
        	r2[1] = seqout
        	r2[3] = qualout
        	r2_good.append("".join(r2))
        	r1_good.append("".join(r1))
        continue
    else:
        r1_bad.append("".join(r1))
        r2_bad.append("".join(r2))
    continue
        

""" Create output files """
out_r1_good = open(outfolder+f1_basename+".fq", "wb")
out_r1_good.write("".join(r1_good))
out_r1_good.close()

out_r2_good = open(outfolder+f2_basename+".fq", "wb")
out_r2_good.write("".join(r2_good))
out_r2_good.close()

out_r1_bad_dbr = gzip.open(outfolder+f1_basename+"_bad_DBRs.fq.gz", "wb")
out_r1_bad_dbr.write("".join(r1_bad))
out_r1_bad_dbr.close()

out_r2_bad_dbr = gzip.open(outfolder+f2_basename+"_bad_DBRs.fq.gz", "wb")
out_r2_bad_dbr.write("".join(r2_bad))
out_r2_bad_dbr.close()

out_r1_dups = gzip.open(outfolder+f1_basename+"_duplicates.fq.gz", "wb")
out_r1_dups.write("".join(r1_dups))
out_r1_dups.close()

out_r2_dups = gzip.open(outfolder+f2_basename+"_duplicates.fq.gz", "wb")
out_r2_dups.write("".join(r2_dups))
out_r2_dups.close()

out_reports = open(outfolder+f2_basename+"_report.txt", "w")
report = ['Accepted reads:\t'+str(len(r1_good))]
report.append('Duplicated reads:\t'+str(len(r1_dups)))
report.append('Length DBR+1:\t'+str(len(pcr_dups_pls1)))
report.append('Length DBR-1:\t'+str(len(pcr_dups_min1)))
report.append('Reads with DBR not matched:\t'+str(len(r1_bad)))
report.append('Note: this script optimizes removal and may slightly split/underestimate DBR counts')
report.append('\nSequence\tDuplicates')
pcr_dups = sorted(pcr_dups.items(), key=lambda x:x[1], reverse=True)
for dup in pcr_dups:
    if dup[1] > 1:
        report.append(dup[0] + "\t" + str(dup[1]))
out_reports.write("\n".join(report))
out_reports.close()

print "\tProcessing of "+ options.f1 + " and " + options.f2 + " complete"
