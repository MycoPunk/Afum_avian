#!/usr/bin/env python3

import os, csv, sys, re
outdir="AF100_8_9_patterns"
if not os.path.exists(outdir):
    os.mkdir(outdir)
input="AF100_8_9.v4.snpEff.tab"
if len (sys.argv) > 1:
    input = sys.argv[1]

base=os.path.basename(input)
stem=os.path.splitext(base)
outfile = "%s.patterns.tsv" % (stem[0])
print("outfile=%s"%(outfile))

global_patterns = {'ref': {}, 'alt': {} }

header = []
with open(input,"r") as snptab:
    tabparser = csv.reader(snptab,delimiter="\t")
    header = next(tabparser)
    
    for i in range(len(header)):
        header[i] = re.sub(r'\[\d+\]','',header[i])
        header[i] = re.sub(r':GT','',header[i])
    total_cols = len(header)

    for row in tabparser:
        ref = row[2]
        pattern = {'ref': [], 'alt': [] }

        for col in range(4,total_cols-1):
            if row[col] is '.' or row[col] is '':
                next
            if row[col] == ref:
                pattern['ref'].append(header[col])
            else:
                pattern['alt'].append(header[col])

        for type in pattern.keys():
            if len(pattern[type]) > 0:
                strpat = ",".join(sorted(pattern[type]))            
                if strpat not in global_patterns[type]:
                    global_patterns[type][strpat] = []
                global_patterns[type][strpat].append(row)

with open(outfile,"w") as patout:
    csvout = csv.writer(patout,delimiter="\t")
    csvout.writerow(["PATTERN","COUNT"])
    for type in ['alt']:
        for pat in global_patterns[type].keys():
            if "1F1SW" not in pat and "AF100-3-A" not in pat and "AFIS_2101" not in pat and "HMR_AF706" not in pat:                
                csvout.writerow([pat,len(global_patterns[type][pat])])
                with open(os.path.join(outdir,"%s.tsv"%(pat)),"w") as outpatlines:
                    outpatcsv = csv.writer(outpatlines,delimiter="\t")
                    outpatcsv.writerow(header)
                    for line in global_patterns[type][pat]:
                        outpatcsv.writerow(line)
            
