#!/usr/bin/env python3
import csv, sys, re
# process a vcf-to-tab result INDEL table file
if len(sys.argv) < 2:
    print("Expecting one filename argument")
    exit()

file = sys.argv[1]
ofile = "indel_matrix.fasaln"
strains = []
aln = {}
with open(file,"r") as tablefh:
    reader = csv.reader(tablefh,delimiter="\t")
    for indel in reader:
        ref = re.sub(r'\/','',indel[2])
        
        if len(strains) == 0:
            aln['REF'] = ""
            for n in range(3,len(indel),1):
                strains.append(indel[n])
                aln[indel[n]] = ""
        else:
            #print(indel[0],indel[1],indel[3],ref,alt)
            ref = indel[2]
            genolookup = {}
            genolookup[ref] = 0
            code = 1
            alts = []
            for allele in range(3,len(indel),1):
                al = re.sub(r'\/','',indel[allele])
                alts.append(al)
            
            alts = set(alts)
            for n in alts:
                if n != ref and n != '.':
                    genolookup[n] = code
                    code += 1
            i = 0
            aln['REF'] += str(0)
            for n in range(3,len(indel),1):
                allele = re.sub(r'\/','',indel[n])
                if allele.startswith('.'):
                    aln[strains[i]] += "-"
                elif allele not in genolookup:
#                    print(strains[i],allele,genolookup[allele])
#                else:
                    print("Cannot find '%s' in %s" % (allele,genolookup))
                    exit()
                else:
                    aln[strains[i]] += str(genolookup[allele])

                i += 1

#with open(ofile,"w") as fh:
#    fh.write("%5d %5d\n" % (len(strains),len(strains[0])))
#    for strain in aln:
#        fh.write("%-10s %s\n" % (strain,aln[strain]))

with open(ofile,"w") as fh:
    for strain in aln:
        fh.write(">%s\n%s\n" % (strain,aln[strain]))
