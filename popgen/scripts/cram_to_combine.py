#!/usr/bin/env python3
import csv, sys, re, os

threads="2"
strains = {}
fmt="cram"
inext="cram"
outext="cram"
outdir="cram"
genome="genome/FungiDB-39_AfumigatusAf293_Genome.fasta"
indir = "aln"
file = "SRA_samples.csv"
with open(file,"r") as tablefh:
    reader = csv.reader(tablefh,delimiter=",")
    header = next(reader)
    for line in reader:
        runAcc = line[0]
        strain = line[1]
        strain = re.sub("\s+","_",strain)
        strain = re.sub("/","-",strain)
        alnfile = "%s.%s"%(os.path.join(indir,runAcc),inext)
        if strain in strains:
            strains[strain].append(alnfile)
        else:
            strains[strain] = [alnfile]

file="Novogene_samples.csv"
with open(file,"r") as tablefh:
    reader = csv.reader(tablefh,delimiter=",")
    header = next(reader)
    for line in reader:
        strain = line[0]
        strain = re.sub("\s+","_",strain)
        strain = re.sub("/","-",strain)
        alnfile = "%s.%s"%(os.path.join(indir,runAcc),inext)
        if strain in strains:
            strains[strain].append(alnfile)
        else:
            strains[strain] = [alnfile]

file="preasm_samples.csv"
with open(file,"r") as tablefh:
    reader = csv.reader(tablefh,delimiter=",")
    header = next(reader)
    for line in reader:
        strain = line[0]
        strain = re.sub("\s+","_",strain)
        strain = re.sub("/","-",strain)
        alnfile = "%s.%s"%(os.path.join(indir,strain),inext)
        if strain in strains:
            strains[strain].append(alnfile)
        else:
            strains[strain] = [alnfile]

file="Af100_samples.csv"
with open(file,"r") as tablefh:
    reader = csv.reader(tablefh,delimiter=",")
    for line in reader:
        strain = line[0]
        runAcc = line[0]
        strain = re.sub("\s+","_",strain)
        strain = re.sub("/","-",strain)
        alnfile = "%s.%s"%(os.path.join(indir,strain),inext)
        if strain in strains:
            strains[strain].append(alnfile)
        else:
            strains[strain] = [alnfile]

for stdfile in [ "UCSF_201902-03_samples.csv", "UK_Imperial_samples.csv", "AF100_DMC_201912.csv", "COH_round2_samples.csv"]:
    with open(stdfile,"r") as tablefh:
        reader = csv.reader(tablefh,delimiter=",")
        header = next(reader)
        for line in reader:
            strain = line[0]
            runAcc = line[0]
            strain = re.sub("\s+","_",strain)
            strain = re.sub("/","-",strain)
            alnfile = "%s.%s"%(os.path.join(indir,runAcc),inext)
            if strain in strains:
                strains[strain].append(alnfile)
            else:
                strains[strain] = [alnfile]

file="NZ_samples.csv"
with open(file,"r") as tablefh:
    reader = csv.reader(tablefh,delimiter=",")
    next(reader)
    for line in reader:
        strain = line[0]
        strain = re.sub("\s+","_",strain)
        strain = re.sub("/","-",strain)
        alnfile = "%s.%s"%(os.path.join(indir,strain),inext)
        if strain in strains:
            strains[strain].append(alnfile)
        else:
            strains[strain] = [alnfile]

for strain in strains:
    outfile = "%s.%s"%(os.path.join(outdir,strain),outext)
    if not os.path.exists(outfile):
        if len(strains[strain]) == 1:
            print("ln -s %s %s" %(os.path.join("..",strains[strain][0]),
                                  outfile))
            print("ln -s %s.crai %s.crai" %(os.path.join("..",
                                                         strains[strain][0]),
                                            outfile))

        else:
            print("samtools merge -O %s --reference %s --threads %s %s %s"
                  %(fmt,genome,threads,outfile,
                    " ".join(strains[strain])))
