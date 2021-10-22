def get_vcf_header(file,compress):
    if compress:
        fh = gzip.open(file,"rt")
    else:
        fh = open(file,"r")
    for line in fh:
        if line.startswith("##"):
            continue
        elif line.startswith("#CHROM"):
            vcf_header = line.strip().split("\t")
            vcf_header[0] = vcf_header[0][1:] # strip leading #
            return vcf_header



vcf = None
if args.vcf.endswith(".gz"):
    vcf_cols = get_vcf_header(args.vcf,1)
    vcf = pd.read_csv(args.vcf,compression="gzip",
                      sep='\t', comment = '#',
                      low_memory = True,
                      header=None,
                      names=vcf_cols)
else:
    vcf_cols = get_vcf_header(args.vcf,0)
    print("cols are ",vcf_cols)
    vcf = pd.read_csv(args.vcf,
                      sep='\t', comment = '#',
                      low_memory = True,
                      header=None,
                      names=vcf_cols)

print(vcf.head())
