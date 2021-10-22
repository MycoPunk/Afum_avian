#!/usr/bin/bash

pushd gvcf
ln -s /bigdata/stajichlab/shared/projects/Population_Genomics/A_fumigatus/align/*/Variants/*.gz .
ln -s /bigdata/stajichlab/shared/projects/Population_Genomics/A_fumigatus/align/*/Variants/*.gz.tbi .

for n in $(cat ../to_remove.txt)
do
	rm $n.g.vcf.gz $n.g.vcf.gz.tbi
done
popd
