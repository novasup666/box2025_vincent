#! /usr/bin/env bash

set -e
set -o pipefail
set -u

curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/734/005/GCA_009734005.2_ASM973400v2/GCA_009734005.2_ASM973400v2_genomic.fna.gz \
	| seqtk seq  \
	> e_faecium.fa

