#! /usr/bin/env bash

set -e
set -o pipefail
set -u

cat ref.fa \
	| gzip \
	> output/ref.fa.gz

cp reads.bfast.fastq.gz output/reads.fq.gz

