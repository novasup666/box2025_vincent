#! /usr/bin/env bash

set -e
set -o pipefail
set -u

dwgsim -C 5 \
	-1 150 -2 0 \
	-e 0.01 -E 0.01 \
	-r 0.01 \
	ref.fa reads
