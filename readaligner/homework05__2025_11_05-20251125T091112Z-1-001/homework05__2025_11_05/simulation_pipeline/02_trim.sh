#! /usr/bin/env bash

set -e
set -o pipefail
set -u

(
head -c 50000 e_faecium.fa
echo
tail -n 2 e_faecium.fa | head -c 50000
) > ref.fa
