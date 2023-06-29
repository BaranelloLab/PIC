#!/bin/bash

python2 extract_genomic_locations.py --gtf ./hg38.ncbiRefSeq.sorted.gtf --TSSup 500 --TSSdown 500 --GBdown 1000 --gsize hg38.genome
