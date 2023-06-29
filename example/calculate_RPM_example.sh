#!/bin/bash
PIC_path=

python2 $PIC_path/tag_RPM.py --chip $bam_path/example.markdup.sorted.bam --input $bam_path/blank_input.bam --min 30 --tss $PIC_path/gene_tss.bed --gb $PIC_path/gene_GB.bed

