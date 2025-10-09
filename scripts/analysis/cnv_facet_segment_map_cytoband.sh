#!/bin/bash

## Download cytoband annotation file from UCSC
wget https://hgdownload.soe.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz && gunzip cytoBand.txt.gz
sort -k1,1 -k2,2n cytoBand.txt > cytoBand.sorted.txt

sort -k1,1 -k2,2n cnv_facet_segment.bed > cnv_facet_segment.sorted.bed

## 50% overlap to assign cytoband, stric, single band focus
bedtools map -a cnv_facet_segment.sorted.bed -b cytoBand.sorted.txt -c 4 -o collapse -f 0.3 > annotated_segments.bed