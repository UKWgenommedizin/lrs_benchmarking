#!/bin/bash

grep 'Real' ../cram/mm2-*/*.hg38.mm2*.cram.log | grep -v '1k' | sed 's/_/\t/g' | sed 's/\.hg38\./\t/g' | sed 's/.*HG0/HG0/g' | sed 's/\.cram.*Real/\tReal/g' > mm2.run_metrics.tsv

grep 'Run Time:\|CPU Time:\|Peak RSS:' ../cram/pbmm2-*/*.hg38.pbmm2*.cram.log
