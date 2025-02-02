#!/bin/bash

total_buckets=${1:-2}
awk -F '\t' '{ print $2 }' design.csv | tail -n+2 > gff_listing.txt
python scheduling.py gff_listing.txt --total_buckets ${total_buckets}
