#!/bin/bash

find $1 -maxdepth 1 -type f -name "*.gff" -exec realpath {} \; > gff_file_list.txt
