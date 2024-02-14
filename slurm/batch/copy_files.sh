#!/bin/bash

G4HUNTER=/storage/group/izg5139/default/akshatha/gquad/raw_data/g4hunter
rm -rf $G4HUNTER
mkdir $G4HUNTER
cd $G4HUNTER

rsync -av /storage/group/izg5139/default/external/quadrupia_database/g4/g4db_everything/archaea_prettified/genomic/ .
rsync -av /storage/group/izg5139/default/external/quadrupia_database/g4/g4db_everything/bacteria_prettified/genomic/ .
rsync -av /storage/group/izg5139/default/external/quadrupia_database/g4/g4db_everything/eukaryota_prettified/genomic/ .
rsync -av /storage/group/izg5139/default/external/quadrupia_database/g4/g4db_everything/viral_prettified/genomic/ .

REGEX=/storage/group/izg5139/default/akshatha/gquad/raw_data/regex
rm -rf $REGEX
mkdir $REGEX
cd $REGEX

rsync -av /storage/group/izg5139/default/external/quadrupia_database/regex/regexdb_everything/archaea/genomic/ .
rsync -av /storage/group/izg5139/default/external/quadrupia_database/regex/regexdb_everything/bacteria/genomic/ .
rsync -av /storage/group/izg5139/default/external/quadrupia_database/regex/regexdb_everything/eukaryota/*/genomic/ . 
rsync -av /storage/group/izg5139/default/external/quadrupia_database/regex/regexdb_everything/viral/genomic/ .
