#!/bin/bash

s=$1

og=/private/groups/treehouse/staging-unrestricted/SeqGroom/primary/original
new=/scratch/jmroger/processing-v2/bam

cp ${og}/${s}/${s}.bam ${new}/${s}/${s}.bam
