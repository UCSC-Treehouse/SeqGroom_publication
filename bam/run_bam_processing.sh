#!/bin/bash

parallel -j 3 bash process_bam.sh :::: bam_sample_list
