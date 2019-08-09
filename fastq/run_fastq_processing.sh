#!/bin/bash

parallel -j 3 bash process_fastq.sh :::: fastq_sample_list
