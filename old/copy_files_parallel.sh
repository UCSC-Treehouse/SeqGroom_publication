#!/bin/bash

parallel -j 6 bash copy_files.sh :::: sample_list
