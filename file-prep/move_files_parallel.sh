#!/bin/bash

parallel -j 6 bash move_files.sh :::: sample_list
