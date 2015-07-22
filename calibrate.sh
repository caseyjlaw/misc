#!/usr/bin/env bash
#
# claw, 14jun13
# script file for 'submit'. most params for the job controlled here
# except for:
# -- <dataarea>/batch_leanpipedt.py, which defines transient pipeline params
# -- code/submit.req, which defines job name and number of nodes

# useful variables
export WORK=/lustre/aoc/projects/fasttransients/14jun13
export CODE=/lustre/aoc/projects/fasttransients/code
export NAME=14jun13

# move to data
cp calibrate.sh ${WORK}/calibrate_${NAME}.sh  # for record keeping
cd ${WORK}

# Run distribute. Default is to do all DMs on one node.
${CODE}/batch_cal2.py 14A-425_sb29272092_1.56821.73800049769 14A-425_${NAME}.ms 1 1 &> calibrate_${NAME}.txt
