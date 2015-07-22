#!/usr/bin/env bash
#
# claw, 14jun13
# script file for 'submit'. most params for the job controlled here
# except for:
# -- <dataarea>/batch_leanpipedt.py, which defines transient pipeline params
# -- code/submit.req, which defines job name and number of nodes

args=("$@")

# useful variables
export NAME=${args[0]}    # directory name as in "14jun13"
export ASDM=${args[1]}    # asdmname given to distribute.pl etc
export PROJCODE='14A-425'
export WORK=/lustre/aoc/projects/fasttransients/${NAME}
export CODE=/lustre/aoc/projects/fasttransients/code

# prepare nodelist for distribute
cut -d ',' -f 1 cluster.conf > ${WORK}/nodelist.txt      # find allocated nodes
#cat ../code/nodelist_fixed.txt >> nodelist.txt    # add in fixed set of nodes

# move to data
cp search.sh ${WORK}/search_${NAME}.sh  # for record keeping
cd ${WORK}

# Run distribute. Default is to do all DMs on one node.
${CODE}/distribute.pl $ASDM ${PROJCODE}_${NAME}.ms nodelist.txt &> distribute_${NAME}.txt
