#!/bin/bash

# Script to run the Forecast Evaluation System Products.
# EDIT THE PATH BELOW:
export HOME=/home/dmelecio
export HOMEEVAL=/home/dmelecio/Evaluation_System/uwrf_evaluation_system
export ASOS=/home/dmelecio/Evaluation_System/uwrf_evaluation_system/obs_station_day_minus_0
export SCRIPTS=/home/dmelecio/Evaluation_System/uwrf_evaluation_system/scripts
source ${HOME}/.bashrc

#######################################################
# Time and date variables
#######################################################

# Date of interest (doi) (Operational)
yyyy_doi="`date -d "yesterday" '+%Y'`"
mm_doi="`date -d "yesterday" '+%m'`"
dd_doi="`date -d "yesterday" '+%d'`"
full_doi="`date -d "yesterday" '+%Y-%m-%d'` 00:00:00"

# Date of Interest (doi) (Testing and Debugging)
yyyy_doi="2020"
mm_doi="01"
dd_doi="02"
full_doi="${yyyy_doi}-${mm_doi}-${dd_doi} 00:00:00"

#######################################################
# Sed Changes
#######################################################

# Change to the home directory of the Eval. System.
cd ${HOMEEVAL}

# Changes in the dl_ny_asos.py script
rm cuerg_asos.py
sed -e '{
       s;YEARSTART;'`date -d "${yyyy_doi}-${mm_doi}-${dd_doi} - 1 day" '+%Y'`';
       s;MONSTART;'`date -d "${yyyy_doi}-${mm_doi}-${dd_doi} - 1 day" '+%m'`';
       s;DAYSTART;'`date -d "${yyyy_doi}-${mm_doi}-${dd_doi} - 1 day" '+%d'`';
       s;YEAREND;'`date -d "${yyyy_doi}-${mm_doi}-${dd_doi} + 1 day" '+%Y'`';
       s;MONEND;'`date -d "${yyyy_doi}-${mm_doi}-${dd_doi} + 1 day" '+%m'`';
       s;DAYEND;'`date -d "${yyyy_doi}-${mm_doi}-${dd_doi} + 1 day" '+%d'`';
       }'  ${ASOS}/dl_ny_asos.py > ./cuerg_asos.py

# Changes in Product01 file
rm ${HOMEEVAL}/cuerg_P01.Rmd
sed 's;CHANGE_DATE_HERE;"${full_doi}";' ${SCRIPTS}/Product01-Forecast-Hour-Based_Evaluation.md > ./cuerg_P01.Rmd
