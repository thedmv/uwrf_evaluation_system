#!/bin/bash

# Script to run the Forecast Evaluation System Products.
# EDIT THE PATH BELOW:
export HOME=/home/dmelecio
export HOMEEVAL=${HOME}/Evaluation_System/uwrf_evaluation_system
export ASOS=${HOMEEVAL}/obs_station_day_minus_0
export WRFD0=${HOMEEVAL}/forecast_day_minus_0_cuerg
export WRFD1=${HOMEEVAL}/forecast_day_minus_1_cuerg
export WRFD2=${HOMEEVAL}/forecast_day_minus_2_cuerg
export SCRIPTS=${HOMEEVAL}/scripts
export PYTHON=${HOME}/miniconda2/bin/python
export RR=/usr/bin/Rscript
export WRFTSDIR=/data/CUERG_PROJECTS/New_York_Forecast/output_prev
source ${HOME}/.bashrc

# Clean up previous ASOS files
rm ${ASOS}/JFK*.txt
rm ${ASOS}/LGA*.txt
rm ${ASOS}/NYC*.txt

#######################################################
# Time and date variables
#######################################################

# Date of interest (doi) (Operational)
yyyy_doi="`date -d "yesterday" '+%Y'`"
mm_doi="`date -d "yesterday" '+%m'`"
dd_doi="`date -d "yesterday" '+%d'`"
full_doi="`date -d "yesterday" '+%Y-%m-%d'` 00:00:00"

# Date of Interest (doi) (Testing and Debugging)
#yyyy_doi="2020"
#mm_doi="01"
#dd_doi="02"
#full_doi="${yyyy_doi}-${mm_doi}-${dd_doi} 00:00:00"

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
rm cuerg_P01.Rmd
sed -e "{
       s;CHANGE_DATE_HERE;${full_doi};
       }" ${SCRIPTS}/Product01-Forecast-Hour-Based_Evaluation.md > ./cuerg_P01.Rmd

########################################################
# Download the ASOS data
########################################################

# Run the python script for downloading ASOS data
rm mystations
ln -sf $ASOS/mystations .
$PYTHON cuerg_asos.py

# The ASOS data file names all have a .txt suffix
mv *.txt ${ASOS}/.
########################################################
# Get WRF Data
########################################################

# Remove previously used WRF data
rm ${WRFD0}/*.TS
rm ${WRFD1}/*.TS
rm ${WRFD2}/*.TS

# WRF data for the forecast initialized the same day as the date of interest
ln -sf ${WRFTSDIR}/${yyyy_doi}/${mm_doi}/${dd_doi}/*.d03.TS ${WRFD0}/.

# WRF data for the forecast initialized the day before the date of interest
wrfd1_doiY="`date -d "${yyyy_doi}-${mm_doi}-${dd_doi} - 1 day" '+%Y'`"
wrfd1_doim="`date -d "${yyyy_doi}-${mm_doi}-${dd_doi} - 1 day" '+%m'`"
wrfd1_doid="`date -d "${yyyy_doi}-${mm_doi}-${dd_doi} - 1 day" '+%d'`"
ln -sf ${WRFTSDIR}/${wrfd1_doiY}/${wrfd1_doim}/${wrfd1_doid}/*.d03.TS ${WRFD1}/.

# WRF data for the forecast initialized two days before the date of interest
wrfd2_doiY="`date -d "${yyyy_doi}-${mm_doi}-${dd_doi} - 2 day" '+%Y'`"
wrfd2_doim="`date -d "${yyyy_doi}-${mm_doi}-${dd_doi} - 2 day" '+%m'`"
wrfd2_doid="`date -d "${yyyy_doi}-${mm_doi}-${dd_doi} - 2 day" '+%d'`"
ln -sf ${WRFTSDIR}/${wrfd2_doiY}/${wrfd2_doim}/${wrfd2_doid}/*.d03.TS ${WRFD2}/.

#########################################################
# Product 01
#########################################################
rm "${yyyy_doi}-${mm_doi}-${dd_doi}_P01.html"
$RR -e "rmarkdown::render(input = 'cuerg_P01.Rmd', output_format = 'html_document', output_file = '${yyyy_doi}-${mm_doi}-${dd_doi}_P01.html')"

# Move the CSV files into the same folder as the plots
rm ./${yyyy_doi}-${mm_doi}-${dd_doi}/*.csv
mv *.csv ./${yyyy_doi}-${mm_doi}-${dd_doi}/.
