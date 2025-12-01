#!/usr/bin/env bash

base_location="/home/gdc/public/prs_methods/data/test/sim_1"
base_location="/home/gdc/public/prs_methods/data/test/sim_2"
anc1_prefix="AFR_simulation"
anc2_prefix="EUR_simulation"


mkdir -p "${base_location}/logs"
mv "${base_location}"/*.log "${base_location}/logs"

mkdir -p "${base_location}/summary_statistics"
mv "${base_location}"/*.txt "${base_location}/summary_statistics"
mv "${base_location}"/*sumstats* "${base_location}/summary_statistics"
mv "${base_location}"/ancestry* "${base_location}/summary_statistics"

mkdir -p "${base_location}/target_population"
mv "${base_location}/${anc1_prefix}"* "${base_location}/target_population"

mkdir -p "${base_location}/training_population"
mv "${base_location}/${anc2_prefix}"* "${base_location}/training_population"
