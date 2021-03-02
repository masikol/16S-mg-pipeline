#!/usr/bin/env bash

# Kernel script for metagenomic reads preprocessing pipeline
# Usage:
#   ./preprocessing-runner.sh INPUT_DIR METADATA_FILE OUTPUT_DIR

# These arguments are already verified in a script
#   that launches this script: the "preprocessing pipeline" one.
# So, we won't check these arguents once again.

# First argument -- input directory. Verify it.
input_dir=`realpath $1`
# Second argument -- metadata file. Verify it.
metadata_file=`realpath $2`
# Third argument -- output directory
outdir=`realpath $3`

# Check if all read file pairs are complete
echo -n "Primary validation..."
for f in ${input_dir}/*_R1*.fastq*; do # check forward-to-reverse mapping
  r="${f/_R1/_R2}"
  if [[ ! -f "${r}" ]]; then
    echo -e "\nError: no reverse (R2) file found for forward (R1) file \`${f}\`."
    exit 1
  fi
done
for r in ${input_dir}/*_R2*.fastq*; do # check reverse-to-forward mapping
  f="${r/_R2/_R1}"
  if [[ ! -f "${f}" ]]; then
    echo -e "\nError: no forward (R1) file found for reverse (R2) file \`${r}\`."
    exit 1
  fi
done
echo "ok"

# Define some variables
global_preprocess_outdir="${outdir}/without_primers"
global_preprocess_log="${global_preprocess_outdir}/preprocess16S_global.log"
curr_preprocess_outdir="${input_dir}/curr_outdir"
curr_preprocess_log="${curr_preprocess_outdir}/preprocess16S.log"

# Get and check "filter_merge_classify" script
filter_merge_classify=`realpath ./src/filter-merge-classify.R`
if [[ ! -f "${filter_merge_classify}" ]]; then
  echo -e "\aError: file \`filter-merge-classify.R\` does not exist at \`${filter_merge_classify}\`!"
  exit 1
fi


# Run preprocess16S over all pairs of reads
echo -e "\nRunning preprocess16S: it will remove crosstalks and trim primers."
mkdir -pv "${global_preprocess_outdir}"

for f in ${input_dir}/*_R1*.fastq*; do
  r="${f/_R1/_R2}"

  # Run preprocess16S
  "${PREPROCESS16S}" -1 "${f}" -2 "${r}" -x 0.7 -o "${curr_preprocess_outdir}"

  # If eror occured -- exit
  if [[ $? != 0 ]]; then
    echo "Error in preprocess16S. Exiting now."
    exit 1
  fi

  # Move valid files to 'outdir' and save current log
  for valid_file in ${curr_preprocess_outdir}/*.16S.fastq*; do
    mv -v ${valid_file} ${global_preprocess_outdir}
    echo -en "\n" >> "${global_preprocess_log}"
    cat "${curr_preprocess_log}" >> "${global_preprocess_log}"
  done

  # Remove temporary directory
  rm -rf "${curr_preprocess_outdir}"
  if [[ $? != 0 ]]; then
    echo -e "\nError: cannot remove temporary directory."
    echo "It is strange and might be fatal. Exitting now."
    exit 1
  fi
done


# Run preprocessing R script: filter, calculate error models, merge, remove chimeras, classify and create phyloseq object
echo -e "\nStarting R preprocessing script. It will filter reads, merge and classify them."
echo -e "What is important, it will create phyloseq object (\`phyloseq_object.rds\`), which you need to import to GenePiper for further analyses.\n"
Rscript --vanilla "${filter_merge_classify}" "${global_preprocess_outdir}" "${metadata_file}"
if [[ $? != 0 ]]; then
  echo "R preprocessing script failed! Exitting now."
  exit 1
fi

exit 0

