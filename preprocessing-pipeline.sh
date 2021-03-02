#!/usr/bin/env bash

# Wrapper script for metagenomic reads preprocessing pipeline
# Usage:
#   ./preprocessing-pipeline.sh INPUT_DIR METADATA_FILE

# First argument -- input directory. Verify it.
input_dir=$1
if [[ -z "${input_dir}" ]]; then
  echo -e "\aError: you must specify input directory as first argument."
  exit 1
fi
if [[ ! -d "${input_dir}" ]]; then
  echo -e "\aError: input directory \`${input_dir}\` does not exist!"
  exit 1
fi
input_dir=`realpath "${input_dir}"`

# Second argument -- metadata file. Verify it.
metadata_file=$2
if [[ -z "${metadata_file}" ]]; then
  echo -e "\aError: you must specify metadata file as second argument."
  exit 1
fi
if [[ ! -f "${metadata_file}" ]]; then
  echo -e "\aError: metadata file \`${metadata_file}\` does not exist!"
  exit 1
fi
metadata_file=`realpath "${metadata_file}"`

# Check if SampleID column is present in metadata file
header=`head -n 1 "${metadata_file}"`
if [[ -z `echo ",${header}," | grep ',SampleID,'` ]]; then
  echo -e "\aError: no \`SampleID\` column in metadata file detected!"
  exit 1
fi

# Load configurations
config_file=`realpath ./src/preprocessing-pipeline.conf`
if [[ ! -f "${config_file}" ]]; then
  echo -e "\aError: cannot find configuration file at \`${config_file}\`."
  echo "Please, create it and configure."
  exit 1
fi
source "${config_file}"

if [[ -z "${PREPROCESS16S}" ]]; then
  echo -e "\aError: path to preprocess16S.py is not specified in cofiguration file \`${config_file}\`"
  echo 'Please, specify it.'
  exit 1
fi
if [[ ! -f "${PREPROCESS16S}" ]]; then
  echo -e "\aError: cannot find preprocess16S.py executable at \`${PREPROCESS16S}\`!".
  echo "Please, make sure that configurations in file \`${config_file}\` are correct."
  exit 1
fi

# Check if reference file for taxonomic classification exists
if [[ -z "${REF_CLASSIF_PATH}" ]]; then
  echo -e "\aError: path to reference file for taxonomic classification is not specified in cofiguration file \`${config_file}\`"
  echo 'Please, specify it.'
  exit 1
fi
if [[ ! -f "${REF_CLASSIF_PATH}" ]]; then
  echo -e "\aError: reference file for taxonomic classification does not exist: \`${REF_CLASSIF_PATH}\`!".
  echo "Please, make sure that configurations in file \`${config_file}\` are correct."
  echo 'You can download it with following commands (if you are connected to the Internet):'
  echo '  cd /home/justme/Metagenomics/Silva-training-set/'
  echo '  wget https://zenodo.org/record/3986799/files/silva_nr99_v138_wSpecies_train_set.fa.gz'
  exit 1
fi


# Check log directory in config file
if [[ -z "${LOG_DIR}" ]]; then
  echo -e "\aError: path to log directory is not specified in cofiguration file \`${config_file}\`"
  echo 'Please, specify it.'
  exit 1
fi


# Configure outdir path
outdir=`realpath "${input_dir}/preprocessed"`
if [[ -d "${outdir}" ]]; then
  echo "Error: output directory exists: \`${outdir}\`"
  echo "Please, rename it or remove."
  exit 1
fi


# Create log file
if [[ ! -d "${LOG_DIR}" ]]; then
  mkdir -vp "${LOG_DIR}"
  if [[ $? != 0 ]]; then
    echo -e "Error: cannot create log directory \`${LOG_DIR}\`."
    exit 1
  fi
fi
# Get number of files in log dir
# There are always one extra line in output of `ls -lA`. E.g. "total 0"
n_files_in_LOG_DIR=`ls -lA "${LOG_DIR}" | wc -l`
log_file="${LOG_DIR}/preprocessing-pipeline-${n_files_in_LOG_DIR}.log"
touch "${log_file}"
if [[ $? != 0 ]]; then
  echo -e "\n\aWarning: cannot create log file \`${log_file}\`!"
  echo "It is not a fatal error. We just won't save log. Pretty sad situation, though.\n"
  log_file='/dev/null'
fi

# Get and check "runner" script
runner=`realpath ./src/preprocessing-runner.sh`
if [[ ! -f "${runner}" ]]; then
  echo -e "\aError: file \`preprocessing-runner.sh\` does not exist at \`${runner}\`!"
  exit 1
fi


# Print some stuff
echo "Starting preprocessing pipeline."
echo "Input directory: \`${input_dir}\`."
echo "Log file: \`${log_file}\`."
echo "Metadata file: \`${metadata_file}\`."
echo -e "Output directory: \`${outdir}\`\n"

# Launch runner and write it's output to log file
bash "${runner}" "${input_dir}" "${metadata_file}" "${outdir}" |& tee "${log_file}"


global_preprocess_log="${outdir}/without_primers/preprocess16S_global.log"
if [[ -f "${global_preprocess_log}" ]]; then
  echo -e "\npreprocess16S logs go below:\n" >> "${log_file}"
  cat "${global_preprocess_log}" >> "${log_file}"
  rm "${global_preprocess_log}"
else
  echo -e "\a\nWarning: cannot save logs from preprocess16S."
  echo "It's log file should be here: \`${global_preprocess_log}\`, but this file does not exist."
  echo "It is not a fatal error, this log just won't be saved."
fi


# Print some stuff again: it worth reminding
echo -e "\nInput directory: \`${input_dir}\`."
echo "Log file: \`${log_file}\`."
echo "Metadata file: \`${metadata_file}\`."
echo -e "Output directory: \`${outdir}\`\n"

exit 0
