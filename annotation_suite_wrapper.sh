#!/bin/bash

# Copyright (c) 2020 Memorial Sloan Kettering Cancer Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# verify python3 and java installations
command -v python3 >/dev/null 2>&1 || { echo "python3 is required to run this program - aborting..." >&2; exit 1; }
command -v java >/dev/null 2>&1 || { echo "java is required to run this program - aborting..." >&2; exit 1; }

function usage {
    echo "annotation_suite_wrapper.sh"
    echo -e "\t-i | --input-directory               input data directory for processing mutation data files [REQUIRED]"
    echo -e "\t-o | --output-directory              output directory to write processed and annotated mutation data files to [REQUIRED]"
    echo -e "\t-m | --merged-mutation-file          path to write the merged mutation file for the center [REQUIRED]"
    echo -e "\t-c | --center-name                   name of the center being processed [REQUIRED]"
    echo -e "\t-s | --sequence-source               name of the sequence source used by the center (i.e., WXS, WGS) [REQUIRED]"
    echo -e "\t-p | --annotation-scripts-home       path to the annotation suite scripts directory [REQUIRED]"
}

# parse input arguments
for i in "$@"; do
case $i in
    -i=*|--input-directory=*)
    INPUT_DATA_DIRECTORY="${i#*=}"
    echo -e "\tINPUT_DATA_DIRECTORY=${INPUT_DATA_DIRECTORY}"
    shift
    ;;
    -o=*|--output-directory=*)
    OUTPUT_DATA_DIRECTORY="${i#*=}"
    echo -e "\tOUTPUT_DATA_DIRECTORY=${OUTPUT_DATA_DIRECTORY}"
    shift
    ;;
    -m=*|--merged-mutation-file=*)
    MERGED_MUTATION_FILENAME="${i#*=}"
    echo -e "\tMERGED_MUTATION_FILENAME=${MERGED_MUTATION_FILENAME}"
    shift
    ;;
    -c=*|--center-name=*)
    CENTER_NAME="${i#*=}"
    echo -e "\tCENTER_NAME=${CENTER_NAME}"
    shift
    ;;
    -s=*|--sequence-source=*)
    SEQUENCE_SOURCE="${i#*=}"
    echo -e "\tSEQUENCE_SOURCE=${SEQUENCE_SOURCE}"
    shift
    ;;
    -p=*|--annotation-scripts-home=*)
    ANNOTATION_SUITE_SCRIPTS_HOME="${i#*=}"
    echo -e "\tANNOTATION_SUITE_SCRIPTS_HOME=${ANNOTATION_SUITE_SCRIPTS_HOME}"
    shift
    ;;
    *)
    ;;
esac
done

if [[ -z "${INPUT_DATA_DIRECTORY}" || -z "${OUTPUT_DATA_DIRECTORY}" || -z "${MERGED_MUTATION_FILENAME}" || -z "${CENTER_NAME}" || -z "${SEQUENCE_SOURCE}" || -z "${ANNOTATION_SUITE_SCRIPTS_HOME}" ]]; then
    usage
    exit 1
fi

# check for presence of SSL cert in ${ANNOTATION_SUITE_SCRIPTS_HOME}
JAVA_SSL_ARGS=""
SSL_CERT_PATH=${ANNOTATION_SUITE_SCRIPTS_HOME}/AwsSsl.truststore
if ! [ -f ${SSL_CERT_PATH} ] ; then
    echo "Could not find SSL certificate: ${SSL_CERT_PATH} - please make sure this certificate exists and is present in ${ANNOTATION_SUITE_SCRIPTS_HOME}. Exiting..."
    exit 1
else
    JAVA_SSL_ARGS="-Djavax.net.ssl.trustStore=${SSL_CERT_PATH}"
fi


PROCESSED_SUB_DIR_NAME="${OUTPUT_DATA_DIRECTORY}/processed"
ANNOTATED_SUB_DIR_NAME="${OUTPUT_DATA_DIRECTORY}/annotated"
FILE_EXTENSIONS_LIST="vcf,maf,txt" # text files are treated as MAFs to handle names like data_mutations_extended.txt

STANDARDIZE_MUTATION_DATA_SCRIPT=${ANNOTATION_SUITE_SCRIPTS_HOME}/standardize_mutation_data.py
GENOME_NEXUS_ANNOTATOR_JAR=${ANNOTATION_SUITE_SCRIPTS_HOME}/annotator.jar
MERGE_MAFS_SCRIPT=${ANNOTATION_SUITE_SCRIPTS_HOME}/merge_mafs.py

GENOME_NEXUS_ANNOTATOR_ISOFORM="uniprot"
GENOME_NEXUS_ANNOTATOR_POST_SIZE=1000

function initAndCleanWorkingDirectory {
    # (1): directory to be cleaned / initialized
    data_dir=$1
    if [ -d "${data_dir}" ] ; then
        rm -f "${data_dir}"/*
    else
        mkdir "${data_dir}"
        if [ $? -gt 0 ] ; then
            echo -e "\n[ERROR], Error encountered while attempting to create directory: ${data_dir}"
            exit 2
        fi
    fi
}

function standardizeMutationFilesFromDirectory {
    # processed files will be written to ${OUTPUT_DATA_DIRECTORY}/processed
    echo -e "\t[INFO] standardizeMutationFilesFromDirectory(), standardized mutation files from ${INPUT_DATA_DIRECTORY} will be written to ${PROCESSED_SUB_DIR_NAME}"
    python3 ${STANDARDIZE_MUTATION_DATA_SCRIPT} --input-directory "${INPUT_DATA_DIRECTORY}" --output-directory "${PROCESSED_SUB_DIR_NAME}" --center ${CENTER_NAME} --sequence-source ${SEQUENCE_SOURCE} --extensions ${FILE_EXTENSIONS_LIST}
    if [ $? -gt 0 ] ; then
        echo -e "\n[ERROR] standardizeMutationFilesFromDirectory(), error encountered while running ${STANDARDIZE_MUTATION_DATA_SCRIPT}"
        exit 1
    fi
}

function annotateMAF {
    # (1): input file to annotate
    input_file="$1"
    output_file=${ANNOTATED_SUB_DIR_NAME}/$(basename "${input_file}").annotated
    echo -e "\t[INFO] annotateMAF(), annotating MAF: ${input_file} --> ${output_file}"
    java -Xmx48g ${JAVA_SSL_ARGS} -jar ${GENOME_NEXUS_ANNOTATOR_JAR} --filename "${input_file}" --output-filename "${output_file}" --isoform-override ${GENOME_NEXUS_ANNOTATOR_ISOFORM} -p ${GENOME_NEXUS_ANNOTATOR_POST_SIZE} -r
    if [ $? -gt 0 ] ; then
        echo -e "\n[ERROR] annotateMAF(), error encountered while running the genome nexus annotation pipeline"
        exit 1
    fi
}

function annotateStandardizedMAFs {
    for f in "${PROCESSED_SUB_DIR_NAME}"/* ; do
        if [ -f "${f}" ] ; then
            annotateMAF "${f}"
        fi
    done
}

function mergeMAFsInDirectory {
    echo -e "\n[INFO] mergeMAFsInDirectory(), merging MAFs from directory ${ANNOTATED_SUB_DIR_NAME} into ${MERGED_MUTATION_FILENAME}"
    python3 ${MERGE_MAFS_SCRIPT} --input-mafs-directory "${ANNOTATED_SUB_DIR_NAME}" --output-maf "${MERGED_MUTATION_FILENAME}"
}

# init and cleanup working directories where processed (standardized) mutation data files are written to
# as well as where the annotated mutation files are written to
initAndCleanWorkingDirectory ${OUTPUT_DATA_DIRECTORY}
initAndCleanWorkingDirectory ${PROCESSED_SUB_DIR_NAME}
initAndCleanWorkingDirectory ${ANNOTATED_SUB_DIR_NAME}

# standardize mutation files from ${INPUT_DATA_DIRECTORY}
standardizeMutationFilesFromDirectory
# annotate MAFs in ${PROCESSED_SUB_DIR_NAME}
annotateStandardizedMAFs
# merge MAFs in ${ANNOTATED_SUB_DIR_NAME}
mergeMAFsInDirectory
