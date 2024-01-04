#!/bin/bash

# Contains helper functions used in annotation_suite_wrapper.sh

# Function creates a directory, clears it if it already exists
# Arguments:
#   $1: data_dir - directory to be cleaned / initialized
function initAndCleanWorkingDirectory {
    data_dir=$1
    if [ -d "${data_dir}" ]; then
        rm -f "${data_dir}"/*
    else
        mkdir "${data_dir}"
        if [ $? -gt 0 ]; then
            echo -e "\n[ERROR], Error encountered while attempting to create directory: ${data_dir}"
            exit 2
        fi
    fi
}

# Function checks if a filepath exists and is a file
# Arguments:
#   $1: path_to_file - filepath of file to be checked
function check_file_existence {
    path_to_file=$1 # Get the file path from the function argument
    if [ -f "${path_to_file}" ]; then
        echo "File: ${path_to_file} exists"
    else
        echo -e "\n[ERROR], File: ${path_to_file} doesn't exist"
        exit 2
    fi
}
