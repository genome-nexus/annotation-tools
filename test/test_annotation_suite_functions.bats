#!/usr/bin/env bats

# Test cases for annotation_suite_functions.sh
load '../annotation_suite_functions.sh'

@test "test_initAndCleanWorkingDirectory_create_dir_when_it_exists" {
    mkdir "test_dir"
    touch "test_dir/file1.txt"
    touch "test_dir/file2.txt"
    run initAndCleanWorkingDirectory "test_dir"
    [ "$status" -eq 0 ]
    [ "$(ls -A "test_dir")" = "" ]
    rm -rf "test_dir"
}

@test "test_initAndCleanWorkingDirectory_create_dir_non_existing" {
    run initAndCleanWorkingDirectory "non_existing_dir"
    [ "$status" -eq 0 ]
    [ -d "non_existing_dir" ]
    rm -rf "non_existing_dir"
}

@test "test_initAndCleanWorkingDirectory_create_dir_error" {
    # Creating a file with the same name as the directory
    touch "error_dir"
    run initAndCleanWorkingDirectory "error_dir"
    [ "$status" -eq 2 ]
    [[ "${lines[0]}" == "\n[ERROR], Error encountered while attempting to create directory: error_dir" ]]
    rm -rf "error_dir"
}

@test "test_check_file_existence_file_does_not_exist" {
    run check_file_existence "filename.txt"
    [ "$status" -eq 2 ]
    [[ "${lines[0]}" == "[ERROR], File: filename.txt doesn't exist" ]]
}

@test "test_check_file_existence_file_exists" {
    mkdir "test_output_dir"
    touch "test_output_dir/file_temp.txt"
    run check_file_existence "test_output_dir/file_temp.txt"
    [ "$status" -eq 0 ]
    [ -f "test_output_dir/file_temp.txt" ]
    rm -rf "test_output_dir"
}
