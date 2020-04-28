#!/usr/bin/env python3

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

import os
import sys
import argparse

# Script for merging MAFs
# Creates one aggregate header taken from column names acorss MAFs
# Writes out rows while mapping to header
# NAs and equivalents replaced with ''

TUMOR_SAMPLE_BARCODE_COLUMN = "Tumor_Sample_Barcode"
HUGO_SYMBOL_COLUMN = "Hugo_Symbol"
ENTREZ_GENE_ID_COLUMN = "Entrez_Gene_Id"

NULL_OR_MISSING_VALUES = ["", "NA", "N/A", None]

def process_datum(val):
    """ Cleans up datum. """
    try:
        vfixed = val.strip()
    except AttributeError:
        vfixed = ""
    if vfixed in NULL_OR_MISSING_VALUES:
        return ""
    return vfixed

def get_file_header(data_file):
    """
        Returns header from MAF
        Assumes file is tab-delimited
        Assumes header is first non-commented line
    """
    header = []
    with open(data_file, "r") as header_source:
        for line in header_source:
            if not line.startswith("#"):
                header = list(map(str.strip, line.rstrip().split("\t")))
                break
    if not header:
        print("Could not extract header from mutation data file: %s - Exiting..." % (filename))
        sys.exit(2)
    return header

def merge_headers(data_filenames):
    """
        Generates a merged header from all input files.
        * Also ensures that "Hugo_Symbol" and "Entrez_Gene_Id"
          are the first 2 columns in the merged MAF header.
    """
    merged_header = [HUGO_SYMBOL_COLUMN, ENTREZ_GENE_ID_COLUMN]
    for fname in data_filenames:
        columns_to_add = [column for column in get_file_header(fname) if column not in merged_header]
        merged_header.extend(columns_to_add)
    return merged_header

def merge_input_mafs(input_mafs, output_maf_filename):
    """ Generates a merged MAF given a list of input MAFs. """
    merged_header = merge_headers(input_mafs)

    rows_to_write = []
    for input_maf in input_mafs:
        print('Loading data from ' + input_maf)
        header_processed = False
        file_header = get_file_header(input_maf)
        with open(input_maf, "r") as maf:
            for line in maf:
                if line.startswith("#"):
                    continue
                if not header_processed:
                    header_processed = True
                    continue
                # map row values to current header columns
                mapped_row = dict(zip(file_header, list(map(lambda x: process_datum(x), line.split("\t")))))
                sample_id = mapped_row[TUMOR_SAMPLE_BARCODE_COLUMN]

                # full record for the merged MAF (if mapped_row does not contain a column in merged_header, it is blank)
                normalized_row = list(map(lambda x: mapped_row.get(x, ""), merged_header))
                rows_to_write.append("\t".join(normalized_row))

    with open(output_maf_filename, "w") as output_maf:
        output_maf.write("\t".join(merged_header))
        output_maf.write("\n")
        output_maf.write("\n".join(rows_to_write))
        output_maf.write("\n")
    print('Merged MAF written to: %s' % (output_maf_filename))

def verify_input_mafs(input_mafs):
    for maf in input_mafs:
        if not os.path.isfile(maf):
            print("Could not find MAF by filename: %s" % (maf))
            exit(1)

def usage(parser):
    parser.print_help()
    sys.exit(2)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input-mafs-list", action = "store", dest = "input_mafs_list", help = "comma-delimited list of MAFs to merge")
    parser.add_argument("-d", "--input-mafs-directory", action = "store", dest = "input_mafs_directory", help = "directory containing all MAFs to merge")
    parser.add_argument("-o", "--output-maf", action = "store", dest = "output_maf", help = "output filename for merged MAF [REQUIRED]")

    args = parser.parse_args()
    input_mafs_list = args.input_mafs_list
    input_mafs_directory = args.input_mafs_directory
    output_maf = args.output_maf

    if not output_maf:
        print("Missing required argument: -o | --output-maf")
        usage(parser)

    if (input_mafs_list and input_mafs_directory) or (not input_mafs_list and not input_mafs_directory):
        print("Please choose only one of the following options when running script:  --input-mafs-list | --input-mafs-directory")
        usage(parser)

    if input_mafs_list:
        input_mafs = list(map(str.strip, input_mafs_list.split(",")))
    else:
        input_mafs = []
        for maf in os.listdir(input_mafs_directory):
            full_maf_path = os.path.join(input_mafs_directory, maf)
            if os.path.isfile(full_maf_path):
                input_mafs.append(full_maf_path)
            else:
                print("Skipping non-file element in directory... '%s'" % (full_maf_path))

    verify_input_mafs(input_mafs)
    merge_input_mafs(input_mafs, output_maf)
