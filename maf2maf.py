#!/usr/bin/env python3

# Copyright (c) 2023 Memorial Sloan Kettering Cancer Center
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

# imports
import sys
import os
import click
import re
from chardet import detect

# constants
OUTPUT_FILE = sys.stdout
ERROR_FILE = sys.stderr

MAF_HEADER = [
	"Hugo_Symbol",
	"Entrez_Gene_Id",
	"Center",
	"NCBI_Build",
	"Chromosome",
	"Start_Position",
	"End_Position",
	"Strand",
	"Variant_Classification",
	"Variant_Type",
	"Reference_Allele",
	"Tumor_Seq_Allele1",
	"Tumor_Seq_Allele2",
	"dbSNP_RS",
	"dbSNP_Val_Status",
	"Tumor_Sample_Barcode",
	"Matched_Norm_Sample_Barcode",
	"Match_Norm_Seq_Allele1",
	"Match_Norm_Seq_Allele2",
	"Tumor_Validation_Allele1",
	"Tumor_Validation_Allele2",
	"Match_Norm_Validation_Allele1",
	"Match_Norm_Validation_Allele2",
	"Verification_Status",
	"Validation_Status",
	"Mutation_Status",
	"Sequencing_Phase",
	"Sequence_Source",
	"Validation_Method",
	"Score",
	"BAM_File",
	"Sequencer",
	"HGVSp_Short",
	"t_ref_count",
	"t_alt_count",
	"n_ref_count",
	"n_alt_count",
	"Protein_position",
	"Codons",
	"SWISSPROT",
	"RefSeq",
	"t_depth",
	"n_depth",
	"FILTER",
	"gnomAD_AF",
	"gnomAD_AFR_AF",
	"gnomAD_AMR_AF",
	"gnomAD_ASJ_AF",
	"gnomAD_EAS_AF",
	"gnomAD_FIN_AF",
	"gnomAD_NFE_AF",
	"gnomAD_OTH_AF",
	"gnomAD_SAS_AF"
]

# DEFAULTS
DEFAULT_NCBI_BUILD = "GRCh37"
DEFAULT_STRAND = "+"
DEFAULT_VALIDATION_STATUS = "Unknown"
DEFAULT_VERIFICATION_STATUS = "Unknown"
DEFAULT_MUTATION_STATUS = "Somatic"

# KEY COLUMN NAMES
TUMOR_SEQ_ALLELE1_COLUMNS = ["Tumor_Seq_Allele1", "TumorSeq_Allele1"]
TUMOR_SEQ_ALLELE2_COLUMNS = ["Tumor_Seq_Allele2", "TumorSeq_Allele2"]
MUTATED_FROM_ALLELE_COLUMN = "mutated_from_allele"
MUTATED_TO_ALLELE_COLUMN = "mutated_to_allele"
VARIANT_TYPE_COLUMNS = ["Variant_Type", "VariantType", "mut_type", "mutation_type"]
VARIANT_CLASSIFICATION_COLUMNS = ["Variant_Classification", "class", "Transcript architecture around variant"]
START_POSITION_COLUMNS = ["Start_Position", "Start_position", "POS", "chromosome_start"]
END_POSITION_COLUMNS = ["End_Position", "End_position", "chromosome_end"]
REFERENCE_ALLELE_COLUMNS = ["Reference_Allele", "reference_genome_allele"]
TUMOR_GENOTYPE_COLUMN = "tumour_genotype"
MATCHED_NORMAL_SAMPLE_BARCODE_COLUMNS = ["Match_Normal_Sample_Barcode", "Matched_Norm_Sample_Barcode", "matched_sample_id"]
CHROMOSOME_COLUMNS = ["Chromosome", "chromosome", "CHROM"]
VERIFICATION_STATUS_COLUMNS = ["Verification_Status", "verification_status"]
VALIDATION_METHOD_COLUMNS = ["Validation_Method", "verification_platform"]
MATCHED_NORMAL_SEQ_ALLELE1_COLUMNS = ["Match_Norm_Seq_Allele1", "mutated_to_allele"]
MATCHED_NORMAL_SEQ_ALLELE2_COLUMNS = ["Match_Norm_Seq_Allele2", "control_genotype"]
TUMOR_SAMPLE_BARCODE_COLUMNS = ["Tumor_Sample_Barcode", "analyzed_sample_id", "submitted_sample_id"]

# NULL VALUES
NULL_OR_MISSING_VALUES = set(["", "N/A", None, ".", "?"])

VALID_VARIANT_TYPES = ["SNP", "DNP", "TNP", "ONP", "DEL", "INS"]
MUTATION_FILTER = ["3'Flank", "3'UTR", "5'Flank", "5'UTR", "IGR", "Intron", "RNA", "Silent"]
MUTATION_STATUS_FILTER = ["LOH", "Wildtype"]

# GLOBALS FOR VALIDATING RESULTS
VALID_VARIANT_CLASSIFICATIONS = [
	"Frame_Shift_Del",
	"Frame_Shift_Ins",
	"In_Frame_Del",
	"In_Frame_Ins",
	"Missense_Mutation",
	"Nonsense_Mutation",
	"Silent",
	"Splice_Site",
	"Translation_Start_Site",
	"Nonstop_Mutation",
	"3'UTR",
	"3'Flank",
	"5'UTR",
	"5'Flank",
	"IGR",
	"Intron",
	"RNA",
	"Targeted_Region"
]

# CGI GLOBALS
CGI_VARIANT_CLASS_FILTER = [
	"INTRON",
	"TSS-UPSTREAM",
	"UTR5",
	"UTR3",
	"UTR",
	"SPAN5",
	"SPAN3",
	"SPAN",
	"SYNONYMOUS",
	"IGR",
	"NO-CHANGE",
	"UPSTREAM"
]

CGI_INDEL_VARIANT_CLASSES = [
	"INSERT",
	"DELETE",
	"INSERT+",
	"DELETE+",
	"FRAMESHIFT"
]

# NO DIRECT MAPPING FOR: UTR, SPAN, FRAMESHIFT - These require additional data either from other variant classes
# if datum contains comma-separated variant classes (i.e., NO-CHANGE,DELETE,FRAMESHIFT for one record)
# or value of variant type
# > If only UTR is present then use IGR as default
# > If only SPAN is present then use IGR as default
# > If only FRAMESHIFT is present and variant type is INS/DEL then use Frame_Shift_Ins, Frame_Shift_Del.
#   Otherwise, if variant type is not INS/DEL then use Missense_Mutation
CGI_VARIANT_CLASS_MAP = {
	"INTRON":"Intron",
	"TSS-UPSTREAM":"5'Flank",
	"UPSTREAM":"5'Flank",
	"UTR5":"5'UTR",
	"UTR3":"3'UTR",
	"SPAN5":"5'UTR",
	"SPAN3":"3'UTR",
	"SYNONYMOUS":"Silent",
	"MISSTART":"Translation_Start_Site",
	"DONOR":"Splice_Site",
	"ACCEPTOR":"Splice_Site",
	"DISRUPT":"Splice_Site",
	"NO-CHANGE":"Silent",
	"MISSENSE":"Missense_Mutation",
	"NONSENSE":"Nonsense_Mutation",
	"NONSTOP":"Nonstop_Mutation",
	"DELETE":"In_Frame_Del",
	"INSERT":"In_Frame_Ins",
	"DELETE+":"Frame_Shift_Del",
	"INSERT+":"Frame_Shift_Ins"
}

# problematic files tracker
PROBLEMATIC_FILES_REPORT = {}

def get_file_header(filename):
	""" Returns the file header as a list. """
	header = []
	raw_header_line = ""
	for line in extract_file_data(filename):
		if line.startswith("#"):
			continue
		raw_header_line = line
		header = list(map(str.strip, line.split("\t")))
		break
	if not header:
		print("Could not extract header from mutation data file: %s - Exiting..." % (filename))
		sys.exit(2)
	return (header, raw_header_line)

def detect_file_encoding(filename):
	"""
		Reads the first million bytes of a file
		to detect the type of encoding.
	"""
	with open(filename, "rb") as data_file:
		encoding = detect(data_file.read(1000000))
		if encoding["encoding"] != "ascii":
			print("[WARNING] detect_file_encoding(), Non-ASCII encoding detected in file %s, encoding type detected: %s - non-ASCII characters will be ignored and removed" % (filename, encoding["encoding"]))
		return encoding["encoding"]

def is_missing_data_value(value):
	return (value in NULL_OR_MISSING_VALUES)

def resolve_ref_allele(data):
	""" Resolves reference allele. """
	ref = ""
	for col in REFERENCE_ALLELE_COLUMNS:
		if col in data.keys():
			ref = data[col]
	return ref

def resolve_tumor_seq_alleles(data, ref_allele):
	""" Resolve the tumor seq allele. """
	tum_seq_allele1 = ""
	tum_seq_allele2 = ""

	if MUTATED_TO_ALLELE_COLUMN in data.keys():
		# use ref allele for tumor seq allele 1 if "mutated_from_allele" not present
		# but "mutated_to_allele" is present
		tum_seq_allele1 = data.get(MUTATED_FROM_ALLELE_COLUMN, ref_allele)
		tum_seq_allele2 = data[MUTATED_TO_ALLELE_COLUMN]
		return (tum_seq_allele1, tum_seq_allele2)

	for column in TUMOR_SEQ_ALLELE1_COLUMNS:
		if column in data.keys():
			tum_seq_allele1 = process_datum(data[column])
			if tum_seq_allele1 != "":
				break
	for column in TUMOR_SEQ_ALLELE2_COLUMNS:
		if column in data.keys():
			tum_seq_allele2 = process_datum(data[column])
			if tum_seq_allele2 != "":
				break

	# if both tumor seq alleles are empty then there might be something wrong with the data
	if tum_seq_allele1 == "" and tum_seq_allele2 == "":
		return ("", "")

	# resolve tumor seq allele 1 from the tumor genotype column if it still has not been resolved
	if TUMOR_GENOTYPE_COLUMN in data.keys() and tum_seq_allele1 == "":
		tum_seq_allele1 = re.split("[\/|]", data[TUMOR_GENOTYPE_COLUMN])[0]

	# the importer determines which tumor seq allele to use as the alt allele
	# so simply return the resolved values as they are
	return (tum_seq_allele1, tum_seq_allele2)

def resolve_variant_classification(data, variant_type, ref_allele, alt_allele):
	""" Resolves the variant classification. """

	variant_class = ""
	for column in VARIANT_CLASSIFICATION_COLUMNS:
		if column in data.keys():
			variant_class = process_datum(data[column])
			break
	# if variant classification is valid then return, else try to resolve value
	if variant_class in VALID_VARIANT_CLASSIFICATIONS:
		return variant_class

	# if empty string then assume missense or indel - let annotator resolve correct variant class
	if variant_class == "":
		# if indel then determine whether in frame or out of frame
		if variant_type in ["INS", "DEL"]:
			if variant_type == "INS":
				in_frame_variant = (len(alt_allele) % 3 == 0)
			else:
				in_frame_variant = (len(ref_allele) % 3 == 0)
			if in_frame_variant:
				variant_class = "In_Frame_" + variant_type.title()
			else:
				variant_class = "Frame_Shift_" + variant_type.title()
		else:
			# let annotator figure it out from the VEP response
			variant_class = "Missense_Mutation"
	else:
		variant_class_list = re.split("[\\|, ]", variant_class)
		variant_class = resolve_complex_variant_classification(data, variant_class_list, variant_type)
	return variant_class

def resolve_start_position(data):
	""" Resolves the start position. """
	start_pos = ""
	for column in START_POSITION_COLUMNS:
		if column in data.keys():
			start_pos = process_datum(data[column])
			break
	return start_pos

def resolve_end_position(data, start_pos, variant_type, ref_allele):
	""" Resolves the end position. """
	end_pos = ""
	for column in END_POSITION_COLUMNS:
		if column in data.keys():
			end_pos = process_datum(data[column])
			break
	# if insertion then end pos is start pos + 1
	if variant_type == "INS" or ref_allele == "-":
		try:
			end_pos = str(int(start_pos)+1)
		except ValueError:
			print(data)
			sys.exit(2)

	# resolve end pos from ref allele length if empty string
	if end_pos == "":
		end_pos = str(int(start_pos) + len(ref_allele) - 1)

	return end_pos

def is_valid_chromosome(chromosome):
	return (
			(is_valid_integer(chromosome) and int(chromosome) > 0 and int(chromosome) <= 24) or
			(chromosome.upper() in ["X", "Y", "MT"])
		)

def is_valid_integer(value):
	try:
		int(value)
	except ValueError:
		return False
	return True

def resolve_complex_variant_classification(data, variant_class_list, variant_type):
	"""
		Resolve a complex variant classification.
		Note: Complex variant classifications may be from data generated from CGI.
	"""

	# MISSTART variant classes take precedence over other CGI variant classes
	if "MISSTART" in variant_class_list:
		return CGI_VARIANT_CLASS_MAP["MISSTART"]

	# if length of variant class list is 1 and variant class is in list of filtered CGI variants then return
	# direct mapping of variant class or IGR by default, as in the cases of SPAN and UTR variant classes
	if len(variant_class_list) == 1 and variant_class_list[0] in CGI_VARIANT_CLASS_FILTER:
		return CGI_VARIANT_CLASS_MAP.get(variant_class_list[0], "IGR")

	# check if any CGI variant classes that we normally filter are present in given variant_class_list
	filtered_cgi_var_classes = [var_class for var_class in variant_class_list if var_class in CGI_VARIANT_CLASS_FILTER]
	if len(filtered_cgi_var_classes) > 0:
		# if filtered CGI variant classes are present then return the first one that can be directly mapped to
		# a standard variant classification
		for var_class in filtered_cgi_var_classes:
			if var_class in CGI_VARIANT_CLASS_MAP.keys():
				return CGI_VARIANT_CLASS_MAP[var_class]

	# Map var classes to standard variant classifications or use orig var class name if not found.
	# Variant classes that do not map directly include FRAMESHIFT, SPAN, UTR
	variant_class_candidates = list(map(lambda x: CGI_VARIANT_CLASS_MAP.get(x, x), variant_class_list))

	# splice sites take precedence over other variant classifications
	if "Splice_Site" in variant_class_candidates:
		return "Splice_Site"

	# if variant type is INS/DEL or INSERT/DELETE/FRAMESHIFT/INSERT+/DELETE+ in input variant_class_list then
	# return either In_Frame_Ins/Del or Frame_Shift_Ins/Del
	indel_variant_classes = [var_class for var_class in variant_class_list if var_class in CGI_INDEL_VARIANT_CLASSES]
	if variant_type in ["INS", "DEL"] and len(indel_variant_classes) > 0:
		for var_class in indel_variant_classes:
			if var_class == "FRAMESHIFT":
				return "Frame_Shift_" + variant_type.title()
			else:
				return CGI_VARIANT_CLASS_MAP[var_class]

	# if indel variant classes present but variant_type is not INS or DEL - variant type will need to be fixed later
	for var_class in indel_variant_classes:
		if var_class != "FRAMESHIFT":
			return CGI_VARIANT_CLASS_MAP[var_class]

	# if variant class is not an indel then variant class is either Missense_Mutation, Nonsense_Mutation, Nonstop_Mutation, or NO-CHANGE/Silent
	# if variant type is INS/DEL then return Frame_Shift_Ins/Del if Nonsense_Mutation or Nonstop_Mutation in variant class candidates
	# or In_Frame_Ins/Del if Missense_Mutation in candidates - Nonsense/Nonstop mutations take precedence over Missense
	if variant_type in ["INS", "DEL"]:
		if "Nonsense_Mutation" in variant_class_candidates or "Nonstop_Mutation" in variant_class_candidates:
			return "Frame_Shift_" + variant_type.title()
		elif "Missense_Mutation" in variant_class_candidates:
			return "In_Frame_" + variant_type.title()

	# if variant type not INS/DEL then must be SNP/SNV/SUB
	# Return Nonsense_Mutation, Nonstop_Mutation, or Missense_Mutation with Nonsense/Nonstop taking precedence
	for var_class in ["Nonstop_Mutation", "Nonsense_Mutation","Missense_Mutation"]:
		if var_class in variant_class_candidates:
			return var_class

	# if variant class is empty but variant type is SNP then return missense mutation
	if variant_type in ["SNP", "DNP", "TNP", "ONP"]:
		return "Missense_Mutation"

	# if variant class can't be resolve by this point then arbitrarily return first in list that can be mapped directly
	return variant_class_candidates[0]

def resolve_variant_type(data, ref_allele, tum_seq_allele):
	""" Resolve variant type.  """
	variant_type = ""
	for column in VARIANT_TYPE_COLUMNS:
		if column in data.keys():
			variant_type = process_datum(data[column].upper())
			break
	if variant_type == "1":
		return "SNP"
	elif variant_type == "2":
		return "INS"
	elif variant_type == "3":
		return "DEL"

	# if variant type is empty or not valid then try
	# resolving it based on ref allele and tum seq allele values
	if not variant_type in VALID_VARIANT_TYPES:
		if len(ref_allele) > len(tum_seq_allele):
			variant_type = "DEL"
		elif len(ref_allele) < len(tum_seq_allele):
			variant_type = "INS"
		elif len(ref_allele) == len(tum_seq_allele):
			if len(ref_allele) == 1:
				variant_type = "SNP"
			elif len(ref_allele) == 2:
				variant_type = "DNP"
			elif len(ref_allele) == 3:
				variant_type = "TNP"
			else:
				variant_type = "ONP"

	# if variant type still not in VALID_VARIANT_TYPES then must be SUB/SNV
	if variant_type == "SNV":
		variant_type = "SNP"
	elif variant_type.upper() == "SUB":
		variant_type = "ONP"

	# if variant type is still an empty string then report error to user
	if variant_type == "":
		message = "Could not salvage variant type from alleles [ ref allele = %s , tumor allele = %s ] " % (ref_allele, tum_seq_allele)
		print_warning(message)
	return variant_type

def extract_file_data(filename):
	"""
		Reads in file data with the appropriate encoding
		type and returns file data as ascii.
	"""
	encoding_type = detect_file_encoding(filename)
	with open(filename, encoding = encoding_type) as data_file:
		filedata = []
		for line in data_file.readlines():
			# TODO: report on a line-by-line basis if there are non-ascii characters encountered
			filedata.append(line.encode("ascii", "ignore").decode("ascii"))
		return filedata

def is_valid_mutation_record_to_process(raw_header_line, header, line):
	"""
		Determines whether the current line is a valid record in the mutation file that should be processed.

		A record is considered invalid if:
			- it is a duplicate of the header line
			- the line begins with a '#'
			- the line is empty
	"""
	if line.strip() == raw_header_line.strip():
		return False
	if not line.strip():
		return False
	if line.startswith("#"):
		return False
	return True

def is_malformed_record(filename, record_index, line, header):
	"""
		Verify whether the number of fields in the line matches the number of fields in the header.
	"""
	split_line = list(map(str.strip, line.split("\t")))
	if len(split_line) != len(header):
		message = "[line %s] Encountered record in file where the number of fields in the line does not match the expected number of fields in the header, skipping processing of file..." % (record_index + 1)
		print_warning(message)
		update_problematic_report_for_file(filename, "ERROR", message)
		return True
	return False

def map_data_values_to_header(header, line):
	"""
		Maps data values in 'line' to the column headers.
	"""
	data = dict(zip(header, list(map(str.strip, line.split("\t")))))
	return data

def init_maf_record():
	""" Creates a new MAF record with default values for every header. """
	maf_data = dict(zip(MAF_HEADER, ["" for column in MAF_HEADER]))

	# set defaults
	maf_data["Strand"] = DEFAULT_STRAND
	maf_data["Validation_Status"] = DEFAULT_VALIDATION_STATUS
	return maf_data

def resolve_match_normal_sample_barcode(data):
	""" Resolves the matched normal sample barcode. """
	barcode = ""
	for column in MATCHED_NORMAL_SAMPLE_BARCODE_COLUMNS:
		if column in data.keys():
			barcode = process_datum(data[column])
			break
	return barcode

def resolve_hugo_symbol(data):
	""" Resolves the hugo symbol. """
	hugo_symbol = ""
	for column in ["Hugo_Symbol", "HugoSymbol", "Gene Symbol", "GENE"]:
		if column in data.keys():
			hugo_symbol = process_datum(data[column].split("|")[0])
			break
	if hugo_symbol == "":
		hugo_symbol = "Unknown"
	return hugo_symbol

def process_datum(value):
	"""
		Returns a cleaned up data value.
		A null or 'empty' value is always returned as an empty string.
	"""
	try:
		vfixed = value.strip().replace('"',"")
	except AttributeError:
		vfixed = ""
	if is_missing_data_value(vfixed):
		return ""
	return vfixed

def resolve_ncbi_build(data):
	""" Resolves the NCBI build. """
	ncbi_build = ""
	for col in ["NCBI_Build", "assembly_version"]:
		if col in data.keys():
			ncbi_build = process_datum(data.get(col, ""))
			break
	if ncbi_build == "":
		ncbi_build = DEFAULT_NCBI_BUILD
	return ncbi_build

def resolve_chromosome(data):
	""" Resolves the chromosome. """
	chromosome = ""
	for column in CHROMOSOME_COLUMNS:
		if column in data.keys():
			chromosome = process_datum(data[column]).replace("chr","")
			break
	return chromosome.split("_")[0]

def resolve_dbsnp_rs(data):
	""" Resolves the dbSNP ID. """
	dbsnp_rs = ""
	for column in ["dbSNP_RS", "dbSNP rsID"]:
		if column in data.keys():
			dbsnp_rs = process_datum(data[column])
			break
	return dbsnp_rs

def resolve_sequence_source(data, sequence_source):
	""" Resolves the sequence source. """
	# convert sequence strategy to sequence source values if column present in data
	if "sequencing_strategy" in data.keys():
		if data["sequencing_strategy"] == "1":
			return "WGS"
		elif data["sequencing_strategy"] == "3":
			return "WXS"
	# fall back on parsing sequence source column or fall on default
	seq_source = process_datum(data.get("Sequence_Source", ""))
	if seq_source == "":
		seq_source = sequence_source
	return seq_source

def resolve_validation_method(data):
	""" Resolves validation method. """
	validation_method = ""
	for col in VALIDATION_METHOD_COLUMNS:
		if col in data.keys():
			validation_method = process_datum(data.get(col, ""))
			break
	if validation_method == "60":
		validation_method = "Illumina HiSeq (RNAseq)"
	elif validation_method == "6":
		validation_method = "454 sequencing"
	elif validation_method == "67":
		validation_method = "Ion Torrent PGM"
	elif validation_method == "-888":
		validation_method = "NA"
	return validation_method

def resolve_center_name(data, center_name):
	""" Resolves the chromosome. """
	center = process_datum(data.get("Center", ""))
	if center == "":
		center = center_name
	return center

def resolve_verification_status(data):
	""" Resolves the verification status. """
	ver_status = ""
	for col in VERIFICATION_STATUS_COLUMNS:
		if col in data.keys():
			ver_status = process_datum(data[col])
	if ver_status == "1":
		ver_status = "Verified"
	elif ver_status == "":
		ver_status = DEFAULT_VERIFICATION_STATUS
	return ver_status

def resolve_mutation_status(data):
	""" Resolves the mutation status. """
	mutation_status = process_datum(data.get("Mutation_Status", DEFAULT_MUTATION_STATUS))
	if mutation_status == "":
		mutation_status = DEFAULT_MUTATION_STATUS
	return mutation_status

def resolve_sequencer(data):
	""" Resolves the sequencer. """
	sequencer = ""
	for col in ["Sequencer", "platform"]:
		if col in data.keys():
			sequencer = process_datum(data.get(col, ""))
			break
	if sequencer == "60":
		sequencer = "Illumina HiSeq 2000"
	return sequencer

def resolve_allele_counts_data(data, maf_data):
	"""
		Resolves allele counts data based on which tumor/normal reads columns are present.
	"""

	# get defaults from expected maf columns first
	t_depth = process_datum(data.get("t_depth", ""))
	t_ref_count = process_datum(data.get("t_ref_count", ""))
	t_alt_count = process_datum(data.get("t_alt_count", ""))
	n_depth = process_datum(data.get("n_depth", ""))
	n_ref_count = process_datum(data.get("n_ref_count", ""))
	n_alt_count = process_datum(data.get("n_alt_count", ""))
	if "Tumor_ReadCount_Total" in data.keys():
		# tumor counts
		t_depth = process_datum(data.get("Tumor_ReadCount_Total", ""))
		t_ref_count = process_datum(data.get("Tumor_ReadCount_Ref", ""))
		t_alt_count = process_datum(data.get("Tumor_ReadCount_Alt", ""))
		# normal counts
		n_depth = process_datum(data.get("Normal_ReadCount_Total", ""))
		n_ref_count = process_datum(data.get("Normal_ReadCount_Ref", ""))
		n_alt_count = process_datum(data.get("Normal_ReadCount_Alt", ""))
	if "TTotCovVal" in data.keys():
		# tumor counts
		t_depth = process_datum(data.get("TTotCovVal", ""))
		t_alt_count = process_datum(data.get("TVarCovVal", ""))
		if not t_depth in ["", "None"] and not t_alt_count in ["", "None"]:
			t_ref_count = str(int(t_depth) - int(t_alt_count))
		# normal counts
		n_depth = process_datum(data.get("NTotCovVal", ""))
		n_alt_count = process_datum(data.get("NVarCovVal", ""))
		if not n_depth in ["", "None"] and not n_alt_count in ["", "None"]:
			n_ref_count = str(int(n_depth)-int(n_alt_count))
	elif "Ref_Allele_Coverage" in data.keys():
		# tumor counts
		t_ref_count = process_datum(data.get("Ref_Allele_Coverage", ""))
		t_alt_count = process_datum(data.get("Tumor_Seq_Allele1_Coverage", ""))
		if t_ref_count != "" and t_alt_count != "":
			t_depth = str(int(t_ref_count) + int(t_alt_count))
		# normal counts
		n_ref_count = process_datum(data.get("Normal_Ref_Allele_Coverage", ""))
		n_alt_count = process_datum(data.get("Normal_Seq_Allele1_Coverage", ""))
		if n_ref_count != "" and n_alt_count != "":
			n_depth = str(int(n_ref_count) + int(n_alt_count))
	if "total_read_count" in data.keys():
		# tumor counts
		t_depth = process_datum(data.get("total_read_count", ""))
		t_alt_count = process_datum(data.get("mutant_allele_read_count", ""))
		if not t_depth in ["", "None"] and not t_alt_count in ["", "None"]:
			t_ref_count = str(int(t_depth) - int(t_alt_count))


	maf_data["t_depth"] = t_depth
	maf_data["t_ref_count"] = t_ref_count
	maf_data["t_alt_count"] = t_alt_count
	maf_data["n_depth"] = n_depth
	maf_data["n_ref_count"] = n_ref_count
	maf_data["n_alt_count"] = n_alt_count

	return maf_data

def resolve_variant_allele_data(data, maf_data):
	"""
		Resolves the variant allele data.
		Updates maf_data dict with values for:
		- Variant_Classification
		- Variant_Type
		- Reference_Allele
		- Tumor_Seq_Allele1
		- Tumor_Seq_Allele2
		- Start_Position
		- End_Position
	"""
	ref_allele = resolve_ref_allele(data)
	(tumor_seq_allele1, tumor_seq_allele2) = resolve_tumor_seq_alleles(data, ref_allele)

	# set the general tumor_seq_allele as the first non-ref allele encountered
	# this will be used to resolve the variant classification and variant type
	# if there are no tumor alleles that do not match the ref allele then use empty string
	# in the event that this happens then there might be something wrong with the data itself
	tumor_seq_allele = ""
	for allele in [tumor_seq_allele1, tumor_seq_allele2]:
		if allele != "" and allele != ref_allele:
			tumor_seq_allele = allele
			break

	# resolve start and end positions
	start_pos = resolve_start_position(data)

	# if the alleles share a common prefix then remove and adjust the start position accordingly
	if not is_missing_data_value(ref_allele) and not is_missing_data_value(tumor_seq_allele) and not is_missing_data_value(start_pos):
		common_prefix = os.path.commonprefix([ref_allele, tumor_seq_allele])
		if common_prefix:
			start_pos = str(int(start_pos) + len(common_prefix))
			ref_allele = ref_allele[len(common_prefix):]
			tumor_seq_allele = tumor_seq_allele[len(common_prefix):]
			if not is_missing_data_value(tumor_seq_allele1):
				tumor_seq_allele1 = tumor_seq_allele1[len(common_prefix):]
			if not is_missing_data_value(tumor_seq_allele2):
				tumor_seq_allele2 = tumor_seq_allele2[len(common_prefix):]

	# ref and tumor seq allele might have been updated to remove common prefixes
	# attempt to resolve the variant type based on the potentially updated allele strings
	variant_type = resolve_variant_type(data, ref_allele, tumor_seq_allele)
	variant_class = resolve_variant_classification(data, variant_type, ref_allele, tumor_seq_allele)
	# fix variant type just in case it was missed before
	if variant_class.endswith("INS") and variant_type != "INS":
		variant_type = "INS"
	elif variant_class.endswith("DEL") and variant_type != "DEL":
		variant_type = "DEL"

	# fix ref allele and tum seq allele for INS or DEL variant types
	if variant_type == "INS" and len(ref_allele) == 0:
		ref_allele = "-"
	elif variant_type == "DEL" and len(tumor_seq_allele) == 0:
		tumor_seq_allele = "-"

	end_pos = resolve_end_position(data, start_pos, variant_type, ref_allele)

	maf_data["Variant_Classification"] = variant_class
	maf_data["Variant_Type"] = variant_type
	maf_data["Reference_Allele"] = ref_allele
	maf_data["Tumor_Seq_Allele1"] = tumor_seq_allele1
	maf_data["Tumor_Seq_Allele2"] = tumor_seq_allele2
	maf_data["Start_Position"] = start_pos
	maf_data["End_Position"] = end_pos

	return maf_data

def update_problematic_report_for_file(filename, problem_type, message):
	# get current errors and warnings for given filename and update accordingly
	datafile_errors_and_warnings = PROBLEMATIC_FILES_REPORT.get(filename, {})

	# append warning/error message to appropriate list (either "ERROR" or "WARNING")
	messages_for_problem_type = datafile_errors_and_warnings.get(problem_type, [])
	messages_for_problem_type.append(message)

	# update global "PROBLEMATIC_FILES_REPORT" with latest info
	datafile_errors_and_warnings[problem_type] = messages_for_problem_type
	PROBLEMATIC_FILES_REPORT[filename] = datafile_errors_and_warnings

def resolve_match_norm_seq_alleles(data, maf_data):
	""" Resolves matched normal seq alleles. """
	norm_allele1 = ""
	for col in MATCHED_NORMAL_SEQ_ALLELE1_COLUMNS:
		if col in data.keys():
			norm_allele1 = process_datum(data.get(col,""))
			break

	norm_allele2 = ""
	for col in MATCHED_NORMAL_SEQ_ALLELE2_COLUMNS:
		if col in data.keys():
			norm_allele2 = process_datum(data.get(col,""))
			if col == "control_genotype":
				norm_allele2 = norm_allele2.split("/")[1]
			break
	maf_data["Match_Norm_Seq_Allele1"] = norm_allele1
	maf_data["Match_Norm_Seq_Allele2"] = norm_allele2
	return maf_data

def capture_warnings_for_extracted_maf_record(filename, record_index, maf_record):
	"""
		Adds warnings to 'PROBLEMATIC_FILES_REPORT' global.
		These are non-critical issues encountered during processing.
	"""

	# check chromosome value
	if not is_valid_chromosome(maf_record["Chromosome"]):
		message = "[line %s], invalid chromosome value encountered (%s)" % ((record_index + 1), maf_record["Chromosome"])
		update_problematic_report_for_file(filename, "WARNING", message)

	# check start/end positions are not empty
	if is_missing_data_value(maf_record["Start_Position"]):
		message = "[line %s], 'Start_Position' is missing" % ((record_index + 1))
		update_problematic_report_for_file(filename, "WARNING", message)
	if is_missing_data_value(maf_record["End_Position"]):
		message = "[line %s], 'End_Position' is missing" % ((record_index + 1))
		update_problematic_report_for_file(filename, "WARNING", message)

	# check allele fields are not empty
	if (
		(is_missing_data_value(maf_record["Reference_Allele"]) or maf_record["Reference_Allele"] == "N") and
		(is_missing_data_value(maf_record["Tumor_Seq_Allele1"]) or maf_record["Tumor_Seq_Allele1"] == "N") and
		(is_missing_data_value(maf_record["Tumor_Seq_Allele2"]) or maf_record["Tumor_Seq_Allele2"] == "N")
		):
		message = "[line %s], all allele fields are missing or invalid values ('Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2'): (%s, %s, %s)" % ((record_index + 1), maf_record["Reference_Allele"], maf_record["Tumor_Seq_Allele1"], maf_record["Tumor_Seq_Allele2"])
		update_problematic_report_for_file(filename, "WARNING", message)

def print_data_loading_summary(filename, records_loaded, rejected_vars_by_varclass, rejected_vars_by_mutstatus, rejected_vars_by_inc_data):
	print("\tTotal records loaded %s " % (records_loaded))
	if rejected_vars_by_varclass > 0:
		print("\tTotal records filtered by variant classification: %s" % (rejected_vars_by_varclass))
	if rejected_vars_by_mutstatus > 0:
		print("\tTotal records rejected due to LOH or Wildtype Mutation Status: %s" % (rejected_vars_by_mutstatus))
	if rejected_vars_by_inc_data > 0:
		print("\tTotal records rejected due to incomplete allele data: %s" % (rejected_vars_by_inc_data))

def print_problematic_files_report():
	""" Reports errors and warnings encountered while processing input files. """
	if not PROBLEMATIC_FILES_REPORT:
		print("No errors encountered during standardize_mutation_data.py run...")
		return

	print("\n ## ERRORS AND WARNINGS ENCOUNTERED ##\nEncountered %s files with errors and/or warnings during standardize_mutation_data.py run...\n" % (len(PROBLEMATIC_FILES_REPORT.keys())))
	print("- Warnings are errors which are specific to record(s) within a mutation data file. Causes of such warning messages may be due to empty data values for certain key fields like allele depth fields, for example. The warning message will detail the field that the issue was encountered for.")
	print("- Critical errors arise if data could not be loaded from a mutation data file, such as with data files containing a header but no records or a file where all of the lines are commented out.\n")
	for filename, datafile_errors_and_warnings in PROBLEMATIC_FILES_REPORT.items():
		print("\n%s report:" % (filename))
		for problem_type, messages in sorted(datafile_errors_and_warnings.items()):
			print("\t%s" % (problem_type))
			for m in messages:
				print("\t\t%s" % (m))
		print("\n")

def create_maf_record_from_maf(filename, data, center_name, sequence_source):
	""" Creates a MAF record from a given input MAF. """
	maf_data = init_maf_record()

	# identify the tumor sample barcode column present in the input MAF
	try:
		for col in TUMOR_SAMPLE_BARCODE_COLUMNS:
			if col in data.keys():
				tumor_sample_barcode_col_name = col
				break
		maf_data["Tumor_Sample_Barcode"] = data[tumor_sample_barcode_col_name]
	except AttributeError:
		message = "[ERROR] create_maf_record_from_maf(), Error enountered while trying to identify the tumor sample barcode column to use from MAF. "
		message += "MAF columns found: %s" % (", ".join(data.keys()))
		PROBLEMATIC_FILES_REPORT[filename] = {"ERROR": [message]}
		print_warning(message)
		return None

	# set values for MAF fields
	maf_data["Matched_Norm_Sample_Barcode"] = resolve_match_normal_sample_barcode(data)
	maf_data["Hugo_Symbol"] = resolve_hugo_symbol(data)
	maf_data["Entrez_Gene_Id"] = process_datum(data.get("Entrez_Gene_Id","").split("|")[0])
	maf_data["NCBI_Build"] = resolve_ncbi_build(data)
	maf_data["Chromosome"] = resolve_chromosome(data)
	maf_data["dbSNP_RS"] = resolve_dbsnp_rs(data)
	maf_data["Sequencing_Phase"]  = process_datum(data.get("Sequencing_Phase", ""))
	maf_data["Sequence_Source"] = resolve_sequence_source(data, sequence_source)
	maf_data["Validation_Method"] = resolve_validation_method(data)
	maf_data["Center"] = resolve_center_name(data, center_name)
	maf_data["Verification_Status"] = resolve_verification_status(data)
	maf_data["Mutation_Status"] = resolve_mutation_status(data)
	maf_data["Sequencer"] = resolve_sequencer(data)
	maf_data["FILTER"] = data.get("FILTER", "")

	# if the verification status if "Verified" then the validation status can be set to "Valid"
	if maf_data["Verification_Status"] == "Verified":
		maf_data["Validation_Status"] = "Valid"

	# resolve allele counts for tumor and normal sample
	resolve_allele_counts_data(data, maf_data)
	resolve_variant_allele_data(data, maf_data)
	resolve_match_norm_seq_alleles(data, maf_data)
	return maf_data

def extract_maf_data_from_file(filename, center_name, sequence_source):
	rejected_vars_by_varclass = 0
	rejected_vars_by_mutstatus = 0
	rejected_vars_by_inc_data = 0
	records_loaded = 0

	maf_data = []
	print("\nLoading data from file: %s" % (filename))
	(header, raw_header_line) = get_file_header(filename)
	for record_index, line in enumerate(extract_file_data(filename)):
		if not is_valid_mutation_record_to_process(raw_header_line, header, line):
			continue
		if is_malformed_record(filename, record_index, line, header):
			return None
		data = map_data_values_to_header(header, line)
		maf_record = create_maf_record_from_maf(filename, data, center_name, sequence_source)
		if not maf_record:
			print_warning("Encountered a critical error, see report at end of standardize_mutation_data.py run - skipping processing of the rest of the mutation data file: %s" % (filename))
			return None

		# capture non-critical data issues as warnings (i.e., data issues that would prevent a successful annotation of single record)
		capture_warnings_for_extracted_maf_record(filename, record_index, maf_record)

		maf_data.append(maf_record)
		records_loaded += 1

	# if file passes all checks but we could not extract any data from the file then report as critical error
	if not maf_data:
		message = "Data could not be extracted from file - output file will not be generated."
		update_problematic_report_for_file(filename, "ERROR", message)
	else:
		print_data_loading_summary(filename, records_loaded, rejected_vars_by_varclass, rejected_vars_by_mutstatus, rejected_vars_by_inc_data)
	return maf_data

def write_standardized_mutation_file(maf_data, output_filename):
	""" Writes standardized MAF data to output file. """
	output_file = open(output_filename, "w")
	output_file.write("\t".join(MAF_HEADER))
	for data in maf_data:
		formatted_data = list(map(lambda x: process_datum(data.get(x, "")), MAF_HEADER))
		output_file.write("\n" + "\t".join(formatted_data))
	output_file.write("\n")
	output_file.close()
	print("\nStandardized MAF written to: %s" % (output_filename))

def generate_maf_from_input_data(datafiles, output_directory, center_name, sequence_source):
	for filename in datafiles:
		head, tail = os.path.split(filename)
		extract_data = False
		if not "normal" in filename:
			extract_data = True
		if extract_data:
			output_filename = os.path.join(output_directory, os.path.splitext(tail)[0]+".maf")
			maf_data = extract_maf_data_from_file(filename, center_name, sequence_source)
			# only generate output file if data successfully extracted from input file
			if maf_data:
				write_standardized_mutation_file(maf_data, output_filename)

@click.command()
@click.option('-i', '--input-data', required=True, help='A list of .maf or .txt files or input data directories, separated by commas')
@click.option('-o', '--output-directory', default='maf2maf_output', help='output data directory [optional]')
@click.option('-c','--center', default='', help = 'name of center (standard MAF field = \'Center\') [optional]')
@click.option('-s', '--sequence-source', default='', help = 'Sequencing source (standard MAF field = \'Sequencing_Source\'), e.g., WXS or WGS [optional]')

def main(input_data, output_directory, center, sequence_source):

	# generate input maf files list
	datafiles = []
	for item in input_data.split(','):
		print("\nLoading data from: %s" % (item))
		if os.path.isfile(item):
			datafiles.append(item)
		elif os.path.isdir(item):
			files = [os.path.join(item,x) for x in os.listdir(item) if x.endswith('.maf') or x.endswith('.txt')]
			datafiles.extend(files)
		else: print("The path \'%s\' is ivalid" % item)

	# Create output directory
	if output_directory=='maf2maf_output': #default folder if output path not specified
		output_directory = os.getcwd()+'/maf2maf_output/'

	if not os.path.isdir(output_directory):
		try:
			os.mkdir(output_directory)
			print("created output folder : '{}' ".format(output_directory))
		except:
			print("The output path \'%s\' is ivalid" % output_directory)
			sys.exit(1)

	# Process .maf or .txt files
	generate_maf_from_input_data(datafiles, output_directory, center, sequence_source.upper())
	print_problematic_files_report()

if __name__ == "__main__":
	main()

