import argparse
import requests

MVI_VARIANT_ENDPOINT = "https://myvariant.info/v1/variant"
MVI_VARIANT_ENDPOINT_FIELDS_PARAM = "gnomad_exome.af"
GNOMAD_EXOME_FIELD = "gnomad_exome"
AF_FIELD = "af"

REQUIRED_HEADERS = {"Content-Type" : "application/x-www-form-urlencoded",
		    "accept" : "application/json"}

GNOMAD_FIELD_MVI_TO_ANNOTATION_MAP = {
	"af" : "gnomAD_AF",
	"af_afr" : "gnomAD_AFR_AF",
	"af_amr" : "gnomAD_AMR_AF",
	"af_asj" : "gnomAD_ASJ_AF",
	"af_eas" : "gnomAD_EAS_AF",
	"af_fin" : "gnomAD_FIN_AF",
	"af_nfe" : "gnomAD_NFE_AF",
	"af_oth" : "gnomAD_OTH_AF",
	"af_sas" : "gnomAD_SAS_AF"
}

def add_gnomad_to_maf(input_maf, output_maf):
	"""
		Main logic, iterate through each line of the input maf
		For each line, construct HGVS query and send to MyVariantInfo
		Parse the response for gnomAD fields
		Write each line (with gnomAD fields filled in) into the output MAF
		Will exit overall script if MAF is missing required fields
		for constructing HGVS query
	"""
	header_indices = {}
	with open(input_maf) as f1, open(output_maf, "w") as f2:
		for line in f1.readlines():
			# line is a comment, just write it to new MAF
			if line.startswith("#"):
				f2.write(line)
				continue
			# first non-comment line is header
			# write header to new MAF
			# get header + indices for later processing
			if not header_indices:
				f2.write(line)
				header = line.strip().split("\t")
				header_indices = dict(zip(header, range(len(header))))
			# rest of non-commented lines are records
			# for each record construct query to MyVariantInfo and get gnomAD values
			# write out record with gnomAD values filled in to new MAF
			else:
				record = line.strip().split('\t')
				query = record[header_indices["Query"]] #line.strip() # construct_hgvs_query(chr, start, end, ...)
				variant = get_query(query)
				gnomad_map = get_gnomad_map(variant)
				record_with_gnomad = add_gnomad_values(record, header_indices, gnomad_map)
				f2.write('\t'.join(record_with_gnomad) + "\n")

def add_gnomad_values(record, header_indices, gnomad_map):
	"""
		Returns a copy of the original record but with gnomad values replaced
	"""
	record_with_gnomad = record
	for gnomad_key in gnomad_map:
		record_with_gnomad[header_indices[gnomad_key]] = str(gnomad_map[gnomad_key])
	return record_with_gnomad
		
def get_query(query):
	"""
		Send request to MyVariantInfo /variant GET endpoint for gnomAD fields
		Logs and returns nothing if non-200 status code returned (this will translate to empty records)
		Else returns JSON (dict) version of response
	"""
	my_variant_info_response = requests.get("%s/%s?fields=%s" % (MVI_VARIANT_ENDPOINT, query, MVI_VARIANT_ENDPOINT_FIELDS_PARAM))
	status_code = my_variant_info_response.status_code
	if status_code != 200:
		print "MyVariantInfo request failed with status code %s for query %s" % (status_code, query)
		return None
	return my_variant_info_response.json()

def post_queries(queries):
	return requests.post("%s?fields=%s" % (MVI_VARIANT_ENDPOINT, MVI_VARIANT_ENDPOINT_FIELDS_PARAM), data = "ids=%s" % (",".join(["chr" + query for query in queries])), headers = REQUIRED_HEADERS)	

def get_gnomad_map(variant):
	"""
		Construct a map of required gnomAD values
		Key is the column header name inside the MAF
		Value is taken from MyVariantInfo response
		Missing gnomAD values set to ""
	"""
	gnomad_map = {}
	for gnomad_mvi_field, gnomad_annotation_field in GNOMAD_FIELD_MVI_TO_ANNOTATION_MAP.items():
		try:
			variant_gnomad_af = variant[GNOMAD_EXOME_FIELD][AF_FIELD]
			gnomad_map[gnomad_annotation_field] = variant_gnomad_af[gnomad_mvi_field]
		except:
			gnomad_map[gnomad_annotation_field] = ""
	return gnomad_map

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input-maf", help="path to input maf", required = True)
	parser.add_argument("-o", "--output-maf", help = "path to output maf", required = True)
	
	args = parser.parse_args()
	input_maf = args.input_maf
	output_maf = args.output_maf
	
	add_gnomad_to_maf(input_maf, output_maf)

main()
