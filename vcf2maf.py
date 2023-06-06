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
import json
from chardet import detect
import warnings

# ---------------------------- GLOBALS ----------------------------

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
    "gnomAD_SAS_AF",
]

VCF_FIXED_HEADER_NON_CASE_IDS = [
    "CHROM",
    "POS",
    "ID",
    "REF",
    "ALT",
    "QUAL",
    "FILTER",
    "INFO",
    "FORMAT",
]

# DEFAULTS
DEFAULT_STRAND = "+"
DEFAULT_VALIDATION_STATUS = "Unknown"
DEFAULT_VERIFICATION_STATUS = "Unknown"
DEFAULT_MUTATION_STATUS = "Somatic"

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
    "Targeted_Region",
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
    "UPSTREAM",
]

CGI_INDEL_VARIANT_CLASSES = ["INSERT", "DELETE", "INSERT+", "DELETE+", "FRAMESHIFT"]
# NO DIRECT MAPPING FOR: UTR, SPAN, FRAMESHIFT - These require additional data either from other variant classes
# if datum contains comma-separated variant classes (i.e., NO-CHANGE,DELETE,FRAMESHIFT for one record)
# or value of variant type
# > If only UTR is present then use IGR as default
# > If only SPAN is present then use IGR as default
# > If only FRAMESHIFT is present and variant type is INS/DEL then use Frame_Shift_Ins, Frame_Shift_Del.
#   Otherwise, if variant type is not INS/DEL then use Missense_Mutation
CGI_VARIANT_CLASS_MAP = {
    "INTRON": "Intron",
    "TSS-UPSTREAM": "5'Flank",
    "UPSTREAM": "5'Flank",
    "UTR5": "5'UTR",
    "UTR3": "3'UTR",
    "SPAN5": "5'UTR",
    "SPAN3": "3'UTR",
    "SYNONYMOUS": "Silent",
    "MISSTART": "Translation_Start_Site",
    "DONOR": "Splice_Site",
    "ACCEPTOR": "Splice_Site",
    "DISRUPT": "Splice_Site",
    "NO-CHANGE": "Silent",
    "MISSENSE": "Missense_Mutation",
    "NONSENSE": "Nonsense_Mutation",
    "NONSTOP": "Nonstop_Mutation",
    "DELETE": "In_Frame_Del",
    "INSERT": "In_Frame_Ins",
    "DELETE+": "Frame_Shift_Del",
    "INSERT+": "Frame_Shift_Ins",
}

# KEY COLUMN NAMES
TUMOR_SEQ_ALLELE1_COLUMNS = ["Tumor_Seq_Allele1", "TumorSeq_Allele1"]
TUMOR_SEQ_ALLELE2_COLUMNS = ["Tumor_Seq_Allele2", "TumorSeq_Allele2"]
MUTATED_FROM_ALLELE_COLUMN = "mutated_from_allele"
MUTATED_TO_ALLELE_COLUMN = "mutated_to_allele"
VARIANT_TYPE_COLUMNS = ["Variant_Type", "VariantType", "mut_type", "mutation_type"]
VARIANT_CLASSIFICATION_COLUMNS = [
    "Variant_Classification",
    "class",
    "Transcript architecture around variant",
]
START_POSITION_COLUMNS = ["Start_Position", "Start_position", "POS", "chromosome_start"]
END_POSITION_COLUMNS = ["End_Position", "End_position", "chromosome_end"]
REFERENCE_ALLELE_COLUMNS = ["Reference_Allele", "reference_genome_allele"]
TUMOR_GENOTYPE_COLUMN = "tumour_genotype"
MATCHED_NORMAL_SAMPLE_BARCODE_COLUMNS = [
    "Match_Normal_Sample_Barcode",
    "Matched_Norm_Sample_Barcode",
    "matched_sample_id",
]
CHROMOSOME_COLUMNS = ["Chromosome", "chromosome", "CHROM"]
VERIFICATION_STATUS_COLUMNS = ["Verification_Status", "verification_status"]
VALIDATION_METHOD_COLUMNS = ["Validation_Method", "verification_platform"]
MATCHED_NORMAL_SEQ_ALLELE1_COLUMNS = ["Match_Norm_Seq_Allele1", "mutated_to_allele"]
MATCHED_NORMAL_SEQ_ALLELE2_COLUMNS = ["Match_Norm_Seq_Allele2", "control_genotype"]
TUMOR_SAMPLE_BARCODE_COLUMNS = [
    "Tumor_Sample_Barcode",
    "analyzed_sample_id",
    "submitted_sample_id",
]

# VCF KEYS FOR RESOLVING WHICH VCF PIPIELINE WAS USED
VCF_STRELKA_KEY_COLUMNS = ["AU", "CU", "GU", "TU"]
VCF_CAVEMAN_KEY_COLUMNS = ["FAZ", "FCZ", "FGZ", "FTZ", "RAZ", "RCZ", "RGZ", "RTZ"]
VCF_ION_TORRENT_KEY_COLUMNS = ["AO", "RO"]
VCF_DELLY_KEY_COLUMNS = ["RR", "RV"]
VCF_CGPPINDEL_KEY_COLUMNS = ["PP", "NP", "PR", "NR"]
VCF_ALT_ALLELE_FRACTION_KEY_COLUMNS = ["FA", "DP"]
VCF_MPILEUP_BCFTOOLS_KEY_COLUMNS = ["DV", "DP"]

VCF_REFERENCE_ALLELE_VALUE = "0"
VCF_VARIANT_ALLELE1_VALUE = "1"
VCF_VARIANT_ALLELE2_VALUE = "2"
VCF_DEFAULT_VARIANT_ALLELE_IDX = 1

NULL_OR_MISSING_VALUES = set(["", "N/A", None, ".", "?"])
NULL_OR_MISSING_VCF_VALUES = set(["", ".", "./."])

NORMAL_IDENTIFIERS = ["NORMAL", "DONOR"]

# problematic files tracker
PROBLEMATIC_FILES_REPORT = {}

# ----------------------------------------------------------------


def is_valid_integer(value):
    try:
        int(value)
    except ValueError:
        return False
    return True


def is_missing_vcf_data_value(value):
    return value in NULL_OR_MISSING_VCF_VALUES


def all_values_in_list(values_to_check, input_list):
    for value in values_to_check:
        if not value in input_list:
            return False
    return True


def is_varscan_vcf(vcf_format_data_keys):
    return "RD" in vcf_format_data_keys


def is_somatic_sniper_vcf(vcf_format_data_keys):
    return not "AD" in vcf_format_data_keys and "BCOUNT" in vcf_format_data_keys


def is_strelka_snp_vcf(vcf_format_data_keys):
    return not "AD" in vcf_format_data_keys and all_values_in_list(
        VCF_STRELKA_KEY_COLUMNS, vcf_format_data_keys
    )


def is_strelka_indel_vcf(vcf_format_data_keys):
    return not "AD" in vcf_format_data_keys and "TIR" in vcf_format_data_keys


def is_caveman_vcf(vcf_format_data_keys):
    return not "AD" in vcf_format_data_keys and all_values_in_list(
        VCF_CAVEMAN_KEY_COLUMNS, vcf_format_data_keys
    )


def is_ion_torrent_vcf(vcf_format_data_keys):
    return not "AD" in vcf_format_data_keys and all_values_in_list(
        VCF_ION_TORRENT_KEY_COLUMNS, vcf_format_data_keys
    )


def is_delly_vcf(vcf_format_data_keys):
    return not "AD" in vcf_format_data_keys and all_values_in_list(
        VCF_DELLY_KEY_COLUMNS, vcf_format_data_keys
    )


def is_cgppindel_vcf(vcf_format_data_keys):
    return not "AD" in vcf_format_data_keys and all_values_in_list(
        VCF_CGPPINDEL_KEY_COLUMNS, vcf_format_data_keys
    )


def is_alt_allele_fraction_vcf(vcf_format_data_keys):
    return not "AD" in vcf_format_data_keys and all_values_in_list(
        VCF_ALT_ALLELE_FRACTION_KEY_COLUMNS, vcf_format_data_keys
    )


def is_mpileup_bcftools_vcf(vcf_format_data_keys):
    return not "AD" in vcf_format_data_keys and all_values_in_list(
        VCF_MPILEUP_BCFTOOLS_KEY_COLUMNS, vcf_format_data_keys
    )


def resolve_vcf_allele_depth_values(
    mapped_sample_format_data, vcf_alleles, variant_allele_idx, vcf_data
):
    """
    Resolves the allele depth values based on the type of VCF pipeline identified.

    Support VCF pipelines/methods for resolving allele counts:
    1. VarScan
    2. SomaticSniper
    3. Strelka (SNPs only)
    4. Strelka (INDELs only)
    5. CaVEMan
    6. Ion Torrent
    7. Delly
    8. cgpPINDEL
    9. VCF ALT allele fractions (derive values from allele fractions)
    10. MPileUp/BCFTools
    11. "AD" only contains a single value (does not contain a comma)
    12. If none of the above criteria are met then allele depths are set to empty strings
    """
    ref_count = ""
    alt_count = ""
    depth = ""

    # get list of keys stored in the VCF 'FORMAT' column - this is used
    # to determine how to resolve the alelle counts
    vcf_format_data_keys = mapped_sample_format_data.keys()

    # init allele depth values to list of empty strings matching length of "vcf_alleles"
    allele_depth_values = [""] * len(vcf_alleles)

    # if AD is defined, then parse out all REF/ALT allele depths, or whatever is in it
    if "AD" in vcf_format_data_keys and not is_missing_vcf_data_value(
        mapped_sample_format_data["AD"]
    ):
        # attempt to parse values as int - if not an int then set value to empty string
        allele_depth_values = []
        for value in mapped_sample_format_data["AD"].split(","):
            if is_valid_integer(value):
                allele_depth_values.append(value)
            else:
                allele_depth_values.append("")

    # 1. VarScan VCF: handle VarScan VCF lines where AD contains only 1 depth, and REF allele depth is in RD
    if len(allele_depth_values) == 1 and is_varscan_vcf(vcf_format_data_keys):
        allele_depth_values = [""] * len(vcf_alleles)
        allele_depth_values[0] = process_datum(mapped_sample_format_data["RD"])
        allele_depth_values[variant_allele_idx] = process_datum(
            mapped_sample_format_data["AD"]
        )

    # 2. SomaticSniper: handle SomaticSniper VCF lines, where allele depths must be extracted from BCOUNT
    elif is_somatic_sniper_vcf(vcf_format_data_keys):
        # bcount values are always reported in the order of "A", "C", "G", "T"
        tumor_bcount_values = mapped_sample_format_data["BCOUNT"].split(",")
        read_depth_bases = ["A", "C", "G", "T"]
        read_depth_bases_to_counts_map = {}
        for i, bcount_val in enumerate(tumor_bcount_values):
            corresponding_base = read_depth_bases[i]
            read_depth_bases_to_counts_map[read_depth_bases[i]] = bcount_val

        allele_depth_values = [""] * len(vcf_alleles)
        for i, allele in enumerate(vcf_alleles):
            allele_depth_values[i] = read_depth_bases_to_counts_map.get(allele, "")

    # 3. Strelka (SNP): handle VCF SNV lines by Strelka, where allele depths are in AU:CU:GU:TU
    elif is_strelka_snp_vcf(vcf_format_data_keys):
        # need to convert the read values to an integer so we can sort
        read_depth_bases_to_counts_map = {}
        for k in VCF_STRELKA_KEY_COLUMNS:
            value = -1
            # If multiple tiers, report tier1 counts by picking the first comma-delimited value
            if is_valid_integer(mapped_sample_format_data[k].split(",")[0]):
                value = int(mapped_sample_format_data[k].split(",")[0])
            read_depth_bases_to_counts_map[k.replace("U", "")] = value

        # if the only alt allele is N then set it to the allele with the highest non-ref readcount
        if len(vcf_alleles) == 2 and vcf_alleles[1] == "N":
            sorted_read_depths = sorted(
                read_depth_bases_to_counts_map.items(), lambda k: k[1], reverse=True
            )
            # find highest non-ref readcount
            if sorted_read_depths[0][0] != vcf_alleles[0]:
                vcf_alleles[variant_allele_idx] = sorted_read_depths[0][0]
            else:
                vcf_alleles[variant_allele_idx] = sorted_read_depths[1][0]

        # change values back to string representations of the numbers or empty strings if default value was used
        allele_depth_values = [""] * len(vcf_alleles)
        for i, allele in enumerate(vcf_alleles):
            read_depth_value = str(read_depth_bases_to_counts_map.get(allele, -1))
            # if value was set to default value of -1 from above during the sorting then set to empty string
            if read_depth_value == "-1":
                read_depth_value = ""
            allele_depth_values[i] = read_depth_value

    # 4. Strelka (INDEL): handle VCF INDEL lines by Strelka, where variant allele depth is in TIR, reference allele depth is in TAR.
    elif is_strelka_indel_vcf(vcf_format_data_keys):
        # pick tier1 counts if a comma-delimited value
        allele_depth_values[variant_allele_idx] = mapped_sample_format_data[
            "TIR"
        ].split(",")[0]
        if variant_allele_idx == 1:
            ref_allele_idx = 0
        else:
            ref_allele_idx = 1
        allele_depth_values[ref_allele_idx] = mapped_sample_format_data["TAR"].split(
            ","
        )[0]

    # 5. CaVEMan: handle VCF lines by CaVEMan, where allele depths are in FAZ:FCZ:FGZ:FTZ:RAZ:RCZ:RGZ:RTZ
    elif is_caveman_vcf(vcf_format_data_keys):
        # allele depths are provided for forward strand and reverse strand, add these numbers together to get depth by nucleotide
        caveman_values_as_integers = {}
        for k in VCF_CAVEMAN_KEY_COLUMNS:
            try:
                caveman_values_as_integers[k] = int(mapped_sample_format_data[k])
            except:
                caveman_values_as_integers[k] = 0
        read_depth_bases_to_counts_map["A"] = str(
            caveman_values_as_integers["FAZ"] + caveman_values_as_integers["RAZ"]
        )
        read_depth_bases_to_counts_map["C"] = str(
            caveman_values_as_integers["FCZ"] + caveman_values_as_integers["RCZ"]
        )
        read_depth_bases_to_counts_map["G"] = str(
            caveman_values_as_integers["FGZ"] + caveman_values_as_integers["RGZ"]
        )
        read_depth_bases_to_counts_map["T"] = str(
            caveman_values_as_integers["FTZ"] + caveman_values_as_integers["RTZ"]
        )
        allele_depth_values = [""] * len(vcf_alleles)
        for i, allele in enumerate(vcf_alleles):
            allele_depth_values[i] = read_depth_bases_to_counts_map.get(allele, "")

    # 6. Ion Torrent: handle VCF lines from the Ion Torrent Suite where ALT depths are in AO and REF depths are in RO
    elif is_ion_torrent_vcf(vcf_format_data_keys):
        allele_depth_values = [process_datum(mapped_sample_format_data["RO"])]
        allele_depth_values.extend(mapped_sample_format_data["AO"].split(","))

    # 7. Delly: handle VCF lines from Delly where REF/ALT SV junction read counts are in RR/RV respectively
    elif is_delly_vcf(vcf_format_data_keys):
        allele_depth_values[0] = process_datum(mapped_sample_format_data["RR"])
        allele_depth_values[variant_allele_idx] = process_datum(
            mapped_sample_format_data["RV"]
        )

    # 8. cgpPINDEL: handle VCF lines from cgpPindel, where ALT depth and total depth are in PP:NP:PR:NR
    elif is_cgppindel_vcf(vcf_format_data_keys):
        # reference allele depth and depths for any other ALT alleles must be left undefined
        allele_depth_values[variant_allele_idx] = str(
            float(mapped_sample_format_data["PP"])
            + float(mapped_sample_format_data["NP"])
        )
        mapped_sample_format_data["DP"] = str(
            float(mapped_sample_format_data["PR"])
            + float(mapped_sample_format_data["NR"])
        )

    # 9. VCF ALT allele fractions: handle VCF lines with ALT allele fraction in FA, which needs to be multiplied by DP to get AD
    elif is_alt_allele_fraction_vcf(
        vcf_format_data_keys
    ) and not is_missing_vcf_data_value(mapped_sample_format_data["DP"]):
        allele_depth_values[variant_allele_idx] = "%.0f" % (
            float(mapped_sample_format_data["FA"])
            * float(mapped_sample_format_data["DP"])
        )

    # 10. MPileUp/BCFTools: handle VCF lines from mpileup/bcftools where DV contains the ALT allele depth
    elif is_mpileup_bcftools_vcf(vcf_format_data_keys):
        allele_depth_values[variant_allele_idx] = process_datum(
            mapped_sample_format_data["DV"]
        )

    # 11. AD single value: handle VCF lines where AD contains only 1 value, that we can assume is the variant allele
    elif "AD" in vcf_format_data_keys and len(allele_depth_values) == 1:
        allele_depth_values[variant_allele_idx] = process_datum(
            mapped_sample_format_data["AD"]
        )

    # 12. For all other lines where number of depths is not equal to number of alleles, blank out the depths
    elif len(allele_depth_values) != len(vcf_alleles):
        allele_depth_values = [""] * len(vcf_alleles)

    ref_count = allele_depth_values[0]
    alt_count = allele_depth_values[variant_allele_idx]
    vcf_dp_value = process_datum(mapped_sample_format_data.get("DP", ""))
    # Sanity check that REF/ALT allele depths are lower than the total depth
    if (
        "DP" in vcf_format_data_keys
        and not is_missing_vcf_data_value(vcf_dp_value)
        and (
            not is_missing_vcf_data_value(ref_count)
            and float(ref_count) > float(vcf_dp_value)
            or (
                not is_missing_vcf_data_value(alt_count)
                and float(alt_count) > float(vcf_dp_value)
            )
            or (
                not is_missing_vcf_data_value(ref_count)
                and not is_missing_vcf_data_value(alt_count)
                and ((float(ref_count) + float(alt_count)) > float(vcf_dp_value))
            )
        )
    ):
        mapped_sample_format_data["DP"] = str(sum(map(float, allele_depth_values)))

    # if we have REF/ALT allele depths, but no DP, then set DP equal to the sum of all ADs
    if (
        not is_missing_vcf_data_value(ref_count)
        and not is_missing_vcf_data_value(alt_count)
    ) and (
        not "DP" in vcf_format_data_keys
        or is_missing_vcf_data_value(mapped_sample_format_data["DP"])
    ):
        mapped_sample_format_data["DP"] = str(sum(map(float, allele_depth_values)))
    mapped_sample_format_data["AD"] = ",".join(allele_depth_values)

    try:
        depth = mapped_sample_format_data["DP"]
    except:
        message = (
            "DP could not be resolved for current record in VCF: %s - using default value of empty string..."
            % (str(vcf_data))
        )
        print_warning(message)

    # if depth has been resolved but not ref_count, alt_count then calculate the counts
    # if allele frequency vcf field "AF" exists
    if (
        not is_missing_vcf_data_value(depth)
        and (
            "AF" in mapped_sample_format_data
            and not is_missing_vcf_data_value(mapped_sample_format_data["AF"])
        )
        and (
            is_missing_vcf_data_value(ref_count) or is_missing_vcf_data_value(alt_count)
        )
    ):
        # check if ref count or alt count are still missing but AF VCF field is available
        if is_missing_vcf_data_value(ref_count):
            ref_count = str(
                round(float(depth) * float(mapped_sample_format_data["AF"]))
            )
        if is_missing_vcf_data_value(alt_count) and not is_missing_vcf_data_value(
            ref_count
        ):
            alt_count = str(round(float(depth) - float(ref_count)))

    return (ref_count, alt_count, depth)


def get_vcf_variant_allele_idx(
    tumor_sample_format_data, normal_sample_format_data, vcf_alleles
):
    """
    Determines the variant allele index to use based on genotype information if available.
    If genotype information is not available or a call could not be made for a sample at
    a given locus then use a variant allele idx of 1 by default.

    Allele values:
    - 0 = reference allele (i.e., what is in the 'REF' field)
    - 1 = the first allele listed in 'ALT'
    - 2 = the second allele listed in 'ALT'
    """

    variant_allele_idx = []
    # if genotype information is available and a call was made for the sample then
    # choose the first non-REF allele seen in the sample genotype
    if "GT" in tumor_sample_format_data.keys() and not is_missing_vcf_data_value(
        tumor_sample_format_data["GT"]
    ):
        tumor_sample_genotype_info = re.split("[\/|]", tumor_sample_format_data["GT"])
        variant_allele_idx = [
            allele
            for allele in tumor_sample_genotype_info
            if allele != VCF_REFERENCE_ALLELE_VALUE
            and not is_missing_vcf_data_value(allele)
        ]
        # if possible, choose the first non-REF tumor allele that is also not in normal genotype
        # if this check results in variant_allele_idx being empty or the index value is larger than
        # the size of the alleles then use the default value "VCF_DEFAULT_VARIANT_ALLELE_IDX"
        if (
            normal_sample_format_data != None
            and "GT" in normal_sample_format_data.keys()
            and not is_missing_vcf_data_value(normal_sample_format_data["GT"])
        ):
            normal_sample_genotype_info = re.split(
                "[\/|]", normal_sample_format_data["GT"]
            )
            variant_allele_idx = [
                allele
                for allele in tumor_sample_genotype_info
                if allele != VCF_REFERENCE_ALLELE_VALUE
                and not is_missing_vcf_data_value(allele)
                and not allele in normal_sample_genotype_info
            ]

    # if the idx found is undefined or is equal to or larger than the size of the vcf alleles list
    # then use the vcf default value for the variant allele index
    if (
        variant_allele_idx == []
        or not is_valid_integer(variant_allele_idx[0])
        or int(variant_allele_idx[0]) >= len(vcf_alleles)
    ):
        return VCF_DEFAULT_VARIANT_ALLELE_IDX
    return int(variant_allele_idx[0])


def resolve_vcf_counts_data(
    vcf_data, maf_data, matched_normal_sample_id, tumor_sample_data_col
):
    """Resolves VCF allele counts data."""
    vcf_alleles = [vcf_data["REF"]]
    vcf_alleles.extend(vcf_data["ALT"].split(","))

    tumor_sample_format_data = vcf_data["MAPPED_TUMOR_FORMAT_DATA"]
    normal_sample_format_data = None
    if matched_normal_sample_id in vcf_data.keys():
        normal_sample_format_data = vcf_data["MAPPED_NORMAL_FORMAT_DATA"]

    variant_allele_idx = get_vcf_variant_allele_idx(
        tumor_sample_format_data, normal_sample_format_data, vcf_alleles
    )
    (t_ref_count, t_alt_count, t_depth) = resolve_vcf_allele_depth_values(
        tumor_sample_format_data, vcf_alleles, variant_allele_idx, vcf_data
    )
    maf_data["t_ref_count"] = t_ref_count
    maf_data["t_alt_count"] = t_alt_count
    maf_data["t_depth"] = t_depth

    # only resolve values for normal allele depths if "NORMAL" data is present in VCF
    if normal_sample_format_data:
        (n_ref_count, n_alt_count, n_depth) = resolve_vcf_allele_depth_values(
            normal_sample_format_data, vcf_alleles, variant_allele_idx, vcf_data
        )
        maf_data["n_ref_count"] = n_ref_count
        maf_data["n_alt_count"] = n_alt_count
        maf_data["n_depth"] = n_depth
    return maf_data


def resolve_sequence_source(data):
    """Resolves the sequence source."""
    # convert sequence strategy to sequence source values if column present in data
    if "sequencing_strategy" in data.keys():
        if data["sequencing_strategy"] == "1":
            return "WGS"
        elif data["sequencing_strategy"] == "3":
            return "WXS"
    # fall back on parsing sequence source column or fall on default
    seq_source = process_datum(data.get("Sequence_Source", ""))
    return seq_source


def resolve_chromosome(data):
    """Resolves the chromosome."""
    chromosome = ""
    for column in CHROMOSOME_COLUMNS:
        if column in data.keys():
            chromosome = process_datum(data[column]).replace("chr", "")
            break
    return chromosome.split("_")[0]


def resolve_hugo_symbol(data):
    """Resolves the hugo symbol."""
    hugo_symbol = ""
    for column in ["Hugo_Symbol", "HugoSymbol", "Gene Symbol", "GENE"]:
        if column in data.keys():
            hugo_symbol = process_datum(data[column].split("|")[0])
            break
    if hugo_symbol == "":
        hugo_symbol = "Unknown"
    return hugo_symbol


def init_maf_record():
    """Creates a new MAF record with default values for every header."""
    maf_data = dict(zip(MAF_HEADER, ["" for column in MAF_HEADER]))

    # set defaults
    maf_data["Strand"] = DEFAULT_STRAND
    maf_data["Validation_Status"] = DEFAULT_VALIDATION_STATUS
    return maf_data


def resolve_vcf_allele(vcf_data):
    """Resolves the VCF alleles."""
    ref_allele = ""
    alt_allele = ""

    if vcf_data["ALT"] in ["<DEL>", "<DUP>", "<INV>", "<TRA>"]:
        if vcf_data["REF"] == "N" or vcf_data["REF"] == "":
            ref_allele = vcf_data["INFO"].get("CONSENSUS", "")
        else:
            ref_allele = vcf_data["REF"]

        if ref_allele != "N" and ref_allele != "":
            if vcf_data["ALT"] == "<DEL>":
                alt_allele = "-"
            if vcf_data["ALT"] == "<INV>":
                alt_allele = ref_allele[::-1]
            elif vcf_data["ALT"] == "<DUP>":
                alt_allele = ref_allele * 2
    else:
        ref_allele = process_datum(vcf_data["REF"].split(",")[0])
        alt_allele = process_datum(vcf_data["ALT"].split(",")[0])
        if ref_allele == "":
            ref_allele = "-"
        if alt_allele == "":
            alt_allele = "-"
    return ref_allele, alt_allele


def resolve_start_position(data):
    """Resolves the start position."""
    start_pos = ""
    for column in START_POSITION_COLUMNS:
        if column in data.keys():
            start_pos = process_datum(data[column])
            break
    return start_pos


def resolve_vcf_variant_type(ref_allele, alt_allele):
    """Resolves the VCF variant type."""
    variant_type = ""

    # first check if indel
    if ref_allele == "-" or len(ref_allele) < len(alt_allele):
        variant_type = "INS"
    elif alt_allele == "-" or len(alt_allele) < len(ref_allele):
        variant_type = "DEL"
    else:
        # check whether variant type is type of polymorphism
        if len(ref_allele) == len(alt_allele):
            if len(ref_allele) == 1:
                variant_type = "SNP"
            elif len(ref_allele) == 2:
                variant_type = "DNP"
            elif len(ref_allele) == 3:
                variant_type = "TNP"
            else:
                variant_type = "ONP"

    # if variant type is still empty then report
    if variant_type == "":
        message = (
            "Could not salvage variant type from alleles [ ref allele = %s , tumor allele = %s ] "
            % (ref_allele, alt_allele)
        )
        print_warning(message)

    return variant_type


def resolve_end_position(data, start_pos, variant_type, ref_allele):
    """Resolves the end position."""
    end_pos = ""
    for column in END_POSITION_COLUMNS:
        if column in data.keys():
            end_pos = process_datum(data[column])
            break
    # if insertion then end pos is start pos + 1
    if variant_type == "INS" or ref_allele == "-":
        try:
            end_pos = str(int(start_pos) + 1)
        except ValueError:
            print(data)
            sys.exit(2)

    # resolve end pos from ref allele length if empty string
    if end_pos == "":
        end_pos = str(int(start_pos) + len(ref_allele) - 1)

    return end_pos


def resolve_center_name(data, center_name):
    """Resolves the chromosome."""
    center = process_datum(data.get("Center", ""))
    if center == "":
        center = center_name
    return center


def resolve_sequence_source(data, sequence_source):
    """Resolves the sequence source."""
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
    if (
        len(variant_class_list) == 1
        and variant_class_list[0] in CGI_VARIANT_CLASS_FILTER
    ):
        return CGI_VARIANT_CLASS_MAP.get(variant_class_list[0], "IGR")

    # check if any CGI variant classes that we normally filter are present in given variant_class_list
    filtered_cgi_var_classes = [
        var_class
        for var_class in variant_class_list
        if var_class in CGI_VARIANT_CLASS_FILTER
    ]
    if len(filtered_cgi_var_classes) > 0:
        # if filtered CGI variant classes are present then return the first one that can be directly mapped to
        # a standard variant classification
        for var_class in filtered_cgi_var_classes:
            if var_class in CGI_VARIANT_CLASS_MAP.keys():
                return CGI_VARIANT_CLASS_MAP[var_class]

    # Map var classes to standard variant classifications or use orig var class name if not found.
    # Variant classes that do not map directly include FRAMESHIFT, SPAN, UTR
    variant_class_candidates = list(
        map(lambda x: CGI_VARIANT_CLASS_MAP.get(x, x), variant_class_list)
    )

    # splice sites take precedence over other variant classifications
    if "Splice_Site" in variant_class_candidates:
        return "Splice_Site"

    # if variant type is INS/DEL or INSERT/DELETE/FRAMESHIFT/INSERT+/DELETE+ in input variant_class_list then
    # return either In_Frame_Ins/Del or Frame_Shift_Ins/Del
    indel_variant_classes = [
        var_class
        for var_class in variant_class_list
        if var_class in CGI_INDEL_VARIANT_CLASSES
    ]
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
        if (
            "Nonsense_Mutation" in variant_class_candidates
            or "Nonstop_Mutation" in variant_class_candidates
        ):
            return "Frame_Shift_" + variant_type.title()
        elif "Missense_Mutation" in variant_class_candidates:
            return "In_Frame_" + variant_type.title()

    # if variant type not INS/DEL then must be SNP/SNV/SUB
    # Return Nonsense_Mutation, Nonstop_Mutation, or Missense_Mutation with Nonsense/Nonstop taking precedence
    for var_class in ["Nonstop_Mutation", "Nonsense_Mutation", "Missense_Mutation"]:
        if var_class in variant_class_candidates:
            return var_class

    # if variant class is empty but variant type is SNP then return missense mutation
    if variant_type in ["SNP", "DNP", "TNP", "ONP"]:
        return "Missense_Mutation"

    # if variant class can't be resolve by this point then arbitrarily return first in list that can be mapped directly
    return variant_class_candidates[0]


def resolve_variant_classification(data, variant_type, ref_allele, alt_allele):
    """Resolves the variant classification."""

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
                in_frame_variant = len(alt_allele) % 3 == 0
            else:
                in_frame_variant = len(ref_allele) % 3 == 0
            if in_frame_variant:
                variant_class = "In_Frame_" + variant_type.title()
            else:
                variant_class = "Frame_Shift_" + variant_type.title()
        else:
            # let annotator figure it out from the VEP response
            variant_class = "Missense_Mutation"
    else:
        variant_class_list = re.split("[\\|, ]", variant_class)
        variant_class = resolve_complex_variant_classification(
            data, variant_class_list, variant_type
        )
    return variant_class


def resolve_vcf_variant_allele_data(vcf_data, maf_data):
    """Resolves the reference allele and tumor allele values."""
    # get ref allele and alt allele values

    ref_allele, alt_allele = resolve_vcf_allele(vcf_data)
    start_pos = resolve_start_position(vcf_data)
    variant_type = ""
    end_pos = ""
    variant_class = ""

    if (ref_allele != "N" and ref_allele != "") and alt_allele != "":
        # indels from vcf need to be shifted by one nucleotide and start position needs to be incremented by one
        if ref_allele[0] == alt_allele[0] and ref_allele != alt_allele:
            # shift ref allele and alt allele by one nucleotide, set as "-" if len == 1
            if len(ref_allele) == 1:
                ref_allele = "-"
            else:
                ref_allele = ref_allele[1:]
            if len(alt_allele) == 1:
                alt_allele = "-"
            else:
                alt_allele = alt_allele[1:]

            # fix start position value
            if start_pos != "":
                start_pos = str(int(start_pos) + 1)

        # resolve variant type, end position, and variant class
        variant_type = resolve_vcf_variant_type(ref_allele, alt_allele)
        end_pos = resolve_end_position(vcf_data, start_pos, variant_type, ref_allele)
        variant_class = resolve_variant_classification(
            vcf_data, variant_type, ref_allele, alt_allele
        )

    maf_data["Variant_Classification"] = variant_class
    maf_data["Variant_Type"] = variant_type
    maf_data["Reference_Allele"] = ref_allele
    maf_data["Tumor_Seq_Allele1"] = ref_allele
    maf_data["Tumor_Seq_Allele2"] = alt_allele
    maf_data["Start_Position"] = start_pos
    maf_data["End_Position"] = end_pos
    return maf_data


def resolve_vcf_matched_normal_allele_data(
    vcf_data, maf_data, matched_normal_sample_id
):
    """Resolves VCF matched normal seq allele data if normal genotype info is available."""

    # if normal genotype info is unavailable then assume normal seq alleles are ref/ref homozygous
    maf_data["Match_Norm_Seq_Allele1"] = maf_data["Reference_Allele"]
    maf_data["Match_Norm_Seq_Allele2"] = maf_data["Reference_Allele"]

    vcf_alleles = [vcf_data["REF"]]
    vcf_alleles.extend(vcf_data["ALT"].split(","))

    if matched_normal_sample_id in vcf_data.keys():
        normal_sample_format_data = vcf_data["MAPPED_NORMAL_FORMAT_DATA"]
        if "GT" in normal_sample_format_data.keys() and not is_missing_vcf_data_value(
            normal_sample_format_data["GT"]
        ):
            normal_sample_genotype_info = re.split(
                "[\/|]", normal_sample_format_data["GT"]
            )
            match_norm_seq_allele1 = ""
            match_norm_seq_allele2 = ""
            if len(normal_sample_genotype_info) == 1 and is_valid_integer(
                normal_sample_genotype_info[0]
            ):
                match_norm_seq_allele1 = vcf_alleles[
                    int(normal_sample_genotype_info[0])
                ]
                match_norm_seq_allele2 = vcf_alleles[
                    int(normal_sample_genotype_info[0])
                ]
            else:
                if is_valid_integer(normal_sample_genotype_info[0]):
                    match_norm_seq_allele1 = vcf_alleles[
                        int(normal_sample_genotype_info[0])
                    ]
                if is_valid_integer(normal_sample_genotype_info[1]):
                    match_norm_seq_allele2 = vcf_alleles[
                        int(normal_sample_genotype_info[1])
                    ]
            maf_data["Match_Norm_Seq_Allele1"] = match_norm_seq_allele1
            maf_data["Match_Norm_Seq_Allele2"] = match_norm_seq_allele2
    return maf_data


def create_maf_record_from_vcf(
    sample_id,
    center_name,
    sequence_source,
    vcf_data,
    is_germline_data,
    matched_normal_sample_id,
    tumor_sample_data_col,
):
    """
    Creates MAF record from VCF data.
    """
    # init maf record
    maf_data = init_maf_record()

    # set easy to resolve values
    maf_data["Tumor_Sample_Barcode"] = sample_id
    maf_data["Matched_Norm_Sample_Barcode"] = matched_normal_sample_id
    maf_data["Center"] = center_name
    maf_data["Hugo_Symbol"] = resolve_hugo_symbol(vcf_data)
    maf_data["Chromosome"] = resolve_chromosome(vcf_data)
    maf_data["Sequence_Source"] = resolve_sequence_source(vcf_data, sequence_source)
    maf_data["Verification_Status"] = DEFAULT_VERIFICATION_STATUS
    maf_data["FILTER"] = vcf_data.get("FILTER", "")
    if is_germline_data:
        maf_data["Mutation_Status"] = "GERMLINE"
    else:
        maf_data["Mutation_Status"] = DEFAULT_MUTATION_STATUS

    # resolve counts and variant allele data
    resolve_vcf_counts_data(
        vcf_data, maf_data, matched_normal_sample_id, tumor_sample_data_col
    )
    resolve_vcf_variant_allele_data(vcf_data, maf_data)
    resolve_vcf_matched_normal_allele_data(vcf_data, maf_data, matched_normal_sample_id)
    return maf_data


def is_missing_data_value(value):
    return value in NULL_OR_MISSING_VALUES


def print_warning(message):
    print(message)


def process_datum(value):
    """
    Returns a cleaned up data value.
    A null or 'empty' value is always returned as an empty string.
    """
    try:
        vfixed = value.strip().replace('"', "")
    except AttributeError:
        vfixed = ""
    if is_missing_data_value(vfixed):
        return ""
    return vfixed


def extract_vcf_format_info_data(
    vcf_data, tumor_sample_data_col, matched_normal_sample_id
):
    """
    Extract data from VCF 'INFO' and 'FORMAT' columns.

    The 'INFO' columns stores data as a string of semi-colon separated key-value pairs.
    This function formats the data in the 'INFO' column into an actual dictionary.

    Example:
            vcf_data["INFO"] = { key1;value1, key2:value2, ... }

    The 'FORMAT' column stores keys corresponding to values in the 'TUMOR' and 'NORMAL' column
    (if 'NORMAL' column is present in the VCF header). This function formats the data in the
    'FORMAT' column as a list of keys.

    Example:
            vcf_data["FORMAT"] = [ key1, key2, ... ]
    """

    parsed_vcf_info = {}
    for key_value_pair in vcf_data["INFO"].split(";"):
        key = process_datum(key_value_pair.split("=")[0])
        try:
            value = process_datum(key_value_pair.split("=")[1])
        except:
            value = ""
        parsed_vcf_info[key] = value

        # if 'FUNC' is a key in the 'INFO' data then promote key-value pairs to primary
        # key-value pairs in vcf_data dictionary - "FUNC" may contain important tumor or
        # normal sample data such as reads, alleles, etc.
        if key.startswith("FUNC"):
            # must replace single quotes with double quotes to avoid
            # errors when loading string as a json
            vcf_func_data = json.loads(str(value.replace("'", '"')))
            for elem in vcf_func_data[0].items():
                vcf_data[elem[0].upper().strip()] = process_datum(elem[1])

    vcf_data["INFO"] = parsed_vcf_info

    # the VCF "FORMAT" column stores keys that correspond to values in the "TUMOR"
    # column and "NORMAL" column (if "NORMAL" is present) and the keys are
    # separated by colons ":"
    parsed_vcf_format_keys = list(
        map(lambda v: process_datum(v), vcf_data["FORMAT"].split(":"))
    )
    vcf_data["FORMAT"] = parsed_vcf_format_keys

    # map the corresponding values in the VCF tumor column to the appropriate "FORMAT" keys
    # which is also separated by colons (":")
    tumor_vcf_format_values = list(
        map(lambda v: process_datum(v), vcf_data[tumor_sample_data_col].split(":"))
    )
    vcf_data["MAPPED_TUMOR_FORMAT_DATA"] = dict(
        zip(parsed_vcf_format_keys, tumor_vcf_format_values)
    )

    # repeat above to map the corresponding values in the VCF normal column (if present)
    if matched_normal_sample_id in vcf_data.keys():
        normal_vcf_format_values = list(
            map(
                lambda v: process_datum(v),
                vcf_data[matched_normal_sample_id].split(":"),
            )
        )
        vcf_data["MAPPED_NORMAL_FORMAT_DATA"] = dict(
            zip(parsed_vcf_format_keys, normal_vcf_format_values)
        )
    return vcf_data


def map_data_values_to_header(header, line):
    """
    Maps data values in 'line' to the column headers.
    """
    data = dict(zip(header, list(map(str.strip, line.split("\t")))))
    return data


def is_malformed_record(filename, record_index, line, header):
    """
    Verify whether the number of fields in the line matches the number of fields in the header.
    """
    split_line = list(map(str.strip, line.split("\t")))
    if len(split_line) != len(header):
        message = (
            "[line %s] Encountered record in file where the number of fields in the line does not match the expected number of fields in the header, skipping processing of file..."
            % (record_index + 1)
        )
        print_warning(message)
        update_problematic_report_for_file(filename, "ERROR", message)
        return True
    return False


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


def get_vcf_file_header(filename):
    vcf_file_header = []
    raw_header_line = ""
    for line in extract_file_data(filename):
        if line.startswith("#CHROM"):
            raw_header_line = line
            vcf_file_header = list(map(str.strip, line.replace("#", "").split("\t")))
            break
    if not vcf_file_header:
        print("Could not extract VCF header from file: %s - Exiting..." % (filename))
        sys.exit(2)
    return (vcf_file_header, raw_header_line)


def detect_file_encoding(filename):
    """
    Reads the first million bytes of a file
    to detect the type of encoding.
    """
    with open(filename, "rb") as data_file:
        encoding = detect(data_file.read(1000000))
        if encoding["encoding"] != "ascii":
            print(
                "[WARNING] detect_file_encoding(), Non-ASCII encoding detected in file %s, encoding type detected: %s - non-ASCII characters will be ignored and removed"
                % (filename, encoding["encoding"])
            )
        return encoding["encoding"]


def extract_file_data(filename):
    """
    Reads in file data with the appropriate encoding
    type and returns file data as ascii.
    """
    encoding_type = detect_file_encoding(filename)
    with open(filename, encoding=encoding_type) as data_file:
        filedata = []
        for line in data_file.readlines():
            # TODO: report on a line-by-line basis if there are non-ascii characters encountered
            filedata.append(line.encode("ascii", "ignore").decode("ascii"))
        return filedata


def get_vcf_sample_and_normal_ids(filename, tumor_id, normal_id):
    """
    Returns the sample name, tumor sample barcode column,
    and normal column/sample id from the VCF file.

    NOTE: The VCF fixed header may not always have "TUMOR" as the column name
    where the tumor sample data is stored. Even so, it can still be assumed that one of
    the columns following the non-case ID fixed VCF headers contains the tumor data
    (if the column header doesn't contain the keywords 'normal' or 'donor').

    The "NORMAL" column is optional and can have a different name other than "NORMAL".
    As with the "TUMOR" column, it can be assumed that if there is an additional column
    present with the "TUMOR" column then this is the column to refer to for extracting
    "NORMAL" sample data.
    """

    vcf_file_header = []
    vcf_meta_header = dict()
    for line in extract_file_data(filename):
        if line.startswith("#CHROM"):
            vcf_file_header = list(map(str.strip, line.replace("#", "").split("\t")))
            break
        elif line.startswith("##"):
            key, *val_splits = list(map(str.strip, line.replace("##", "").split("=")))
            vcf_meta_header[key] = "=".join(val_splits)
    # get the case id columns based on which columns in the header are not part of the fixed VCF header
    case_ids_cols = [
        col for col in vcf_file_header if col not in VCF_FIXED_HEADER_NON_CASE_IDS
    ]

    tumor_sample_data_col_name = None
    tumor_sample_id = None
    matched_normal_sample_id = None

    # if the tumor_id and normal_id are provided, then use them to find the sample data columns
    if tumor_id != "" and normal_id != "":
        if tumor_id in case_ids_cols:
            tumor_sample_data_col_name = tumor_sample_id = tumor_id
        if normal_id in case_ids_cols:
            matched_normal_sample_id = normal_id
        if (tumor_sample_data_col_name is None) or (matched_normal_sample_id is None):
            print(
                f"\nWARNING: The provided Tumor ID: '{tumor_id}' or Normal ID: '{normal_id}' were not found in the VCF file. \nParsing the names and sample ids in the old way (from column headers)..."
            )
        else:
            return (
                tumor_sample_id,
                tumor_sample_data_col_name,
                matched_normal_sample_id,
            )

    if (
        "normal_sample" in vcf_meta_header
        and vcf_meta_header["normal_sample"] not in case_ids_cols
    ):
        raise Exception(
            f"There is normal_sample={vcf_meta_header['normal_sample']} in the header, but no respective column found."
        )
    if (
        "tumor_sample" in vcf_meta_header
        and vcf_meta_header["tumor_sample"] not in case_ids_cols
    ):
        raise Exception(
            f"There is tumor_sample={vcf_meta_header['tumor_sample']} in the header, but no respective column found."
        )
    if (
        "tumor_sample" in vcf_meta_header and "normal_sample" not in vcf_meta_header
    ) or ("normal_sample" in vcf_meta_header and "tumor_sample" not in vcf_meta_header):
        raise Exception(
            f"The tumor_sample and normal_sample are expected together in the header. Found only one of them."
        )
    sample_columns_num = len(case_ids_cols)
    if sample_columns_num == 0:
        raise Exception("No sample column found")
    if sample_columns_num == 1 and any(
        id in case_ids_cols[0].upper() for id in NORMAL_IDENTIFIERS
    ):
        raise Exception(
            "There is only one sample column and it has NORMAL label. No tumor sample column present."
        )
    if sample_columns_num > 2:
        raise Exception(
            f"Expected max 2 sample columns for tumor and normal sample. But found {sample_columns_num} columns."
        )

    # if 'tumor_sample' is in vcf_meta_header (meaning 'normal_sample' will also be), use that:
    if "tumor_sample" in vcf_meta_header:
        tumor_sample_data_col_name = tumor_sample_id = vcf_meta_header["tumor_sample"]
        matched_normal_sample_id = vcf_meta_header["normal_sample"]
        return (tumor_sample_id, tumor_sample_data_col_name, matched_normal_sample_id)

    # otherwise, parse the names and sample ids in the old way (from column headers):
    if len(case_ids_cols) == 1 and case_ids_cols:
        tumor_sample_data_col_name = case_ids_cols[0]
        matched_normal_sample_id = "NORMAL"
    elif len(case_ids_cols) == 2:
        if not any(id in case_ids_cols[0].upper() for id in NORMAL_IDENTIFIERS):
            tumor_sample_data_col_name = case_ids_cols[0]
            matched_normal_sample_id = case_ids_cols[1]
        else:
            tumor_sample_data_col_name = case_ids_cols[1]
            matched_normal_sample_id = case_ids_cols[0]

    if tumor_sample_data_col_name.upper() == "TUMOR":
        tumor_sample_id = os.path.basename(filename).replace(".vcf", "")
    else:
        tumor_sample_id = tumor_sample_data_col_name

    return (tumor_sample_id, tumor_sample_data_col_name, matched_normal_sample_id)


def is_valid_chromosome(chromosome):
    return (
        is_valid_integer(chromosome) and int(chromosome) > 0 and int(chromosome) <= 24
    ) or (chromosome.upper() in ["X", "Y", "MT"])


def print_problematic_files_report():
    """Reports errors and warnings encountered while processing input files."""
    if not PROBLEMATIC_FILES_REPORT:
        print("No errors encountered during vcf2maf.py run...")
        return

    print(
        "\n ## ERRORS AND WARNINGS ENCOUNTERED ##\nEncountered %s files with errors and/or warnings during vcf2maf.py run...\n"
        % (len(PROBLEMATIC_FILES_REPORT.keys()))
    )
    print(
        "- Warnings are errors which are specific to record(s) within a mutation data file. Causes of such warning messages may be due to empty data values for certain key fields like allele depth fields, for example. The warning message will detail the field that the issue was encountered for."
    )
    print(
        "- Critical errors arise if data could not be loaded from the file, such as with data files containing a header but no records or a file where all of the lines are commented out.\n"
    )
    for filename, datafile_errors_and_warnings in PROBLEMATIC_FILES_REPORT.items():
        print("\n%s report:" % (filename))
        for problem_type, messages in sorted(datafile_errors_and_warnings.items()):
            print("\t%s" % (problem_type))
            for m in messages:
                print("\t\t%s" % (m))
        print("\n")


def update_problematic_report_for_file(filename, problem_type, message):
    # get current errors and warnings for given filename and update accordingly
    datafile_errors_and_warnings = PROBLEMATIC_FILES_REPORT.get(filename, {})

    # append warning/error message to appropriate list (either "ERROR" or "WARNING")
    messages_for_problem_type = datafile_errors_and_warnings.get(problem_type, [])
    messages_for_problem_type.append(message)

    # update global "PROBLEMATIC_FILES_REPORT" with latest info
    datafile_errors_and_warnings[problem_type] = messages_for_problem_type
    PROBLEMATIC_FILES_REPORT[filename] = datafile_errors_and_warnings


def capture_warnings_for_extracted_maf_record(filename, record_index, maf_record):
    """
    Adds warnings to 'PROBLEMATIC_FILES_REPORT' global.
    These are non-critical issues encountered during processing.
    """

    # check chromosome value
    if not is_valid_chromosome(maf_record["Chromosome"]):
        message = "[line %s], invalid chromosome value encountered (%s)" % (
            (record_index + 1),
            maf_record["Chromosome"],
        )
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
        (
            is_missing_data_value(maf_record["Reference_Allele"])
            or maf_record["Reference_Allele"] == "N"
        )
        and (
            is_missing_data_value(maf_record["Tumor_Seq_Allele1"])
            or maf_record["Tumor_Seq_Allele1"] == "N"
        )
        and (
            is_missing_data_value(maf_record["Tumor_Seq_Allele2"])
            or maf_record["Tumor_Seq_Allele2"] == "N"
        )
    ):
        message = (
            "[line %s], all allele fields are missing or invalid values ('Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2'): (%s, %s, %s)"
            % (
                (record_index + 1),
                maf_record["Reference_Allele"],
                maf_record["Tumor_Seq_Allele1"],
                maf_record["Tumor_Seq_Allele2"],
            )
        )
        update_problematic_report_for_file(filename, "WARNING", message)


def print_data_loading_summary(
    filename,
    records_loaded,
    rejected_vars_by_varclass,
    rejected_vars_by_mutstatus,
    rejected_vars_by_inc_data,
):
    print("\tTotal records loaded %s " % (records_loaded))
    if rejected_vars_by_varclass > 0:
        print(
            "\tTotal records filtered by variant classification: %s"
            % (rejected_vars_by_varclass)
        )
    if rejected_vars_by_mutstatus > 0:
        print(
            "\tTotal records rejected due to LOH or Wildtype Mutation Status: %s"
            % (rejected_vars_by_mutstatus)
        )
    if rejected_vars_by_inc_data > 0:
        print(
            "\tTotal records rejected due to incomplete allele data: %s"
            % (rejected_vars_by_inc_data)
        )


def extract_vcf_data_from_file(
    filename, center_name, sequence_source, tumor_id, normal_id
):
    """Load data from a VCF file."""
    print("Loading data from file: %s" % (filename))

    rejected_vars_by_varclass = 0
    rejected_vars_by_mutstatus = 0
    rejected_vars_by_inc_data = 0

    (
        sample_id,
        tumor_sample_data_col,
        matched_normal_sample_id,
    ) = get_vcf_sample_and_normal_ids(filename, tumor_id, normal_id)
    is_germline_data = "germline" in filename
    maf_data = []
    (vcf_file_header, raw_header_line) = get_vcf_file_header(filename)
    for record_index, line in enumerate(extract_file_data(filename)):
        if not is_valid_mutation_record_to_process(
            raw_header_line, vcf_file_header, line
        ):
            continue
        if is_malformed_record(filename, record_index, line, vcf_file_header):
            return None

        # map vcf data values to columns in file header
        vcf_data = map_data_values_to_header(vcf_file_header, line)
        # do some additional processing on the VCF columns "INFO", "FORMAT", the identified tumor_sample_data_col, and the NORMAL column (if NORMAL is present)
        extract_vcf_format_info_data(
            vcf_data, tumor_sample_data_col, matched_normal_sample_id
        )
        maf_record = create_maf_record_from_vcf(
            sample_id,
            center_name,
            sequence_source,
            vcf_data,
            is_germline_data,
            matched_normal_sample_id,
            tumor_sample_data_col,
        )

        # capture non-critical data issues as warnings (i.e., data issues that would prevent a successful annotation of single record)
        capture_warnings_for_extracted_maf_record(filename, record_index, maf_record)

        maf_data.append(maf_record)

    # if file passes all checks but we could not extract any data from the file then report as critical error
    if not maf_data:
        message = (
            "Data could not be extracted from file - output file will not be generated."
        )
        update_problematic_report_for_file(filename, "ERROR", message)
    else:
        print_data_loading_summary(
            filename,
            len(maf_data),
            rejected_vars_by_varclass,
            rejected_vars_by_mutstatus,
            rejected_vars_by_inc_data,
        )
    return maf_data


def write_standardized_mutation_file(maf_data, output_filename):
    """Writes standardized MAF data to output file."""
    output_file = open(output_filename, "w")
    output_file.write("\t".join(MAF_HEADER))
    for data in maf_data:
        formatted_data = list(map(lambda x: process_datum(data.get(x, "")), MAF_HEADER))
        output_file.write("\n" + "\t".join(formatted_data))
    output_file.write("\n")
    output_file.close()
    print("Standardized MAF written to: %s\n" % (output_filename))


def generate_maf_from_input_data(
    datafiles, output_directory, center_name, sequence_source, tumor_id, normal_id
):
    for filename in datafiles:
        head, tail = os.path.split(filename)
        extract_data = False
        if not "normal" in filename:
            extract_data = True
        if extract_data:
            output_filename = os.path.join(
                output_directory, os.path.splitext(tail)[0] + ".maf"
            )
            maf_data = extract_vcf_data_from_file(
                filename, center_name, sequence_source, tumor_id, normal_id
            )
            # only generate output file if data successfully extracted from input file
            if maf_data:
                write_standardized_mutation_file(maf_data, output_filename)


@click.command()
@click.option(
    "-i",
    "--input-data",
    required=True,
    help="A list of .vcf files or input data directories, separated by commas",
)
@click.option(
    "-o",
    "--output-directory",
    default="vcf2maf_output",
    help="output data directory [optional]",
)
@click.option(
    "-c",
    "--center",
    default="",
    help="name of center (standard MAF field = 'Center') [optional]",
)
@click.option(
    "-s",
    "--sequence-source",
    default="",
    help="Sequencing source (standard MAF field = 'Sequencing_Source'), e.g., WXS or WGS [optional]",
)
@click.option(
    "-t",
    "--tumor-id",
    default="",
    help="The ID of the tumor sample utilized in the genotype columns of the VCF file. [optional]",
)
@click.option(
    "-n",
    "--normal-id",
    default="",
    help="The ID of the normal sample utilized in the genotype columns of the VCF file. [optional]",
)

def main(input_data, output_directory, center, sequence_source, tumor_id, normal_id):

    # generate vcf files list
    datafiles = []
    for item in input_data.split(","):
        if os.path.isfile(item):
            datafiles.append(item)
        elif os.path.isdir(item):
            files = [
                os.path.join(item, x) for x in os.listdir(item) if x.endswith(".vcf")
            ]
            datafiles.extend(files)
        else:
            print("The path '%s' is ivalid" % item)

    # Create output directory
    if (
        output_directory == "vcf2maf_output"
    ):  # default folder if output path not specified
        output_directory = os.getcwd() + "/vcf2maf_output/"

    if not os.path.isdir(output_directory):
        try:
            os.mkdir(output_directory)
            print("created output folder : '{}' ".format(output_directory))
        except:
            print("The output path '%s' is ivalid" % output_directory)
            sys.exit(1)

    # Process .vcf files
    generate_maf_from_input_data(
        datafiles,
        output_directory,
        center,
        sequence_source.upper(),
        tumor_id,
        normal_id,
    )
    print_problematic_files_report()


if __name__ == "__main__":
    main()
