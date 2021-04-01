#!/usr/bin/env python3

import re

VALID_SPECIAL_CHROMOSOMES = {"X", "Y"} # MT chromosome references not handled
MINIMAL_GENOMIC_OFFSET = 1
SPECIES_TO_SEX_CHROMOSOME_OFFSET = {"Homo sapiens" : "23", "Mus musculus" : "20", "Rattus norvegicus" : "21", "Bos taurus" : "30", "Danio rerio" : None, "Caenorhabditis elegans" : "6", "Saccharomyces cerevisiae" : None}
NUCLEOTIDE_SEQUENCE_RE_PATTERN  = "^([ATGC]+)$"

def is_string(value):
    if value == None:
        return False
    return isinstance(value, str)

def is_positive_integer(string_value):
    if not is_string(string_value):
        return False
    try:
        i = int(string_value)
        return i > 0
    except ValueError:
        return False

def offsets_are_valid(start, end):
    if not is_string(start) or not is_string(end):
        return False
    try:
        s = int(start)
        e = int(end)
        return s >= MINIMAL_GENOMIC_OFFSET and e >= MINIMAL_GENOMIC_OFFSET and s <= e
    except ValueError:
        return False

def chromosome_is_valid(chromosome):
    if not is_string(chromosome):
        return False
    trimmed_chromosome = chromosome.strip()
    if trimmed_chromosome in VALID_SPECIAL_CHROMOSOMES:
        return True
    return is_positive_integer(trimmed_chromosome)

def allele_is_valid(allele):
    return is_string(allele) # note : any string (even empty) might be acceptable in one of the provided tumor seq allele positions

def allele_is_empty(allele):
    stripped_allele = allele.strip()
    return stripped_allele == "NA" or stripped_allele == ""

def genomic_variant_has_proper_format(genomic_variant):
    if genomic_variant == None:
        return False
    if not isinstance(genomic_variant, list):
        return False
    if not len(genomic_variant) == 6:
        return False
    [ chromosome, start, end, ref_allele, seq_allele1, seq_allele2 ] = genomic_variant
    if not chromosome_is_valid(chromosome):
        return False
    if not offsets_are_valid(start, end):
        return False
    if not allele_is_valid(ref_allele) or not allele_is_valid(seq_allele1) or not allele_is_valid(seq_allele2):
        return False
    if allele_is_empty(seq_allele1) and allele_is_empty(seq_allele2):
        return False
    return True

def species_code_is_valid(species_code):
    if not is_string(species_code):
        return False
    return species_code in SPECIES_TO_SEX_CHROMOSOME_OFFSET

def increment_integer_string(string_value):
    return str(int(string_value) + 1) 

def normalize_chromosome(chromosome, sex_chromosome_offset):
    normalized_chromosome = chromosome.strip()
    if normalized_chromosome[0:3] == "chr":
        normalized_chromosome = normalized_chromosome[3:]
    if sex_chromosome_offset == None:
        return normalized_chromosome
    x_chromosome_offset = sex_chromosome_offset
    if normalized_chromosome == x_chromosome_offset:
        return "X"
    y_chromosome_offset = increment_integer_string(x_chromosome_offset)
    if normalized_chromosome == y_chromosome_offset:
        return "Y"
    return normalized_chromosome

def normalize_offset(offset):
    return offset.strip()
        
def normalize_allele(allele):
    return allele.strip().upper()
        
def normalize_genomic_variant(genomic_variant, sex_chromosome_offset):
    [ chromosome, start, end, ref_allele, seq_allele1, seq_allele2 ] = genomic_variant
    normalized_genomic_variant = []
    normalized_genomic_variant.append(normalize_chromosome(chromosome, sex_chromosome_offset))
    normalized_genomic_variant.append(normalize_offset(start))
    normalized_genomic_variant.append(normalize_offset(end))
    normalized_genomic_variant.append(normalize_allele(ref_allele))
    normalized_genomic_variant.append(normalize_allele(seq_allele1))
    normalized_genomic_variant.append(normalize_allele(seq_allele2))
    return normalized_genomic_variant

def seq_allele_ambiguity_exists(ref_allele, seq_allele1, seq_allele2):
    if allele_is_empty(seq_allele1) or seq_allele1 == ref_allele:
        return False
    if allele_is_empty(seq_allele2) or seq_allele2 == ref_allele:
        return False
    if seq_allele1 == "-" and re.match(NUCLEOTIDE_SEQUENCE_RE_PATTERN, seq_allele2):
        return True
    if seq_allele2 == "-" and re.match(NUCLEOTIDE_SEQUENCE_RE_PATTERN, seq_allele1):
        return True
    return False

def select_variant_allele(ref_allele, seq_allele1, seq_allele2):
    if allele_is_empty(seq_allele1) or seq_allele1 == ref_allele:
        return seq_allele2
    if allele_is_empty(seq_allele2) or seq_allele2 == ref_allele:
        return seq_allele1
    if seq_allele_ambiguity_exists(ref_allele, seq_allele1, seq_allele2):
        if seq_allele2 == "-":
            return seq_allele1
        return seq_allele2
    return seq_allele1

def genomic_variant_to_hgvs(genomic_variant, species_name):
    if not genomic_variant_has_proper_format(genomic_variant):
        return None
    if not species_code_is_valid(species_name):
        return None
    # normalize genomic location
    chromosome_offset = SPECIES_TO_SEX_CHROMOSOME_OFFSET[species_name];
    normalized_genomic_variant = normalize_genomic_variant(genomic_variant, chromosome_offset)
    [ chromosome, start, end, ref_allele, seq_allele1, seq_allele2 ] = normalized_genomic_variant
    # detect type and do conversoin
    selected_variant_allele = select_variant_allele(ref_allele, seq_allele1, seq_allele2)
    if ref_allele == "-" or len(ref_allele) == 0 or ref_allele == "NA" or ref_allele.find("--") != -1:
        return chromosome + ":g." + start + "_" + increment_integer_string(start) + "ins" + selected_variant_allele
    if selected_variant_allele == '-' or len(selected_variant_allele) == 0:
        return chromosome + ":g." + start + "_" + end + "del" + ref_allele;
    if len(ref_allele) > 1 or len(selected_variant_allele) > 1:
        return chromosome + ":g." + start + "_" + end + "del" + ref_allele + "ins" + selected_variant_allele;
    else:
        return chromosome + ":g." + start + ref_allele + ">" + selected_variant_allele;
    print("not handled yet : received " + str(genomic_variant))
