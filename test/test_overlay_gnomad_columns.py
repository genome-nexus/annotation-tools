#!/usr/bin/env python3
# run tests from annotation-tools directory with:
#     python -m unittest discover


import unittest
import variant_notation_converter

class TestOverlayGnomadColumns(unittest.TestCase):

    def test_all_cases(self):
        test_case = []
        # SNP tests
        #    template
        test_case.append([["6", "123456", "123456", "A", "A", "G"], "6:g.123456A>G"])
        #    improper input value count
        test_case.append([["6", "123456", "123456", "A", "A", "G", "extra"], None])
        test_case.append([["6", "123456", "123456", "A", "A"], None])
        test_case.append([["6"], None])
        test_case.append([[], None])
        test_case.append([None, None])
        #    invalid field values
        test_case.append([["", "123456", "123456", "A", "A", "G"], None])
        test_case.append([["NA", "123456", "123456", "A", "A", "G"], None])
        test_case.append([[None, "123456", "123456", "A", "A", "G"], None])
        test_case.append([["-", "123456", "123456", "A", "A", "G"], None])
        test_case.append([[" ", "123456", "123456", "A", "A", "G"], None])
        test_case.append([["Z", "123456", "123456", "A", "A", "G"], None])
                        # note : for non-human genomes, chromosome numbers go above 22. accept all positive natural numbers as valid chromosome identifiers
        test_case.append([["0", "123456", "123456", "A", "A", "G"], None])
        test_case.append([["-1", "123456", "123456", "A", "A", "G"], None])
        test_case.append([["6", "", "123456", "A", "A", "G"], None])
        test_case.append([["6", "NA", "123456", "A", "A", "G"], None])
        test_case.append([["6", None, "123456", "A", "A", "G"], None])
        test_case.append([["6", "-", "123456", "A", "A", "G"], None])
        test_case.append([["6", " ", "123456", "A", "A", "G"], None])
        test_case.append([["6", "A", "123456", "A", "A", "G"], None])
        test_case.append([["6", "-3", "123456", "A", "A", "G"], None])
        test_case.append([["6", "0", "123456", "A", "A", "G"], None])
        test_case.append([["6", "123456", "", "A", "A", "G"], None])
        test_case.append([["6", "123456", "NA", "A", "A", "G"], None])
        test_case.append([["6", "123456", None, "A", "A", "G"], None])
        test_case.append([["6", "123456", "-", "A", "A", "G"], None])
        test_case.append([["6", "123456", " ", "A", "A", "G"], None])
        test_case.append([["6", "123456", "A", "A", "A", "G"], None])
        test_case.append([["6", "123456", "-3", "A", "A", "G"], None])
        test_case.append([["6", "123456", "0", "A", "A", "G"], None])
        test_case.append([["", "123456", "123456", "", "A", "G"], None])
        test_case.append([["", "123456", "123456", "NA", "A", "G"], None])
        test_case.append([["", "123456", "123456", None, "A", "G"], None])
        test_case.append([["", "123456", "123456", "-", "A", "G"], None])
        test_case.append([["", "123456", "123456", " ", "A", "G"], None])
        test_case.append([["", "123456", "123456", "A", "", ""], None])
        test_case.append([["", "123456", "123456", "A", "InvalidValue", "InvalidValue"], None]) # 'NA' is a valid nucleotide string
        test_case.append([["", "123456", "123456", "A", None, None], None])
        #    valid cases with sex chromosomes
        test_case.append([["X", "22125504", "22125504", "A", "A", "C"], "X:g.22125504A>C"])
        test_case.append([["23", "22125504", "22125504", "A", "A", "C"], "X:g.22125504A>C"])
        test_case.append([["Y", "22125504", "22125504", "A", "A", "C"], "Y:g.22125504A>C"])
        test_case.append([["24", "22125504", "22125504", "A", "A", "C"], "Y:g.22125504A>C"])
        #    valid case with mitochondrial chromosome (not implemented)
        #test_case.append([["MT", "8993", "8993", "T", "C", "C"], "m.8993T>C"])
        #    valid insertions
        test_case.append([["17", "7579470", "7579471", "-", "C", "C"], "17:g.7579470_7579471insC"])
        test_case.append([["17", "7579590", "7579591", "-", "CC", "CC"], "17:g.7579590_7579591insCC"])
        test_case.append([["17", "7577111", "7577112", "-", "CCCCCC", "CCCCCC"], "17:g.7577111_7577112insCCCCCC"])
        #    invalid insertions (commented out ... in our current logic the end offset is ignored for constructing insert variants
        #### test_case.append([["17", "7579470", "7579472", "-", "C", "C"], None]) # excessive offset
        #### test_case.append([["17", "7579471", "7579470", "-", "C", "C"], None]) # descending offset
        #### test_case.append([["17", "7579471", "7579471", "-", "C", "C"], None]) # empty offset
        #### #    valid deletions
        test_case.append([["17", "7578276", "7578276", "A", "-", "-"], "17:g.7578276_7578276delA"])
        test_case.append([["17", "7578276", "7578277", "AG", "-", "-"], "17:g.7578276_7578277delAG"])
        #    invalid deletions
        test_case.append([["17", "7578276", "7578275", "A", "-", "-"], None]) # descending offset
        ####test_case.append([["17", "7578276", "7578277", "AGGGGGGGGGGGGGGGGGGG", "-", "-"], None]) # non-matching ref allele ... not an error -- our code does not consider the ref allele string for deletions .. only the offsets matter
        #    valid indels (include dnp and tnp)
        test_case.append([["17", "7579547", "7579551", "GGGGA", "CCC", "CCC"], "17:g.7579547_7579551delGGGGAinsCCC"]) # delete 5, insert 3
        test_case.append([["17", "7579547", "7579548", "GG", "CC", "CC"], "17:g.7579547_7579548delGGinsCC"]) # delete 2, insert 2
        test_case.append([["17", "7579547", "7579549", "GGG", "CCC", "CCC"], "17:g.7579547_7579549delGGGinsCCC"]) # delete 3, insert 3
        #    invalid indels
        test_case.append([["17", "7579548", "7579547", "GG", "CC", "CC"], None]) # descending offset
        ####test_case.append([["17", "7579548", "7579548", "GG", "CC", "CC"], None]) # non-matching ref allele ... not an error -- our code does not consider the ref allele string for deletions .. only the offsets matter
        #    valid snps
        test_case.append([["17", "7577121", "7577121", "G", "C", "C"], "17:g.7577121G>C"])
        test_case.append([["17", "7577121", "7577121", "AGGGGGGGGGGGGGGGGGGG", "C", "C"], "17:g.7577121_7577121delAGGGGGGGGGGGGGGGGGGGinsC"]) # surprisingly, ensembl VEP acutally does the right thing for this HGVS string
        #    invalid snps
        #    (no examples, because our logic will ignore elements of the variant when necessary to construct an hgvs query which is valid (such as ignoring the 'end' or 'ref_allele' fields when that removes inconsistency - see the previous example where "1 NT" ref_allele 'AGGGGGGGGGGGGGGGGGGG' is ignored to make a valid SNP
        #    test correct selection of ambiguous tumor seq allele
        test_case.append([["17", "7577121", "7577121", "G", "G", "C"], "17:g.7577121G>C"]) # seq allele 1 matches ref allele (SNP)
        test_case.append([["17", "7577121", "7577121", "G", "C", "G"], "17:g.7577121G>C"]) # seq allele 2 matches ref allele (SNP)
        test_case.append([["17", "7577121", "7577121", "G", "-", "C"], "17:g.7577121G>C"]) # seq allele 1 is missing (SNP)
        test_case.append([["17", "7577121", "7577121", "G", "C", "-"], "17:g.7577121G>C"]) # seq allele 2 is missing (SNP)
        test_case.append([["17", "7577121", "7577121", "G", "C", "T"], "17:g.7577121G>C"]) # two alternatives (SNP) NOTE : According to our current code .. this is not ambiguous (even though it actually is) .. so pick #1
        test_case.append([["17", "7577121", "7577121", "G", "T", "C"], "17:g.7577121G>T"]) # two alternatives (SNP) NOTE : According to our current code .. this is not ambiguous (even though it actually is) .. so pick #1
        test_case.append([["17", "7579470", "7579471", "-", "-", "C"], "17:g.7579470_7579471insC"]) # seq allele 1 matches ref allele (INS)
        test_case.append([["17", "7579470", "7579471", "-", "C", "-"], "17:g.7579470_7579471insC"]) # seq allele 2 matches ref allele (INS)
        test_case.append([["17", "7579470", "7579471", "-", "C", "T"], "17:g.7579470_7579471insC"]) # two alternatives (INS) NOTE : According to our current code .. this is not ambiguous (even though it actually is) .. so pick #1
        test_case.append([["17", "7579470", "7579471", "-", "T", "C"], "17:g.7579470_7579471insT"]) # two alternatives (INS) NOTE : According to our current code .. this is not ambiguous (even though it actually is) .. so pick #1
        test_case.append([["17", "7578276", "7578277", "AG", "AG", "-"], "17:g.7578276_7578277delAG"]) # seq allele 1 matches ref allele (DEL)
        test_case.append([["17", "7578276", "7578277", "AG", "-", "AG"], "17:g.7578276_7578277delAG"]) # seq allele 2 matches ref allele (DEL)
                                                                                                    # note : if alternatives are delete v. other(snp/indel/ins) the other is preferred
        test_case.append([["17", "7579547", "7579549", "GGG", "GGG", "CCC"], "17:g.7579547_7579549delGGGinsCCC"]) # seq allele 1 matches ref allele (INDEL)
        test_case.append([["17", "7579547", "7579549", "GGG", "CCC", "GGG"], "17:g.7579547_7579549delGGGinsCCC"]) # seq allele 2 matches ref allele (INDEL)
        test_case.append([["17", "7579547", "7579549", "GGG", "-", "CCC"], "17:g.7579547_7579549delGGGinsCCC"]) # seq allele 1 is missing (INDEL)
        test_case.append([["17", "7579547", "7579549", "GGG", "CCC", "-"], "17:g.7579547_7579549delGGGinsCCC"]) # seq allele 2 is missing (INDEL)
        test_case.append([["17", "7579547", "7579549", "GGG", "CCC", "TTT"], "17:g.7579547_7579549delGGGinsCCC"]) # two alternatives (INDEL)NOTE : According to our current code .. this is not ambiguous (even though it actually is) .. so pick #1
        test_case.append([["17", "7579547", "7579549", "GGG", "TTT", "CCC"], "17:g.7579547_7579549delGGGinsTTT"]) # two alternatives (INDEL)NOTE : According to our current code .. this is not ambiguous (even though it actually is) .. so pick #1
        for [case, expected_response] in test_case:
            response = str(variant_notation_converter.genomic_variant_to_hgvs(case, "Homo sapiens"))
            self.assertEqual(str(response), str(expected_response))

if __name__ == '__main__':
    unittest.main()
