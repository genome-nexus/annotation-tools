#!/usr/bin/env python3

import unittest
import tempfile
from standardize_mutation_data import extract_vcf_data_from_file

class StandardizeMutationDataTests(unittest.TestCase):

    def test_extract_vcf_data_from_file_no_samples(self):
        _, vcf = tempfile.mkstemp()
        with open(vcf, 'w') as f:
           f.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n"
            "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\n"
           )
        with self.assertRaises(Exception) as exc:
            extract_vcf_data_from_file(vcf, 'center name 1', 'sequence source 1')
        self.assertEqual("No sample column found", str(exc.exception))

    def test_extract_vcf_data_from_file_1_sample(self):
        _, vcf = tempfile.mkstemp()
        with open(vcf, 'w') as f:
           f.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\n"
           ) 
        maf_data = extract_vcf_data_from_file(vcf, 'center name 1', 'sequence source 1')
        self.assertEqual(1, len(maf_data)) 
        maf_row = maf_data[0]
        self.assertEqual('S1', maf_row['Tumor_Sample_Barcode']) 
        self.assertEqual('NORMAL', maf_row['Matched_Norm_Sample_Barcode']) 

    def test_extract_vcf_data_from_file_2_samples(self):
        _, vcf = tempfile.mkstemp()
        with open(vcf, 'w') as f:
           f.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
            "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\n"
           ) 
        maf_data = extract_vcf_data_from_file(vcf, 'center name 1', 'sequence source 1')
        self.assertEqual(1, len(maf_data)) 
        maf_row = maf_data[0]
        self.assertEqual('S1', maf_row['Tumor_Sample_Barcode']) 
        self.assertEqual('S2', maf_row['Matched_Norm_Sample_Barcode']) 

    def test_extract_vcf_data_from_file_tumor(self):
        _, vcf = tempfile.mkstemp()
        with open(vcf, 'w') as f:
           f.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\n"
            "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\n"
           )
        maf_data = extract_vcf_data_from_file(vcf, 'center name 1', 'sequence source 1')
        self.assertEqual(1, len(maf_data)) 
        maf_row = maf_data[0]
        self.assertTrue(maf_row['Tumor_Sample_Barcode']) 
        self.assertNotEqual('S1', maf_row['Tumor_Sample_Barcode']) 
        self.assertEqual('NORMAL', maf_row['Matched_Norm_Sample_Barcode']) 

    def test_extract_vcf_data_from_file_3_samples(self):
        _, vcf = tempfile.mkstemp()
        with open(vcf, 'w') as f:
           f.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n"
            "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1|0:48:8:51,51\n"
           )
        with self.assertRaises(Exception) as exc:
            extract_vcf_data_from_file(vcf, 'center name 1', 'sequence source 1')
        self.assertEqual("Expected max 2 sample columns for tumor and normal sample. But found 3 columns.", str(exc.exception))

    def test_extract_vcf_data_from_file_tumor_normal_swap(self):
        # this test is just showing that the name of the columns does not matter...it will be the order
        # that matters in this case: first column is seen as Tumor_Sample_Barcode, second column as Matched_Norm_Sample_Barcode:
        _, vcf = tempfile.mkstemp()
        with open(vcf, 'w') as f:
           f.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR\n"
            "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\n"
           )
        maf_data = extract_vcf_data_from_file(vcf, 'center name 1', 'sequence source 1')
        self.assertEqual(1, len(maf_data))
        maf_row = maf_data[0]
        self.assertEqual('NORMAL', maf_row['Tumor_Sample_Barcode'])
        self.assertEqual('TUMOR', maf_row['Matched_Norm_Sample_Barcode'])

    def test_extract_vcf_data_from_file_1_normal_sample(self):
        _, vcf = tempfile.mkstemp()
        with open(vcf, 'w') as f:
           f.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\n"
            "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\n"
           )
        with self.assertRaises(Exception) as exc:
            extract_vcf_data_from_file(vcf, 'center name 1', 'sequence source 1')
        self.assertEqual("There is only one sample column and it has NORMAL label. No tumor sample column present.", str(exc.exception))

    def test_extract_vcf_data_from_file_2_samples_specified_in_header(self):
        _, vcf = tempfile.mkstemp()
        with open(vcf, 'w') as f:
           f.write(
            '##normal_sample=S1\n'
            '##tumor_sample=S2\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
            "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\n"
           )
        maf_data = extract_vcf_data_from_file(vcf, 'center name 1', 'sequence source 1')
        self.assertEqual(1, len(maf_data))
        maf_row = maf_data[0]
        self.assertEqual('S2', maf_row['Tumor_Sample_Barcode'])
        self.assertEqual('S1', maf_row['Matched_Norm_Sample_Barcode'])

    def test_extract_vcf_data_from_file_2_samples_specified_in_header_reversed_order(self):
        # same as above, but with column order reversed in CHROM header line...should still give the same result, just
        # to show that it is really reading from header and not from column order:
        _, vcf = tempfile.mkstemp()
        with open(vcf, 'w') as f:
           f.write(
            '##normal_sample=S1\n'
            '##tumor_sample=S2\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS2\tS1\n"
            "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\n"
           )
        maf_data = extract_vcf_data_from_file(vcf, 'center name 1', 'sequence source 1')
        self.assertEqual(1, len(maf_data))
        maf_row = maf_data[0]
        self.assertEqual('S2', maf_row['Tumor_Sample_Barcode'])
        self.assertEqual('S1', maf_row['Matched_Norm_Sample_Barcode'])

    def test_extract_vcf_data_from_file_normal_sample_refers_to_non_existing_column(self):
        _, vcf = tempfile.mkstemp()
        with open(vcf, 'w') as f:
           f.write(
            '##normal_sample=S1\n'
            '##tumor_sample=S2\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS2\n"
            "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\n"
           )
        with self.assertRaises(Exception) as exc:
            extract_vcf_data_from_file(vcf, 'center name 1', 'sequence source 1')
        self.assertEqual("There is normal_sample=S1 in the header, but no respective column found.", str(exc.exception))

    def test_extract_vcf_data_from_file_tumor_sample_refers_to_non_existing_column(self):
        _, vcf = tempfile.mkstemp()
        with open(vcf, 'w') as f:
           f.write(
            '##normal_sample=S1\n'
            '##tumor_sample=S2\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\n"
           )
        with self.assertRaises(Exception) as exc:
            extract_vcf_data_from_file(vcf, 'center name 1', 'sequence source 1')
        self.assertEqual("There is tumor_sample=S2 in the header, but no respective column found.", str(exc.exception))

    def test_extract_vcf_data_from_file_header_samples_have_precedence_over_the_rest_of_strategies(self):
        _, vcf = tempfile.mkstemp()
        with open(vcf, 'w') as f:
           f.write(
            '##normal_sample=TUMOR\n'
            '##tumor_sample=NORMAL\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL\n"
            "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\n"
           )
        maf_data = extract_vcf_data_from_file(vcf, 'center name 1', 'sequence source 1')
        self.assertEqual(1, len(maf_data))
        maf_row = maf_data[0]
        self.assertEqual('NORMAL', maf_row['Tumor_Sample_Barcode'])
        self.assertEqual('TUMOR', maf_row['Matched_Norm_Sample_Barcode'])

    def test_extract_vcf_data_from_file_multiple_equals_in_header(self):
        _, vcf = tempfile.mkstemp()
        with open(vcf, 'w') as f:
           f.write(
            '##tumor_sample=A=B\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA=B\n"
            "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\n"
           )
        with self.assertRaises(Exception) as exc:
            extract_vcf_data_from_file(vcf, 'center name 1', 'sequence source 1')
        self.assertIn("The tumor_sample and normal_sample are expected together", str(exc.exception))

if __name__=='__main__':
    unittest.main()
