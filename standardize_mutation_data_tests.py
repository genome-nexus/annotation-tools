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
        self.assertEqual("Can't detect the tumor sample column among 3 columns", str(exc.exception))

if __name__=='__main__':
    unittest.main()
