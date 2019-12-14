#! /usr/bin/env python

import unittest
import sys, os, gzip, tempfile, shutil, filecmp, tarfile
import nanomonsv 
from .check_download import check_download

class TestMain(unittest.TestCase):

    def setUp(self):

        def extract_tar_gz(input_tar_gz_file, out_path):
            tar = tarfile.open(input_tar_gz_file)
            tar.extractall(out_path)
            tar.close()

        # prepare reference genome
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        check_download("https://storage.googleapis.com/friend1ws_package_data/common/GRCh37.fa", \
                       cur_dir + "/resource/reference_genome/GRCh37.fa")

        self.parser = nanomonsv.arg_parser.create_parser()


    def test_parse(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        
        input_bam = cur_dir + "/resource/bam/test_tumor.bam"
        output_prefix = tmp_dir + "/test_tumor"

        nanomonsv_parse_args = ["parse", input_bam, output_prefix]
        print("nanomonsv " + ' '.join(nanomonsv_parse_args))

        args = self.parser.parse_args(nanomonsv_parse_args)
        nanomonsv.run.parse_main(args)

        with gzip.open(tmp_dir + "/test_tumor.rearrangement.sorted.bedpe.gz", 'rt') as hin: record_num = len(hin.readlines())
        self.assertTrue(record_num == 26)

        with gzip.open(tmp_dir + "/test_tumor.deletion.sorted.bed.gz", 'rt') as hin: record_num = len(hin.readlines())
        self.assertTrue(record_num == 205)

        with gzip.open(tmp_dir + "/test_tumor.insertion.sorted.bed.gz", 'rt') as hin: record_num = len(hin.readlines())
        self.assertTrue(record_num == 35)

        shutil.rmtree(tmp_dir)


    def test_get(self):
        
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        tumor_bam = cur_dir + "/resource/bam/test_tumor.bam"
        ctrl_bam = cur_dir + "/resource/bam/test_ctrl.bam"
        ref_genome = cur_dir + "/resource/reference_genome/GRCh37.fa"
        ctrl_prefix = cur_dir + "/data/test_ctrl/test_ctrl"
        
        shutil.copytree(cur_dir + "/data/test_tumor", tmp_dir + "/test_tumor") 
        tumor_prefix_dst = tmp_dir + "/test_tumor/test_tumor"
        tumor_prefix_src = cur_dir + "/data/test_tumor/test_tumor"

        nanomonsv_get_args = ["get", tumor_prefix_dst, tumor_bam, ref_genome, 
                              "--control_prefix", ctrl_prefix, "--control_bam", ctrl_bam]

        print("nanomonsv " + ' '.join(nanomonsv_get_args))

        args = self.parser.parse_args(nanomonsv_get_args)
        nanomonsv.run.get_main(args)

        with open(tumor_prefix_dst + ".nanomonsv.result.txt", 'r') as hin: record_num = len(hin.readlines())
        self.assertTrue(record_num == 4)

        with open(tumor_prefix_dst + ".nanomonsv.supporting_read.txt", 'r') as hin: record_num = len(hin.readlines()) 
        self.assertTrue(record_num == 32) 
   
        shutil.rmtree(tmp_dir)


    def test_validate(self):
        
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        sv_list = cur_dir + "/data/test_tumor/test_tumor.nanomonsv.result.txt"
        tumor_bam = cur_dir + "/resource/bam/test_tumor.bam"
        output_file = tmp_dir + "/test_tumor.validate.txt"
        ctrl_bam = cur_dir + "/resource/bam/test_ctrl.bam"
        ref_genome = cur_dir + "/resource/reference_genome/GRCh37.fa"

        nanomonsv_validate_args = ["validate", sv_list, tumor_bam, output_file, ref_genome, "--control_bam", ctrl_bam]
        print("nanomonsv " + ' '.join(nanomonsv_validate_args))
 
        args = self.parser.parse_args(nanomonsv_validate_args)
        nanomonsv.run.validate_main(args)

        self.assertTrue(filecmp.cmp(output_file, sv_list, shallow=False))
        shutil.rmtree(tmp_dir)


if __name__ == "__main__":
   unittest.main()
