#! /usr/bin/env python

import unittest
import sys, os, gzip, tempfile, shutil, filecmp, tarfile
import nanomonsv 
from .check_download import check_download

class TestMain(unittest.TestCase):

    def setUp(self):

        # prepare reference genome
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        check_download("https://1000genomes.s3.us-east-1.amazonaws.com/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa",
                          cur_dir + "/resource/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa")
         
        self.parser = nanomonsv.arg_parser.create_parser()

    def test_parse_tumor(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        
        input_bam = cur_dir + "/resource/bam/test_tumor.bam"
        output_prefix = tmp_dir + "/test_tumor"

        nanomonsv_parse_args = ["parse", input_bam, output_prefix]
        print("nanomonsv " + ' '.join(nanomonsv_parse_args))

        args = self.parser.parse_args(nanomonsv_parse_args)
        nanomonsv.run.parse_main(args)

        with gzip.open(tmp_dir + "/test_tumor.rearrangement.sorted.bedpe.gz", 'rt') as hin: record_num1 = len(hin.readlines())
        with gzip.open(tmp_dir + "/test_tumor.deletion.sorted.bed.gz", 'rt') as hin: record_num2 = len(hin.readlines())
        with gzip.open(tmp_dir + "/test_tumor.insertion.sorted.bed.gz", 'rt') as hin: record_num3 = len(hin.readlines())

        print([record_num1, record_num2, record_num3])
        self.assertTrue(record_num1 == 96)
        self.assertTrue(record_num2 == 179)
        self.assertTrue(record_num3 == 165)

        shutil.rmtree(tmp_dir)
        #shutil.rmtree(cur_dir + "/data/test_tumor")
        #shutil.copytree(tmp_dir, cur_dir + "/data/test_tumor")

    """
    def test_parse_control(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        
        input_bam = cur_dir + "/resource/bam/test_ctrl.bam"
        output_prefix = tmp_dir + "/test_ctrl"

        nanomonsv_parse_args = ["parse", input_bam, output_prefix]
        print("nanomonsv " + ' '.join(nanomonsv_parse_args))

        args = self.parser.parse_args(nanomonsv_parse_args)
        nanomonsv.run.parse_main(args)

        with gzip.open(tmp_dir + "/test_ctrl.rearrangement.sorted.bedpe.gz", 'rt') as hin: record_num1 = len(hin.readlines())
        with gzip.open(tmp_dir + "/test_ctrl.deletion.sorted.bed.gz", 'rt') as hin: record_num2 = len(hin.readlines())
        with gzip.open(tmp_dir + "/test_ctrl.insertion.sorted.bed.gz", 'rt') as hin: record_num3 = len(hin.readlines())

        print([record_num1, record_num2, record_num3])
        self.assertTrue(record_num1 == 0)
        self.assertTrue(record_num2 == 153)
        self.assertTrue(record_num3 == 158)

        shutil.rmtree(tmp_dir)
        #shutil.rmtree(cur_dir + "/data/test_ctrl")
        #shutil.copytree(tmp_dir, cur_dir + "/data/test_ctrl")
    """

    def test_parse_tumor_cram(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        
        input_bam = cur_dir + "/resource/bam/test_tumor.cram"
        output_prefix = tmp_dir + "/test_tumor"
        ref_genome = cur_dir + "/resource/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"

        nanomonsv_parse_args = ["parse", input_bam, output_prefix, "--reference", ref_genome]
        print("nanomonsv " + ' '.join(nanomonsv_parse_args))

        args = self.parser.parse_args(nanomonsv_parse_args)
        nanomonsv.run.parse_main(args)

        with gzip.open(tmp_dir + "/test_tumor.rearrangement.sorted.bedpe.gz", 'rt') as hin: record_num1 = len(hin.readlines())
        with gzip.open(tmp_dir + "/test_tumor.deletion.sorted.bed.gz", 'rt') as hin: record_num2 = len(hin.readlines())
        with gzip.open(tmp_dir + "/test_tumor.insertion.sorted.bed.gz", 'rt') as hin: record_num3 = len(hin.readlines())

        print([record_num1, record_num2, record_num3])
        self.assertTrue(record_num1 == 96)
        self.assertTrue(record_num2 == 179)
        self.assertTrue(record_num3 == 165)

        shutil.rmtree(tmp_dir)

    def test_get1_1(self):
        # cram
        # with control
        # with --single_bnd --use_racon
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        tumor_bam = cur_dir + "/resource/bam/test_tumor.cram"
        ctrl_bam = cur_dir + "/resource/bam/test_ctrl.cram"
        ref_genome = cur_dir + "/resource/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
        ctrl_prefix = cur_dir + "/data/test_ctrl/test_ctrl"
        
        shutil.copytree(cur_dir + "/data/test_tumor", tmp_dir + "/test_tumor") 
        tumor_prefix = tmp_dir + "/test_tumor/test_tumor"

        nanomonsv_get_args = ["get", tumor_prefix, tumor_bam, ref_genome, 
                              "--control_prefix", ctrl_prefix, "--control_bam", ctrl_bam, "--single_bnd", "--use_racon"]

        print("nanomonsv " + ' '.join(nanomonsv_get_args))

        args = self.parser.parse_args(nanomonsv_get_args)
        nanomonsv.run.get_main(args)

        with open(tumor_prefix + ".nanomonsv.result.txt", 'r') as hin: record_num1 = len(hin.readlines())
        with open(tumor_prefix + ".nanomonsv.supporting_read.txt", 'r') as hin: record_num2 = len(hin.readlines()) 

        print([record_num1, record_num2])
        self.assertTrue(record_num1 == 12)
        self.assertTrue(record_num2 == 93) 

        shutil.rmtree(tmp_dir)

    def test_get1_2(self):
        # cram
        # with control
        # without --single_bnd
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        tumor_bam = cur_dir + "/resource/bam/test_tumor.cram"
        ctrl_bam = cur_dir + "/resource/bam/test_ctrl.cram"
        ref_genome = cur_dir + "/resource/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
        ctrl_prefix = cur_dir + "/data/test_ctrl/test_ctrl"
        
        shutil.copytree(cur_dir + "/data/test_tumor", tmp_dir + "/test_tumor") 
        tumor_prefix = tmp_dir + "/test_tumor/test_tumor"

        nanomonsv_get_args = ["get", tumor_prefix, tumor_bam, ref_genome, 
                              "--control_prefix", ctrl_prefix, "--control_bam", ctrl_bam]

        print("nanomonsv " + ' '.join(nanomonsv_get_args))

        args = self.parser.parse_args(nanomonsv_get_args)
        nanomonsv.run.get_main(args)

        with open(tumor_prefix + ".nanomonsv.result.txt", 'r') as hin: record_num1 = len(hin.readlines())
        with open(tumor_prefix + ".nanomonsv.supporting_read.txt", 'r') as hin: record_num2 = len(hin.readlines()) 

        print([record_num1, record_num2])
        self.assertTrue(record_num1 == 12)
        self.assertTrue(record_num2 == 92) 

        shutil.rmtree(tmp_dir)

    def test_get1_3(self):
        # cram
        # without control
        # with --single_bnd --use_racon 
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        tumor_bam = cur_dir + "/resource/bam/test_tumor.cram"
        ref_genome = cur_dir + "/resource/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
        
        shutil.copytree(cur_dir + "/data/test_tumor", tmp_dir + "/test_tumor") 
        tumor_prefix = tmp_dir + "/test_tumor/test_tumor"

        nanomonsv_get_args = ["get", tumor_prefix, tumor_bam, ref_genome, 
                              "--single_bnd", "--use_racon"]

        print("nanomonsv " + ' '.join(nanomonsv_get_args))

        args = self.parser.parse_args(nanomonsv_get_args)
        nanomonsv.run.get_main(args)

        with open(tumor_prefix + ".nanomonsv.result.txt", 'r') as hin: record_num1 = len(hin.readlines())
        with open(tumor_prefix + ".nanomonsv.supporting_read.txt", 'r') as hin: record_num2 = len(hin.readlines()) 

        print([record_num1, record_num2])
        self.assertTrue(record_num1 == 20)
        self.assertTrue(record_num2 == 167) 

        shutil.rmtree(tmp_dir)

    def test_get1_4(self):
        # cram
        # without control
        # without --single_bnd
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        tumor_bam = cur_dir + "/resource/bam/test_tumor.cram"
        ref_genome = cur_dir + "/resource/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
        
        shutil.copytree(cur_dir + "/data/test_tumor", tmp_dir + "/test_tumor") 
        tumor_prefix = tmp_dir + "/test_tumor/test_tumor"

        nanomonsv_get_args = ["get", tumor_prefix, tumor_bam, ref_genome, ]

        print("nanomonsv " + ' '.join(nanomonsv_get_args))

        args = self.parser.parse_args(nanomonsv_get_args)
        nanomonsv.run.get_main(args)

        with open(tumor_prefix + ".nanomonsv.result.txt", 'r') as hin: record_num1 = len(hin.readlines())
        with open(tumor_prefix + ".nanomonsv.supporting_read.txt", 'r') as hin: record_num2 = len(hin.readlines()) 

        print([record_num1, record_num2])
        self.assertTrue(record_num1 == 20)
        self.assertTrue(record_num2 == 168) 

        shutil.rmtree(tmp_dir)

    def test_get1_5(self):
        # bam
        # with control
        # with --single_bnd --use_racon
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        tumor_bam = cur_dir + "/resource/bam/test_tumor.bam"
        ctrl_bam = cur_dir + "/resource/bam/test_ctrl.bam"
        ref_genome = cur_dir + "/resource/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
        ctrl_prefix = cur_dir + "/data/test_ctrl/test_ctrl"
        
        shutil.copytree(cur_dir + "/data/test_tumor", tmp_dir + "/test_tumor") 
        tumor_prefix = tmp_dir + "/test_tumor/test_tumor"

        nanomonsv_get_args = ["get", tumor_prefix, tumor_bam, ref_genome, 
                              "--control_prefix", ctrl_prefix, "--control_bam", ctrl_bam, "--single_bnd", "--use_racon"]

        print("nanomonsv " + ' '.join(nanomonsv_get_args))

        args = self.parser.parse_args(nanomonsv_get_args)
        nanomonsv.run.get_main(args)

        with open(tumor_prefix + ".nanomonsv.result.txt", 'r') as hin: record_num1 = len(hin.readlines())
        with open(tumor_prefix + ".nanomonsv.supporting_read.txt", 'r') as hin: record_num2 = len(hin.readlines()) 

        print([record_num1, record_num2])
        self.assertTrue(record_num1 == 12)
        self.assertTrue(record_num2 == 93) 

        shutil.rmtree(tmp_dir)

    def test_get1_6(self):
        # bam
        # with control
        # without --single_bnd
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        tumor_bam = cur_dir + "/resource/bam/test_tumor.bam"
        ctrl_bam = cur_dir + "/resource/bam/test_ctrl.bam"
        ref_genome = cur_dir + "/resource/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
        ctrl_prefix = cur_dir + "/data/test_ctrl/test_ctrl"
        
        shutil.copytree(cur_dir + "/data/test_tumor", tmp_dir + "/test_tumor") 
        tumor_prefix = tmp_dir + "/test_tumor/test_tumor"

        nanomonsv_get_args = ["get", tumor_prefix, tumor_bam, ref_genome, 
                              "--control_prefix", ctrl_prefix, "--control_bam", ctrl_bam]

        print("nanomonsv " + ' '.join(nanomonsv_get_args))

        args = self.parser.parse_args(nanomonsv_get_args)
        nanomonsv.run.get_main(args)

        with open(tumor_prefix + ".nanomonsv.result.txt", 'r') as hin: record_num1 = len(hin.readlines())
        with open(tumor_prefix + ".nanomonsv.supporting_read.txt", 'r') as hin: record_num2 = len(hin.readlines()) 

        print([record_num1, record_num2])
        self.assertTrue(record_num1 == 12)
        self.assertTrue(record_num2 == 92) 

        shutil.rmtree(tmp_dir)

    def test_get1_7(self):
        # bam
        # without control
        # with --single_bnd --use_racon
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        tumor_bam = cur_dir + "/resource/bam/test_tumor.bam"
        ref_genome = cur_dir + "/resource/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
        
        shutil.copytree(cur_dir + "/data/test_tumor", tmp_dir + "/test_tumor") 
        tumor_prefix = tmp_dir + "/test_tumor/test_tumor"

        nanomonsv_get_args = ["get", tumor_prefix, tumor_bam, ref_genome, 
                              "--single_bnd", "--use_racon"]

        print("nanomonsv " + ' '.join(nanomonsv_get_args))

        args = self.parser.parse_args(nanomonsv_get_args)
        nanomonsv.run.get_main(args)

        with open(tumor_prefix + ".nanomonsv.result.txt", 'r') as hin: record_num1 = len(hin.readlines())
        with open(tumor_prefix + ".nanomonsv.supporting_read.txt", 'r') as hin: record_num2 = len(hin.readlines()) 

        print([record_num1, record_num2])
        self.assertTrue(record_num1 == 20)
        self.assertTrue(record_num2 == 167) 

        shutil.rmtree(tmp_dir)

    def test_get1_8(self):
        # bam
        # without control
        # without --single_bnd
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        tumor_bam = cur_dir + "/resource/bam/test_tumor.bam"
        ref_genome = cur_dir + "/resource/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
        
        shutil.copytree(cur_dir + "/data/test_tumor", tmp_dir + "/test_tumor") 
        tumor_prefix = tmp_dir + "/test_tumor/test_tumor"

        nanomonsv_get_args = ["get", tumor_prefix, tumor_bam, ref_genome, ]

        print("nanomonsv " + ' '.join(nanomonsv_get_args))

        args = self.parser.parse_args(nanomonsv_get_args)
        nanomonsv.run.get_main(args)

        with open(tumor_prefix + ".nanomonsv.result.txt", 'r') as hin: record_num1 = len(hin.readlines())
        with open(tumor_prefix + ".nanomonsv.supporting_read.txt", 'r') as hin: record_num2 = len(hin.readlines()) 

        print([record_num1, record_num2])
        self.assertTrue(record_num1 == 20)
        self.assertTrue(record_num2 == 168) 

        shutil.rmtree(tmp_dir)

    def test_get2(self):
        # with --control_panel_prefix
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        tumor_bam = cur_dir + "/resource/bam/test_tumor.bam"
        ctrl_bam = cur_dir + "/resource/bam/test_ctrl.bam"
        ref_genome = cur_dir + "/resource/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
        ctrl_prefix = cur_dir + "/data/test_ctrl/test_ctrl"
        #ctrl_panel_prefix = "/home/aiokada/sandbox/nanomonsv-tutorial/control_panel/hprc_year1_data_freeze_nanopore_minimap2_2_24_merge_control/hprc_year1_data_freeze_nanopore_minimap2_2_24_merge_control"
        ctrl_panel_prefix = cur_dir + "/resource/control_panel/control_panel"

        shutil.copytree(cur_dir + "/data/test_tumor", tmp_dir + "/test_tumor") 
        tumor_prefix = tmp_dir + "/test_tumor/test_tumor"

        nanomonsv_get_args = ["get", tumor_prefix, tumor_bam, ref_genome, 
                              "--control_prefix", ctrl_prefix, "--control_bam", ctrl_bam,
                              "--control_panel_prefix", ctrl_panel_prefix]

        print("nanomonsv " + ' '.join(nanomonsv_get_args))

        args = self.parser.parse_args(nanomonsv_get_args)
        nanomonsv.run.get_main(args)

        with open(tumor_prefix + ".nanomonsv.result.txt", 'r') as hin: record_num1 = len(hin.readlines())
        with open(tumor_prefix + ".nanomonsv.supporting_read.txt", 'r') as hin: record_num2 = len(hin.readlines()) 

        print([record_num1, record_num2])
        self.assertTrue(record_num1 == 11)
        self.assertTrue(record_num2 == 84) 

        shutil.rmtree(tmp_dir)

    def test_get3(self):
        # with --use_racon
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        tumor_bam = cur_dir + "/resource/bam/test_tumor.bam"
        ctrl_bam = cur_dir + "/resource/bam/test_ctrl.bam"
        ref_genome = cur_dir + "/resource/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
        ctrl_prefix = cur_dir + "/data/test_ctrl/test_ctrl"
        
        shutil.copytree(cur_dir + "/data/test_tumor", tmp_dir + "/test_tumor") 
        tumor_prefix = tmp_dir + "/test_tumor/test_tumor"

        nanomonsv_get_args = ["get", tumor_prefix, tumor_bam, ref_genome, 
                              "--control_prefix", ctrl_prefix, "--control_bam", ctrl_bam, "--use_racon"]

        print("nanomonsv " + ' '.join(nanomonsv_get_args))

        args = self.parser.parse_args(nanomonsv_get_args)
        nanomonsv.run.get_main(args)

        with open(tumor_prefix + ".nanomonsv.result.txt", 'r') as hin: record_num1 = len(hin.readlines())
        with open(tumor_prefix + ".nanomonsv.supporting_read.txt", 'r') as hin: record_num2 = len(hin.readlines()) 

        print([record_num1, record_num2])
        self.assertTrue(record_num1 == 12)
        self.assertTrue(record_num2 == 93) 

        shutil.rmtree(tmp_dir)

    def test_get4_1(self):
        # with --processes 4
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        tumor_bam = cur_dir + "/resource/bam/test_tumor.bam"
        ctrl_bam = cur_dir + "/resource/bam/test_ctrl.bam"
        ref_genome = cur_dir + "/resource/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
        ctrl_prefix = cur_dir + "/data/test_ctrl/test_ctrl"

        shutil.copytree(cur_dir + "/data/test_tumor", tmp_dir + "/test_tumor") 
        tumor_prefix = tmp_dir + "/test_tumor/test_tumor"

        nanomonsv_get_args = ["get", tumor_prefix, tumor_bam, ref_genome, 
                              "--control_prefix", ctrl_prefix, "--control_bam", ctrl_bam,
                              "--processes", "4"]

        print("nanomonsv " + ' '.join(nanomonsv_get_args))

        args = self.parser.parse_args(nanomonsv_get_args)
        nanomonsv.run.get_main(args)

        with open(tumor_prefix + ".nanomonsv.result.txt", 'r') as hin: record_num1 = len(hin.readlines())
        with open(tumor_prefix + ".nanomonsv.supporting_read.txt", 'r') as hin: record_num2 = len(hin.readlines()) 

        print([record_num1, record_num2])
        self.assertTrue(record_num1 == 12)
        self.assertTrue(record_num2 == 92) 

        shutil.rmtree(tmp_dir)

    def test_get4_2(self):
        # with --processes 8, empty file
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        tumor_bam = cur_dir + "/resource/bam/test_tumor.bam"
        ctrl_bam = cur_dir + "/resource/bam/test_ctrl.bam"
        ref_genome = cur_dir + "/resource/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
        ctrl_prefix = cur_dir + "/data/test_ctrl/test_ctrl"

        shutil.copytree(cur_dir + "/data/test_tumor_n1", tmp_dir + "/test_tumor") 
        tumor_prefix = tmp_dir + "/test_tumor/test_tumor"

        nanomonsv_get_args = ["get", tumor_prefix, tumor_bam, ref_genome, 
                              "--control_prefix", ctrl_prefix, "--control_bam", ctrl_bam,
                              "--processes", "4"]

        print("nanomonsv " + ' '.join(nanomonsv_get_args))

        args = self.parser.parse_args(nanomonsv_get_args)
        nanomonsv.run.get_main(args)

        with open(tumor_prefix + ".nanomonsv.result.txt", 'r') as hin: record_num1 = len(hin.readlines())
        with open(tumor_prefix + ".nanomonsv.supporting_read.txt", 'r') as hin: record_num2 = len(hin.readlines()) 

        print([record_num1, record_num2])
        self.assertTrue(record_num1 == 1)
        self.assertTrue(record_num2 == 0) 

        shutil.rmtree(tmp_dir)

    def test_get4_3(self):
        # with --processes 4
        # switch tumor - control
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        tumor_bam = cur_dir + "/resource/bam/test_ctrl.bam"
        ctrl_bam = cur_dir + "/resource/bam/test_tumor.bam"
        ref_genome = cur_dir + "/resource/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
        ctrl_prefix = cur_dir + "/data/test_tumor/test_tumor"

        shutil.copytree(cur_dir + "/data/test_tumor", tmp_dir + "/test_tumor") 
        tumor_prefix = tmp_dir + "/test_tumor/test_tumor"

        nanomonsv_get_args = ["get", tumor_prefix, tumor_bam, ref_genome, 
                              "--control_prefix", ctrl_prefix, "--control_bam", ctrl_bam,
                              "--processes", "4"]

        print("nanomonsv " + ' '.join(nanomonsv_get_args))

        args = self.parser.parse_args(nanomonsv_get_args)
        nanomonsv.run.get_main(args)

        with open(tumor_prefix + ".nanomonsv.result.txt", 'r') as hin: record_num1 = len(hin.readlines())
        with open(tumor_prefix + ".nanomonsv.supporting_read.txt", 'r') as hin: record_num2 = len(hin.readlines()) 

        print([record_num1, record_num2])
        self.assertTrue(record_num1 == 1)
        self.assertTrue(record_num2 == 0) 

        shutil.rmtree(tmp_dir)

    def test_validate(self):
        
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        sv_list = cur_dir + "/resource/validate/test_tumor_sv_list.txt"
        tumor_bam = cur_dir + "/resource/bam/test_tumor.bam"
        output_file = tmp_dir + "/test_tumor.validate.txt"
        ctrl_bam = cur_dir + "/resource/bam/test_ctrl.bam"
        ref_genome = cur_dir + "/resource/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"

        nanomonsv_validate_args = ["validate", sv_list, tumor_bam, output_file, ref_genome, "--control_bam", ctrl_bam]
        print("nanomonsv " + ' '.join(nanomonsv_validate_args))
 
        args = self.parser.parse_args(nanomonsv_validate_args)
        nanomonsv.run.validate_main(args)

        self.assertTrue(filecmp.cmp(output_file, sv_list, shallow=False))
        shutil.rmtree(tmp_dir)
        #shutil.copyfile(output_file, sv_list)

if __name__ == "__main__":
   unittest.main()

