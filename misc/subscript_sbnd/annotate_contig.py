#! /usr/bin/env python3

# Usage: python integrate_sbnd.py <output_prefix> <reference_fasta>
# Generates:
#   {output_prefix}.nanomonsv.bwa.txt   — BWA alignment of sbnd contigs
#   {output_prefix}.nanomonsv.rmsk.txt  — RepeatMasker annotation of sbnd contigs
#   {output_prefix}.class.txt           — Contig classification (Simple_repeat, Satellite, Plain_SV, L1_Mediated_Del, Complex)

import sys, subprocess, os, shutil, csv
from collections import namedtuple
import pysam


def proc_rmsk(input_file, output_file):

    with open(input_file, 'r') as hin, open(output_file, 'w') as hout:
        print("Contig_ID\tContig_Len\tQuery_Align_Start\tQuery_Align_End\tQuery_Align_Strand\tRepeat_NAME\tRepeat_Class\tRepeat_Align_Start\tRepeat_Align_End", file = hout)
        for line in hin:
            F = line.rstrip('\n').split()
            if len(F) <= 1 or F[0] in ["SW", "score"]: continue
            if F[0] == "There": continue

            contigs = F[4].split(',')
            contig_id, contig_len = ','.join(contigs[:-1]), contigs[-1]
            if F[8] == 'C': F[8] = '-'
            if F[8] == '+':
                print(f"{contig_id}\t{contig_len}\t{F[5]}\t{F[6]}\t{F[8]}\t{F[9]}\t{F[10]}\t{F[11]}\t{F[12]}", file = hout)
            else:
                print(f"{contig_id}\t{contig_len}\t{F[5]}\t{F[6]}\t{F[8]}\t{F[9]}\t{F[10]}\t{F[13]}\t{F[12]}", file = hout)


def proc_sam(input_sam, output_file):

    samfile = pysam.AlignmentFile(input_sam, 'r')

    hout = open(output_file, 'w')
    print("Contig_ID\tContig_Len\tQuery_Align_Start\tQuery_Align_End\tQuery_Align_Strand\tTarget_Align_Chromosome\tTarget_Align_Start\tTarget_Align_End\tMapping_Quality", file = hout)
    for read in samfile.fetch():

        if read.is_unmapped: continue
        if read.is_secondary: continue

        query_name = read.query_name
        query_strand = '-' if read.is_reverse else '+'
        reference_name = read.reference_name
        reference_start = read.reference_start + 1
        reference_end = read.reference_end
        mapping_quality = read.mapping_quality
        query_length = read.infer_read_length()

        cigartuples = read.cigartuples
        left_hard_clipping_size, right_hard_clipping_size = 0, 0
        if cigartuples[0][0] == 5: left_hard_clipping_size = cigartuples[0][1]
        if cigartuples[-1][0] == 5: right_hard_clipping_size = cigartuples[-1][1]

        if not read.is_supplementary:
            if query_strand == '+':
                query_start = read.query_alignment_start + 1
                query_end = read.query_alignment_end
            else:
                query_start = query_length - read.query_alignment_end + 1
                query_end = query_length - read.query_alignment_start
        else:
            if query_strand == '+':
                query_start = left_hard_clipping_size + 1
                query_end = query_length - right_hard_clipping_size
            else:
                query_start = right_hard_clipping_size + 1
                query_end = query_length - left_hard_clipping_size

        contigs = query_name.split(',')
        contig_id, contig_len = ','.join(contigs[:-1]), contigs[-1]
        print(f"{contig_id}\t{contig_len}\t{query_start}\t{query_end}\t{query_strand}\t{reference_name}\t{reference_start}\t{reference_end}\t{mapping_quality}", file = hout)

    hout.close()


Rmsk_record = namedtuple("Rmsk_record", "QA_Start QA_End QA_Strand RName RClass RA_Start RA_End")

Bwa_record = namedtuple("Bwa_record", "QA_Start QA_End QA_Strand TA_Chr TA_Start TA_End MQ")

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                  'W': 'W', 'S': 'S', 'M': 'K', 'K': 'M', 'R': 'Y', 'Y': 'R',
                  'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N'}

    return("".join(complement.get(base, base) for base in reversed(seq)))

class Contig_info:

    def __init__(self, contig_id, contig):

        bchr, bpos, bstrand, bid = contig_id.split(',')

        self.bp_chr = bchr
        self.bp_pos = int(bpos)
        self.bp_strand = bstrand
        self.bp_id = bid
        self.contig = contig
        self.contig_len = len(contig)
        self.contig_class = "None"
        self.rmsk_info = []
        self.bwa_info = []

        self.chr1 = bchr
        self.chr2 = None
        self.pos1 = int(bpos)
        self.pos2 = None
        self.dir1 = bstrand
        self.dir2 = None
        self.inseq = None

        self.repeat_size_ratio_thres = 0.8

    def sv_key(self):

        if self.chr2 is not None:

            if self.inseq == '': self.inseq = '---'

            if self.chr1 > self.chr2 or (self.chr1 == self.chr2 and self.pos1 > self.pos2):
                self.chr1, self.chr2 = self.chr2, self.chr1
                self.pos1, self.pos2 = self.pos2, self.pos1
                self.dir1, self.dir2 = self.dir2, self.dir1
                self.inseq = reverse_complement(self.inseq)

            return(f'{self.chr1},{self.pos1},{self.dir1},{self.chr2},{self.pos2},{self.dir2},{self.inseq}')
        else:
            return("---")


    def check_classification(self):

        self.check_simple_satellite()
        if self.contig_class == "None":
            self.check_plain_SV()
        if self.contig_class == "None":
            self.check_L1_mediated_deletion()
        if self.contig_class == "None":
            self.contig_class = "Complex"

    def check_simple_satellite(self):

        simple_repeat_size = 0
        satellite_size = 0
        satellite_alpha_size = 0
        for rec in self.rmsk_info:
            if rec.RClass == "Simple_repeat":
                simple_repeat_size = simple_repeat_size + rec.QA_End - rec.QA_Start
            if rec.RClass == "Satellite":
                satellite_size = satellite_size + rec.QA_End - rec.QA_Start
            if rec.RClass == "Satellite/centr":
                satellite_alpha_size = satellite_alpha_size + rec.QA_End - rec.QA_Start

        if simple_repeat_size / self.contig_len > self.repeat_size_ratio_thres:
            self.contig_class = "Simple_repeat"
        if satellite_size / self.contig_len > self.repeat_size_ratio_thres:
            self.contig_class = "Satellite"
        if satellite_alpha_size / self.contig_len > self.repeat_size_ratio_thres:
            self.contig_class = "Satellite/centr"


    def check_plain_SV(self):

        for rec in self.bwa_info:
            if rec.QA_Start < 100 and rec.QA_End - rec.QA_Start >= 2000 and rec.MQ >= 40:
                self.contig_class = "Plain_SV"

                self.chr2 = rec.TA_Chr
                if rec.QA_Strand == '+':
                    self.pos2 = rec.TA_Start
                    self.dir2 = '-'
                else:
                    self.pos2 = rec.TA_End
                    self.dir2 = '+'

                self.inseq = self.contig[:(rec.QA_Start - 1)]

    def check_L1_mediated_deletion(self):

        self.bwa_info = sorted(self.bwa_info, key = lambda x: x.QA_Start)

        L1HS_segments = []
        for rec in self.rmsk_info:
            if rec.RName in ["L1HS", "L1P1", "L1PA2"]: L1HS_segments.append((rec.QA_Start, rec.QA_End))

        L1HS_intersection_size = 0
        query_size = 0
        for i in range(0, min(len(self.bwa_info) - 1, 2)):
            rec = self.bwa_info[i]
            query_size = query_size + rec.QA_End - rec.QA_Start
            for lseg in L1HS_segments:
                if lseg[1] >= rec.QA_Start and lseg[0] <= rec.QA_End:
                    L1HS_intersection_size = L1HS_intersection_size + min(lseg[1], rec.QA_End) - max(lseg[0], rec.QA_Start)

            if float(L1HS_intersection_size) / query_size > 0.8:
                if self.bwa_info[i + 1].QA_End - self.bwa_info[i + 1].QA_Start >= 2000 and self.bwa_info[i + 1].MQ >= 40:
                    self.contig_class = "L1_Mediated_Del"

                    self.chr2 = self.bwa_info[i + 1].TA_Chr
                    if self.bwa_info[i + 1].QA_Strand == '+':
                        self.pos2 = self.bwa_info[i + 1].TA_Start
                        self.dir2 = '-'
                    else:
                        self.pos2 = self.bwa_info[i + 1].TA_End
                        self.dir2 = '+'

                    self.inseq = self.contig[:(self.bwa_info[i + 1].QA_Start - 1)]

                elif self.bwa_info[i + 1].QA_End - self.bwa_info[i + 1].QA_Start >= 2000 and self.bwa_info[i + 1].TA_Chr == self.bp_chr:
                    if (self.bwa_info[i + 1].QA_Strand == self.bp_strand == '+' and -100 < self.bwa_info[i + 1].TA_Start - self.bp_pos < 10000) or \
                        (self.bwa_info[i + 1].QA_Strand == self.bp_strand == '-' and -10000 < self.bwa_info[i + 1].TA_End - self.bp_pos < 100):
                        self.contig_class = "L1_Mediated_Del"

                        self.chr2 = self.bwa_info[i + 1].TA_Chr
                        if self.bwa_info[i + 1].QA_Strand == '+':
                            self.pos2 = self.bwa_info[i + 1].TA_Start
                            self.dir2 = '-'
                        else:
                            self.pos2 = self.bwa_info[i + 1].TA_End
                            self.dir2 = '+'

                        self.inseq = self.contig[:(self.bwa_info[i + 1].QA_Start - 1)]


if __name__ == "__main__":

    output_prefix = sys.argv[1]
    reference_fasta = sys.argv[2]

    # Generate FASTA from sbnd contigs
    with open(output_prefix + ".nanomonsv.sbnd.result.txt", 'r') as hin, open(output_prefix + ".tmp.fasta", 'w') as hout:
        for F in csv.DictReader(hin, delimiter = "\t"):
            print(f'>{F["Chr_1"]},{F["Pos_1"]},{F["Dir_1"]},{F["SV_ID"]},{len(F["Contig"])}\n{F["Contig"]}', file = hout)

    # Run RepeatMasker
    os.makedirs(output_prefix + ".tmp.rmsk", exist_ok = True)
    subprocess.check_call(["RepeatMasker", "-species", "human", output_prefix + ".tmp.fasta", "-dir", output_prefix + ".tmp.rmsk"])

    boutput_prefix = os.path.basename(output_prefix)
    proc_rmsk(output_prefix + ".tmp.rmsk/" + boutput_prefix + ".tmp.fasta" + ".out", output_prefix + ".nanomonsv.rmsk.txt")

    # Run BWA
    with open(output_prefix + ".tmp.bwa.sam", 'w') as hout:
        subprocess.check_call(["bwa", "mem", "-h", "200", reference_fasta, output_prefix + ".tmp.fasta"], stdout = hout)

    proc_sam(output_prefix + ".tmp.bwa.sam", output_prefix + ".nanomonsv.bwa.txt")

    # Clean up temporary files
    shutil.rmtree(output_prefix + ".tmp.rmsk")
    os.remove(output_prefix + ".tmp.fasta")
    os.remove(output_prefix + ".tmp.bwa.sam")

    # Load contig sequences
    contig_id2contig = {}
    with open(output_prefix + ".nanomonsv.sbnd.result.txt", 'r') as hin:
        for F in csv.DictReader(hin, delimiter = '\t'):
            contig_id = f'{F["Chr_1"]},{F["Pos_1"]},{F["Dir_1"]},{F["SV_ID"]}'
            contig_id2contig[contig_id] = F["Contig"]

    # Load rmsk annotations
    contig_id2contig_info = {}
    with open(output_prefix + ".nanomonsv.rmsk.txt", 'r') as hin:
        for F in csv.DictReader(hin, delimiter = '\t'):
            contig_id = F["Contig_ID"]
            if contig_id not in contig_id2contig_info:
                contig_id2contig_info[contig_id] = Contig_info(contig_id, contig_id2contig[contig_id])
            contig_id2contig_info[contig_id].rmsk_info.append(
                Rmsk_record(int(F["Query_Align_Start"]), int(F["Query_Align_End"]), F["Query_Align_Strand"], F["Repeat_NAME"], F["Repeat_Class"], int(F["Repeat_Align_Start"]), int(F["Repeat_Align_End"]))
            )

    # Load bwa annotations
    with open(output_prefix + ".nanomonsv.bwa.txt", 'r') as hin:
        for F in csv.DictReader(hin, delimiter = '\t'):
            contig_id = F["Contig_ID"]
            if contig_id not in contig_id2contig_info:
                contig_id2contig_info[contig_id] = Contig_info(contig_id, contig_id2contig[contig_id])
            contig_id2contig_info[contig_id].bwa_info.append(
                Bwa_record(int(F["Query_Align_Start"]), int(F["Query_Align_End"]), F["Query_Align_Strand"], F["Target_Align_Chromosome"], int(F["Target_Align_Start"]), int(F["Target_Align_End"]), int(F["Mapping_Quality"]))
            )

    # Classify contigs and output class.txt
    with open(output_prefix + ".class.txt", 'w') as hout:
        for contig_id in contig_id2contig_info:
            contig_info = contig_id2contig_info[contig_id]
            contig_info.check_classification()
            print(f'{contig_id}\t{contig_info.contig_class}\t{contig_info.sv_key()}', file = hout)
