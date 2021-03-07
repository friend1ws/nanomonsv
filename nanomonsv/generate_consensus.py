#! /usr/bin/env python3

import sys, os, subprocess, shutil, statistics, logging
from collections import Counter
import pysam
import parasail

from .logger import get_logger

logger = get_logger(__name__)


class Consensus_generator(object):

    def __init__(self, output_file, use_racon, debug):

        self.tmp_dir = output_file + ".tmp_dir"
        os.makedirs(self.tmp_dir, exist_ok = True)

        self.hout = open(output_file, 'w')
        self.temp_key = None 
        self.readid2is_bp = {}
        self.temp_support_read_file_h = None
        self.use_racon = False if use_racon == False else True
        self.debug = debug
        self.parasail_error = []

        self.start_margin = 120
        self.min_inclusion_ratio = 0.9
        self.min_inclusion_count = 3

    def __del__(self):
        self.hout.close()
        if len(self.parasail_error) > 0:
            logger.debug(f"Alignment by parasail failed in {self.parasail_error}")

        if not self.debug:
            shutil.rmtree(self.tmp_dir)

    def initialize(self, temp_key):
        if self.temp_support_read_file_h is not None: self.temp_support_read_file_h.close()
        self.temp_key = temp_key
        self.temp_support_read_file_h = open(self.tmp_dir + '/' + self.temp_key + ".supporting_read.fa", 'w')
        self.readid2is_bp = {}

    def add_support_read_seq(self, readid, seq, size):
        # skip if the seq corresponding to the input readid is already printed 
        if readid not in self.readid2is_bp:
            print(f">{readid}\n{seq}", file = self.temp_support_read_file_h)
            if size in ['-', '+']: # indication of break point type supporting read
                self.readid2is_bp[readid] = True
            else:
                self.readid2is_bp[readid] = False

    # function used for obtaining consensus sequence from multiple alignment result
    def get_consensus_from_mafft_result(self, input_file):

        id2seq = {}
        with open(input_file, 'r') as hin:
            for line in hin:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    tid = line
                    id2seq[tid] = ''
                else:
                    id2seq[tid] = id2seq[tid] + line

        ind2bases = {}
        for tid in id2seq:
            seq = id2seq[tid]
            for i in range(len(seq)):
                if i not in ind2bases: ind2bases[i] = []
                ind2bases[i].append(seq[i])

        seq_len = len(list(ind2bases))
        consensus = ''
        for i in range(seq_len):
            mycounter = Counter(ind2bases[i] )
            consensus = consensus + mycounter.most_common()[0][0]

        consensus = consensus.replace('-', '').upper()

        return(consensus)


    def print_consensus_mafft(self):

        # if self.read_num >= 3: return

        hout_m = open(self.tmp_dir + '/' + self.temp_key + ".mafft_result.fa", 'w')
        subprocess.check_call(["mafft", self.tmp_dir + '/' + self.temp_key + ".supporting_read.fa"], stdout = hout_m, stderr = subprocess.DEVNULL)
        hout_m.close()
                
        tconsensus = self.get_consensus_from_mafft_result(self.tmp_dir + '/' + self.temp_key + ".mafft_result.fa")
        print(self.temp_key + '\t' + tconsensus, file = self.hout)

        
    def generate_paf_file(self, query_fasta, target_fasta, output_file):

        user_matrix = parasail.matrix_create("ACGT", 2, -2)
        paf_rec_count = 0

        with open(target_fasta, 'r') as hin:
            for line in hin:
                if line.startswith('>'): 
                    tid = line.rstrip('\n').split(' ')[0].lstrip('>')
                else:
                    tseq = line.rstrip('\n')

        with open(query_fasta, 'r') as hin, open(output_file, 'w') as hout:
            for line in hin:
                if line.startswith('>'):
                    qid = line.rstrip('\n').lstrip('>')
                else:
                    qseq = line.rstrip('\n')
                    
                    res = parasail.ssw(qseq, tseq, 3, 1, user_matrix)
                    if res is not None:
                        print(f"{qid}\t{len(qseq)}\t{res.read_begin1}\t{res.read_end1}\t+\t" +
                            f"{tid}\t{len(tseq)}\t{res.ref_begin1}\t{res.ref_end1}\t*\t*\t60", file = hout)
                        paf_rec_count = paf_rec_count + 1
                    else:
                        self.parasail_error.append((qid, tid))

        return(paf_rec_count)


    def print_consensus_racon(self):

        target_flag = False
        with open(self.tmp_dir + '/' + self.temp_key + ".supporting_read.fa", 'r') as hin, \
            open(self.tmp_dir + '/' + self.temp_key + ".tmp.seg.first.fa", 'w') as hout: 
            for line in hin:
                tid = line.rstrip('\n').lstrip('>')
                tseq = hin.readline().rstrip('\n')

                if not self.readid2is_bp[tid]:
                    print(f'>{tid}\n{tseq}', file = hout)
                    target_flag = True
                    break

        if target_flag == False:
            logger.debug(f"Template sequence could not be found for for {self.temp_key}")
            return

        paf_rec_count = self.generate_paf_file(self.tmp_dir + '/' + self.temp_key + ".supporting_read.fa",
            self.tmp_dir + '/' + self.temp_key + ".tmp.seg.first.fa",
            self.tmp_dir + '/' + self.temp_key + ".parasail.paf")
                        
        if paf_rec_count < 3: 
            logger.debug(f"Not enough PAF records for the first round consensus generation for {self.temp_key}")
            return

        try: 
            with open(self.tmp_dir + '/' + self.temp_key + ".racon1.fa", 'w') as hout:
                subprocess.check_call(["racon", "-u", 
                    self.tmp_dir + '/' + self.temp_key + ".supporting_read.fa",
                    self.tmp_dir + '/' + self.temp_key + ".parasail.paf",
                    self.tmp_dir + '/' + self.temp_key + ".tmp.seg.first.fa"],
                    stdout = hout, stderr = subprocess.DEVNULL)
        except Exception as e:
            logger.warning(f'{e}')
            return

        with open(self.tmp_dir + '/' + self.temp_key + ".racon1.fa", 'r') as hin, \
            open(self.tmp_dir + "/" + self.temp_key + ".racon1.mod.fa", 'w') as hout:
            tid = hin.readline().lstrip('>').split(' ')[0]
            print(f">{tid}_2nd", file = hout)
            tseq = hin.readline().rstrip('\n')
            print(tseq, file = hout)

        paf_rec_count = self.generate_paf_file(self.tmp_dir + '/' + self.temp_key + ".supporting_read.fa",
            self.tmp_dir + '/' + self.temp_key + ".racon1.mod.fa",
            self.tmp_dir + '/' + self.temp_key + ".parasail2.paf")
        
        if paf_rec_count < 3: 
            logger.debug(f"Not enough PAF records for the second round consensus generation for {self.temp_key}")
            return

        try:
            with open(self.tmp_dir + '/' + self.temp_key + ".racon2.fa", 'w') as hout:
                subprocess.check_call(["racon", "-u",
                    self.tmp_dir + '/' + self.temp_key + ".supporting_read.fa",
                    self.tmp_dir + '/' + self.temp_key + ".parasail2.paf",
                    self.tmp_dir + '/' + self.temp_key + ".racon1.mod.fa"],
                    stdout = hout, stderr = subprocess.DEVNULL)
        except Exception as e: 
            logger.warning(f'{e}')
            return

        with open(self.tmp_dir + "/" + self.temp_key + ".racon2.fa") as hin:
            header = hin.readline()
            tconsensus = hin.readline().rstrip('\n')
            print(f"{self.temp_key}\t{tconsensus}", file = self.hout)


    def print_consensus(self):

        if self.temp_support_read_file_h is not None: self.temp_support_read_file_h.close()
        if self.use_racon:
            self.print_consensus_racon()
        else:
            self.print_consensus_mafft()


    def print_consensus_sbnd(self):

        if self.temp_support_read_file_h is not None: self.temp_support_read_file_h.close()

        with open(self.tmp_dir + '/' + self.temp_key + "_ava_minimap2.paf", 'w') as hout:
            subprocess.check_call(["minimap2", "-x", "ava-ont", 
                self.tmp_dir + '/' + self.temp_key + ".supporting_read.fa",
                self.tmp_dir + '/' + self.temp_key + ".supporting_read.fa"],
                stderr = subprocess.DEVNULL, stdout = hout)

        readid2inclusion_count = {}
        with open(self.tmp_dir + '/' + self.temp_key + "_ava_minimap2.paf") as hin:
            for line in hin:
                row = line.rstrip('\n').split('\t')
                if int(row[2]) <= self.start_margin and (float(row[3]) - float(row[2])) / float(row[1]) >= self.min_inclusion_ratio:
                    if row[5] not in readid2inclusion_count: readid2inclusion_count[row[5]] = 0
                    readid2inclusion_count[row[5]] = readid2inclusion_count[row[5]] + 1
                if int(row[7]) <= self.start_margin and (float(row[8]) - float(row[7])) / float(row[6]) >= self.min_inclusion_ratio:
                    if row[0] not in readid2inclusion_count: readid2inclusion_count[row[0]] = 0
                    readid2inclusion_count[row[0]] = readid2inclusion_count[row[0]] + 1

        if len(readid2inclusion_count) == 0: return
        readid_max = max(readid2inclusion_count, key = readid2inclusion_count.get)
        if readid2inclusion_count[readid_max] < self.min_inclusion_count: return

        # first round racon
        with open(self.tmp_dir + '/' + self.temp_key + ".supporting_read.fa", 'r') as hin, \
            open(self.tmp_dir + '/' + self.temp_key + "_ref.fa", 'w') as hout:
            seq_read_count = 0
            while True:
                tid = hin.readline().rstrip('\n').lstrip('>')
                tseq = hin.readline().rstrip('\n')
                if tid == readid_max:
                    print(f">{tid}\n{tseq}", file = hout)
                    break
                seq_read_count = seq_read_count + 1
                if tid == '' or seq_read_count > 10000: 
                    logger.warning(f"Something inconsistent happend when choosing template reads for {self.temp_key}")
                    return

        with open(self.tmp_dir + '/' + self.temp_key + "_ova_minimap2.paf", 'w') as hout:
            subprocess.check_call(["minimap2", "-x", "map-ont", self.tmp_dir + '/' + self.temp_key + "_ref.fa",
                self.tmp_dir + '/' + self.temp_key + ".supporting_read.fa"], 
                stderr = subprocess.DEVNULL, stdout = hout)

        paf_rec_count = 0
        with open(self.tmp_dir + '/' + self.temp_key + "_ova_minimap2.paf", 'r') as hin:
            for line in hin:
                paf_rec_count = paf_rec_count + 1

        if paf_rec_count < 3:
            logger.debug(f"Not enough PAF records for the first round consensus generation for {self.temp_key}")
            return

        consensus = ''
        try:
            with open(self.tmp_dir + '/' + self.temp_key + "_ref_polished.fa", 'w') as hout:
                subprocess.check_call(["racon", self.tmp_dir + '/' + self.temp_key + ".supporting_read.fa",
                    self.tmp_dir + '/' + self.temp_key + "_ova_minimap2.paf",
                    self.tmp_dir + '/' + self.temp_key + "_ref.fa"], 
                    stderr = subprocess.DEVNULL, stdout = hout)
        except Exception as e:
            logger.warning(f'{e}')
            return

        with open(self.tmp_dir + '/' + self.temp_key + "_ref_polished.fa", 'r') as hin:
            for line in hin:
                if line.startswith('>'): continue
                consensus = consensus + line.rstrip('\n')

        if len(consensus) < 1000: return
                        
        # second round racon
        with open(self.tmp_dir + '/' + self.temp_key + "_ref_2nd.fa", 'w') as hout:
            print(f">{self.temp_key}_1st_polished_seq\n{consensus}", file = hout)

        with open(self.tmp_dir + '/' + self.temp_key + "_ova_minimap2_2nd.paf", 'w') as hout:
            subprocess.check_call(["minimap2", "-x", "map-ont", self.tmp_dir + '/' + self.temp_key + "_ref_2nd.fa",
                self.tmp_dir + '/' + self.temp_key + ".supporting_read.fa"],
                stderr = subprocess.DEVNULL, stdout = hout)

        paf_rec_count = 0
        with open(self.tmp_dir + '/' + self.temp_key + "_ova_minimap2_2nd.paf", 'r') as hin:
            for line in hin:
                paf_rec_count = paf_rec_count + 1
            
        if paf_rec_count < 3:
            logger.debug(f"Not enough PAF records for the first round consensus generation for {self.temp_key}")
            return

        consensus = ''
        try:
            with open(self.tmp_dir + '/' + self.temp_key + "_ref_polished_2nd.fa", 'w') as hout:
                subprocess.check_call(["racon", self.tmp_dir + '/' + self.temp_key + ".supporting_read.fa",
                    self.tmp_dir + '/' + self.temp_key + "_ova_minimap2_2nd.paf",
                    self.tmp_dir + '/' + self.temp_key + "_ref_2nd.fa"],
                    stderr = subprocess.DEVNULL, stdout = hout)
        except Exception as e:
            logger.warning(f'{e}')
            return

        with open(self.tmp_dir + '/' + self.temp_key + "_ref_polished_2nd.fa", 'r') as hin:
            for line in hin:
                if line.startswith('>'): continue
                consensus = consensus + line.rstrip('\n')

        print(f"{self.temp_key}\t{consensus}", file = self.hout)



def generate_consensus(input_file, output_file, use_racon = False, debug = False): 

    if debug: logger.setLevel(logging.DEBUG)

    consensus_generator = Consensus_generator(output_file, use_racon, debug)
        
    with open(input_file, 'r') as hin:
        for line in hin:
            tkey, treadid, tsize, tseq = line.rstrip('\n').split('\t')
            # when new key appears transaction            
            if tkey != consensus_generator.temp_key:
                if consensus_generator.temp_key is not None:
                    consensus_generator.print_consensus()
                consensus_generator.initialize(tkey)

            consensus_generator.add_support_read_seq(treadid, tseq, tsize)

        consensus_generator.print_consensus()

    del consensus_generator 


def generate_consensus_sbnd(input_file, output_file, use_racon = False, debug = False):

    if debug: logger.setLevel(logging.DEBUG)

    # generate contig for single breakend
    consensus_generator_sbnd = Consensus_generator(output_file, use_racon, debug)
    with open(input_file, 'r') as hin:
        for line in hin:
            tkey, treadid, tsize, tseq = line.rstrip('\n').split('\t')
            # when new key appears transaction            
            if tkey != consensus_generator_sbnd.temp_key:
                if consensus_generator_sbnd.temp_key is not None:
                    consensus_generator_sbnd.print_consensus_sbnd()
                consensus_generator_sbnd.initialize(tkey)

            consensus_generator_sbnd.add_support_read_seq(treadid, tseq, tsize)

        if consensus_generator_sbnd.temp_key is not None:
            consensus_generator_sbnd.print_consensus_sbnd()

    del consensus_generator_sbnd

    
