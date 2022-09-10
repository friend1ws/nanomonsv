#! /usr/bin/env python3

import csv, itertools
import pysam, parasail

from nanomonsv.my_seq import reverse_complement
from nanomonsv.logger import get_logger

logger = get_logger(__name__)


def filter_indel_in_simple_repeat(tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq, simple_repeat_dist_margin = 30):

    if tchr1 == tchr2 and tdir1 == '+' and tdir2 == '-':
        sv_size = tpos2 - tpos1 + len(tinseq) - 1

        tabix_error_flag = False
        try:
            records = simple_repeat_tb.fetch(tchr1, max(tpos1 - simple_repeat_dist_margin + 1, 0), 
                tpos1 + simple_repeat_dist_margin)
        except Exception as inst:
            logger.warning("%s: %s" % (type(inst), inst.args))
            tabix_error_flag = True

        if tabix_error_flag == False:
            for record_line in records:
                record = record_line.split('\t')
                if tpos1 >= int(record[1]) - simple_repeat_dist_margin and \
                    int(tpos2) <= int(record[2]) + simple_repeat_dist_margin:
                    return True

        return False




def post_filter_main(args):

    with open(args.sv_list_file, 'r') as hin:
        dreader = csv.DictReader(hin, delimiter = '\t')
        header = dreader.fieldnames

        for F in dreader:
            tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq = F["Chr_1"], int(F["Pos_1"]), F["Dir_1"], F["Chr_2"], int(F["Pos_2"]), F["Dir_2"], F["Inserted_Seq"]

            simple_repeat_flag = filter_indel_in_simple_repeat(tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq)

            if F["Is_Filter"] == "PASS": 
                F["Is_Filter"] = "Simple_repeat"
            else:
                F["Is_Filter"] = F["Is_Filter"] + ';' + "Simple_repeat"

            print('\t'.join(F.values()))
 

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(prog = "nanomonsv_post_filter",
        description = "Post filtering script for the result of nanomonsv")

    parser.add_argument("sv_list_file", type = str,
                        help = "Path to the nanomonsv result file")

    parser.add_argument("output_file", type = str,
                        help = "Path to the output file")

    parser.add_argument("--simple_repeat_bed", metavar = "simpleRepeat.bed.gz", type = str, default = None,
                        help = "Path to the tabix indexed simple repeat bed file")

    args = parser.parse_args()

    post_filter_main(args)


