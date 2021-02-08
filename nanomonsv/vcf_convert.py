#! /usr/bin/en python3

import sys, csv, datetime
import pysam
from .version import __version__
from .my_seq import reverse_complement


def genomesv2vcf_convert(result_file, output_vcf, reference):

    today_str = datetime.datetime.today().strftime("%Y%m%d")

    header = '##fileformat=VCFv4.3\n'\
             f'##fileDate={today_str}\n'\
             f'##source=nanomonsv-{__version__}\n'\
             f'##reference={reference}'

    ref_tb = pysam.FastaFile(reference)

    for (tchr, tlen) in zip(ref_tb.references, ref_tb.lengths):
        header = header + '\n' + f"##contig=<ID={tchr},length={tlen}>"

    header = header + '\n' + \
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n'\
            '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">\n'\
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\n'\
            '##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakend">\n'\
            '##INFO=<ID=SVINSLEN,Number=.,Type=Integer,Description="Length of insertion">\n'\
            '##INFO=<ID=SVINSSEQ,Number=.,Type=String,Description="Sequence of insertion">\n'\
            '##ALT=<ID=DEL,Description="Deletion">\n'\
            '##ALT=<ID=INS,Description="Insertion">\n'\
            '##ALT=<ID=DUP,Description="Duplication">\n'\
            '##FORMAT=<ID=TR,Number=.,Type=Integer,Description="The number of reads around the breakpoints">\n'\
            '##FORMAT=<ID=VR,Number=.,Type=Integer,Description="The number of variant supporting reads determined in the validation realignment step">\n'\
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tCONTROL'

    with open(result_file, 'r') as hin, open(output_vcf, 'w') as hout:

        print(header,  file = hout)
        dreader = csv.DictReader(hin, delimiter = '\t')
        fieldname_list = dreader.fieldnames
        is_control = True if "Checked_Read_Num_Control" in fieldname_list and "Supporting_Read_Num_Control" in fieldname_list else False

        for F in dreader:

            tchrom = F["Chr_1"]
            tid = F["SV_ID"]
            tqual = '.'
            tfilter = F["Is_Filter"]

            if F["Inserted_Seq"] != "---":
                tsvinsseq = F["Inserted_Seq"]
                tsvinslen = len(F["Inserted_Seq"])
            else:
                tsvinsseq = ''
                tsvinslen = 0
            
            tformat_sample = f'TR:VR\t{F["Checked_Read_Num_Tumor"]}:{F["Supporting_Read_Num_Tumor"]}'
            if is_control:
                tformat_sample = tformat_sample + f'\t{F["Checked_Read_Num_Control"]}:{F["Supporting_Read_Num_Control"]}' 

            if F["Chr_1"] == F["Chr_2"] and F["Dir_1"] == '+' and F["Dir_2"] == '-':

                tpos = int(F["Pos_1"])
                tref = ref_tb.fetch(tchrom, tpos - 1, tpos)
                tsvlen = int(F["Pos_2"]) - int(F["Pos_1"]) - 1
                tend = int(F["Pos_2"]) - 1

                # Deletion
                if tsvlen > tsvinslen:
                    talt = "<DEL>"
                    tsvlen = int(F["Pos_2"]) - int(F["Pos_1"]) - 1
                    tinfo = f"END={tend};SVTYPE=DEL;SVLEN=-{tsvlen}"

                    if tsvinslen != 0:
                        tinfo = tinfo + f";SVINSLEN={tsvinslen};SVINSSEQ={tsvinsseq}"

                # Insertion
                else:
                    talt = "<INS>"
                    tinfo = f"END={tend};SVTYPE=INS;SVINSLEN={tsvinslen};SVINSSEQ={tsvinsseq}"

                print(f"{tchrom}\t{tpos}\t{tid}\t{tref}\t{talt}\t{tqual}\t{tfilter}\t{tinfo}\t{tformat_sample}", file = hout)

            # Duplication
            elif F["Chr_1"] == F["Chr_2"] and F["Dir_1"] == '-' and F["Dir_2"] == '+': 

                tpos = int(F["Pos_1"]) - 1
                tref = ref_tb.fetch(tchrom, tpos - 1, tpos)
                talt = "<DUP>"
                tend = int(F["Pos_2"]) 
                tsvlen = int(F["Pos_2"]) - int(F["Pos_1"]) + 1
                tinfo = f"END={tend};SVTYPE=DUP;SVLEN={tsvlen}"
                if tsvinslen != 0:
                    tinfo = tinfo + f";SVINSLEN={tsvinslen};SVINSSEQ={tsvinsseq}"

                print(f"{tchrom}\t{tpos}\t{tid}\t{tref}\t{talt}\t{tqual}\t{tfilter}\t{tinfo}\t{tformat_sample}", file = hout)

            # Breakend 
            else:
                tchrom1 = F["Chr_1"]
                tpos = int(F["Pos_1"])
                tref = ref_tb.fetch(tchrom1, tpos - 1, tpos)
                tbracket = ']' if F["Dir_2"] == '+' else '['
                if F["Dir_1"] == '+':
                    talt = f'{tref}{tsvinsseq}{tbracket}{F["Chr_2"]}:{F["Pos_2"]}{tbracket}'
                else:
                    talt = f'{tbracket}{F["Chr_2"]}:{F["Pos_2"]}{tbracket}{tsvinsseq}{tref}' 

                tinfo = f"SVTYPE=BND;MATEID={tid}_0"
                if tsvinslen != 0: tinfo = tinfo + f";SVINSLEN={tsvinslen};SVINSSEQ={tsvinsseq}"

                print(f"{tchrom1}\t{tpos}\t{tid}_0\t{tref}\t{talt}\t{tqual}\t{tfilter}\t{tinfo}\t{tformat_sample}", file = hout)

                tchrom2 = F["Chr_2"]
                tpos = int(F["Pos_2"])
                tref = ref_tb.fetch(tchrom2, tpos - 1, tpos)
                tbracket = ']' if F["Dir_1"] == '+' else '['
                tsvinsseq = reverse_complement(tsvinsseq)
                if F["Dir_2"] == '+':
                    talt = f'{tref}{tsvinsseq}{tbracket}{F["Chr_1"]}:{F["Pos_1"]}{tbracket}'
                else:
                    talt = f'{tbracket}{F["Chr_1"]}:{F["Pos_1"]}{tbracket}{tsvinsseq}{tref}'

                tinfo = f"SVTYPE=BND;MATEID={tid}_1"
                if tsvinslen != 0: tinfo = tinfo + f";SVINSLEN={tsvinslen};SVINSSEQ={tsvinsseq}"

                print(f"{tchrom2}\t{tpos}\t{tid}_1\t{tref}\t{talt}\t{tqual}\t{tfilter}\t{tinfo}\t{tformat_sample}", file = hout)


if __name__ == "__main__":

    import sys
    result_file = sys.argv[1]
    output_vcf = sys.argv[2]
    reference = sys.argv[3]
    genomesv2vcf_convert(result_file, output_vcf, reference)
