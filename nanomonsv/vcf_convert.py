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
            '##FILTER=<ID=Duplicate_with_close_SV,Description="When multiple SVs that share breakpoints in close proximity are detected, all but one SVs are filtered.">\n'\
            '##FILTER=<ID=Duplicate_with_insertion,Description="Breakend SVs that are inferred to be the same as any of detected insertions">\n'\
            '##FILTER=<ID=Duplicate_with_close_insertion,Description="When multiple insertions in close proximity are detected, all but one insertions are filtered.">\n'\
            '##FILTER=<ID=SV_with_decoy,Description="SVs involving decoy contigs">\n'\
            '##FILTER=<ID=Too_small_size,Description="Insertions whose size is below the threshould (currently 100bp)">\n'\
            '##FILTER=<ID=Too_low_VAF,Description="SVs whose variant allele frequencies are inferred to be low">\n'\
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n'\
            '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n'\
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\n'\
            '##INFO=<ID=MATEID,Number=1,Type=String,Description="ID of mate breakend">\n'\
            '##INFO=<ID=SVINSLEN,Number=1,Type=Integer,Description="Length of insertion">\n'\
            '##INFO=<ID=SVINSSEQ,Number=1,Type=String,Description="Sequence of insertion">\n'\
            '##ALT=<ID=DEL,Description="Deletion">\n'\
            '##ALT=<ID=INS,Description="Insertion">\n'\
            '##ALT=<ID=DUP,Description="Duplication">\n'\
            '##FORMAT=<ID=TR,Number=1,Type=Integer,Description="The number of reads around the breakpoints">\n'\
            '##FORMAT=<ID=VR,Number=1,Type=Integer,Description="The number of variant supporting reads determined in the validation realignment step">'

    with open(result_file, 'r') as hin, open(output_vcf, 'w') as hout:

        dreader = csv.DictReader(hin, delimiter = '\t')
        fieldname_list = dreader.fieldnames
        is_control = True if "Checked_Read_Num_Control" in fieldname_list and "Supporting_Read_Num_Control" in fieldname_list else False

        if is_control:
            header = header + '\n' + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tCONTROL"
        else:
            header = header + '\n' + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR"

        print(header,  file = hout)

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
                if tref == '' or tref is None: continue
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
                elif tsvlen > 0:
                    talt = "<INS>"
                    tinfo = f"END={tend};SVTYPE=INS;SVINSLEN={tsvinslen};SVINSSEQ={tsvinsseq}"

                else:
                    continue
                print(f"{tchrom}\t{tpos}\t{tid}\t{tref}\t{talt}\t{tqual}\t{tfilter}\t{tinfo}\t{tformat_sample}", file = hout)

            # Duplication
            elif F["Chr_1"] == F["Chr_2"] and F["Dir_1"] == '-' and F["Dir_2"] == '+' and F["Pos_1"] != '1': 

                tpos = int(F["Pos_1"]) - 1
                tref = ref_tb.fetch(tchrom, tpos - 1, tpos)
                if tref == '' or tref is None: continue
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
                tpos1 = int(F["Pos_1"])
                tref1 = ref_tb.fetch(tchrom1, tpos1 - 1, tpos1)
                if tref1 == '' or tref1 is None: continue

                tchrom2 = F["Chr_2"]
                tpos2 = int(F["Pos_2"])
                tref2 = ref_tb.fetch(tchrom2, tpos2 - 1, tpos2)
                if tref2 == '' or tref2 is None: continue

                tbracket = ']' if F["Dir_2"] == '+' else '['
                if F["Dir_1"] == '+':
                    talt1 = f'{tref1}{tsvinsseq}{tbracket}{tchrom2}:{tpos2}{tbracket}'
                else:
                    talt1 = f'{tbracket}{tchrom2}:{tpos2}{tbracket}{tsvinsseq}{tref2}' 

                tinfo1 = f"SVTYPE=BND;MATEID={tid}_1"
                if tsvinslen != 0: tinfo1 = tinfo1 + f";SVINSLEN={tsvinslen};SVINSSEQ={tsvinsseq}"

                print(f"{tchrom1}\t{tpos1}\t{tid}_0\t{tref1}\t{talt1}\t{tqual}\t{tfilter}\t{tinfo1}\t{tformat_sample}", file = hout)

                # tchrom2 = F["Chr_2"]
                # tpos = int(F["Pos_2"])
                # tref = ref_tb.fetch(tchrom2, tpos - 1, tpos)
                # if tref == '' or tref is None: continue
                tbracket = ']' if F["Dir_1"] == '+' else '['
                tsvinsseq = reverse_complement(tsvinsseq)
                if F["Dir_2"] == '+':
                    talt2 = f'{tref2}{tsvinsseq}{tbracket}{tchrom1}:{tpos1}{tbracket}'
                else:
                    talt2 = f'{tbracket}{tchrom1}:{tpos1}{tbracket}{tsvinsseq}{tref2}'

                tinfo2 = f"SVTYPE=BND;MATEID={tid}_0"
                if tsvinslen != 0: tinfo2 = tinfo2 + f";SVINSLEN={tsvinslen};SVINSSEQ={tsvinsseq}"

                print(f"{tchrom2}\t{tpos2}\t{tid}_1\t{tref2}\t{talt2}\t{tqual}\t{tfilter}\t{tinfo2}\t{tformat_sample}", file = hout)


if __name__ == "__main__":

    import sys
    result_file = sys.argv[1]
    output_vcf = sys.argv[2]
    reference = sys.argv[3]
    genomesv2vcf_convert(result_file, output_vcf, reference)
