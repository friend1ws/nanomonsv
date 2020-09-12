#! /usr/bin/env python

import csv, sys

sv_get_file = sys.argv[1]
output_file = sys.argv[2]

with open(sv_get_file, 'r') as hin, open(output_file, 'w') as hout:

    dreader = csv.DictReader(hin, delimiter = '\t')

    # print header
    print('\t'.join(dreader.fieldnames) + '\t' + "SV_Type", file = hout)

    for F in dreader:

        sv_type = None
        inseq = '' if F["Inserted_Seq"] == '---' else F["Inserted_Seq"]
        if F["Chr_1"] != F["Chr_2"]:
            sv_type = "Translocation"
        elif F["Dir_1"] == F["Dir_2"]:
            sv_type = "Inversion"
        elif F["Dir_1"] == '-' and F["Dir_2"] == '+':
            sv_type = "Duplication"
        elif F["Dir_1"] == '+' and F["Dir_2"] == '-':
            if len(inseq) > int(F["Pos_2"]) - int(F["Pos_1"]) + 1:
                sv_type = "Insertion"
            else:
                sv_type = "Deletion"

        print('\t'.join(F.values()) + '\t' + sv_type, file = hout)



