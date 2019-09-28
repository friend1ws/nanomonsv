#! /usr/bin/env python

import itertools
def bp_check(input_file, output_file, check_margin = 50, minimum_ambiguity = 20):
    
    hout = open(output_file, 'w')

    temp_read_name = ''
    query2target = {}
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[11] == "True": continue

            if F[0] != temp_read_name:
                query_list = list(query2target)
                if temp_read_name != '' and len(query_list) > 1:

                    bp_flag = False
                    for qpos_comb in list(itertools.combinations(query_list, 2)):

                        
                        if qpos_comb[0][1] < qpos_comb[1][1]:
                            qpos1, qpos2 = qpos_comb[0], qpos_comb[1]
                        else:
                            qpos1, qpos2 = qpos_comb[1], qpos_comb[0]

                        # if the first region completely covers the second region
                        if qpos1[2] >= qpos2[2]: continue

                        # if there is a significant overlap         
                        if (qpos1[2] - qpos2[1]) / (qpos2[2] - qpos1[2]) >= 0.2: continue

                        if abs(qpos2[1] - qpos1[2]) <= check_margin:
                            bp_flag = True
                            tchr1, tstart1, tend1, tmapQ1, tnumM1, tnumI1, tnumD1, tis_supp1 = query2target[qpos1]
                            tchr2, tstart2, tend2, tmapQ2, tnumM2, tnumI2, tnumD2, tis_supp2 = query2target[qpos2]      

                            if qpos2[1] - qpos1[2] > 0:
                                outward_ambiguity, inward_ambiguity = minimum_ambiguity, max(qpos2[1] - qpos1[2], minimum_ambiguity)
                            else:
                                outward_ambiguity, inward_ambiguity = max(qpos1[2] - qpos2[1], minimum_ambiguity), minimum_ambiguity

                            bchr1, bchr2 = tchr1, tchr2
                            if qpos1[4] == '+': 
                                bstart1, bend1, bstrand1 = int(tend1) - inward_ambiguity, int(tend1) + outward_ambiguity, '+'
                            else:
                                bstart1, bend1, bstrand1 = int(tstart1) - outward_ambiguity, int(tstart1) + inward_ambiguity, '-'
                    
                            if qpos2[4] == '+':
                                bstart2, bend2, bstrand2 = int(tstart2) - outward_ambiguity, int(tstart2) + inward_ambiguity, '-'
                            else:
                                bstart2, bend2, bstrand2 = int(tend2) - inward_ambiguity, int(tend2) + outward_ambiguity, '+'

                            bread_name = qpos1[0]

                            binfo1 = ','.join([str(qpos1[1]), str(qpos1[2]), tmapQ1, tnumM1, tnumI1, tnumD1, tis_supp1])
                            binfo2 = ','.join([str(qpos2[1]), str(qpos2[2]), tmapQ2, tnumM2, tnumI2, tnumD2, tis_supp2])

                            if bchr1 > bchr2 or (bchr1 == bchr2 and bstart1 > bstart2):
                                bchr1, bstart1, bend1, bstrand1, binfo1, bchr2, bstart2, bend2, bstrand2, binfo2 = \
                                    bchr2, bstart2, bend2, bstrand2, binfo2, bchr1, bstart1, bend1, bstrand1, binfo1

                            print('\t'.join([bchr1, str(bstart1), str(bend1), bchr2, str(bstart2), str(bend2), bread_name, "0",
                                             bstrand1, bstrand2, binfo1, binfo2]))

                            
                        # 00000538-979f-4e25-9be5-247d605273e2    138     9341    9357    -       5       97828126        97837738        60      8703    501     910     False   False

                        # 0000058b-a141-43fd-8001-0253a1d18ec4    39      3009    5124    +       12      47383760        47386790        60      2884    87      147     False   False
                    """
                    if bp_flag:
                        for query in query2target:
                            query_list = '\t'.join([str(x) for x in query])
                            target_list = '\t'.join([str(x) for x in query2target[query]])

                            print(query_list + '\t' + target_list)
                    """

                temp_read_name = F[0]
                query2target = {}

            query2target[(F[0], int(F[1]), int(F[2]), int(F[3]), F[4])] = (F[5], int(F[6]), int(F[7]), F[8], F[9], F[10], F[11], F[12])




if __name__ == "__main__":

    import sys
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    bp_check(input_file, output_file)


