#! /usr/bin/env python

import gzip, statistics
import pysam

def cluster_rearrangement(input_file, output_file, cluster_margin_size = 100):

    hout = open(output_file, 'w')

    merged_bedpe = {}
    with gzip.open(input_file, 'rt') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
    
            for key in sorted(merged_bedpe):

                bchr1, bstart1, bend1, bdir1, bchr2, bstart2, bend2, bdir2 = key
                breadids, binfo1, binfo2 = merged_bedpe[key]

                if F[0] != bchr1 or int(F[1]) > int(bend1) + cluster_margin_size: 

                    print('\t'.join([bchr1, str(bstart1), str(bend1), bchr2, str(bstart2), str(bend2), \
                       breadids, '0', bdir1, bdir2, binfo1, binfo2]), file = hout)
                    del merged_bedpe[key]


            match_flag = False
            for key in sorted(merged_bedpe):

                bchr1, bstart1, bend1, bdir1, bchr2, bstart2, bend2, bdir2 = key 
                breadids, binfo1, binfo2 = merged_bedpe[key]

                if F[0] == bchr1 and F[3] == bchr2 and F[8] == bdir1 and F[9] == bdir2 and \
                  int(F[2]) >= bstart1 and int(F[1]) <= bend1 and int(F[5]) >= bstart2 and int(F[4]) <= bend2:
                        
                    new_start1 = min(int(F[1]), bstart1)
                    new_end1 = max(int(F[2]), bend1)
                    new_start2 = min(int(F[4]), bstart2)
                    new_end2 = max(int(F[5]), bend2)
                            
                    new_key = (bchr1, new_start1, new_end1, bdir1, bchr2, new_start2, new_end2, bdir2)
                    new_readids = breadids + ';' + F[6]
                    new_info1 = binfo1 + ';' + F[10]
                    new_info2 = binfo2 + ';' + F[11]

                    del merged_bedpe[key]

                    merged_bedpe[new_key] = (new_readids, new_info1, new_info2)
                    match_flag = True
                    break

                    
            if not match_flag:
                new_key = (F[0], int(F[1]), int(F[2]), F[8], F[3], int(F[4]), int(F[5]), F[9])
                merged_bedpe[new_key] = (F[6], F[10], F[11])


    for key in sorted(merged_bedpe):
        bchr1, bstart1, bend1, bdir1, bchr2, bstart2, bend2, bdir2 = key
        breadids, binfo1, binfo2 = merged_bedpe[key]

        print('\t'.join([bchr1, str(bstart1), str(bend1), bchr2, str(bstart2), str(bend2), \
          breadids, '0', bdir1, bdir2, binfo1, binfo2]), file = hout)

    hout.close()


def filt_clustered_rearrangement1(input_file, output_file, read_num_thres = 3, median_mapQ_thres = 40, max_overhang_size_thres = 300):

    hout = open(output_file, 'w') 
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            read_ids = F[6].split(';')
            if len(read_ids) < read_num_thres: continue

            info1 = F[10].split(';')
            info2 = F[11].split(';')

            median_mapQ1 = statistics.median([int(x.split(',')[5]) for x in info1])
            median_mapQ2 = statistics.median([int(x.split(',')[5]) for x in info2])
        
            max_overhang_size1 = max([abs(int(x.split(',')[2]) - int(x.split(',')[0])) for x in info1])
            max_overhang_size2 = max([abs(int(x.split(',')[2]) - int(x.split(',')[0])) for x in info2]) 
        
            if median_mapQ1 < median_mapQ_thres or median_mapQ2 < median_mapQ_thres: continue
            if max_overhang_size1 < max_overhang_size_thres or max_overhang_size2 < max_overhang_size_thres: continue

            print('\t'.join(F), file = hout)


    hout.close()



def filt_clustered_rearrangement2(input_file, output_file, control_junction_bedpe, control_read_num_thres = 0, control_check_margin = 50):

    hout = open(output_file, 'w')
    if control_junction_bedpe is not None: control_junction_db = pysam.TabixFile(control_junction_bedpe)
    
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            tchr1, tstart1, tend1, tchr2, tstart2, tend2 = F[0], F[1], F[2], F[3], F[4], F[5]
            tdir1, tdir2 = F[8], F[9]

            is_short_deletion = True if tchr1 == tchr2 and tdir1 == '+' and tdir2 == '-' and int(tend2) - int(tstart1) <= 500 else False

            control_flag = False
            if control_junction_bedpe is not None:
                tabix_error_flag = False
                try:
                    records = control_junction_db.fetch(F[0], int(tstart1) - 200, int(tend1) + 200)
                except:
                    tabix_error_flag = True

                if not tabix_error_flag:
                    for record_line in records:
                        record = record_line.split('\t')

                        if is_short_deletion:
                            if tchr1 == record[0] and tchr2 == record[3] and tdir1 == record[8] and tdir2 == record[9] and \
                                int(tstart1) <= int(record[5]) and int(tend2) >= int(record[1]):   
                                control_flag = True
                        else:
                            if tchr1 == record[0] and tdir1 == record[8] and int(tend1) >= int(record[1]) - control_check_margin and int(tstart1) <= int(record[2]) + control_check_margin and \
                                tchr2 == record[3] and tdir2 == record[9] and int(tend2) >= int(record[4]) - control_check_margin and int(tstart2) <= int(record[5]) + control_check_margin:
                                control_flag = True

            if not control_flag:
                print('\t'.join(F), file = hout)


    hout.close()
    if control_junction_bedpe is not None: control_junction_db.close()



def cluster_deletion(input_file, output_file, deletion_cluster_margin_size = 10, check_margin_size = 50):

    dcms = deletion_cluster_margin_size
    hout = open(output_file, 'w')
    
    merged_bedpe = {}
    with gzip.open(input_file, 'rt') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if int(F[4]) < 90: continue
 
            for key in sorted(merged_bedpe):

                # bchr, bstart, bend = key
                bchr1, bstart1, bend1, bdir1, bchr2, bstart2, bend2, bdir2 = key
                breadids, bsize, binfo = merged_bedpe[key]

                if F[0] != bchr1 or int(F[1]) > int(bend1) + check_margin_size: 

                    print('\t'.join([bchr1, str(bstart1), str(bend1), bchr2, str(bstart2), str(bend2), \
                       breadids, '0', bdir1, bdir2, bsize, binfo]), file = hout)
                    del merged_bedpe[key]

                    # print('\t'.join([bchr, str(bstart), str(bend), breadids, bsize, '+', binfo]), file = hout)
                    # del merged_bedpe[key]


            match_flag = False
            for key in sorted(merged_bedpe):

                # bchr, bstart, bend = key 
                bchr1, bstart1, bend1, bdir1, bchr2, bstart2, bend2, bdir2 = key 
                breadids, bsize, binfo = merged_bedpe[key]

                # if F[0] == bchr and abs(int(F[1]) - int(bstart)) <= deletion_cluster_margin_size and \
                #     abs(int(F[2]) - int(bend)) <= deletion_cluster_margin_size:
                # if F[0] == bchr1 and F[3] == bchr2 and F[8] == bdir1 and F[9] == bdir2 and \
                #   int(F[2]) >= bstart1 and int(F[1]) <= bend1 and int(F[5]) >= bstart2 and int(F[4]) <= bend2:
                if F[0] == bchr1 and int(F[1]) - dcms <= bend1 and int(F[1]) + dcms >= bstart1 and \
                  int(F[2]) - dcms <= bend2 and int(F[2]) + dcms >= bstart2:

                    new_start1 = min(int(F[1]) - dcms, bstart1)
                    new_end1 = max(int(F[1]) + dcms, bend1)
                    new_start2 = min(int(F[2]) - dcms, bstart2)
                    new_end2 = max(int(F[2]) + dcms, bend2)
                            
                    new_key = (bchr1, new_start1, new_end1, bdir1, bchr2, new_start2, new_end2, bdir2)

                    # ney_key = key
                    new_size = bsize + ';' + F[4]
                    new_readids = breadids + ';' + F[3]
                    new_info = binfo + ';' + F[6]

                    del merged_bedpe[key]

                    merged_bedpe[new_key] = (new_readids, new_size, new_info)
                    match_flag = True
                    break

                    
            if not match_flag:
                new_key = (F[0], int(F[1]) - dcms, int(F[1]) + dcms, '+', F[0], int(F[2]) - dcms, int(F[2]) + dcms, '-') 
                # new_key = (F[0], F[1], F[2])
                merged_bedpe[new_key] = (F[3], F[4], F[6])


        for key in sorted(merged_bedpe):
            # bchr, bstart, bene = key
            bchr1, bstart1, bend1, bdir1, bchr2, bstart2, bend2, bdir2 = key
            breadids, bsize, binfo = merged_bedpe[key]

            print('\t'.join([bchr1, str(bstart1), str(bend1), bchr2, str(bstart2), str(bend2), \
              breadids, '0', bdir1, bdir2, bsize, binfo]), file = hout)

            # if F[0] != bchr or int(F[1]) > int(bend) + check_margin_size:  
            #     print('\t'.join([bchr, str(bstart), str(bend), breadids, bsize, '+', binfo]), file = hout)

    hout.close()



def filt_clustered_deletion1(input_file, output_file, read_num_thres = 3, median_mapQ_thres = 40, max_overhang_size_thres = 300):

    hout = open(output_file, 'w')
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            read_ids = list(set(F[6].split(';')))
            if len(read_ids) < read_num_thres: continue

            info = F[11].split(';')

            median_mapQ = statistics.median([int(x.split(',')[5]) for x in info])
            non_secondary_readnum = len([x.split(',')[10] for x in info if x.split(',')[10] == "False"])
            # max_overhang_size = max([abs(int(x.split(',')[1]) - int(x.split(',')[0])) for x in info2])

            if median_mapQ < median_mapQ_thres: continue
            if non_secondary_readnum < read_num_thres: continue
            # if max_overhang_size1 < max_overhang_size_thres or max_overhang_size2 < max_overhang_size_thres: continue

            print('\t'.join(F), file = hout)


    hout.close()


def filt_clustered_deletion2(input_file, output_file, control_junction_bedpe, control_read_num_thres = 0, control_check_margin = 50):

    hout = open(output_file, 'w')
    if control_junction_bedpe is not None: control_junction_db = pysam.TabixFile(control_junction_bedpe)
    
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            tchr1, tstart1, tend1, tchr2, tstart2, tend2 = F[0], F[1], F[2], F[3], F[4], F[5]
            median_size = statistics.median([int(x) for x in F[10].split(';')])

            control_flag = False
            if control_junction_bedpe is not None:
                tabix_error_flag = False
                try:
                    records = control_junction_db.fetch(F[0], int(tstart1) - 50, int(tend2) + 50)
                except:
                    tabix_error_flag = True

                if not tabix_error_flag:
                    for record_line in records:
                        record = record_line.split('\t')

                        if tchr1 == record[0] and int(tstart1) - control_check_margin <= int(record[2]) and int(tend2) + control_check_margin >= int(record[1]) \
                            and int(record[2]) - int(record[1]) >= median_size * 0.5:
                            control_flag = True


            if not control_flag:
                print('\t'.join(F), file = hout)


    hout.close()
    if control_junction_bedpe is not None: control_junction_db.close()

"""
if __name__ == "__main__":

    import sys

    tumor_prefix = sys.argv[1]
    control_prefix = sys.argv[2]

    cluster_junction(tumor_prefix + ".junction.sorted.bedpe.gz", tumor_prefix + ".junction.sorted.clustered.bedpe")

    filt_clustered_junction1(tumor_prefix + ".junction.sorted.clustered.bedpe", tumor_prefix + ".junction.sorted.clustered.filt1.bedpe")

    filt_clustered_junction2(tumor_prefix + ".junction.sorted.clustered.filt1.bedpe", tumor_prefix + ".junction.sorted.clustered.filt2.bedpe", 
                             control_prefix + ".junction.sorted.bedpe.gz")
"""

