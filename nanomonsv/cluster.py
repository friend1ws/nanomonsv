#! /usr/bin/env python

import sys, gzip, statistics
import pysam

class Sv_cluster(object):

    def __init__(self, tchr1, tstart1, tend1, tdir1, tchr2, tstart2, tend2, tdir2,
        treadid, tsize, tinfo1, tinfo2):

        self.chr1 = tchr1
        self.start1 = tstart1
        self.end1 = tend1
        self.dir1 = tdir1
        self.chr2 = tchr2
        self.start2 = tstart2
        self.end2 = tend2
        self.dir2 = tdir2
        self.readids = [treadid]
        self.size = [tsize]
        self.info1 = [tinfo1]
        self.info2 = [tinfo2]


class Sv_clusterer(object):

    def __init__(self, svtype, output_file, control_junction_bedpe = None, bp_bed = None,
        cluster_margin_size = None, size_margin_ratio = None, maximum_local_variant_num = None, skip_margin = None,
        read_num_thres = None, median_mapQ_thres = None, 
        max_overhang_size_thres = None, control_read_num_thres = None, control_check_margin = None):

        self.hout = open(output_file, 'w')
        self.control_tb = None
        if control_junction_bedpe is not None:
            self.control_tb = pysam.TabixFile(control_junction_bedpe)

        self.sv_cluster_list = []
        self.svtype = svtype
        self.control_junction_bedpe = control_junction_bedpe
        self.bp_bed = bp_bed
        self.cluster_margin_size = cluster_margin_size
        self.size_margin_ratio = size_margin_ratio
        self.maximum_local_variant_num = maximum_local_variant_num
        self.skip_margin = skip_margin
        self.read_num_thres = read_num_thres
        self.median_mapQ_thres = median_mapQ_thres
        self.max_overhang_size_thres = max_overhang_size_thres
        self.control_read_num_thres = control_read_num_thres
        self.control_check_margin = control_check_margin

        self.next_pos_after_skip = 0
        self.temp_chr = None


    def check_mergeability(self, tchr1, tstart1, tend1, tdir1, tchr2, tstart2, tend2, tdir2, 
        treadid, tsize, tinfo1, tinfo2):
    
        for cluster in self.sv_cluster_list:

            if treadid in cluster.readids: continue

            # pairs of chromosome and direction should be the same.
            # two pairs of start and end should overlap.   
            is_mergeable = False
            if self.svtype == "rearrangement": 
                if tchr1 == cluster.chr1 and tchr2 == cluster.chr2 and tdir1 == cluster.dir1 and tdir2 == cluster.dir2 and \
                    tend1 >= cluster.start1 and tstart1 <= cluster.end1 and tend2 >= cluster.start2 and tstart2 <= cluster.end2: 
                    is_mergeable = True

            elif self.svtype in ["insertion", "deletion"]:

                if tchr1 == cluster.chr1 and tstart1 <= cluster.end1 and tend1 >= cluster.start1 and \
                    tstart2 <= cluster.end2 and tend2 >= cluster.start1 and \
                    tsize > (1.0 - self.size_margin_ratio) * float(min(cluster.size)) and \
                    tsize < (1.0 + self.size_margin_ratio) * float(min(cluster.size)):
                    is_mergeable = True
                   
            if not is_mergeable: continue

            cluster.start1 = min(tstart1, cluster.start1)
            cluster.end1 = max(tend1, cluster.end1)
            cluster.start2 = min(tstart2, cluster.start2)
            cluster.end2 = max(tend2, cluster.end2)

            cluster.readids.append(treadid)
            cluster.size.append(tsize)
            cluster.info1.append(tinfo1)
            cluster.info2.append(tinfo2)
                        
            return cluster
            
        return None


    def add_new_cluster(self, tchr1, tstart1, tend1, tdir1, tchr2, tstart2, tend2, tdir2, 
        treadid, tsize, tinfo1, tinfo2):

        self.sv_cluster_list.append(
            Sv_cluster(tchr1, tstart1, tend1, tdir1, tchr2, tstart2, tend2, tdir2,
                treadid, tsize, tinfo1, tinfo2))


    def filter_rearrangement_cluster(self, cl):
        
        is_filter = False
        if len(cl.readids) < self.read_num_thres: is_filter = True

        median_mapQ1 = statistics.median([int(x.split(',')[5]) for x in cl.info1])
        median_mapQ2 = statistics.median([int(x.split(',')[5]) for x in cl.info2])
        if median_mapQ1 < self.median_mapQ_thres or median_mapQ2 < self.median_mapQ_thres: is_filter = True

        max_overhang_size1 = max([abs(int(x.split(',')[2]) - int(x.split(',')[0])) for x in cl.info1])
        max_overhang_size2 = max([abs(int(x.split(',')[2]) - int(x.split(',')[0])) for x in cl.info2]) 
        if max_overhang_size1 < self.max_overhang_size_thres or max_overhang_size2 < self.max_overhang_size_thres: is_filter = True

        if is_filter == True: return(True)

        if self.control_tb is not None:

            is_short_deletion = False
            if cl.chr1 == cl.chr2 and cl.dir1 == '+' and cl.dir2 == '-' and cl.end2 - cl.start1 <= 500:
                is_short_deletion = True

            control_flag = False
            tabix_error_flag = False
            try:
                records = self.control_tb.fetch(cl.chr1, max(0, cl.start1 - 200), cl.end1 + 200)
            except:
                # logger.warning(f"Tabix fetch error at {cl.chr1}:{max(0, cl.start1 - 200)}-{cl.end1 + 200}")
                tabix_error_flag = True

            if not tabix_error_flag:
                for record_line in records:
                    rec = record_line.split('\t')

                    if cl.chr1 != rec[0] or cl.chr2 != rec[3] or cl.dir1 != rec[8] or cl.dir2 != rec[9]: continue
        
                    if is_short_deletion:
                        if cl.start1 <= int(rec[5]) and cl.end2 >= int(rec[1]):
                            control_flag = True
                    else:
                        if cl.end1 >= int(rec[1]) - self.control_check_margin and \
                            cl.start1 <= int(rec[2]) + self.control_check_margin and \
                            cl.end2 >= int(rec[4]) - self.control_check_margin and \
                            cl.start2 <= int(rec[5]) + self.control_check_margin:
                            control_flag = True

            if control_flag == True: is_filter = True 

        return(is_filter)

        
    def flush_sv_cluster_list(self, current_chr, current_pos):

        remove_cluster = []
        for i in range(len(self.sv_cluster_list)):
            cl = self.sv_cluster_list[i]
            if current_chr != cl.chr1 or current_pos > cl.end1 + self.cluster_margin_size:
                remove_cluster.append(cl)

                if self.svtype == "rearrangement":

                    if self.filter_rearrangement_cluster(cl): continue
                    print_line_readids = ';'.join(cl.readids)
                    print_line_info1 = ';'.join(cl.info1)
                    print_line_info2 = ';'.join(cl.info2)
        
                    print(f'{cl.chr1}\t{cl.start1}\t{cl.end1}\t{cl.chr2}\t{cl.start2}\t{cl.end2}\t' +
                        f'{print_line_readids}\t0\t{cl.dir1}\t{cl.dir2}\t{print_line_info1}\t{print_line_info2}',
                        file = self.hout)

        for cl in remove_cluster:
            self.sv_cluster_list.remove(cl)


    def check_exceeding_local_variant_num(self, current_chr, current_pos):

        if len(self.sv_cluster_list) > self.maximum_local_variant_num:
            print(f'Exceeded maximum number of local variants at {current_chr}:{current_pos}', file = sys.stderr)
            print(f'Skip {current_chr}:{current_pos + self.skip_margin}', file = sys.stderr)
            self.sv_cluster_list = []
            self.next_pos_after_skip = current_pos + self.skip_margin


def cluster_supporting_reads(input_file, output_file, svtype, control_junction_bedpe = None, bp_bed = None, 
    cluster_margin_size = 100, indel_cluster_margin_size = 10, size_margin_ratio = 0.2, min_indel_size = 90,
    maximum_local_variant_num = 1000, skip_margin = 5000, 
    read_num_thres = 3, median_mapQ_thres = 40, max_overhang_size_thres = 300,
    control_read_num_thres = 0, control_check_margin = 50):

    sv_clusterer = Sv_clusterer(svtype, output_file, control_junction_bedpe = control_junction_bedpe, bp_bed = bp_bed,
        cluster_margin_size = cluster_margin_size, size_margin_ratio = size_margin_ratio,
        maximum_local_variant_num = maximum_local_variant_num, skip_margin = skip_margin,
        read_num_thres = read_num_thres, median_mapQ_thres = median_mapQ_thres, 
        max_overhang_size_thres = max_overhang_size_thres, control_read_num_thres = control_read_num_thres,
        control_check_margin = control_check_margin)
    
    with gzip.open(input_file, 'rt') as hin: # , open(output_file, 'w') as hout:
        for line in hin:
            if sv_clusterer.svtype in ["insertion", "deletion"]: # insertion or deletion

                tchr1, tstart, tend, treadid, tsize, tdir, tinfo1 = ne.rstrip('\n').split('\t')
                tstart1 = int(tstart) - indel_cluster_margin_size
                tend1 = int(start) + indel_cluster_margin_size
                tstart2 = int(tend) - indel_cluster_margin_size
                tend2 = int(tend) + indel_cluster_margin_size
                tchr2, tdir1, tdir2 = tchr1, '+', '-'
                tsize, tinfo2 = int(tsize), None

                # this should be parametrized    
                if tsize < min_indel_size: continue
                if tend1 < next_pos_after_skip: continue
                
            elif sv_clusterer.svtype == "rearrangement": # rearrangement
                tchr1, tstart1, tend1, tchr2, tstart2, tend2, treadid, _, tdir1, tdir2, tinfo1, tinfo2 = \
                    line.rstrip('\n').split('\t')
                tstart1, tend1, tstart2, tend2 = int(tstart1), int(tend1), int(tstart2), int(tend2)
                tsize = None

            else:
                pass
                # logger.error("The svtype argument should be either of rearrangement, insertion or deletion.")
                # sys.exit(1)
 
            if sv_clusterer.temp_chr != tchr1:   
                sv_clusterer.next_pos_after_skip = 0
                sv_clusterer.temp_chr = tchr1

            # flush out existing cluster SVs wholse supporting reads have parsed alreadly 
            #   (considering the chromosome and coordinates).
            sv_clusterer.flush_sv_cluster_list(tchr1, tend1)

            # check the new SV can be merged into any existing SV clusters
            cret = sv_clusterer.check_mergeability(tchr1, tstart1, tend1, tdir1, tchr2, tstart2, tend2, tdir2, 
                treadid, tsize, tinfo1, tinfo2)           

            # create new SV cluster when new SV key cannot be merged into any existing SV clusteres
            if cret is None:
                sv_clusterer.add_new_cluster(tchr1, tstart1, tend1, tdir1, tchr2, tstart2, tend2, tdir2,
                    treadid, tsize, tinfo1, tinfo2)

            # check whether the number of locally clustered SVs and if it exceeds the threshould, 
            # then we ignore those SV clusteres and skip these regions.
            sv_clusterer.check_exceeding_local_variant_num(tchr1, tend1)


        sv_clusterer.flush_sv_cluster_list("EOF", 0)


if __name__ == "__main__":

    import sys
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    svtype = sys.argv[3]
    control = sys.argv[4]

    cluster_supporting_reads(input_file, output_file, svtype, control_junction_bedpe = control)
             
