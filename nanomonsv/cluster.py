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


    def __init__(self, cluster_margin_size = None, size_margin_ratio = None):
        self.sv_cluster_list = []
        self.cluster_margin_size = cluster_margin_size
        self.size_margin_ratio = size_margin_ratio


    def check_meageability(self, tchr1, tstart1, tend1, tdir1, tchr2, tstart2, tend2, tdir2, 
        treadid, tsize, tinfo1, tinfo2):
    
        for cluster in self.sv_cluster_list:

            if treadid in cluster.readids: continue

            # pairs of chromosome and direction should be the same.
            # two pairs of start and end should overlap.    
            if tchr1 == cluster.chr1 and tchr2 == cluster.chr2 and tdir1 == cluster.dir1 and tdir2 == cluster.dir2 and \
                tend1 >= cluster.start1 and tstart1 <= cluster.end1 and tend2 >= cluster.start2 and tstart2 <= cluster.end2: 

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

    
    def flush_sv_cluster_list(self, hout, current_chr, current_pos):

        remove_cluster = []
        for i in range(len(self.sv_cluster_list)):
            cl = self.sv_cluster_list[i]
            if cl.chr1 != current_chr or cl.end1 > current_pos + self.cluster_margin_size:
                remove_cluster.append(cl)

                print_line_readids = ';'.join(cl.readids)
                print_line_info1 = ';'.join(cl.info1)
                print_line_info2 = ';'.join(cl.info2)
 
                print(f'{cl.chr1}\t{cl.start1}\t{cl.end1}\t{cl.chr2}\t{cl.start2}\t{cl.end2}\t' +
                    f'{print_line_readids}\t0\t{cl.dir1}\t{cl.dir2}\t{print_line_info1}\t{print_line_info2}',
                    file = hout)

        for cl in remove_cluster:
            self.sv_cluster_list.remove(cl)

                

def cluster_rearrangement(input_file, output_file, is_indel = False, cluster_margin_size = 100,
    indel_cluster_margin_size = 10, size_margin_ratio = 0.2, maximum_local_variant_num = 1000, skip_margin = 5000):

    sv_clusterer = Sv_clusterer(cluster_margin_size, size_margin_ratio = size_margin_ratio)
    local_variant_num = 0
    skip_pos = 0

    with gzip.open(input_file, 'rt') as hin, open(output_file, 'w') as hout:
        for line in hin:
            if is_indel: # insertion or deletion
                tchr1, tstart, tend, treadid, tsize, tdir, tinfo1 = ne.rstrip('\n').split('\t')
                tstart1 = int(tstart) - indel_cluster_margin_size
                tend1 = int(start) + indel_cluster_margin_size
                tstart2 = int(tend) - indel_cluster_margin_size
                tend2 = int(tend) + indel_cluster_margin_size
                tchr2, tdir1, tdir2 = tchr1, '+', '-'
                tsize, tinfo2 = int(tsize), None
            else: # rearrangement
                tchr1, tstart1, tend1, tchr2, tstart2, tend2, treadid, _, tdir1, tdir2, tinfo1, tinfo2 = \
                    line.rstrip('\n').split('\t')
                tstart1, tend1, tstart2, tend2 = int(tstart1), int(tend1), int(tstart2), int(tend2)
                tsize = None
 
            if tend2 < skip_pos: continue

            # flush out existing cluster SVs wholse supporting reads have parsed alreadly 
            #   (considering the chromosome and coordinates).
            sv_clusterer.flush_sv_cluster_list(hout, tchr1, tend2)

            # check the new SV can be merged into any existing SV clusters
            cret = sv_clusterer.check_meageability(tchr1, tstart1, tend1, tdir1, tchr2, tstart2, tend2, tdir2, 
                treadid, tsize, tinfo1, tinfo2)           

            # create new SV cluster when new SV key cannot be merged into any existing SV clusteres
            if cret is None:
                sv_clusterer.add_new_cluster(tchr1, tstart1, tend1, tdir1, tchr2, tstart2, tend2, tdir2,
                    treadid, tsize, tinfo1, tinfo2)

            local_variant_num = local_variant_num + 1
            if local_variant_num > maximum_local_variant_num:
                print(f'Exceeded maximum number of local variants at {tchr1}:{tend}', file = sys.stderr)
                print(f'Skip {tchr1}:{tend + skip_margin}', file = sys.stderr)
                sv_clusterer.sv_cluster_list = []
                local_variant_num = 0
                skip_pos = tend1 + skip_margin

        sv_clusterer.flush_sv_cluster_list(hout, "EOF", 0)

"""
def cluster_insertion_deletion(input_file, output_file, deletion_cluster_margin_size = 10, check_margin_size = 50, 
    size_margin_ratio = 0.2, maximum_unique_pairs = 100, maximum_local_variant_num = 1000):

    dcms = deletion_cluster_margin_size
    sv_clusterer = Sv_clusterer(cluster_margin_size, size_margin_ratio = size_margin_ratio)
    local_variant_num = 0
    skip_pos = 0

    with gzip.open(input_file, 'rt') as hin, open(output_file, 'w') as hout:
        for line in hin:
            tchr1, tstart, tend, tchr2, treadid, tsize, tdir1, tdir2, tinfo1 = \
                line.rstrip('\n').split('\t')
            tstart1, tend1, tstart2, tend2 = int(tstart1), int(tend1), int(tstart2), int(tend2)

            if tend
            tsize = None

            # flush out existing cluster SVs wholse supporting reads have parsed alreadly 
            #   (considering the chromosome and coordinates).
            sv_clusterer.flush_sv_cluster_list(hout, tchr1, tend2)

            # check the new SV can be merged into any existing SV clusters
            cret = sv_clusterer.check_meageability(tchr1, tstart - dcms, tstart + dcms, '+', 
                tchr2, tend - dcms, tend + dcms, '-', treadid, tsize, tinfo1)

            # create new SV cluster when new SV key cannot be merged into any existing SV clusteres
            if cret is None:
                sv_clusterer.add_new_cluster(tchr1, tstart1- dcms, tend1 + dcms, tdir1, 
                    tchr2, tstart2 - dcms, tend2 + dcms, tdir2, treadid, tsize, tinfo1, None)

            sv_clusterer.variant_num_without_any_flush = sv_clusterer.variant_num_without_any_flush + 1
            if sv_clusterer.variant_num_without_any_flush > maximum_local_variant_num:
                print("Exceeded maximum number of local variants at %s:%s" % (F[0], F[1]), file = sys.stderr)
                print("Skip %s:%s" % (F[0], str(int(F[1]) + 10 * check_margin_size)), file = sys.stderr)
                sv_clusterer.sv_cluster_list = []
                local_variant_num = 0
                skip_pos = int(F[1]) + 10 * check_margin_size
"""

if __name__ == "__main__":

    import sys
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    cluster_rearrangement(input_file, output_file)
             
