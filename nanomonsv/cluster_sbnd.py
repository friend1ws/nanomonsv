#! /usr/bin/env python3

import sys, gzip, statistics, logging
import pysam

from nanomonsv.logger import get_logger

logger = get_logger(__name__)

class Sbnd_cluster(object):

    def __init__(self, tchr, tstart, tend, tdir, treadid, tinfo):

        self.chr = tchr
        self.start = tstart
        self.end = tend
        self.dir = tdir
        self.readids = [treadid]
        self.info = [tinfo]


class Sbnd_clusterer(object):

    def __init__(self, output_file, control_bed = None, control_panel_bed = None,
        cluster_margin_size = None, read_num_thres = None, median_mapQ_thres = None):

        self.hout = open(output_file, 'w')
        self.control_tb = None
        if control_bed is not None:
            self.control_tb = pysam.TabixFile(control_bed)

        self.control_panel_tb = None
        if control_panel_bed is not None:
            self.control_panel_tb = pysam.TabixFile(control_panel_bed)

        self.sbnd_cluster_list = []
        self.cluster_margin_size = cluster_margin_size
        self.read_num_thres = read_num_thres
        self.median_mapQ_thres = median_mapQ_thres


    def __del__(self):
        self.hout.close()
        if self.control_tb is not None: self.control_tb.close()
        if self.control_panel_tb is not None: self.control_panel_tb.close()

 
    def check_mergeability(self, tchr, tstart, tend, tdir, treadid, tinfo):
   
        for cluster in self.sbnd_cluster_list:

            if treadid in cluster.readids: continue

            # pairs of chromosome and direction should be the same.
            # regions from start to end should overlap.   
            is_mergeable = False
            if tchr == cluster.chr and tdir == cluster.dir and tend >= cluster.start and tstart <= cluster.end:
                is_mergeable = True

            if not is_mergeable: continue

            cluster.start = min(tstart, cluster.start)
            cluster.end = max(tend, cluster.end)

            cluster.readids.append(treadid)
            cluster.info.append(tinfo)
                        
            return cluster
            
        return None


    def add_new_cluster(self, tchr, tstart, tend, tdir, treadid, tinfo):

        self.sbnd_cluster_list.append(
            Sbnd_cluster(tchr, tstart, tend, tdir, treadid, tinfo))


    def filter_single_breakend_cluster(self, cl):
        
        is_filter = False
        if len(cl.readids) < self.read_num_thres: is_filter = True

        median_mapQ = statistics.median([int(x.split(',')[5]) for x in cl.info])
        if median_mapQ < self.median_mapQ_thres: is_filter = True

        if is_filter == True: return(True)

        control_flag = False
        if self.control_tb is not None:

            tabix_error_flag = False
            try:
                records = self.control_tb.fetch(cl.chr, max(0, cl.start - 10), cl.end + 10)
            except Exception as e:
                logger.debug(f'{e}')
                tabix_error_flag = True

            if not tabix_error_flag:
                for record_line in records:
                    rec = record_line.split('\t')

                    if cl.chr != rec[0] or cl.dir != rec[5]: continue
        
                    if cl.end >= int(rec[2]) and cl.start <= int(rec[2]):
                        control_flag = True

        control_panel_flag = False
        if self.control_panel_tb is not None:

            tabix_error_flag = False
            try:
                records = self.control_panel_tb.fetch(cl.chr, max(0, cl.start - 10), cl.end + 10)
            except Exception as e:
                logger.debug(f'{e}')
                tabix_error_flag = True
                
            if not tabix_error_flag:
                for record_line in records:
                    rec = record_line.split('\t')
                    
                    if cl.chr != rec[0] or cl.dir != rec[5]: continue
                    
                    if cl.end >= int(rec[2]) and cl.start <= int(rec[2]):
                        control_panel_flag = True
            

        if control_flag == True or control_panel_flag == True: is_filter = True 

        return(is_filter)


    def flush_sv_cluster_list(self, current_chr, current_pos):

        remove_cluster = []
        for i in range(len(self.sbnd_cluster_list)):
            cl = self.sbnd_cluster_list[i]
            if current_chr != cl.chr or current_pos > cl.end + self.cluster_margin_size:
                remove_cluster.append(cl)

                if self.filter_single_breakend_cluster(cl): continue
                print_line_readids = ';'.join(cl.readids)
                print_line_info = ';'.join(cl.info)
        
                print(f'{cl.chr}\t{cl.start}\t{cl.end}\t{print_line_readids}\t0\t{cl.dir}\t{print_line_info}',
                    file = self.hout)
                        
        for cl in remove_cluster:
            self.sbnd_cluster_list.remove(cl)


def cluster_supporting_reads_sbnd(input_file, output_file, control_bed = None, control_panel_bed = None, 
    cluster_margin_size = 100, sbnd_cluster_margin_size = 20, 
    read_num_thres = 3, median_mapQ_thres = 20, debug = False):

    if debug: logger.setLevel(logging.DEBUG)


    sbnd_clusterer = Sbnd_clusterer(output_file, control_bed = control_bed, control_panel_bed = control_panel_bed,
        cluster_margin_size = cluster_margin_size, read_num_thres = read_num_thres, 
        median_mapQ_thres = median_mapQ_thres)
    
    with gzip.open(input_file, 'rt') as hin: 
        for line in hin:

            tchr, _, tpos, treadid, _, tdir, tinfo = line.rstrip('\n').split('\t')
            tstart = max(int(tpos) - sbnd_cluster_margin_size, 0)
            tend = int(tpos) + sbnd_cluster_margin_size

            # skip single breakend from secondary alignment
            if tinfo.split(',')[9] == "True": continue

            # flush out existing cluster whose supporting reads have parsed alreadly 
            #   (considering the chromosome and coordinates).
            sbnd_clusterer.flush_sv_cluster_list(tchr, int(tpos))

            # check the new single breakend can be merged into any existing clusters
            cret = sbnd_clusterer.check_mergeability(tchr, tstart, tend, tdir, treadid, tinfo)

            # create new single breakend cluster when new key cannot be merged into any existing clusteres
            if cret is None:
                sbnd_clusterer.add_new_cluster(tchr, tstart, tend, tdir, treadid, tinfo)
 

        sbnd_clusterer.flush_sv_cluster_list("EOF", 0)


if __name__ == "__main__":

    import sys
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    control = sys.argv[3]

    cluster_supporting_reads_single_breakend(input_file, output_file, control_bed = control)
             
