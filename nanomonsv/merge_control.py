#! /usr/bin/env python3

import sys, gzip, os, subprocess
from .logger import get_logger

logger = get_logger(__name__)


def merge_bed_by_sv(bedpe_file, hout, prefix_label, is_bedpe = False):

    temp_key, temp_count = None, 0
    with gzip.open(bedpe_file, 'rt') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if is_bedpe:
                key = '\t'.join(F[:6] + ['*'] + F[7:10])
            else:
                key = '\t'.join(F[:3] + ['*'] + F[4:6])
            if key != temp_key:
                if temp_key is not None:
                    print(f'{temp_key}\t{temp_count}\t{prefix_label}', file = hout)
                temp_key = key 
                temp_count = 0
            temp_count = temp_count + 1

        if temp_key is not None:
            print(f'{temp_key}\t{temp_count}\t{prefix_label}', file = hout)


def merge_bed_by_sample(input_bedpe, output_bedpe, is_bedpe = False):

    temp_key, temp_counts, temp_samples = None, [], []
    with open(input_bedpe, 'r') as hin, open(output_bedpe, 'w') as hout:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if is_bedpe:
                key = '\t'.join(F[:10])
            else:
                key = '\t'.join(F[:6])
            if key != temp_key:
                if temp_key is not None:
                    temp_counts_print_line = ','.join(temp_counts)
                    temp_samples_print_line = ','.join(temp_samples)
                    print(f'{temp_key}\t{temp_counts_print_line}\t{temp_samples_print_line}', file = hout)
                temp_key, temp_counts, temp_samples = key, [], []

            if is_bedpe:
                temp_counts.append(F[10])
                temp_samples.append(F[11])
            else:
                temp_counts.append(F[6])
                temp_samples.append(F[7])

        if temp_key is not None:
            temp_counts_print_line = ','.join(temp_counts)
            temp_samples_print_line = ','.join(temp_samples)
            print(f'{temp_key}\t{temp_counts_print_line}\t{temp_samples_print_line}', file = hout)


def sort_bed(input_bed, output_bed, is_bedpe = False):

    with open(output_bed, 'w') as hout:
        if is_bedpe:
            subprocess.check_call(["sort", "-k1,1", "-k2,2n", "-k3,3n", "-k4,4", "-k5,5n", "-k6,6n", input_bed], stdout = hout)
        else:
            subprocess.check_call(["sort", "-k1,1", "-k2,2n", "-k3,3n", input_bed], stdout = hout)


def bgzip_tabix(input_bed, output_bedgz):

    with open(output_bedgz, 'w') as hout: 
        subprocess.check_call(["bgzip", "-f", "-c", input_bed], stdout = hout)
     
    subprocess.check_call(["tabix", "-p", "bed", output_bedgz])


def merge_control_from_parse_files(prefix_list_file, output_prefix):

    output_dir = os.path.dirname(output_prefix)
    if output_dir != '' and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open(prefix_list_file, 'r') as hin:
        for line in hin:
            prefix = line.rstrip('\n')
    
            if not os.path.exists(prefix + ".rearrangement.sorted.bedpe.gz"):
                logger.error(f'No file: {prefix}/.rearrangement.sorted.bedpe.gz')
                sys.exit(1)
            if not os.path.exists(prefix + ".insertion.sorted.bed.gz"):
                logger.error(f'No file: {prefix}/.insertion.sorted.bed.gz')
                sys.exit(1)
            if not os.path.exists(prefix + ".deletion.sorted.bed.gz"):
                logger.error(f'No file: {prefix}/.deletion.sorted.bed.gz')
                sys.exit(1)
            if not os.path.exists(prefix + ".bp_info.sorted.bed.gz"):
                logger.error(f'No file: {prefix}/.bp_info.sorted.bed.gz')
                sys.exit(1)

    hout_r = open(output_prefix + ".rearrangement.sorted.bedpe.gz.tmp1", 'w')
    hout_i = open(output_prefix + ".insertion.sorted.bed.gz.tmp1", 'w')
    hout_d = open(output_prefix + ".deletion.sorted.bed.gz.tmp1", 'w')
    hout_b = open(output_prefix + ".bp_info.sorted.bed.gz.tmp1", 'w')
    prefix_ind = 0
    with open(prefix_list_file, 'r') as hin:
        for line in hin:
            prefix = line.rstrip('\n')
            prefix_label = os.path.basename(prefix) + '_' + str(prefix_ind)

            merge_bed_by_sv(prefix + ".rearrangement.sorted.bedpe.gz", hout_r, prefix_label, is_bedpe = True)
            merge_bed_by_sv(prefix + ".insertion.sorted.bed.gz", hout_i, prefix_label)
            merge_bed_by_sv(prefix + ".deletion.sorted.bed.gz", hout_d, prefix_label)
            merge_bed_by_sv(prefix + ".bp_info.sorted.bed.gz", hout_b, prefix_label)

            prefix_ind = prefix_ind + 1


    hout_r.close()
    hout_i.close()
    hout_d.close()
    hout_b.close()

    sort_bed(output_prefix + ".rearrangement.sorted.bedpe.gz.tmp1", output_prefix + ".rearrangement.sorted.bedpe.gz.tmp2", is_bedpe = True)
    sort_bed(output_prefix + ".insertion.sorted.bed.gz.tmp1", output_prefix + ".insertion.sorted.bed.gz.tmp2")
    sort_bed(output_prefix + ".deletion.sorted.bed.gz.tmp1", output_prefix + ".deletion.sorted.bed.gz.tmp2")
    sort_bed(output_prefix + ".bp_info.sorted.bed.gz.tmp1", output_prefix + ".bp_info.sorted.bed.gz.tmp2")

    merge_bed_by_sample(output_prefix + ".rearrangement.sorted.bedpe.gz.tmp2", output_prefix + ".rearrangement.sorted.bedpe.gz.tmp3", is_bedpe = True)
    merge_bed_by_sample(output_prefix + ".insertion.sorted.bed.gz.tmp2", output_prefix + ".insertion.sorted.bed.gz.tmp3")
    merge_bed_by_sample(output_prefix + ".deletion.sorted.bed.gz.tmp2", output_prefix + ".deletion.sorted.bed.gz.tmp3")
    merge_bed_by_sample(output_prefix + ".bp_info.sorted.bed.gz.tmp2", output_prefix + ".bp_info.sorted.bed.gz.tmp3")

    bgzip_tabix(output_prefix + ".rearrangement.sorted.bedpe.gz.tmp3", output_prefix + ".rearrangement.sorted.bedpe.gz")
    bgzip_tabix(output_prefix + ".insertion.sorted.bed.gz.tmp3", output_prefix + ".insertion.sorted.bed.gz")
    bgzip_tabix(output_prefix + ".deletion.sorted.bed.gz.tmp3", output_prefix + ".deletion.sorted.bed.gz")
    bgzip_tabix(output_prefix + ".bp_info.sorted.bed.gz.tmp3", output_prefix + ".bp_info.sorted.bed.gz")

    os.remove(output_prefix + ".rearrangement.sorted.bedpe.gz.tmp1")
    os.remove(output_prefix + ".rearrangement.sorted.bedpe.gz.tmp2")
    os.remove(output_prefix + ".rearrangement.sorted.bedpe.gz.tmp3")
    os.remove(output_prefix + ".insertion.sorted.bed.gz.tmp1")
    os.remove(output_prefix + ".insertion.sorted.bed.gz.tmp2")
    os.remove(output_prefix + ".insertion.sorted.bed.gz.tmp3")
    os.remove(output_prefix + ".deletion.sorted.bed.gz.tmp1")
    os.remove(output_prefix + ".deletion.sorted.bed.gz.tmp2")
    os.remove(output_prefix + ".deletion.sorted.bed.gz.tmp3")
    os.remove(output_prefix + ".bp_info.sorted.bed.gz.tmp1")
    os.remove(output_prefix + ".bp_info.sorted.bed.gz.tmp2")
    os.remove(output_prefix + ".bp_info.sorted.bed.gz.tmp3")


