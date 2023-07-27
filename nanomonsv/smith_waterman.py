#! /usr/bin/env python3

import numpy as np

def sw_jump(contig, region1_seq, region2_seq, match_score = 1, mismatch_penalty = 2, gap_cost = 2, jump_cost = 2, minimum_score = 10):

    H1 = np.zeros((len(contig) + 1, len(region1_seq) + 1), np.int32)
    H1_path = np.zeros((len(contig) + 1, len(region1_seq) + 1), np.int32)
    H2 = np.zeros((len(contig) + 1, len(region2_seq) + 1), np.int32)
    H2_path = np.zeros((len(contig) + 1, len(region2_seq) + 1), np.int32)
    H2_origin_i = np.zeros((len(contig) + 1, len(region2_seq) + 1), np.int32)
    H2_origin_j = np.zeros((len(contig) + 1, len(region2_seq) + 1), np.int32)
    H2_origin_score = np.zeros((len(contig) + 1, len(region2_seq) + 1), np.int32)

    
    for i in range(1, H1.shape[0]):
        for j in range(1, H1.shape[1]):
            match = H1[i - 1, j - 1] + (match_score if contig[i - 1] == region1_seq[j - 1] else - mismatch_penalty)
            delete = H1[i - 1, j] - gap_cost
            insert = H1[i, j - 1] - gap_cost
            tscores = (float("-inf"), match, delete, insert)
            H1[i, j] = max(tscores)
            H1_path[i,j] = tscores.index(max(tscores))

    H1_max = np.amax(H1, axis = 1)
    H1_argmax = np.argmax(H1, axis = 1)

    for i in range(1, H2.shape[0]):
        for j in range(1, H2.shape[1]):
            match = H2[i - 1, j - 1] + (match_score if contig[i - 1] == region2_seq[j - 1] else - mismatch_penalty)
            delete = H2[i - 1, j] - gap_cost
            insert = H2[i, j - 1] - gap_cost

            jump = float("-inf")
            org_i = np.argmax(H1_max[:i])
            if H1_max[org_i] >= minimum_score:

                jump = H1_max[org_i] + (match_score if contig[i - 1] == region2_seq[j - 1] else - mismatch_penalty)
                jump_diff = i - org_i
                jump = jump - jump_cost * np.log(jump_diff)
 
            tscores = (float("-inf"), match, delete, insert, jump)
            H2[i, j] = max(tscores)
            H2_path[i, j] = tscores.index(max(tscores))

            if H2_path[i, j] == 4:
                H2_origin_i[i, j] = np.argmax(H1_max[:i])
                H2_origin_j[i, j]  = H1_argmax[H2_origin_i[i, j]]
                H2_origin_score[i, j] = H1_max[org_i]

    # i_cur, j2_cur = np.unravel_index(H2.argmax(), H2.shape)
    # a_string, b_string = a[i_end - 1], b[j_end - 1]
    if np.max(H2[H2.shape[0] - 1, :]) >= np.max(H2[:, H2.shape[1] - 1]):
        j2_cur = np.unravel_index(H2[H2.shape[0] - 1, :].argmax(), (1, H2.shape[1]))[1]
        i_cur = H2.shape[0] - 1
    else:
        j2_cur = H2.shape[1] - 1 
        i_cur = np.unravel_index(H2[:, H2.shape[1] - 1].argmax(), (H2.shape[0], 1))[0]


    contig_match, region1_seq_match, region2_seq_match = '', '', ''
    i_end, i2_end, j2_end = i_cur, i_cur, j2_cur # one-based position

    match_count, deletion_count, insertion_count, jump_count = 0, 0, 0, 0
    H1_score = 0
    scores = []

    while i_cur > 0 and j2_cur > 0:
        scores.append(H2[i_cur, j2_cur])
        if H2_path[i_cur, j2_cur] == 1:
            contig_match, region2_seq_match = contig[i_cur - 1] + contig_match, region2_seq[j2_cur - 1] + region2_seq_match
            i_cur, j2_cur = i_cur - 1, j2_cur - 1
            match_count = match_count + 1
        elif H2_path[i_cur, j2_cur] == 2:
            contig_match, region2_seq_match = contig[i_cur - 1] + contig_match, '-' + region2_seq_match
            i_cur, j2_cur = i_cur - 1, j2_cur
            deletion_count = deletion_count + 1
        elif H2_path[i_cur, j2_cur] == 3: 
            contig_match, region2_seq_match = '-' + contig_match, region2_seq[j2_cur - 1] + region2_seq_match
            i_cur, j2_cur = i_cur, j2_cur - 1
            insertion_count = insertion_count + 1
        elif H2_path[i_cur, j2_cur] ==  4:
            contig_match, region2_seq_match = contig[i_cur - 1] + contig_match, region2_seq[j2_cur - 1] + region2_seq_match
            i1_end = H2_origin_i[i_cur, j2_cur]
            j1_end = H2_origin_j[i_cur, j2_cur]
            H1_score = H2_origin_score[i_cur, j2_cur]
            i2_start = i_cur
            j2_start = j2_cur # one-based alignment start of second fragment
            """
            i2_start, i1_end = i_cur, i_cur - 1
            j2_start = j2_cur # one-based alignment start of second fragment
            i_cur = i_cur - 1
            """
            jump_count = jump_count + 1
            break
        elif H2_path[i_cur, j2_cur] == 0:
            break

    if H1_score < minimum_score or H2.max() - H1_score < minimum_score or jump_count == 0:
        return(None)

    contig_match = ' ' + contig_match
    # j1_end = H1_argmax[i_cur]
    j1_cur = j1_end
    i_cur = i1_end

    while i_cur > 0 and j1_cur > 0:
        scores.append(H1[i_cur, j1_cur])
        if H1_path[i_cur, j1_cur] == 1:
            contig_match, region1_seq_match = contig[i_cur - 1] + contig_match, region1_seq[j1_cur - 1] + region1_seq_match 
            i_cur, j1_cur = i_cur - 1, j1_cur - 1
            match_count = match_count + 1
        elif H1_path[i_cur, j1_cur] == 2:
            contig_match, region1_seq_match = contig[i_cur - 1] + contig_match, '-' + region1_seq_match
            i_cur, j1_cur = i_cur - 1, j1_cur
            deletion_count = deletion_count + 1
        elif H1_path[i_cur, j1_cur] == 3:
            contig_match, region1_seq_match = '-' + contig_match, region1_seq[j1_cur - 1] + region1_seq_match
            i_cur, j1_cur = i_cur, j1_cur - 1
            insertion_count = insertion_count + 1
        elif H1_path[i_cur, j1_cur] == 0:
            break

    i1_start = i_cur
    j1_start = j1_cur + 1

    return(H2.max(), (i1_start, i1_end, i2_start, i2_end), (j1_start, j1_end), (j2_start, j2_end), contig_match, region1_seq_match + ' ' + region2_seq_match)



if __name__ == "__main__":

    A = "AGGAGGGGGTTTCAGATTCCCCCATCCAAAGGAGGTTTTGGCCCCATGGGGAATGAAACAAAGACACATGAACACAGCTGCAAATATGCCCCACCTCTCCCTGCCCCCCCAACTTGCCGGGCCAATGCACCTGCAGTCCCTCACCCTTCTCAGTTCATACCTACTCCCCAGCTCTCTCACTCTTGTATTCTAGATCCCTCAGTACACCAGCGCCCGGGAAGTCAGGGCCGGGCCAGGACTTTGACGTGCTGGTGGCAGACTGAAGGCCAACTCCCATGGGGAGTCCCCAAGGAGAAGGTGGAATCAGGGAAAAAAGGTGAATGGATGACTCAGTCACAGATGGAGCCACAAAGATGAAGCCATTCATTCCACCCACCACCACCACCACAAAAAAAGAACATCTGAGCCCCCCATGAGTCAGGGATGAGCTGGTGCACCTGCCCTCAGGGGCAGCCAATAGAAAAGTGCACCCCTAATGATCTGAAGTAACTTGAAATCACCAGTCATCCTACCCAAAATGCAAACTCAGTGGTTTCCCTCACCTACTCAAACCCTTTAAGGATGCAAAGCACATCTGGTTTTGCTTGCATT"

    B1 = "CAAAGGGGAGTCTCCCAAGGAGAAGAGCCCAGGGCGCAAGGAGCAA"
    B2 = "TTACAATATTATCTGATGCTGAAAGGTGGGGATCAGGGAAAAG"

    score, region1, region2, region3, A_match, B_match = sw_jump(A, B1, B2)

    print(score)
    print(region1)
    print(region2)
    print(region3)
    print(A_match)
    print(B_match)


