import re
import sys
import json
import argparse
import multiprocessing
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from collections import defaultdict
from util import _fa_or_fq, _load_fa, _load_fq, _load_full_fa, _load_full_fq

def phased_ng50(fname, column, genome_size):
    lens = [int(line.split()[column]) for line in open(fname)]
    gs = genome_size if genome_size > 0 else sum(lens)
    lens.sort(key=lambda x: -x)
    acuu = 0
    for l in lens:
        acuu += l
        if acuu >= gs / 2:
            print('phased block NG50 {} total blcok size {}'.format(l, sum(lens)), file = sys.stderr)
            break

def NX_stats(fname):
    length = []
    ftype = _fa_or_fq(fname)
    if ftype == None:
        print('Unrecognized file type: {}'.format(fname), file = sys.stderr)
        exit(1)
    for _, seq in eval('_load_' + ftype)(fname):
        length.append(len(seq))
    length.sort(reverse = True)
    total = sum(length)
    accu = 0
    n = [10, 25, 50, 75, 90]
    nx, lx = [], []
    index = 0
    for i in range(len(n)):
        while accu < total * n[i] / 100:
            accu += length[index]
            index += 1
        nx.append(length[index])
        lx.append(index)
    print('Total base {}, number of contig {}'.format(total, len(length)), file = sys.stderr)
    print('max length {} min length {}'.format(max(length), min(length)), file = sys.stderr)
    for i in range(len(n)):
        print('N{} {} L{} {}'.format(n[i], nx[i], n[i], lx[i]), file = sys.stderr)

def hamming_error_eval(fname):
    count = defaultdict(lambda: [0, 0])
    for i, line in enumerate(open(fname)):
        if i == 0: continue
        items = line.split()
        ss = (int(items[2]), int(items[3]))
        count[items[0]][0] += max(ss)
        count[items[0]][1] += min(ss)
    for k, v in count.items():
        print('{} {:>.6f}%'.format(k, 100 * v[1] / (v[0] + v[1])), file = sys.stderr)

def switch_error_eval(fname, short_num, short_limit):
    scaffold, prev_scaffold, name, prev_block = '', '', '', ''
    start, end, block_start, block_end = -1, -1, -1, -1
    short_start, short_end = -1, -1
    num_markers, num_switch, num_short_switch = 0, 0, 0
    distance_from_switch = -1
    is_short_switch = True
    hap1, hap2 = '', ''

    total_switch, total_marker = 0, 0
    short_switch = [0] * 10
    switches = []
    blocks = []
    blocks_range = [0] * 4

    for line in open(fname):
        items = line.strip().split()
        name = items[3]
        if name != 'gap':
            if hap1 == '':
                hap1 = name
                continue
            elif name != hap1:
                hap2 = name
                break
    print('hap1 {} hap2 {}'.format(hap1, hap2), file = sys.stderr)
    for line in open(fname):
        items = line.strip().split()
        scaffold = items[0]
        start = int(items[1])
        end = int(items[2])
        name = items[3]
        if prev_scaffold == '' or scaffold != prev_scaffold:
            if prev_scaffold != '':
                num_markers -= num_short_switch
                if num_markers > 1:
                    # print('{}\t{}\t{}\t{}\t{:>.3f}\t{}\t{}'.format(prev_scaffold, block_start, block_end, prev_block, num_switch/num_markers, num_switch, num_markers))
                    total_marker += num_markers
                    total_switch += num_switch
                if num_short_switch > 0:
                    block_start = short_start
                    block_end = short_end
                    num_markers = num_short_switch
                    # switches.append(num_short_switch)
                    prev_block = hap1 if name == hap1 else hap2
                    if num_markers > 1:
                        # print('{}\t{}\t{}\t{}\t{:>.3f}\t{}\t{}'.format(prev_scaffold, block_start, block_end, prev_block, num_switch/num_markers, num_switch, num_markers))
                        total_marker += num_markers
                        # total_switch += num_switch
            prev_block = name
            num_switch = 0
            num_short_switch = 0
            distance_from_switch = 0
            num_markers = 0
            block_start = start
            block_end = end
            if name == hap1 or name == hap2:
                num_markers += 1
            is_short_switch = False
        else:
            num_markers += 1
            if name == prev_block or prev_block == 'unknown':
                block_end = end
                if prev_block == 'unknown':
                    prev_block = name
                if is_short_switch:
                    num_switch += num_short_switch
                    switches.append(num_short_switch)
                    short_switch[num_short_switch // 10] += 1
                    block_size = short_end - short_start
                    blocks.append((block_size - 20) // num_short_switch)
                    block_size = (block_size - 20) // num_short_switch
                    if block_size == 1:
                        blocks_range[0] += 1
                    elif block_size < 50:
                        blocks_range[1] += 1
                    elif block_size < 100:
                        blocks_range[2] += 1
                    else:
                        blocks_range[3] += 1
                    num_short_switch = 0
                    is_short_switch = False
            elif name != prev_block:
                short_end = end
                num_short_switch += 1
                if not is_short_switch:
                    short_start = start
                else:
                    distance_from_switch = short_end - short_start
                    if num_short_switch >= short_num or distance_from_switch > short_limit:
                        num_markers -= num_short_switch
                        if num_markers > 1:
                            # print('{}\t{}\t{}\t{}\t{:>.3f}\t{}\t{}'.format(scaffold, block_start, block_end, prev_block, num_switch/num_markers, num_switch, num_markers))
                            total_marker += num_markers
                            total_switch += num_switch
                        block_start = short_start
                        block_end = end
                        num_markers = num_short_switch
                        num_switch = 0
                        num_short_switch = 0
                        prev_block = name
                        is_short_switch = False
                        continue
                is_short_switch = True
        prev_scaffold = scaffold
    if num_markers > 1:
        # print('{}\t{}\t{}\t{}\t{:>.3f}\t{}\t{}'.format(scaffold, block_start, block_end, prev_block, num_switch/num_markers, num_switch, num_markers))
        total_marker += num_markers
        total_switch += num_switch
    print('total_switch {} total_marker {} switch_error {:>.6f}'.format(total_switch, total_marker, total_switch * 100/total_marker), file = sys.stderr)
    print('max_switch {} min_switch {} total_switch {}'.format(max(switches), min(switches), sum(switches)), file = sys.stderr)
    print('max_block {} min_block {} total_block {}'.format(max(blocks), min(blocks), sum(blocks)), file = sys.stderr)
    for i in range(4):
        print('{} {}'.format(i * 50, blocks_range[i]), file = sys.stderr)
    # for i in range(10):
    #     print('{} {}'.format(i * 10, short_switch[i]), file = sys.stderr)

def phase_eval(rfname, b1fname, b2fname):
    hap1, hap2 = {}, {}
    contig = {}
    for line in open(rfname):
        items = line.strip().split()
        hap1.setdefault(items[0], []).append(items[1])
        hap2.setdefault(items[0], []).append(items[2])
    h1, h2 = '', ''
    for file in [b1fname, b2fname]:
        if h1 != '' and h2 != '':
            break
        for line in open(file):
            items = line.strip().split()
            if h1 == '':
                h1 = items[3]
            elif h1 != items[3]:
                h2 = items[3]
                break
    for file in [b1fname, b2fname]:
        prev = ''
        paternal, maternal = 0, 0
        for line in open(file):
            items = line.strip().split()
            if prev == '' or prev != items[0]:
                if prev != '':
                    contig[prev] = 0 if paternal > maternal else 1
                paternal, maternal = 0, 0
                prev = items[0]
                if items[3].startswith('paternal'):
                    maternal += int(items[5])
                    paternal += int(items[6])
                else:
                    paternal += int(items[5])
                    maternal += int(items[6])
        contig[prev] = 0 if paternal > maternal else 1
    phase = [0, 0]
    for hap in [hap1, hap2]:
        for ctg, blocks in hap.items():
            paternal, maternal = 0, 0
            for b in blocks:
                if b not in contig:
                    continue
                if contig[b] == 0:
                    paternal += 1
                else:
                    maternal += 1
            for b in blocks:
                if b not in contig:
                    continue
                if paternal > maternal and contig[b] == 1:
                    print('{} {} {} {}'.format(ctg, b, paternal, maternal), file = sys.stderr)
                elif paternal < maternal and contig[b] == 0:
                    print('{} {} {} {}'.format(ctg, b, paternal, maternal), file = sys.stderr)
            if paternal > maternal:
                phase[0] += paternal
                phase[1] += maternal
            else:
                phase[0] += maternal
                phase[1] += paternal
    print('correct phased {} incorrect phased {} phase error {:>.6f}'.format(phase[0], phase[1], phase[1] * 100/(phase[0] + phase[1])), file = sys.stderr)

def phase_eval1(rfname, cfname):
    contig = {}
    for line in open(cfname):
        if line.startswith('Assembly'): continue
        items = line.strip().split()
        contig[items[1]] = [int(items[2]), int(items[3])]
    hap1, hap2 = {}, {}
    for line in open(rfname):
        items = line.strip().split()
        hap1.setdefault(items[0], []).append(items[1])
        hap2.setdefault(items[0], []).append(items[2])
    for hap in [hap1, hap2]:
        phase = [0, 0]
        for ctg, blocks in hap.items():
            paternal, maternal = 0, 0
            for b in blocks:
                if b not in contig:
                    continue
                if contig[b][0] > contig[b][1]:
                    paternal += 1
                elif contig[b][0] < contig[b][1]:
                    maternal += 1
            if paternal > maternal:
                phase[0] += paternal
                phase[1] += maternal
            elif paternal < maternal:
                phase[0] += maternal
                phase[1] += paternal
        print('correct phased {} incorrect phased {} phase error {:>.6f}%'.format(phase[0], phase[1], phase[1] * 100/(phase[0] + phase[1])), file = sys.stderr)

def phase_shasta_eval(cfname, pfname):
    contig = {}
    for line in open(cfname):
        if not line.startswith('Assembly-Phased'): continue
        items = line.strip().split()
        contig[items[1]] = [int(items[2]), int(items[3])]
    hap1, hap2 = {}, {}
    for line in open(pfname):
        if not line.startswith('PR'): continue
        items = line.strip().split(',')
        c = '.'.join(items[0].split('.')[:2])
        if items[1] == '1':
            hap1.setdefault(c, []).append(items[0])
        elif items[1] == '-1':
            hap2.setdefault(c, []).append(items[0])
    
    for hap in [hap1, hap2]:
        phase = [0, 0]
        for ctg, blocks in hap.items():
            paternal, maternal = 0, 0
            for b in blocks:
                if b not in contig:
                    continue
                if contig[b][0] > contig[b][1]:
                    paternal += 1
                elif contig[b][0] < contig[b][1]:
                    maternal += 1
            if paternal > maternal:
                phase[0] += paternal
                phase[1] += maternal
            elif paternal < maternal:
                phase[0] += maternal
                phase[1] += paternal
        print('correct phased {} incorrect phased {} phase error {:>.6f}%'.format(phase[0], phase[1], phase[1] * 100/(phase[0] + phase[1] + 1)), file = sys.stderr)

def _load_matrix(mfname):
    header, matrix = [], []
    for line in open(mfname):
        if 'F' in line:
            if len(header) > 0 and len(matrix) > 0:
                yield header, matrix
            header = line.strip().split()
            matrix.clear()
        else:
            items = line.strip().split()
            matrix.append(items)
    yield header, matrix

def total_contacts(mfname):
    total, t_diag = 0, 0
    for h, m in _load_matrix(mfname):
        assert len(h) == len(m)
        for i in range(len(m)):
            for j in range(len(m)):
                if m[i][j] != m[j][i]:
                    print(h[i], h[j], m[i][j], m[j][i], file = sys.stderr)
                if i == j:
                    t_diag += int(m[i][j])
                else:
                    total += int(m[i][j])
    print('total contacts {} total diagonal {}'.format(total, t_diag), file = sys.stderr)

import pysam
def hic_contact_eval(cfname, h1fname, afname, mapq, filter):
    pair = {}
    for line in open(afname):
        items = line.strip().split()
        pair[items[0]] = items[1]
        pair[items[1]] = items[0]
    haplotigs = {}
    for line in open(cfname, 'r'):
        if line.startswith('Assembly'): continue
        items = line.strip().split()
        if int(items[2]) + int(items[3]) == 0:
            haplotigs[items[1]] = 0
        elif int(items[2]) > int(items[3]):
            haplotigs[items[1]] = 1
        else:
            haplotigs[items[1]] = -1
    total = 0
    prev = ''
    ctg1, ctg2 = '', ''
    mapq1, mapq2 = 0, 0
    mapping = []
    nm1, nm2 = 0, 0
    bf = pysam.AlignmentFile(h1fname, 'rb')
    for line in bf:
        if filter and line.mapping_quality < mapq: continue
        if line.query_name != prev:
            prev = line.query_name
            ctg1 = line.reference_name
            mapq1 = line.mapping_quality
            nm1 = line.get_tag('NM')
        else:
            total += 1
            ctg2 = line.reference_name
            mapq2 = line.mapping_quality
            nm2 = line.get_tag('NM')
            if ctg1[0: 7] == ctg2[0: 7] and ctg1 in pair and ctg2 in pair and ctg1 != ctg2 and ctg2 != pair[ctg1]:
                nm = 0
                if nm1 >= 5 or nm2 >= 5:
                    nm = 1
                mapping.append((ctg1, ctg2, mapq1, mapq2, nm))
    
    res = [0, 0]
    count = 0
    snp = [0, 0]
    zero = 0
    for ctg1, ctg2, mapq1, mapq2, nm in mapping:
        if ctg1 not in haplotigs or ctg2 not in haplotigs:
            res[0] += 1
            count += 1
            continue
        if haplotigs[ctg1] * haplotigs[ctg2] >= 0:
            res[0] += 1
            if haplotigs[ctg1] + haplotigs[ctg2] == 0:
                zero += 1
            if mapq1 < mapq or mapq2 < mapq or nm:
                snp[0] += 1
        else:
            res[1] += 1
            if mapq1 < mapq or mapq2 < mapq or nm:
                snp[1] += 1
    print('total {} correct {} incorrect {} count {} snp correct {} snp incorrect {}'.format(total, res[0], res[1], count, snp[0], snp[1]), file = sys.stderr)
    print('zero {}'.format(zero), file = sys.stderr)
    print('{}\t{}\t{}\t{}\t{}'.format(mapq, res[0], res[1], snp[0], snp[1]), file = sys.stdout)

def process_one_contig(name, ctg, length, block):
    contig = [0] * (length + 1)
    pattern = re.compile(r'(\d+)([MIDNSHP=X])')
    for pos, cigar_str in ctg:
        ref_pos = pos
        for match in pattern.finditer(cigar_str):
            l, op = int(match.group(1)), match.group(2)
            if op == 'S' or op == 'H':
                # ref_pos += l
                pass
            elif op == 'D':
                for i in range(l):
                    if ref_pos + i >= len(contig):
                        print('name {} pos {} ref_pos {} l {} contig length {}'.format(name, pos, ref_pos, l, len(contig)), file = sys.stderr)
                    if contig[ref_pos + i] < 0:
                        contig[ref_pos + i] = -1
                ref_pos += l
            elif op == 'X':
                for i in range(l):
                    if ref_pos + i >= len(contig):
                        print('name {} pos {} ref_pos {} l {} contig length {}'.format(name, pos, ref_pos, l, len(contig)), file = sys.stderr)
                    if contig[ref_pos + i] < 1:
                        contig[ref_pos + i] = 0
                    contig[ref_pos + i] = 0
                ref_pos += l
            elif op == '=':
                for i in range(l):
                    if ref_pos + i >= len(contig):
                        print('name {} pos {} ref_pos {} l {} contig length {}'.format(name, pos, ref_pos, l, len(contig)), file = sys.stderr)
                    contig[ref_pos + i] = 1
                ref_pos += l

    accu = []
    accu_match = 0
    # for i in range(block, length):
    #     accu_match = contig[i - block: i].count(1)
    #     accu.append(accu_match / block)
    for i in range(0, length, block):
        end = i + block if i + block < length else length
        accu_match = contig[i: end].count(1)
        accu.append(accu_match / (end - i))
    # for i in range(length):
    #     accu_len = i % block + 1
    #     if accu_len - 1 == 0:
    #         if i != 0:
    #             accu.append(accu_match / block)
    #         accu_match = 0
    #     if contig[i] == 1: # match
    #         accu_match += 1
            # accu.append(accu_match / accu_len)
        # elif contig[i] == 0: # mismatch
        #     accu.append(accu_match / accu_len)
        # else: # deletion
        #     accu.append(0)
        
    return name, accu

def identity_curve(fname, threads, thres, block, prefix):
    contigs, ctg_length = {}, {}
    for line in open(fname, 'r'):
        if line.startswith('@SQ'): 
            items = line.strip().split()
            name = items[1][3:]
            length = int(items[2][3:])
            ctg_length[name] = length
        elif line.startswith('@'): continue
        else:
            items = line.strip().split()
            if len(items) < 6: continue
            if items[2] == '*': continue
            pos = int(items[3]) - 1
            assert pos >= 0
            contigs.setdefault(items[2], []).append((pos, items[5]))
    
    results = []
    pool = multiprocessing.Pool(threads)
    for name, ctg in contigs.items():
        results.append(pool.apply_async(process_one_contig, (name, ctg, ctg_length[name], block,)))
    pool.close()
    pool.join()
    identity = {}
    with open(prefix + '.identity.txt', 'w') as fw:
        for res in results:
            n, accu = res.get()
            if len(accu) == 0: continue
            stat = [x < thres for x in accu]
            prev, prev_pos = stat[0], 0
            for i in range(1, len(stat)):
                if stat[i] != prev:
                    fw.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(n, prev_pos, i, prev, i - prev_pos, sum(accu[prev_pos:i]) / (i - prev_pos)))
                    identity.setdefault(n, []).append([prev, i - prev_pos])
                    prev = stat[i]
                    prev_pos = i
            identity.setdefault(n, []).append([prev, len(stat) - prev_pos])
            fw.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(n, prev_pos, len(stat), prev, len(stat) - prev_pos, sum(accu[prev_pos:]) / (len(stat) - prev_pos)))
    results.clear()

    region = {}
    for ctg, stat in identity.items():
        acc_true, acc_false = 0, 0
        for i in range(len(stat)):
            if stat[i][1] < 10:
                if stat[i][0] == True:
                    acc_true += stat[i][1]
                else:
                    acc_false += stat[i][1]
            else:
                tmp, tmp_len = stat[i][0], stat[i][1]
                if acc_true > acc_false and (acc_true > 10 or acc_false > 10):
                    region.setdefault(ctg, []).append((True, acc_true + acc_false))
                elif acc_false > acc_true and (acc_true > 10 or acc_false > 10):
                    region.setdefault(ctg, []).append((False, acc_true + acc_false))
                else:
                    tmp_len += acc_true + acc_false
                acc_true, acc_false = 0, 0
                region.setdefault(ctg, []).append((tmp, tmp_len))
        if acc_true > acc_false:
            region.setdefault(ctg, []).append((True, acc_true + acc_false))
        elif acc_false > acc_true:
            region.setdefault(ctg, []).append((False, acc_true + acc_false))
        else:
            region.setdefault(ctg, []).append((None, acc_true + acc_false))
    # print(region['contig_10_1'], file = sys.stderr)
    compact = {}
    for ctg, stat in region.items():
        prev = 0
        suml = stat[prev][1]
        if len(stat) < 3: continue
        for i in range(1, len(stat)):
            if stat[i][0] != stat[prev][0]:
                # if i == len(stat) - 1 and stat[i][1] < 10:
                #     suml += stat[i][1]
                compact.setdefault(ctg, []).append([stat[prev][0], suml])
                suml = stat[i][1]
                prev = i
            else:
                suml += stat[i][1]
        compact.setdefault(ctg, []).append([stat[prev][0], suml])
    for ctg, stat in compact.items():
        if len(stat) < 2: continue
        if stat[-1][1] < 10:
            stat[-2][1] += stat[-1][1]
            stat.pop()
    region.clear()
    # print(compact['contig_10_1'], file = sys.stderr)
    with open(prefix + '.region', 'w') as fw:
        for ctg, stat in compact.items():
            if len(stat) < 3: continue
            acc = 0
            for i in range(len(stat)):
                fw.write('{}\t{}\t{}\t{}\t{}\n'.format(ctg, acc, acc + stat[i][1], stat[i][0], stat[i][1]))
                acc += stat[i][1]
    # with open(prefix + '.region', 'w') as fw:
    #     for ctg, stat in identity.items():
    #         first_pos, sum_s, sum_e, tmp = 0, 0, 0, 0
    #         if stat[0][1] > 10:
    #             if stat[0][0] == True:
    #                 first = 'T'
    #             else:
    #                 first = 'F'
    #         else:
    #             first = 'N'
    #         for i in range(1, len(stat)):
    #             if stat[i][1] > 10 or i == len(stat) - 1:
    #                 if i != len(stat) - 1 or stat[i][1] > 10:
    #                     second = 'T' if stat[i][0] == True else 'F'
    #                 else:
    #                     second = 'N'
    #                 accu_false, accu_true = 0, 0
    #                 if stat[first_pos][1] < 10:
    #                     if stat[first_pos][0] == False:
    #                         accu_false += stat[first_pos][1]
    #                     else:
    #                         accu_true += stat[first_pos][1]
    #                 for j in range(first_pos + 1, i):
    #                     if stat[j][0] == False:
    #                         accu_false += stat[j][1]
    #                     else:
    #                         accu_true += stat[j][1]
    #                 middle = 'T' if accu_true > accu_false else 'F'
    #                 mid_len = 0
    #                 for j in range(first_pos + 1, i):
    #                     mid_len += stat[j][1]
    #                 if ctg == 'contig_10_1':
    #                     print('first {} second {} middle {} accu_true {} accu_false {} first_pos {} i {} mid_len {} tmp {}'.format(first, second, middle, accu_true, accu_false, first_pos, i, mid_len, tmp), file = sys.stderr)
    #                 if first == 'N':
    #                     if middle != second:
    #                         sum_e = sum_s + stat[first_pos][1] + mid_len + tmp
    #                         tmp = 0
    #                         if second == 'N':
    #                             sum_e += stat[i][1]
    #                         fw.write('{}\t{}\t{}\t{}\t{}\n'.format(ctg, sum_s, sum_e, middle, sum_e - sum_s))
    #                         first_pos = i
    #                         first = second
    #                         sum_s = sum_e
    #                     else:
    #                         tmp += mid_len + stat[first_pos][1]
    #                         first = second
    #                         first_pos = i
    #                 else:
    #                     if mid_len < 10:
    #                         if first == second:
    #                             tmp += mid_len + stat[first_pos][1]
    #                             first_pos = i
    #                             continue
    #                         else:
    #                             if first != middle and second != 'N':
    #                                 sum_e = sum_s + stat[first_pos][1] + tmp
    #                                 tmp = 0
    #                                 fw.write('{}\t{}\t{}\t{}\t{}\n'.format(ctg, sum_s, sum_e, first, sum_e - sum_s))
    #                                 first = second
    #                                 first_pos = i
    #                                 sum_s = sum_e
    #                             else:
    #                                 if second != 'N':
    #                                     sum_e = sum_s + stat[first_pos][1] + mid_len
    #                                 else:
    #                                     sum_e = sum_s + stat[first_pos][1] + mid_len + stat[i][1]
    #                                 sum_e += tmp
    #                                 tmp = 0
    #                                 fw.write('{}\t{}\t{}\t{}\t{}\n'.format(ctg, sum_s, sum_e, first, sum_e - sum_s))
    #                                 first = second
    #                                 first_pos = i
    #                                 sum_s = sum_e
    # plt.figure()
    # name = prefix + '.identity.'
    # for res in results:
    #     n, accu = res.get()
    #     ofname = name + n + '.png'
    #     plt.plot(accu)
    #     plt.xlabel('position')
    #     plt.ylabel('identity')
    #     # plt.ylim(0.99, 1)
    #     plt.grid()
    #     plt.savefig(ofname)
    #     plt.clf()

class PafRecord:
    __slots__ = ['qname', 'qlen', 'qstart', 'qend', 'strand', 'tname', 'tlen', 'tstart', 'tend', 'cigar', 'match', 'best']
    def __init__(self, line):
        pattern = re.compile(r'(\d+)([MIDNSHP=X])')
        items = line.strip().split()
        self.qname = items[0]
        self.qlen = int(items[1])
        self.qstart = int(items[2])
        self.qend = int(items[3])
        self.strand = items[4]
        self.tname = items[5]
        self.tlen = int(items[6])
        self.tstart = int(items[7])
        self.tend = int(items[8])
        for i in range(12, len(items)):
            if items[i].startswith('cg:Z:'):
                self.cigar = items[i][5:]
                break
        total = 0
        for match in pattern.finditer(self.cigar):
            l, op = int(match.group(1)), match.group(2)
            if op == '=':
                total += l
        self.match = total
        self.best = False

    def __key(self):
        return (self.qname, self.tname, self.qstart, self.qend, self.tstart, self.tend)
    
    def __hash__(self):
        return hash(self.__key())
    
    def __eq__(self, other):
        if isinstance(other, PafRecord):
            return self.__key() == other.__key()
        return NotImplemented
    
    def __ne__(self, other):
        if isinstance(other, PafRecord):
            return self.__key() != other.__key()
        return NotImplemented
    
    def IsBest(self):
        return self.best

class Node:
    __slots__ = ['start', 'end']
    def __init__(self):
        self.start = []
        self.end = []

class Link:
    __slots__ = ['fro', 'fori', 'to', 'tori', 'cigar']
    def __init__(self, line):
        items = line.strip().split()
        # if len(items) < 6: pass
        self.fro = items[1]
        self.fori = items[2]
        self.to = items[3]
        self.tori = items[4]
        self.cigar = int(items[5][:-1])

def validate_paf(paf, links):
    if paf.qname not in links:
        return False
    pattern = re.compile(r'(\d+)([MIDNSHP=X])')
    max_match = 0
    for match in pattern.finditer(paf.cigar):
        l, op = int(match.group(1)), match.group(2)
        if op == '=':
            max_match = max(max_match, l)
    # if paf.qname == 'utg055031l':
    #     for link in links[paf.qname]:
    #         print('{} {} {} {} {} {} {}'.format(paf.qlen - max_match, paf.qname, link.fro, link.fori, link.to, link.tori, link.cigar), file = sys.stderr)
    for link in links[paf.qname]:
        if link.cigar < paf.qlen - max_match:
            return False
    return True

def search_start(paf, cands):
    if len(cands) == 0:
        return None
    else:
        return cands[0]

def search_paf(pafs, qname):
    for paf in pafs:
        if paf.qname == qname:
            return paf
    return None

from collections import deque
from collections import Counter
def bfs(nodes, source, sink, dir, size, cands, debug = False):
    cand_c = Counter(cands)
    # print(cand_c, file = sys.stderr)
    queue = deque([(source, dir, [source])])
    # if debug:
    #     print(cands, file = sys.stderr)
    while queue:
        (vertex, d, path) = queue.popleft()
        if len(path) > 30:
            return None
        path_c = Counter(path)
        if len(path) > size + 5: 
            return None
        if d == '+':
            neighbors = nodes[vertex].end
        else:
            neighbors = nodes[vertex].start
        nghs = []
        for n in neighbors:
            if path_c[n[0]] < cand_c[n[0]]:
            # if n[0] not in path:
                nghs.append(n)
        # if debug:
        #     print('vertex {} dir {} path {} nghs {} neighbors {}'.format(vertex, d, path, nghs, neighbors), file = sys.stderr)
        for n in nghs:
            if n[0] == sink:
                yield path + [n[0]]
            elif n[0] in cands:
                queue.append((n[0], '+' if n[1] == 0 else '-', path + [n[0]]))
    return None

def search_best_path(source, sink, candidate, nodes):
    paths = []
    tmp = [c.qname for c in candidate]
    if source is not None:
        tmp.append(source.qname)
    if sink is not None:
        tmp.append(sink.qname)
    if len(candidate) == 0:
        return [source.qname if source is not None else None, sink.qname if sink is not None else None]
    if source is None and sink is not None:
        start = []
        # print(tmp, file = sys.stderr)
        for c in candidate:
            if c.tstart == 0:
                start.append(c)
        if len(start) == 0:
            return None
        for s in start:
            # print('source {} sink {} dir {} len {}'.format(s.qname, sink.qname, s.strand, len(candidate)), file = sys.stderr)
            paths.extend(list(bfs(nodes, s.qname, sink.qname, s.strand, len(candidate), tmp, True)))
        mismatch = {}
        if len(paths) == 0:
            return None
        for i in range(len(paths)):
            mismatch[i] = 0
            paf = None
            for p in paths[i][:-1]:
                paf = search_paf(candidate, p)
                if paf is None or (paf is not None and paf.tend < paf.tlen):
                    mismatch[i] += 1000000
                else:
                    mismatch[i] += paf.qlen - paf.match
                if paf is not None:
                    prev = paf
        min_mismatch = min(mismatch.values())
        for i in range(len(paths)):
            if mismatch[i] == min_mismatch:
                return paths[i]
    elif source is not None and sink is None:
        # print(tmp, file = sys.stderr)
        end = []
        for c in candidate:
            if c.tend == c.tlen:
                end.append(c)
        if len(end) == 0:
            return None
        for e in end:
            # print('source {} sink {} dir {} len {}'.format(source.qname, e.qname, source.strand, len(candidate)), file = sys.stderr)
            paths.extend(list(bfs(nodes, source.qname, e.qname, source.strand, len(candidate), tmp)))
        mismatch = {}
        if len(paths) == 0:
            return None
        for i in range(len(paths)):
            mismatch[i] = 0
            prev = None
            for p in paths[i][1:]:
                paf = search_paf(candidate, p)
                if paf is None or (prev is not None and paf.tstart > prev.tend):
                    mismatch[i] += 1000000
                else:
                    mismatch[i] += paf.qlen - paf.match
                if paf is not None:
                    prev = paf
        min_mismatch = min(mismatch.values())
        for i in range(len(paths)):
            if mismatch[i] == min_mismatch:
                return paths[i]
    elif source is not None and sink is not None:
        debug = False
        paths.extend(list(bfs(nodes, source.qname, sink.qname, source.strand, len(candidate), tmp, debug)))
        mismatch = {}
        if len(paths) == 0:
            return None
        for i in range(len(paths)):
            mismatch[i] = 0
            prev = None
            for p in paths[i][1: -1]:
                paf = search_paf(candidate, p)
                if paf is None or (prev is not None and paf.tstart > prev.tend):
                    mismatch[i] += 1000000
                else:
                    mismatch[i] += paf.qlen - paf.match
                if paf is not None:
                    prev = paf
        min_mismatch = min(mismatch.values())
        for i in range(len(paths)):
            if mismatch[i] == min_mismatch:
                return paths[i]
    else:
        # print(tmp, file = sys.stderr)
        start, end = [], []
        for c in candidate:
            if c.tstart == 0:
                start.append(c)
            if c.tend == c.tlen:
                end.append(c)
        if len(start) == 0 or len(end) == 0:
            return None
        for s in start:
            for e in end:
                # print('source {} sink {} dir {} len {}'.format(s.qname, e.qname, s.strand, len(candidate)), file = sys.stderr)
                paths.extend(list(bfs(nodes, s.qname, e.qname, s.strand, len(candidate), tmp)))
        mismatch = {}
        if len(paths) == 0:
            return None
        for i in range(len(paths)):
            mismatch[i] = 0
            prev = None
            for p in paths[i]:
                paf = search_paf(candidate, p)
                if paf is None or (prev is not None and paf.tstart > prev.tend):
                    mismatch[i] += 1000000
                else:
                    mismatch[i] += paf.qlen - paf.match
                if paf is not None:
                    prev = paf
        min_mismatch = min(mismatch.values())
        for i in range(len(paths)):
            if mismatch[i] == min_mismatch:
                return paths[i]
        
def dup_remove(paths):
    path = []
    for p in paths:
        if p is not None and p not in path:
            path.append(p)
    return path

def generate_path_paf(path, pafs):
    ret = []
    index = 0
    for p in pafs:
        if p.qname in path[index]:
            if ret == [] and p.tstart != 0:
                print('first paf {} {} {}'.format(p.qname, p.tstart, p.tend), file = sys.stderr)
            ret.append(p)
            index += 1
            assert index <= len(path)
            if ret[-1].tend == ret[-1].tlen or index == len(path):
                break
    return ret

def dual_to_block(gfname, pfname, jfname, prefix):
    nodes, edges = {}, []
    for line in open(gfname):
        if line.startswith('S'):
            line = line.strip().split('\t')
            nodes[line[1]] = Node()
        elif line.startswith('L'):
            edges.append(line.strip())
    for e in edges:
        items = e.strip().split()
        if items[1] not in nodes:
            print('an edge between {} and {} but {} does not exist in the file'.format(items[1], items[3], items[1]), file = sys.stderr)
            continue
        if items[3] not in nodes:
            print('an edge between {} and {} but {} does not exist in the file'.format(items[1], items[3], items[3]), file = sys.stderr)
            continue
        # if items[2] == '+' and items[4] == '+':
        #     if items[3] not in nodes[items[1]].end:
        #         nodes[items[1]].end.append(items[3])
        #     if items[1] not in nodes[items[3]].start:
        #         nodes[items[3]].start.append(items[1])
        # elif items[2] == '+' and items[4] == '-':
        #     if items[3] not in nodes[items[1]].end:
        #         nodes[items[1]].end.append(items[3])
        #     if items[1] not in nodes[items[3]].end:
        #         nodes[items[3]].end.append(items[1])
        # elif items[2] == '-' and items[4] == '+':
        #     if items[3] not in nodes[items[1]].start:
        #         nodes[items[1]].start.append(items[3])
        #     if items[1] not in nodes[items[3]].start:
        #         nodes[items[3]].start.append(items[1])
        # elif items[2] == '-' and items[4] == '-':
        #     if items[3] not in nodes[items[1]].start:
        #         nodes[items[1]].start.append(items[3])
        #     if items[1] not in nodes[items[3]].end:
        #         nodes[items[3]].end.append(items[1])
        if items[2] == '+' and items[4] == '+':
            if (items[3], 0) not in nodes[items[1]].end:
                nodes[items[1]].end.append((items[3], 0))
            if (items[1], 1) not in nodes[items[3]].start:
                nodes[items[3]].start.append((items[1], 1))
        elif items[2] == '+' and items[4] == '-':
            if (items[3], 1) not in nodes[items[1]].end:
                nodes[items[1]].end.append((items[3], 1))
            if (items[1], 1) not in nodes[items[3]].end:
                nodes[items[3]].end.append((items[1], 1))
        elif items[2] == '-' and items[4] == '+':
            if (items[3], 0) not in nodes[items[1]].start:
                nodes[items[1]].start.append((items[3], 0))
            if (items[1], 0) not in nodes[items[3]].start:
                nodes[items[3]].start.append((items[1], 0))
        elif items[2] == '-' and items[4] == '-':
            if (items[3], 1) not in nodes[items[1]].start:
                nodes[items[1]].start.append((items[3], 1))
            if (items[1], 0) not in nodes[items[3]].end:
                nodes[items[3]].end.append((items[1], 0))
    edges.clear()

    with open(jfname, 'r') as fr:
        bubble = json.load(fr)
    ends, inside = set(), set()
    for k, v in bubble.items():
        for items in v['bubbles']:
            for b in items['ends']:
                ends.add(b)
            for b in items['inside']:
                inside.add(b)
    ends.difference_update(inside)
    
    ref, query = {}, {}
    for line in open(pfname):
        paf = PafRecord(line)
        if paf.qstart != 0 or paf.qend != paf.qlen: continue
        ref.setdefault(paf.tname, []).append(paf)
        query.setdefault(paf.qname, []).append(paf)

    for n, pafs in query.items():
        pafs.sort(key = lambda x: x.match, reverse = True)
        pafs[0].best = True
        for p in pafs[1:]:
            if p.match == p.qlen or p.qname in ends:
                p.best = True
    
    path = {}
    for ctg, pafs in ref.items():
        path.setdefault(ctg, [])
        pafs.sort(key = lambda x: x.tstart)
        anchors = []
        cands, c = [], []
        print(ctg, file = sys.stderr)
        # if ctg != 'tg000301l_1': continue
        for paf in pafs:
            if paf.qname in ends:
                anchors.append(paf)
                cands.append(c.copy())
                c.clear()
            else:
                c.append(paf)
        
        cands.append(c)
        assert len(anchors) + 1 == len(cands)
        if len(anchors) == 0: continue
        index = 0
        # print('anchors {} cands {}'.format(len(anchors), len(cands)), file = sys.stderr)
        # print('start {} {}'.format(anchors[0].qname, cands[0][0].qname), file = sys.stderr)
        if anchors[index].tstart != 0:
            # print('start {} {}'.format(anchors[0].qname, len(cands[0])), file = sys.stderr)
            p = search_best_path(None, anchors[index], cands[index], nodes)
            while p is None and index < len(anchors) - 1:
                c = []
                index += 1
                for k in range(index + 1):
                    c.extend(cands[k])
                # print('index {} len {} {}'.format(index, len(anchors), p is None), file = sys.stderr)
                p = search_best_path(None, anchors[index], c, nodes)
            if p is not None:
                path[ctg].extend(p)

        while index < len(anchors) - 1:
            p = search_best_path(anchors[index], anchors[index + 1], cands[index + 1], nodes)
            # print('index {} len {} {}'.format(index, len(anchors), p is None), file = sys.stderr)
            # print('source {} sink {} cand {}'.format(anchors[index].qname, anchors[index + 1].qname, len(cands[index + 1])), file = sys.stderr)
            i = index + 1
            while p is None and i < len(anchors) - 1:
                c = []
                for j in range(index + 1, i + 2):
                    c.extend(cands[j])
                p = search_best_path(anchors[index], anchors[i + 1], c, nodes)
                i += 1
            if p is not None:
                path[ctg].extend(p)
            index = i - 1
            index += 1
        # print('last {} {}'.format(anchors[index].qname, anchors[-1].qname), file = sys.stderr)
        if anchors[index].tend < anchors[index].tlen:
            c = []
            for i in range(index + 1, len(cands)):
                c.extend(cands[i])
            # print('{} {}'.format(len(c), c[-1].qname), file = sys.stderr)
            p = search_best_path(anchors[index], None, c, nodes)
            if p is not None:
                path[ctg].extend(p)
        if len(path[ctg]) == 0:
            c = []
            c.extend(anchors)
            for i in range(len(cands)):
                c.extend(cands[i])
            p = search_best_path(None, None, c, nodes)
            if p is not None:
                path[ctg].extend(p)

        if len(path[ctg]) == 0: continue
        path[ctg] = dup_remove(path[ctg])
        path[ctg] = generate_path_paf(path[ctg], pafs)
    
    print('path length {}'.format(len(path['tg000001l_2'])), file = sys.stderr)
    # with open(prefix + '.test', 'w') as fw:
    #     for ctg, pafs in ref.items():
    #         for paf in pafs:
    #             fw.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(paf.qname, paf.qlen, paf.qstart, paf.qend, paf.strand, paf.tname, paf.tlen, paf.tstart, paf.tend, paf.cigar, paf.IsBest()))
    with open(prefix + '.path', 'w') as fw:
        for ctg, ps in path.items():
            for p in ps:
                fw.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(p.qname, p.qlen, p.qstart, p.qend, p.strand, p.tname, p.tlen, p.tstart, p.tend, p.cigar))

def homologous_block(fname, prefix):
    pattern = re.compile(r'(\d+)([MIDNSHP=X])')
    match_c = []
    seq_length = {}
    for line in open(fname, 'r'):
        if line.startswith('@SQ'):
            items = line.strip().split()
            seq_length[items[1][3:]] = int(items[2][3:])
        if line.startswith('@'): continue
        items = line.strip().split()
        if items[5] == '*': continue
        for match in pattern.finditer(items[5]):
            l, op = int(match.group(1)), match.group(2)
            if op == '=' and l < seq_length[items[2]]:
                # if l >= 30000:
                #     print('{} {} {}'.format(items[0], l, op), file = sys.stderr)
                match_c.append(l)
    match_c.sort(reverse = True)
    count = [0] * (match_c[0] // 1000 + 1)
    for c in match_c:
        count[c // 1000] += 1
    for i in range(min(500, len(count))):
        print('{} {}'.format(i * 1000, count[i]), file = sys.stderr)
    print('max {} min {} median {}'.format(match_c[0], match_c[-1], match_c[len(match_c) // 2]), file = sys.stderr)

def eval_hic_link(cfname, lfname, pfname, thres):
    block = {}
    for line in open(cfname):
        if line.startswith('Assembly'): continue
        items = line.strip().split()
        if int(items[2]) > int(items[3]):
            block[items[1]] = 1
        elif int(items[2]) < int(items[3]):
            block[items[1]] = -1
        else:
            block[items[1]] = 0

    i = 0
    head = []
    pair = {}
    for line in open(pfname):
        items = line.strip().split()
        pair[items[0]] = items[1]
        pair[items[1]] = items[0]
        if i < 100:
            head.append(items[0])
            head.append(items[1])
            i += 2
        # if items[1].startswith('h1tg000026l'):
        #     if items[0] not in block or block[items[0]] != -1:
        #         head.append(items[0])
        #         head.append(items[1])
        #     else:
        #         head.append(items[1])
        #         head.append(items[0])

    matrix = [] * len(head)
    for i in range(len(head)):
        matrix.append([0] * len(head))
    count = [0, 0, 0]
    index = 0
    fw = open(lfname + '.stat', 'w')
    for line in open(lfname):
        items = line.strip().split()
        try:
            if int(items[4]) < thres or int(items[8]) < thres: continue
        except IndexError:
            print(line.strip())
        if items[1] == items[5]: continue
        if items[1] not in pair or items[5] not in pair: continue
        if items[1] == pair[items[5]]: continue
        if items[1] not in block or items[5] not in block: 
            count[0] += 1
            count[2] += 1
            continue
        if block[items[1]] * block[items[5]] >= 0:
            count[0] += 1
            fw.write('{}\t1\n'.format(line.strip()))
            # print(line.strip())
        else:
            count[1] += 1
            fw.write('{}\t0\n'.format(line.strip()))
            # print(items[0])
            # if index < 100:
                # print('{} {} {} {}'.format(items[1], items[5], block[items[1]], block[items[5]]), file = sys.stderr)
                # print('{} {}'.format(index, line.strip()), file = sys.stderr)
                # print(items[0])
                # index += 1
        count[2] += 1

        if items[1] in head and items[5] in head:
            matrix[head.index(items[1])][head.index(items[5])] += 1
            matrix[head.index(items[5])][head.index(items[1])] += 1
    fw.close()
    print('correct {} incorrect {} total {} rate {}'.format(count[0], count[1], count[2], count[0] / count[2]), file = sys.stderr)
    # for h in head:
    #     print(h, end = ' ', file = sys.stderr)
    # print(file = sys.stderr)
    # for mtx in matrix:
    #     for m in mtx:
    #         print(m, end = ' ', file = sys.stderr)
    #     print(file = sys.stderr)

def hic_link(args):
    eval_hic_link(args.cfname, args.lfname, args.pfname, args.thres)

def homologous(args):
    homologous_block(args.fname, args.prefix)

def identity(args):
    identity_curve(args.fname, args.threads, args.thres, args.block, args.prefix)

def dual(args):
    dual_to_block(args.gfname, args.pfname, args.jfname, args.prefix)

def hic_contact(args):
    hic_contact_eval(args.cfname, args.hfname, args.afname, args.mapq, args.filter)

def contact(args):
    total_contacts(args.mfname)

def phasedn50(args):
    phased_ng50(args.fname, args.column, args.size)

def nx(args):
    NX_stats(args.fname)

def hamming_error(args):
    hamming_error_eval(args.count)

def switch_error(args):
    switch_error_eval(args.bfname, args.snum, args.slimit)

def phase_error(args):
    phase_eval1(args.rfname, args.cfname)

def phase_shasta_error(args):
    phase_shasta_eval(args.cfname, args.pfname)

def main():
    parser = argparse.ArgumentParser('eval')
    subparser = parser.add_subparsers(title = 'subcommand')

    parser_pn = subparser.add_parser('pn50', help = 'calculate phased N50')
    parser_pn.add_argument('-f', '--fname', dest = 'fname', default = None, metavar = 'path', help = 'phased block size')
    parser_pn.add_argument('-c', '--column', dest = 'column', type = int, default = 2, metavar = 'int', help = 'column')
    parser_pn.add_argument('-g', '--genome_size', dest = 'size', type = int, default = 0, metavar = 'int', help = 'genome size')
    parser_pn.set_defaults(func = phasedn50)

    parser_nx = subparser.add_parser('n50', help = 'N50')
    parser_nx.add_argument('-f', '--fname', dest = 'fname', default = None, metavar = 'path', help = 'genome file name')
    parser_nx.set_defaults(func = nx)

    parser_he = subparser.add_parser('hamming', help = 'calculate hamming error rate in merqury hapmer counts')
    parser_he.add_argument('-c', '--count', dest = 'count', default = None, metavar = 'path', help = '*.hapmer.count in merqury result')
    parser_he.set_defaults(func = hamming_error)

    parser_se = subparser.add_parser('switch', help = 'calculate switch error rate')
    parser_se.add_argument('-b', '--bfname', dest = 'bfname', default = None, metavar = 'path', help = 'block file name')
    parser_se.add_argument('-n', '--snum', dest = 'snum', type = int, default = 100, metavar = 'int', help = 'number of markers in a short switch')
    parser_se.add_argument('-l', '--slimit', dest = 'slimit', type = int, default = 20000, metavar = 'int', help = 'limit of distance in a short switch')
    parser_se.set_defaults(func = switch_error)

    parser_pe = subparser.add_parser('phase', help = 'calculate phase error rate')
    parser_pe.add_argument('-r', '--rfname', dest = 'rfname', default = None, metavar = 'path', help = 'phasing result file name')
    parser_pe.add_argument('-c', '--cfname', dest = 'cfname', default = None, metavar = 'path', help = 'hap-mer counts file name')
    parser_pe.set_defaults(func = phase_error)

    parser_gfase = subparser.add_parser('phase_gfase', help = 'calculate phase error rate in GFAse')
    parser_gfase.add_argument('-c', '--cfname', dest = 'cfname', default = None, metavar = 'path', help = 'hap-mer counts file name')
    parser_gfase.add_argument('-p', '--pfname', dest = 'pfname', default = None, metavar = 'path', help = 'phasing result file name')
    parser_gfase.set_defaults(func = phase_shasta_error)

    parser_contact = subparser.add_parser('contact', help = 'calculate total contacts')
    parser_contact.add_argument('-m', '--mfname', dest = 'mfname', default = None, metavar = 'path', help = 'matrix file name')
    parser_contact.set_defaults(func = contact)

    parser_hic = subparser.add_parser('hic', help = 'calculate hic contacts')
    parser_hic.add_argument('-c', '--cfname', dest = 'cfname', default = None, metavar = 'path', help = 'count file name')
    parser_hic.add_argument('-i', '--hfname', dest = 'hfname', default = None, metavar = 'path', help = 'hic file name=')
    parser_hic.add_argument('-a', '--afname', dest = 'afname', default = None, metavar = 'path', help = 'pair file name')
    parser_hic.add_argument('-q', '--mapq', dest = 'mapq', type = int, default = 1, metavar = 'int', help = 'mapq threshold')
    parser_hic.add_argument('-f', '--filter', dest = 'filter', action = 'store_true', help = 'filter map quality')
    parser_hic.set_defaults(func = hic_contact)

    parser_identity = subparser.add_parser('identity', help = 'calculate identity curve')
    parser_identity.add_argument('-f', '--fname', dest = 'fname', default = None, metavar = 'path', help = 'sam file name')
    parser_identity.add_argument('-t', '--threads', dest = 'threads', type = int, default = 1, metavar = 'int', help = 'number of threads')
    parser_identity.add_argument('-s', '--thres', dest = 'thres', type = float, default = 0.99, metavar = 'float', help = 'threshold of identity')
    parser_identity.add_argument('-b', '--block', dest = 'block', type = int, default = 1000, metavar = 'int', help = 'block size')
    parser_identity.add_argument('-p', '--prefix', dest = 'prefix', default = 'identity', metavar = 'str', help = 'prefix of output file')
    parser_identity.set_defaults(func = identity)

    parser_dual = subparser.add_parser('dual', help = 'convert dual to block')
    parser_dual.add_argument('-g', '--gfname', dest = 'gfname', default = None, metavar = 'path', help = 'genome file name')
    parser_dual.add_argument('-m', '--pfname', dest = 'pfname', default = None, metavar = 'path', help = 'paf file name')
    parser_dual.add_argument('-j', '--jfname', dest = 'jfname', default = None, metavar = 'path', help = 'bubble file name')
    parser_dual.add_argument('-p', '--ofname', dest = 'prefix', default = 'block', metavar = 'str', help = 'output file name')
    parser_dual.set_defaults(func = dual)

    parser_homologyzous = subparser.add_parser('homologous', help = 'calculate homologyzous block')
    parser_homologyzous.add_argument('-f', '--fname', dest = 'fname', default = None, metavar = 'path', help = 'sam file name')
    parser_homologyzous.add_argument('-p', '--prefix', dest = 'prefix', default = 'homologyzous', metavar = 'str', help = 'output file name')
    parser_homologyzous.set_defaults(func = homologous)

    parser_hic_link = subparser.add_parser('hic_link', help = 'calculate hic link')
    parser_hic_link.add_argument('-c', '--cfname', dest = 'cfname', default = None, metavar = 'path', help = 'count file name')
    parser_hic_link.add_argument('-l', '--lfname', dest = 'lfname', default = None, metavar = 'path', help = 'link file name')
    parser_hic_link.add_argument('-p', '--pfname', dest = 'pfname', default = None, metavar = 'path', help = 'paf file name')
    parser_hic_link.add_argument('-s', '--thres', dest = 'thres', type = int, default = 100, metavar = 'int', help = 'threshold of contact')
    parser_hic_link.set_defaults(func = hic_link)

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()