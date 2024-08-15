#!/usr/bin/env python3

import os
import re
import sys
import math
import logging
import argparse
import subprocess
import multiprocessing
from collections import OrderedDict
from collections import defaultdict

from util import _fa_or_fq, _load_fa, _load_fq, _load_full_fa, _load_full_fq

logger = logging.getLogger()

def setLogger(logfile = None):
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('[%(asctime)s] [%(levelname)s] ==> %(message)s', '%Y-%m-%d %H:%M:%S')
    if logfile is not None:
        fh = logging.FileHandler(logfile)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

def _confirm(inconsistent, cov):
    if cov < 10: return 0
    return inconsistent >= math.ceil(cov * 0.75)

def _get_variants(varfile):
    # ctg     pos     base    #A  #C  #G  #T  cov isVar
    # ctg1    1       0       2   1   23  1   27  0
    # ctg2    103     2       20  1   24  2   50  1
    logger.info('Parse variants file {} begin'.format(varfile))
    variants = {}
    ctg, last_line = '', ''
    for line in open(varfile):
        if ctg == '' or not line.startswith(ctg):
            last_line = line
            ctg = last_line.split('\t')[0]
        if line.strip().endswith('1'):
            items = line.strip().split('\t')
            variants.setdefault(items[0], {})[items[1]] = [items[2], 0, 0]
    logger.info('Parse variants file {} end'.format(varfile))
    return variants

def _get_variants_and_index(fname, ctg):
    # ctg     pos     base    #A  #C  #G  #T  cov isVar
    # ctg1    1       0       2   1   23  1   27  0
    # ctg2    103     2       20  1   24  2   50  1
    logger.info('Parse variants file {} and index begin'.format(fname))
    ref, index = [], {}
    for line in open(fname):
        if line.strip().startswith(ctg) and line.strip().endswith('1'):
            items = line.strip().split('\t')
            ref.append('ACGT'[int(items[2])])
            index[items[1]] = len(ref) - 1
    logger.info('Parse variants file {} and index end'.format(fname))
    return ref, index

def _get_ref_length(fname):
    length = {}
    ftype = _fa_or_fq(fname)
    if ftype == None:
        logger.error('Unrecognized file type {}'.format(fname))
        exit(1)
    for head, seq in eval('_load_' + ftype)(fname):
        length[head] = len(seq.strip())
    return length

def _get_unzip_name(name):
    pattern = re.compile(r'([0-9]{6}F[_0-9]*)')
    return re.search(pattern, name).group().rstrip('_')

def pecat_to_unzip(pafname, aafname, ppfname, apfname, prefix):
    # pafname: primary assembly fasta file
    # aafname: alternate assembly fasta file
    # ppfname: primary polish fasta file
    # apfname: alternate polish fasta file
    # prefix:  output prefix

    pri_clean_fname = prefix + '.primary.clean.fasta'
    alt_clean_fname = prefix + '.alt.clean.fasta'

    # convert the name of polish sequences to the format of assembly
    if pafname == None or aafname == None:
        pri_unzip_fname = prefix + '.primary.unzip.fasta'
        fsa_unzip1(ppfname, pri_unzip_fname)
        alt_unzip_fname = prefix + '.alt.unzip.fasta'
        fsa_unzip1(apfname, alt_unzip_fname)
    else:
        pri_rename_fname = prefix + '.primary.rename.fasta'
        fsa_rename(pafname, ppfname, pri_rename_fname)
        alt_rename_fname = prefix + '.alt.rename.fasta'
        fsa_rename(aafname, apfname, alt_rename_fname)

        # convert the name of PECAT to unzip
        pri_unzip_fname = prefix + '.primary.unzip.fasta'
        fsa_unzip(pri_rename_fname, pri_unzip_fname)
        alt_unzip_fname = prefix + '.alt.unzip.fasta'
        fsa_unzip(alt_rename_fname, alt_unzip_fname)

    cmd_pri = 'awk \'{ if(NR%2==1) tmp=$1; if(NR%2==0 && !match(tmp,"_")) print tmp"\\n"$0 }\' '
    logger.info('Command for primary: {}'.format(cmd_pri))
    try:
        subprocess.call(cmd_pri, shell = True, 
                        stdin = open(pri_unzip_fname, 'r'), 
                        stdout = open(pri_clean_fname, 'w'), 
                        stderr = sys.stderr)
    except (subprocess.CalledProcessError, OSError) as e:
        logger.error('Error running awk for primary {}'.format(e.cmd), exc_info = 1)
        exit(1)
    cmd_alt = 'awk \'{ if(NR%2==1) tmp=$1; if(NR%2==0 && match(tmp,"_")) print tmp"\\n"$0 }\' '
    logger.info('Command for alt: {}'.format(cmd_alt))
    try:
        subprocess.call(cmd_alt, shell = True, 
                        stdin = open(alt_unzip_fname, 'r'), 
                        stdout = open(alt_clean_fname, 'w'), 
                        stderr = sys.stderr)
    except (subprocess.CalledProcessError, OSError) as e:
        logger.error('Error running awk for alt {}'.format(e.cmd), exc_info = 1)
        exit(1)

def parse_separate(fname):
    tigs = {}
    records = set()
    logger.info('Parse separated file {}'.format(fname))
    pa = ''
    ftype = _fa_or_fq(fname)
    if ftype == None:
        logger.error('Unrecognized file type {}'.format(fname))
        exit(1)
    for head, seq in eval('_load_full_' + ftype)(fname):
        head = _get_unzip_name(head)
        if head in records:
            logger.error('Tig {} appears multiple times'.format(head))
            exit(1)
        records.add(head)
        tigs[head] = seq
        if len(head.split('_')) == 1:
            if pa == 'a':
                logger.error('Separated file contains both p and a tigs. Try -c instead.')
                exit(1)
            pa = 'p'
        elif len(head.split('_')) == 2:
            if pa == 'p':
                logger.error('Separated file contains both p and a tigs. Try -c instead.')
                exit(1)
            pa = 'a'
        else:
            logger.error('Neither p or a contig {}'.format(head))
            exit(1)
    logger.info('Found {} {} tigs'.format(len(tigs), pa))
    return tigs

def parse_combined(fname):
    primary, alt = {}, {}
    records = set()
    ftype = _fa_or_fq(fname)
    if ftype == None:
        logger.error('Unrecognized file type {}'.format(fname))
        exit(1)
    for head, seq in eval('_load_full_' + ftype)(fname):
        head = _get_unzip_name(head)
        if head in records:
            logger.error('Tig {} appears multiple times'.format(head))
            exit(1)
        records.add(head)
        if len(head.split('_')) == 1:
            primary[head] = seq
        elif len(head.split('_')) == 2:
            alt[head] = seq
        else:
            logger.error('Neither p or a contig {}'.format(head))
            exit(1)
    logger.info('Found {} p tigs and {} a tigs'.format(len(primary), len(alt)))
    return primary, alt

def make_name_mapping(primary, alt):
    mapping = {}
    omitted = []
    for key in primary:
        mapping[key] = []
    for a in alt:
        p = a.split('_')[0]
        if p not in mapping.keys():
            logger.info('P tig {} does not exist, omit corresponsed a tig {}'.format(p, a))
            omitted.append(a)
            continue
        mapping[p].append(a)

    with open('name_mapping.txt', 'w') as fw:
        for p in sorted(mapping.keys()):
            for a in sorted(mapping[p]):
                fw.write('{}\t{}\n'.format(p, a))
    logger.info('Omit {} a tigs without corresponding p tigs'.format(len(omitted)))
    return omitted

def fsa_unzip1(fname, prefix):
    re_main = re.compile("ctg([0-9]+)")
    re_sub = re.compile("ctg([0-9]+)-([0-9]+)-([0-9]+)$")
    index = defaultdict(int)
    namemap = {}
    ftype = _fa_or_fq(fname)
    if ftype == None:
        logger.error('Unrecognized file type {}'.format(fname))
        exit(1)
    for head, seq in eval('_load_' + ftype)(fname):
        if '_' in head: continue
        r = re_main.match(head.split()[0])
        assert r != None
        main, homo = r.group(1), ''
        r = re_sub.match(head.split()[0])
        if r:
            homo = main
        # else:
        #     print(head)
        #     homo = head.split()[1].split(':')[0][3:]
        if homo != '':
            name = '{:>06d}F_{:>03d}'.format(int(homo), index[homo] + 1)
            index[homo] += 1
        else:
            name = '{:>06d}F'.format(int(main))
        namemap[head] = name
    
    with open(prefix, 'w') as fw:
        for head, seq in eval('_load_' + ftype)(fname):
            if head not in namemap: continue
            fw.write('>' + namemap[head] + '\n')
            fw.write(seq + '\n')

def fsa_unzip(fname, prefix):
    re_main = re.compile("ctg([0-9]+)")
    re_sub = re.compile("ctg([0-9]+)-([0-9]+)-([0-9]+)$")
    index = defaultdict(int)
    namemap = {}
    ftype = _fa_or_fq(fname)
    if ftype == None:
        logger.error('Unrecognized file type {}'.format(fname))
        exit(1)
    for head, seq in eval('_load_full_' + ftype)(fname):
        r = re_main.match(head.split()[0])
        assert r != None
        main, homo = r.group(1), ''
        r = re_sub.match(head.split()[0])
        if r:
            homo = main
        else:
            homo = head.split()[1].split(':')[0][3:]
        if homo != '':
            name = '{:>06d}F_{:>03d}'.format(int(homo), index[homo] + 1)
            index[homo] += 1
        else:
            name = '{:>06d}F'.format(int(main))
        namemap[head] = name
    
    with open(prefix, 'w') as fw:
        for head, seq in eval('_load_full_' + ftype)(fname):
            fw.write('>' + namemap[head] + '\n')
            fw.write(seq + '\n')

def get_breakpoints(varfile, rifname, ctgfname, prefix):
    # variants[ctg][pos] = [base, inconsistent, coverage]
    variants = _get_variants(varfile)
    # vars_pos[ctg] = [pos1, pos2, ...]
    vars_pos = {}
    for ctg, var in variants.items():
        vars_pos[ctg] = list(var.keys())
        vars_pos[ctg].sort(key = lambda x: int(x))
    
    # rifname
    # ctg rd  #   ctgpos|rdbase|rdpos|rdbase  ctgpos|rdbase|rdpos|rdbase  ...
    for line in open(rifname):
        debug = False
        items = line.strip().split('\t')
        if len(items) < 5: continue
        name = items[0]
        consistent = []
        prev = '0'
        for sites in items[3:]:
            rec = sites.split('|')
            if rec[1] == '-1' or rec[2] == '-1': continue
            ref_pos = rec[0]
            if variants[name][ref_pos][0] == rec[3]:
                consistent.append([ref_pos, True])
            else:
                consistent.append([ref_pos, False])
            prev = ref_pos
        if len(consistent) < 2: continue
        if debug:
            print(items[1], file = sys.stderr)
        # if one SNP is missing in the read, consistent near it is unknown
        index = vars_pos[name].index(consistent[0][0])
        for i in range(1, len(consistent)):
            if debug and items[1] == '20b62918-d921-4302-9626-ca04d83442a1':
                print('{} {}'.format(vars_pos[name][index + 1], consistent[i][0]), file = sys.stderr)
            if index + 1 >= len(vars_pos[name]):
                logger.error('Index out of range {} {} {}'.format(name, index, len(vars_pos[name])))
            if vars_pos[name][index + 1] != consistent[i][0]:
                index = vars_pos[name].index(consistent[i][0]) - 1
            assert(vars_pos[name][index + 1] == consistent[i][0]), print('{} {} {}'.format(vars_pos[name][index + 1], consistent[i][0], items[1]))
            if vars_pos[name][index] == consistent[i - 1][0]:
                variants[name][consistent[i][0]][2] += 1
                if consistent[i][1] != consistent[i - 1][1]:
                    variants[name][consistent[i][0]][1] += 1
            index += 1

    with open(prefix + '.variants.info', 'w') as f:
        for name, var in variants.items():
            for [pos, [_, inconsistent, cov]] in sorted(var.items(), key = lambda x: int(x[0])):
                f.write(name + '\t' + pos + '\t' + str(inconsistent) + '\t' + str(cov) + '\n')

    bps = {}
    for name, var in variants.items():
        prev = '0'
        for [pos, [_, inconsistent, cov]] in sorted(var.items(), key = lambda x: int(x[0])):
            if _confirm(inconsistent, cov):# or (name in interval and (prev, pos) in interval[name]):
                bps.setdefault(name, []).append([int((int(prev) + int(pos)) / 2), inconsistent, cov])
            prev = pos
    
    ctg_length = _get_ref_length(ctgfname)
    for ctg, length in sorted(ctg_length.items(), key = lambda x: x[0]):
        if ctg not in bps:
            print('{}\t{}\t{}'.format(ctg, 0, 1))
            print('{}\t{}\t{}'.format(ctg, length, 1))
        else:
            print('{}\t{}\t{}'.format(ctg, 0, 1))
            for cand in bps[ctg]:
                if cand[2] != 0:
                    print('{}\t{}\t{}'.format(ctg, cand[0], 1))
                else:
                    print('{}\t{}\t{}'.format(ctg, cand[0], 0))
            print('{}\t{}\t{}'.format(ctg, length, 1))

def view_contig(vfname, rifname, region):
    ## region: chr | chr:beg | chr:beg-end
    items = region.split(':')
    ctg_name = items[0]
    ref, index = _get_variants_and_index(vfname, ctg_name)
    # ensure region is legal
    beg, end = 0, 0
    for pos in index.keys():
        end = max(end, int(pos))
    if len(items) == 2:
        r = items[1].split('-')
        if len(r) != 2:
            logger.error('Invalid region format [chr:beg-end]')
            exit(1)
        if int(r[0]) > int(r[1]):
            logger.error('Invalid region: {}-{}'.format(r[0], r[1]))
        beg, end = int(r[0]), int(r[1])

    new_index, new_ref = {}, [] 
    for pos in sorted(index.keys(), key = lambda x: int(x)):
        base = ref.pop(0)
        if int(pos) >= beg and int(pos) < end:
            new_index[pos] = len(new_ref)
            new_ref.append(base)
    table = []
    table.append(new_ref)

    # rifname
    # ctg rd  #   ctgpos|rdbase|rdpos|rdbase  ctgpos|rdbase|rdpos|rdbase  ...
    for line in open(rifname):
        items = line.strip().split('\t')
        if items[0] != ctg_name: continue
        read = ['N'] * len(new_ref)
        for sites in items[3:]:
            rec = sites.split('|')
            if rec[1] == '-1' or rec[2] == '-1': continue
            ref_pos = rec[0]
            if ref_pos not in new_index.keys(): continue
            rd_base = 'ACGT'[int(rec[3])]
            if rd_base == new_ref[new_index[ref_pos]]:
                read[new_index[ref_pos]] = rd_base
            else:
                read[new_index[ref_pos]] = 'acgt'[int(rec[3])]
        if read.count('N') == len(new_ref): continue
        table.append(read)

    for pos in sorted(new_index.keys(), key = lambda x: int(x)):
        print(pos, end = ',')
    print()

    for read in table[0:]:
        for c in read :
            if c == 'N':
                print(' ', end = ',')
            else:
                print(c, end = ',')
        print()

def _nested(a, b):
    ast, ae = int(a[6]), int(a[7])
    bst, be = int(b[6]), int(b[7])
    if ast < be and ast > bst and ae > bst and ae < be: return True
    return False

def _duplicated(a, b):
    ast, ae = int(a[6]), int(a[7])
    bst, be = int(b[6]), int(b[7])
    if ast == bst and ae == be: return True
    return False

def filter_alt(fname, output):
    # filter nested and duplicated contigs
    primary_dict = {}
    for line in open(fname):
        items = line.strip().split('\t')
        primary_dict.setdefault(items[4], []).append(items)
    
    res = []
    for _, items in primary_dict.items():
        items.sort(key = lambda x: int(x[6]))
        filter = [0] * len(items)
        for i in range(len(items)):
            for j in range(i):
                if _nested(items[i], items[j]):
                    filter[i] = 1
                if _nested(items[j], items[i]):
                    filter[j] = 1
                if _duplicated(items[i], items[j]):
                    filter[j] = 1
        for i in range(len(items)):
            if filter[i] == 0:
                res.append(items[i])

    res.sort(key = lambda x: (x[0], int(x[2])))
    with open(output, 'w') as fw:
        for item in res:
            item[2], item[3] = '0', item[1]
            fw.write('\t'.join(i for i in item) + '\n')

def gap(lf, ll):
    itemsf = lf.split()
    itemsl = ll.split()
    if itemsf[4] == '+':
        return abs((int(itemsl[7]) - int(itemsf[8])))
    elif itemsf[4] == '-':
        return abs((int(itemsf[7]) - int(itemsl[8])))
    return 0

def place_alt(pafname, prefix):
    pafs = {}
    for line in open(pafname):
        items = line.strip().split()
        # if len(items) < 9:
        #     continue
        if int(items[9]) < 5000: continue
        pafs.setdefault(items[0], []).append(line.strip())
    pfw = open(prefix + '.filtered.alt2pri.paf', 'w')
    bfw = open(prefix + '.unfiltered.alt2pri.txt', 'w')
    for ctg, records in pafs.items():
        records.sort(key = lambda x: int(x.split()[2]))
        result = [{'start': -1, 'end': -1, 'score': float('-inf')} for _ in range(len(records))]
        for i in range(len(records)):
            itemsi = records[i].strip().split()
            if result[i]['score'] < int(itemsi[9]):
                result[i]['score'] = int(itemsi[9])
                result[i]['start'] = -1
                result[i]['end'] = -1
            
            for j in range(i + 1, len(records)):
                score = result[i]['score']
                itemsj = records[j].strip().split()
                if itemsi[4] != itemsj[4]: continue
                
                score += int(itemsj[9]) - gap(records[i], records[j])
                if result[j]['score'] < score:
                    result[j]['start'] = i
                    result[j]['end'] = j
                    result[j]['score'] = score
                    result[i]['end'] = j
        max = float('-inf')
        best = -1
        for i in range(len(result)):
            if result[i]['score'] > max:
                best = i
                max = result[best]['score']
        index = [best]
        prev = best
        while result[prev]['start'] != -1:
            index.append(result[prev]['start'])
            prev = result[prev]['start']
        index.reverse()
        itemsf = records[index[0]].split()
        itemsl = records[index[-1]].split()
        if itemsf[4] == '+':
            bfw.write(itemsf[0] + '\t' + itemsf[1] + '\t' + itemsf[2] + '\t' + itemsl[3] + '\t')
            bfw.write(itemsf[5] + '\t' + itemsf[6] + '\t' + itemsf[7] + '\t' + itemsl[8] + '\t' + itemsf[4] + '\n')
        elif itemsf[4] == '-':
            bfw.write(itemsf[0] + '\t' + itemsf[1] + '\t' + itemsf[2] + '\t' + itemsl[3] + '\t')
            bfw.write(itemsf[5] + '\t' + itemsf[6] + '\t' + itemsl[7] + '\t' + itemsf[8] + '\t' + itemsf[4] + '\n')
        for i in index:
            pfw.write(records[i] + '\n')

def _reverse(seq):
    table = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    seq_list = list(seq)
    seq_list.reverse()
    for i in range(len(seq_list)):
        seq_list[i] = table[seq_list[i]]
    return ''.join(seq_list)

def cut_contigs(pfname, afname, cfname, prefix):
    primary, alt, c_pri, l_pri = {}, {}, {}, {}
    pair = []
    # interval in primary and alt contigs, pair information
    for line in open(cfname):
        items = line.strip().split('\t')
        alt.setdefault(items[0], []).append([int(items[2]), int(items[3]), items[8]])
        primary.setdefault(items[4], []).append([int(items[6]), int(items[7])])
        l_pri[items[4]] = int(items[5])
        pair.append([items[0], items[2], items[3], items[4], items[6], items[7]])

    for pri, interval in primary.items():
        interval.sort(key = lambda x: x[0])
        prev = 0
        for beg, end in interval:
            if beg > prev:
                c_pri.setdefault(pri, []).append([prev, beg])
            prev = end
        if l_pri[pri] - prev > 0:
            c_pri.setdefault(pri, []).append([prev, l_pri[pri]])

    ftype = _fa_or_fq(afname)
    if ftype == None:
        logger.error('Unrecognized file type: {}'.format(afname))
        exit(1)
    with open(prefix + '.hap1.cut.fasta', 'w') as falt:
        for head, seq in eval('_load_' + ftype)(afname):
            if head not in alt: continue
            interval = alt[head]
            interval.sort(key = lambda x: x[0])
            for [beg, end, strand] in interval:
                if beg >= end or beg < 0 or end < 0:
                    if beg >= end:
                        logger.error('Invalid interval of {}:{}-{}'.format(head, beg, end))
                    if beg < 0:
                        logger.error('Invalid begining of interval of {}: {}'.format(head, beg))
                    if end < 0:
                        logger.error('Invalid ending of interval of {}: {}'.format(head, end))
                    exit(1)
                falt.write('>' + head + ':' + str(beg) + '-' + str(end) + '\n')
                subseq = seq[beg: end]
                if strand == '-':
                    subseq = _reverse(subseq)
                falt.write(''.join(subseq) + '\n')

    ftype = _fa_or_fq(pfname)
    if ftype == None:
        logger.error('Unrecognized file type: {}'.format(pfname))
        exit(1)
    with open(prefix + '.hap2.cut.fasta', 'w') as f, open(prefix + '.collapsed.cut.fasta', 'w') as fc:
        for head, seq in eval('_load_' + ftype)(pfname):
            if head not in primary:
                length = len(seq)
                c_pri.setdefault(head, []).append([0, length])
                fc.write('>' + head + ':0-' + str(length) + '\n')
                fc.write(seq + '\n')
                continue
            interval = primary[head]
            interval.sort(key = lambda x: x[0])
            for [beg, end] in interval:
                if beg >= end or beg < 0 or end < 0:
                    if beg >= end:
                        logger.error('Invalid interval of {}:{}-{}'.format(head, beg, end))
                    if beg < 0:
                        logger.error('Invalid begining of interval of {}: {}'.format(head, beg))
                    if end < 0:
                        logger.error('Invalid ending of interval of {}: {}'.format(head, end))
                    exit(1)
                f.write('>' + head + ':' + str(beg) + '-' + str(end) + '\n')
                f.write(seq[beg: end] + '\n')
            if head not in c_pri: continue
            interval = c_pri[head]
            interval.sort(key = lambda x: x[0])
            for [beg, end] in interval:
                if beg >= end or beg < 0 or end < 0:
                    if beg >= end:
                        logger.error('Invalid interval of {}:{}-{}'.format(head, beg, end))
                    if beg < 0:
                        logger.error('Invalid begining of interval of {}: {}'.format(head, beg))
                    if end < 0:
                        logger.error('Invalid ending of interval of {}: {}'.format(head, end))
                    exit(1)
                fc.write('>' + head + ':' + str(beg) + '-' + str(end) + '\n')
                fc.write(seq[beg: end] + '\n')

    with open(prefix + '.collapsed.cut.bed', 'w') as f:
        for pri, interval in sorted(c_pri.items(), key = lambda x: x[0]):
            interval.sort(key = lambda x: x[0])
            for [beg, end] in interval:
                f.write(pri + '\t' + str(beg) + '\t' + str(end) + '\n')

    with open(prefix + '.hap1.cut.bed', 'w') as f:
        for ctg, interval in sorted(alt.items(), key = lambda x: x[0]):
            interval.sort(key = lambda x: x[0])
            for [beg, end, _] in interval:
                f.write(ctg + '\t' + str(beg) + '\t' + str(end) + '\n')

    with open(prefix + '.hap2.cut.bed', 'w') as f:
        for pri, interval in sorted(primary.items(), key = lambda x: x[0]):
            interval.sort(key = lambda x: x[0])
            for [beg, end] in interval:
                f.write(pri + '\t' + str(beg) + '\t' + str(end) + '\n')

    with open(prefix + '.pair.cut.txt', 'w') as f:
        pair.sort(key = lambda x: (x[3], int(x[4])))
        for [ctg, cb, ce, pri, pb, pe] in pair:
            f.write(pri + ':' + pb + '-' + pe + '\t' + ctg + ':' + cb + '-' + ce + '\n')

def asm_unzip(pfname, afname, ofname, prefix):
    num_p = 0
    pri_rn = {}
    nmap = {}
    ftype = _fa_or_fq(pfname)
    if ftype == None:
        logger.error('Unrecognized file type {}'.format(pfname))
        exit(1)
    with open(prefix + '.primary.rename.fasta', 'w') as fw:
        for head, seq in eval('_load_' + ftype)(pfname):
            nh = '{:>06d}F'.format(num_p)
            fw.write('>' + nh + '\n')
            fw.write(seq + '\n')
            pri_rn[head] = [nh, 1]
            num_p += 1
            nmap[head] = nh

    olp = {}
    alt = set()
    overlaps = []
    for line in open(ofname):
        overlaps.append(line.strip())
        items = line.strip().split('\t')
        if len(items) < 12: continue
        if items[0].strip() in alt: continue
        alt.add(items[0].strip())
        olp.setdefault(items[5].strip(), []).append([items[0].strip(), items[7].strip()])

    alt_rn = {}
    for pri, alts in olp.items():
        alts.sort(key = lambda x: int(x[1]))
        for alt, _ in alts:
            nh = '{}_{:>03d}'.format(pri_rn[pri][0], pri_rn[pri][1])
            pri_rn[pri][1] += 1
            alt_rn[alt] = nh
            nmap[alt] = nh

    num_a = num_p + 10
    reads = {}
    ftype = _fa_or_fq(afname)
    if ftype == None:
        logger.error('Unrecognized file type {}'.format(afname))
        exit(1)
    with open(prefix + '.alt.rename.fasta', 'w') as fw:
        for head, seq in eval('_load_' + ftype)(afname):
            if head in alt_rn:
                reads[alt_rn[head]] = seq
            else:
                logger.info('Unplaced contig {}'.format(head))
        for name, seq in sorted(reads.items(), key = lambda x: x[0]):
            fw.write('>' + name + '\n')
            fw.write(seq + '\n')

    with open(prefix + '.origin2unzip.mapping', 'w') as fw:
        for h, nh in nmap.items():
            fw.write(h + '\t' + nh + '\n')

    nofname = os.path.basename(ofname)
    overlaps.sort(key = lambda x: (x.split('\t')[5], int(x.split('\t')[7])))
    with open(prefix + '.alt2pri.rename.paf', 'w') as fw:
        for line in overlaps:
            items = line.strip().split('\t')
            items[0] = nmap[items[0]]
            items[5] = nmap[items[5]]
            fw.write('\t'.join(items) + '\n')

def phase_contigs(rfname, pbfname, cbfname, mfname, sfname, prefix):
    # sfname: ctg   pos F|S
    switch_error = {}
    if sfname != None:
        for line in open(sfname):
            items = line.strip().split()
            if items[2] != 'S' and items[2] != 'F' and items[2] != 'C': continue
            switch_error[items[0]] = [items[1], items[2]]
    pair0, pair1 = {}, {}
    for line in open(rfname):
        items = line.strip().split()
        if items[1] in switch_error and items[2] in switch_error and switch_error[items[1]][1] == 'C' and switch_error[items[2]][1] == 'C':
            pair0[items[2]] = items[1]
            pair1[items[1]] = items[2]
        else:
            pair0[items[1]] = items[2]
            pair1[items[2]] = items[1]
    
    bed = []
    for line in open(pbfname):
        items = line.strip().split('\t')
        bed.append(items)
    for line in open(cbfname):
        items = line.strip().split('\t')
        bed.append(items)
    bed.sort(key = lambda x: (x[0], int(x[1])))

    bed0, bed1 = OrderedDict(), OrderedDict()
    for [ctg, start, end] in bed:
        c = ctg + ':' + start + '-' + end
        if c in pair0.keys():
            bed0.setdefault(ctg, []).append(c)
            bed1.setdefault(ctg, []).append(pair0[c])
        elif c in pair1.keys():
            bed0.setdefault(ctg, []).append(pair1[c])
            bed1.setdefault(ctg, []).append(c)
        else:
            bed0.setdefault(ctg, []).append(c)
            bed1.setdefault(ctg, []).append(c)

    with open(prefix + '.phased0.bed', 'w') as fw:
        for c, cbs in bed0.items():
            for cb in cbs:
                if cb in switch_error:
                    items1 = cb.split(':')
                    region1 = items1[1].split('-')
                    if cb in pair0:
                        p = pair0[cb]
                    else:
                        p = pair1[cb]
                    items2 = p.split(':')
                    region2 = items2[1].split('-')
                    if switch_error[cb][1] == 'F':
                        fw.write('{}\t{}\t{}\t{}_0\n'.format(items1[0], region1[0], int(region1[0]) + int(switch_error[cb][0]), c))
                        fw.write('{}\t{}\t{}\t{}_0\n'.format(items2[0], int(region2[0]) + int(switch_error[p][0]), region2[1], c))
                    elif switch_error[cb][1] == 'S':
                        fw.write('{}\t{}\t{}\t{}_0\n'.format(items2[0], region2[0], int(region2[0]) + int(switch_error[p][0]), c))
                        fw.write('{}\t{}\t{}\t{}_0\n'.format(items1[0], int(region1[0]) + int(switch_error[cb][0]), region1[1], c))
                    elif switch_error[cb][1] == 'C':
                        cb = cb.replace(':', '\t')
                        cb = cb.replace('-', '\t')
                        fw.write(cb + '\t' + c + '_0\n')
                else:
                    cb = cb.replace(':', '\t')
                    cb = cb.replace('-', '\t')
                    fw.write(cb + '\t' + c + '_0\n')
            # for cb in cbs:
            #     cb = cb.replace(':', '\t')
            #     cb = cb.replace('-', '\t')
            #     fw.write(cb + '\t' + c + '\n')
    with open(prefix + '.phased1.bed', 'w') as fw:
        for c, cbs in bed1.items():
            for cb in cbs:
                if cb in switch_error:
                    items1 = cb.split(':')
                    region1 = items1[1].split('-')
                    if cb in pair0:
                        p = pair0[cb]
                    else:
                        p = pair1[cb]
                    items2 = p.split(':')
                    region2 = items2[1].split('-')
                    if switch_error[cb][1] == 'F':
                        fw.write('{}\t{}\t{}\t{}_1\n'.format(items1[0], region1[0], int(region1[0]) + int(switch_error[cb][0]), c))
                        fw.write('{}\t{}\t{}\t{}_1\n'.format(items2[0], int(region2[0]) + int(switch_error[p][0]), region2[1], c))
                    elif switch_error[cb][1] == 'S':
                        fw.write('{}\t{}\t{}\t{}_1\n'.format(items2[0], region2[0], int(region2[0]) + int(switch_error[p][0]), c))
                        fw.write('{}\t{}\t{}\t{}_1\n'.format(items1[0], int(region1[0]) + int(switch_error[cb][0]), region1[1], c))
                    elif switch_error[cb][1] == 'C':
                        cb = cb.replace(':', '\t')
                        cb = cb.replace('-', '\t')
                        fw.write(cb + '\t' + c + '_1\n')
                else:
                    cb = cb.replace(':', '\t')
                    cb = cb.replace('-', '\t')
                    fw.write(cb + '\t' + c + '_1\n')
            # for cb in cbs:
            #     cb = cb.replace(':', '\t')
            #     cb = cb.replace('-', '\t')
            #     fw.write(cb + '\t' + c + '\n')

    ctgs = {}
    ftype = _fa_or_fq(mfname)
    if ftype == None:
        logger.error('Unrecognized file type {}'.format(mfname))
        exit(1)
    for head, seq in eval('_load_' + ftype)(mfname):
        ctgs[head] = seq
    
    with open(prefix + '.phased0.fasta', 'w') as fw:
        for c, cbs in bed0.items():
            fw.write('>' + c + '_0\n')
            for cb in cbs:
                if cb in switch_error:
                    if cb in pair0:
                        p = pair0[cb]
                    else:
                        p = pair1[cb]
                    items1 = cb.split(':')
                    region1 = items1[1].split('-')
                    items2 = p.split(':')
                    region2 = items2[1].split('-')
                    if switch_error[cb][1] == 'F':
                        st1 = 0
                        en1 = int(switch_error[cb][0])
                        st2 = int(switch_error[p][0])
                        en2 = int(region2[1]) - int(region2[0])
                        fw.write(ctgs[cb][st1: en1])
                        fw.write(ctgs[p][st2: en2])
                    elif switch_error[cb][1] == 'S':
                        st1 = 0
                        en1 = int(switch_error[p][0])
                        st2 = int(switch_error[cb][0])
                        en2 = int(region1[1]) - int(region1[0])
                        fw.write(ctgs[p][st1: en1])
                        fw.write(ctgs[cb][st2: en2])
                    elif switch_error[cb][1] == 'C':
                        fw.write(ctgs[cb])
                else:
                    fw.write(ctgs[cb])
                # fw.write(ctgs[cb])
            fw.write('\n')

    with open(prefix + '.phased1.fasta', 'w') as fw:
        for c, cbs in bed1.items():
            fw.write('>' + c + '_1\n')
            for cb in cbs:
                if cb in switch_error:
                    if cb in pair0:
                        p = pair0[cb]
                    else:
                        p = pair1[cb]
                    items1 = cb.split(':')
                    region1 = items1[1].split('-')
                    items2 = p.split(':')
                    region2 = items2[1].split('-')
                    if switch_error[cb][1] == 'F':
                        st1 = 0
                        en1 = int(switch_error[cb][0])
                        st2 = int(switch_error[p][0])
                        en2 = int(region2[1]) - int(region2[0])
                        fw.write(ctgs[cb][st1: en1])
                        fw.write(ctgs[p][st2: en2])
                    elif switch_error[cb][1] == 'S':
                        st1 = 0
                        en1 = int(switch_error[p][0])
                        st2 = int(switch_error[cb][0])
                        en2 = int(region1[1]) - int(region1[0])
                        fw.write(ctgs[p][st1: en1])
                        fw.write(ctgs[cb][st2: en2])
                    elif switch_error[cb][1] == 'C':
                        fw.write(ctgs[cb])
                else:
                    fw.write(ctgs[cb])
                # fw.write(ctgs[cb])
            fw.write('\n')

def phase_contigs_dual(rfname, h1bfname, h2bfname, c1bfname, c2bfname, mfname, sfname, prefix):
    # sfname: ctg   pos F|S
    switch_error = {}
    if sfname != None:
        for line in open(sfname):
            items = line.strip().split()
            if items[2] != 'S' and items[2] != 'F' and items[2] != 'C': continue
            switch_error[items[0]] = [items[1], items[2]]
    pair0, pair1 = {}, {}
    for line in open(rfname):
        items = line.strip().split()
        if items[1] in switch_error and items[2] in switch_error and switch_error[items[1]][1] == 'C' and switch_error[items[2]][1] == 'C':
            pair0[items[2]] = items[1]
            pair1[items[1]] = items[2]
        else:
            pair0[items[1]] = items[2]
            pair1[items[2]] = items[1]
    
    bed = []
    for line in open(h1bfname):
        items = line.strip().split('\t')
        bed.append(items)
    for line in open(c1bfname):
        items = line.strip().split('\t')
        bed.append(items)
    bed.sort(key = lambda x: (x[0], int(x[1])))

    bed0, bed1 = OrderedDict(), OrderedDict()
    for [ctg, start, end] in bed:
        c = ctg + ':' + start + '-' + end
        if c in pair0.keys():
            bed0.setdefault(ctg, []).append(c)
        elif c in pair1.keys():
            bed0.setdefault(ctg, []).append(pair1[c])
        else:
            bed0.setdefault(ctg, []).append(c)
    bed = []
    for line in open(h2bfname):
        items = line.strip().split('\t')
        bed.append(items)
    for line in open(c2bfname):
        items = line.strip().split('\t')
        bed.append(items)
    bed.sort(key = lambda x: (x[0], int(x[1])))

    for [ctg, start, end] in bed:
        c = ctg + ':' + start + '-' + end
        if c in pair0.keys():
            bed1.setdefault(ctg, []).append(pair0[c])
        elif c in pair1.keys():
            bed1.setdefault(ctg, []).append(c)
        else:
            bed1.setdefault(ctg, []).append(c)

    with open(prefix + '.phased0.bed', 'w') as fw:
        for c, cbs in bed0.items():
            for cb in cbs:
                if cb in switch_error:
                    items1 = cb.split(':')
                    region1 = items1[1].split('-')
                    if cb in pair0:
                        p = pair0[cb]
                    else:
                        p = pair1[cb]
                    items2 = p.split(':')
                    region2 = items2[1].split('-')
                    if switch_error[cb][1] == 'F':
                        fw.write('{}\t{}\t{}\t{}_0\n'.format(items1[0], region1[0], int(region1[0]) + int(switch_error[cb][0]), c))
                        fw.write('{}\t{}\t{}\t{}_0\n'.format(items2[0], int(region2[0]) + int(switch_error[p][0]), region2[1], c))
                    elif switch_error[cb][1] == 'S':
                        fw.write('{}\t{}\t{}\t{}_0\n'.format(items2[0], region2[0], int(region2[0]) + int(switch_error[p][0]), c))
                        fw.write('{}\t{}\t{}\t{}_0\n'.format(items1[0], int(region1[0]) + int(switch_error[cb][0]), region1[1], c))
                    elif switch_error[cb][1] == 'C':
                        cb = cb.replace(':', '\t')
                        cb = cb.replace('-', '\t')
                        fw.write(cb + '\t' + c + '_0\n')
                else:
                    cb = cb.replace(':', '\t')
                    cb = cb.replace('-', '\t')
                    fw.write(cb + '\t' + c + '_0\n')
    with open(prefix + '.phased1.bed', 'w') as fw:
        for c, cbs in bed1.items():
            for cb in cbs:
                if cb in switch_error:
                    items1 = cb.split(':')
                    region1 = items1[1].split('-')
                    if cb in pair0:
                        p = pair0[cb]
                    else:
                        p = pair1[cb]
                    items2 = p.split(':')
                    region2 = items2[1].split('-')
                    if switch_error[cb][1] == 'F':
                        fw.write('{}\t{}\t{}\t{}_1\n'.format(items1[0], region1[0], int(region1[0]) + int(switch_error[cb][0]), c))
                        fw.write('{}\t{}\t{}\t{}_1\n'.format(items2[0], int(region2[0]) + int(switch_error[p][0]), region2[1], c))
                    elif switch_error[cb][1] == 'S':
                        fw.write('{}\t{}\t{}\t{}_1\n'.format(items2[0], region2[0], int(region2[0]) + int(switch_error[p][0]), c))
                        fw.write('{}\t{}\t{}\t{}_1\n'.format(items1[0], int(region1[0]) + int(switch_error[cb][0]), region1[1], c))
                    elif switch_error[cb][1] == 'C':
                        cb = cb.replace(':', '\t')
                        cb = cb.replace('-', '\t')
                        fw.write(cb + '\t' + c + '_1\n')
                else:
                    cb = cb.replace(':', '\t')
                    cb = cb.replace('-', '\t')
                    fw.write(cb + '\t' + c + '_1\n')

    ctgs = {}
    ftype = _fa_or_fq(mfname)
    if ftype == None:
        logger.error('Unrecognized file type {}'.format(mfname))
        exit(1)
    for head, seq in eval('_load_' + ftype)(mfname):
        ctgs[head] = seq
    
    with open(prefix + '.phased0.fasta', 'w') as fw:
        for c, cbs in bed0.items():
            fw.write('>' + c + '_0\n')
            for cb in cbs:
                if cb in switch_error:
                    if cb in pair0:
                        p = pair0[cb]
                    else:
                        p = pair1[cb]
                    items1 = cb.split(':')
                    region1 = items1[1].split('-')
                    items2 = p.split(':')
                    region2 = items2[1].split('-')
                    if switch_error[cb][1] == 'F':
                        st1 = int(region1[0])
                        en1 = int(region1[0]) + int(switch_error[cb][0])
                        st2 = int(region2[0]) + int(switch_error[p][0])
                        en2 = int(region2[1])
                        fw.write(ctgs[cb][st1: en1])
                        fw.write(ctgs[p][st2: en2])
                    elif switch_error[cb][1] == 'S':
                        st1 = int(region2[0])
                        en1 = int(region2[0]) + int(switch_error[p][0])
                        st2 = int(region1[0]) + int(switch_error[cb][0])
                        en2 = int(region1[1])
                        fw.write(ctgs[p][st1: en1])
                        fw.write(ctgs[cb][st2: en2])
                    elif switch_error[cb][1] == 'C':
                        fw.write(ctgs[cb])
                else:
                    fw.write(ctgs[cb])
            fw.write('\n')

    with open(prefix + '.phased1.fasta', 'w') as fw:
        for c, cbs in bed1.items():
            fw.write('>' + c + '_1\n')
            for cb in cbs:
                if cb in switch_error:
                    if cb in pair0:
                        p = pair0[cb]
                    else:
                        p = pair1[cb]
                    items1 = cb.split(':')
                    region1 = items1[1].split('-')
                    items2 = p.split(':')
                    region2 = items2[1].split('-')
                    if switch_error[cb][1] == 'F':
                        st1 = int(region1[0])
                        en1 = int(region1[0]) + int(switch_error[cb][0])
                        st2 = int(region2[0]) + int(switch_error[p][0])
                        en2 = int(region2[1])
                        fw.write(ctgs[cb][st1: en1])
                        fw.write(ctgs[p][st2: en2])
                    elif switch_error[cb][1] == 'S':
                        st1 = int(region2[0])
                        en1 = int(region2[0]) + int(switch_error[p][0])
                        st2 = int(region1[0]) + int(switch_error[cb][0])
                        en2 = int(region1[1])
                        fw.write(ctgs[p][st1: en1])
                        fw.write(ctgs[cb][st2: en2])
                    elif switch_error[cb][1] == 'C':
                        fw.write(ctgs[cb])
                else:
                    fw.write(ctgs[cb])
            fw.write('\n')

def map_separate(pfname, afname, fmt):
    errfname = pfname + '.err'
    cmd = ['minimap2', '-x', 'asm20', '--secondary=no', '--eqx', '--end-bonus', '100', '-U', '0,2000', '-f', '0.00001']
    if fmt.lower() == 'sam':
        cmd.extend(['-a'])
    elif fmt.lower() == 'paf':
        cmd.extend(['-c'])
    cmd.extend(["'" + pfname + "'", "'" + afname + "'"])
    try:
        result = subprocess.check_output(['/bin/bash', '-c', ' '.join(cmd)], 
                                        stderr = open(errfname, 'w'), encoding = 'utf-8')[:-1]
        result = result.split('\n')
    except subprocess.CalledProcessError as e:
        result = []
        logger.error('Error when mapping {} and {}'.format(os.path.basename(pfname), os.path.basename(afname)))
    finally:
        os.remove(pfname)
        os.remove(afname)
        os.remove(errfname)
        return result

def map_unzip(pfname, afname, mfname, ofname, threads, dir, fmt):
    primary, alt = {}, {}
    # load primary contigs
    ftype = _fa_or_fq(pfname)
    if ftype == None:
        logger.error('Unrecognized file type {}'.format(pfname))
        exit(1)
    for head, seq in eval('_load_' + ftype)(pfname):
        primary[head] = seq
    # laod alt contigs
    ftype = _fa_or_fq(afname)
    if ftype == None:
        logger.error('Unrecognized file type {}'.format(afname))
        exit(1)
    for head, seq in eval('_load_' + ftype)(afname):
        alt[head] = seq
    # load name mapping
    name_mapping = {}
    for line in open(mfname):
        items = line.strip().split('\t')
        name_mapping.setdefault(items[0], []).append(items[1])
    # temporary directory to store temporary files
    os.makedirs(dir, exist_ok = True)
    # each file ends with '.p.fasta' stores one primary contig
    # each file ends with '.a.fasta' stores alt contigs mapped to one contig
    for ref in name_mapping.keys():
        file = os.path.join(dir, ref + '.p.fasta')
        with open(file, 'w') as fw:
            fw.write('>' + ref + '\n')
            fw.write(primary[ref] + '\n')
        file = os.path.join(dir, ref + '.a.fasta')
        with open(file, 'w') as fw:
            for a in name_mapping[ref]:
                fw.write('>' + a + '\n')
                fw.write(alt[a] + '\n')

    # each threads process one contig
    pool = multiprocessing.Pool(threads)
    results = []
    for ref in name_mapping.keys():
        results.append(pool.apply_async(map_separate, (os.path.join(dir, ref + '.p.fasta'), os.path.join(dir, ref + '.a.fasta'), fmt, )))
    pool.close()
    pool.join()

    # merge all results to output file
    if fmt == 'paf':
        ofname += '.paf'
        with open(ofname, 'w') as fw:
            for rs in results:
                for r in rs.get():
                    if len(r) == 0: continue
                    fw.write(r + '\n')
    elif fmt == 'sam':
        ofname += '.sam'
        sq = []
        pg = []
        record = []
        for rs in results:
            for r in rs.get():
                if r.startswith('@SQ'):
                    sq.append(r)
                elif r.startswith('@PG'):
                    pg.append(r)
                else:
                    record.append(r)
        with open(ofname, 'w') as fw:
            for s in sq:
                fw.write(s + '\n')
            for p in pg:
                fw.write(p + '\n')
                break
            for r in record:
                fw.write(r + '\n')
    # remove temporary directory
    os.removedirs(dir)

def fsa_rename(ofname, nfname, output):
    namemap = {}
    ftype = _fa_or_fq(ofname)
    if ftype == None:
        logger.error('Unrecognized file type {}'.format(ofname))
        exit(1)
    for head, seq in eval('_load_full_' + ftype)(ofname):
        namemap[head.split()[0]] = head
    reads = {}
    ftype = _fa_or_fq(nfname)
    if ftype == None:
        logger.error('Unrecognized file type {}'.format(nfname))
        exit(1)
    for head, seq in eval('_load_' + ftype)(nfname):
        reads[head] = seq
    with open(output, 'w') as fw:
        for head, seq in reads.items():
            if head in namemap:
                fw.write('>' + namemap[head] + '\n')
                fw.write(seq + '\n')

def create_scaffolding_lachesis(fname, dir):
    unknown_gap_size = 1000
    clusters_file = os.path.join(dir, 'main_results/clusters.by_name.txt')
    if not os.path.exists(clusters_file):
        logger.error('count not find clusters file \'clusters.by_name.txt\' in {}/main_results/'.format(clusters_file))
        exit(FileNotFoundError)
    N_ordering = 0
    while True:
        f = os.path.join(dir, 'main_results/group{}.ordering'.format(N_ordering))
        if not os.path.exists(f): break
        N_ordering += 1
    if N_ordering == 0:
        logger.error('could not find any ordering files \'group*.ordering\' in {}/main_results/'.format(dir))
    logger.info('find {} ordering files \'group*.ordering\' in {}/main_results/'.format(N_ordering, dir))
    ofname = os.path.join(dir, 'lachesis_scaffold.fasta')
    fw = open(ofname, 'w')
    contigs = {}
    ftype = _fa_or_fq(fname)
    if ftype == None:
        logger.error('Unrecognized file type {}'.format(fname))
        exit(1)
    for head, seq in eval('_load_' + ftype)(fname):
        contigs[head] = seq
    logger.info('find {} contigs/scaffolds in assembly'.format(len(contigs)))
    contigs_used = {}
    for head in contigs.keys():
        contigs_used[head] = 0
    for i in range(N_ordering):
        ordering_file = os.path.join(dir, 'main_results/group{}.ordering'.format(i))
        logger.info('create a scaffold from file {}...'.format(ordering_file))
        if not os.path.exists(ordering_file):
            logger.error('can not find file {}'.format(ordering_file))
            exit(FileNotFoundError)
        scaffold, N_contigs = '', 0
        for line in open(ordering_file, 'r'):
            if line.startswith('#'): continue
            items = line.strip().split()
            ID, name, rc, q = items[0], items[1], items[2], items[3]
            gap_size = 0
            if items[4] == '.':
                gap_size = unknown_gap_size
            else:
                gap_size = int(items[4])
            gap = 'N' * gap_size
            if name not in contigs:
                logger.error('ordering file {} includes contig named {}, not found in fasta file {}'.format(ordering_file, name, fname))
                exit(1)
            contigs_used[name] += 1
            if rc == '1':
                scaffold += _reverse(contigs[name])
            else:
                scaffold += contigs[name]
            scaffold += gap
            N_contigs += 1
        index = len(scaffold) - 1
        while index >= 0:
            if scaffold[index] != 'N': break
            index -= 1
        # assert index > 0
        scaffold = scaffold[0: index + 1]
        length = len(scaffold)
        logger.info('length {}'.format(length))
        if N_contigs == 0: continue
        scaffold_name = 'Lachesis_group{}__{}_contigs__length_{}'.format(i, N_contigs, length)
        fw.write('>{}\n'.format(scaffold_name))
        for j in range(0, length, 80):
            if j + 80 < length:
                fw.write(scaffold[j: j + 80] + '\n')
            else:
                fw.write(scaffold[j: length] + '\n')

    ID, N_unordered_contigs = 0, 0
    for line in open(clusters_file, 'r'):
        if line.startswith('#'): continue
        items = line.strip().split()
        for item in items:
            if contigs_used[item] > 0: continue
            contigs_used[item] += 1
            if item not in contigs:
                logger.error('cluster file {} includes contig named {}, not found in fasta file {}'.format(clusters_file, item, fname))
            fw.write('>{}__unordered_in_group{}\n'.format(item, ID))
            for j in range(0, len(contigs[item]), 80):
                if j + 80 < len(contigs[item]):
                    fw.write(contigs[item][j: j + 80] + '\n')
                else:
                    fw.write(contigs[item][j: len(contigs[item])] + '\n')
            N_unordered_contigs += 1
        ID += 1
    logger.info('write {} contigs that are clustered but not ordered by lachesis'.format(N_unordered_contigs))
    N_unordered_contigs = 0
    for name, seq in contigs.items():
        if contigs_used[name] > 0: continue
        fw.write('>{}__unclustered\n'.format(name))
        for j in range(0, len(seq), 80):
            if j + 80 < len(seq):
                fw.write(seq[j: j + 80] + '\n')
            else:
                fw.write(seq[j: len(seq)] + '\n')
    fw.close()

def create_scaffolding(fname, dir, rfname, prefix):
    unknown_gap_size = 1000
    clusters_file = os.path.join(dir, 'main_results/clusters.by_name.txt')
    if not os.path.exists(clusters_file):
        logger.error('count not find clusters file \'clusters.by_name.txt\' in {}/main_results/'.format(clusters_file))
        exit(FileNotFoundError)
    N_ordering = 0
    while True:
        f = os.path.join(dir, 'main_results/group{}.ordering'.format(N_ordering))
        if not os.path.exists(f): break
        N_ordering += 1
    if N_ordering == 0:
        logger.error('could not find any ordering files \'group*.ordering\' in {}/main_results/'.format(dir))
    logger.info('find {} ordering files \'group*.ordering\' in {}/main_results/'.format(N_ordering, dir))
    if not os.path.exists(rfname):
        logger.error('can not find phasing result file {}'.format(rfname))
    phasing, pair = {}, {}
    for line in open(rfname, 'r'):
        items = line.strip().split()
        phasing[items[1]] = 0
        phasing[items[2]] = 1
        pair[items[1]] = items[2]
        pair[items[2]] = items[1]
    bed, names, names1 = [], [], []
    contigs_used = set()
    for i in range(N_ordering):
        ordering_file = os.path.join(dir, 'main_results/group{}.ordering'.format(i))
        if not os.path.exists(ordering_file):
            logger.error('can not find file {}'.format(ordering_file))
            exit(FileNotFoundError)
        N_contigs = 0
        group = []
        for line in open(ordering_file, 'r'):
            if line.startswith('#'): continue
            items = line.strip().split()
            name, rc = items[1], items[2]
            gap_size = 0
            if items[4] == '.':
                gap_size = unknown_gap_size
            else:
                gap_size = int(items[4])
            contigs_used.add(name)
            if phasing[name] == 0:
                group.append([name, rc, gap_size])
            else:
                group.append([pair[name], rc, gap_size])
            N_contigs += 1
        if N_contigs == 0: continue
        group[-1][2] = 0
        bed.append(group)
        scaffold_name = 'Lachesis_group{}__{}_contigs_'.format(i, N_contigs)
        names.append(scaffold_name)
        names1.append(scaffold_name)

    ID = 0
    for line in open(clusters_file, 'r'):
        if line.startswith('#'): continue
        items = line.strip().split()
        for item in items:
            if item in contigs_used: continue
            contigs_used.add(item)
            group = []
            if phasing[item] == 0:
                group.append([item, '0', 0])
            else:
                group.append([pair[item], '0', 0])
            bed.append(group)
            name = '{}__unordered_in_group{}_'.format(item, ID)
            name1 = '{}__unordered_in_group{}_'.format(pair[item], ID)
            names.append(name)
            names1.append(name1)
        ID += 1
    
    contigs = {}
    ftype = _fa_or_fq(fname)
    if ftype == None:
        logger.error('Unrecognized file type {}'.format(fname))
        exit(1)
    for head, seq in eval('_load_' + ftype)(fname):
        contigs[head] = seq
    logger.info('find {} contigs/scaffolds in assembly'.format(len(contigs)))
    
    # for name, seq in contigs.items():
    #     if name in contigs_used or phasing[name] == 1: continue
    #     contigs_used.add(name)
    #     group = [[name, '0', 0]]
    #     bed.append(group)
    #     name = '{}__unclustered_'.format(name)
    #     names.append(name)
    
    assert len(bed) == len(names), print('{}\t{}'.format(len(bed), len(names)), file = sys.stderr)
    with open(prefix + '.scaffold0.bed', 'w') as fw:
        index = 0
        for items in bed:
            for name, rc, size in items:
                fw.write('{}\t{}\t{}\t{}0\n'.format(name, 0, len(contigs[name]), names[index]))
            index += 1
    with open(prefix + '.scaffold0.fasta', 'w') as fw:
        index = 0
        for items in bed:
            fw.write('>{}0\n'.format(names[index]))
            for name, rc, size in items:
                seq = contigs[name]
                if rc == '1':
                    seq = _reverse(seq)
                fw.write(seq)
                fw.write('N' * size)
            fw.write('\n')
            index += 1
    with open(prefix + '.scaffold1.bed', 'w') as fw:
        index = 0
        for items in bed:
            for name, rc, size in items:
                fw.write('{}\t{}\t{}\t{}1\n'.format(pair[name], 0, len(contigs[pair[name]]), names1[index]))
            index += 1
    with open(prefix + '.scaffold1.fasta', 'w') as fw:
        index = 0
        for items in bed:
            fw.write('>{}1\n'.format(names1[index]))
            for name, rc, size in items:
                seq = contigs[pair[name]]
                if rc == '1':
                    seq = _reverse(seq)
                fw.write(seq)
                fw.write('N' * size)
            fw.write('\n')
            index += 1

def merge_var(vfname, bfname1, bfname2):
    vars = OrderedDict()
    for line in open(vfname, 'r'):
        items = line.strip().split()
        vars.setdefault(items[0], []).append(int(items[1]))

    vars_new = OrderedDict()
    c = 0
    for bfname in [bfname1, bfname2]:
        prev_ctg, prev_end = '', 0
        index = 0
        for line in open(bfname, 'r'):
            items = line.strip().split()
            st, en = int(items[1]), int(items[2])
            if prev_ctg == '' or items[3] != prev_ctg:
                prev_ctg = items[3]
                prev_end = 0
                index = 0
            if items[0] not in vars: continue
            while index < len(vars[items[0]]):
                while index < len(vars[items[0]]) and vars[items[0]][index] < st:
                    index += 1
                if index >= len(vars[items[0]]) or vars[items[0]][index] >= en: break
                ctg_new = items[3] # + '_' + str(c)
                vars_new.setdefault(ctg_new, []).append(prev_end + vars[items[0]][index] - st)
                index += 1
            prev_end += en
        c += 1

    with open(vfname + '.chrom', 'w') as fw:
        for ctg, position in vars_new.items():
            for pos in position:
                fw.write('{}\t{}\n'.format(ctg, pos))

def shasta2unzip_and_mapping(ifname, prefix, ofname, min_length, threads, dir):
    contigs = {}
    ftype = _fa_or_fq(ifname)
    if ftype == None:
        logger.error('Unrecognized file type {}'.format(ifname))
        exit(1)
    for head, seq in eval('_load_' + ftype)(ifname):
        contigs[head] = seq
    logger.info('load {} segments in {}'.format(len(contigs), ifname))
    chain = defaultdict(list)
    optional = {}
    
    for head, seq in contigs.items():
        if head.startswith('PR'):
            items = head.split('.')
            bc, pos, _, h = int(items[1]), int(items[2]), items[3], int(items[4])
            chain[bc].append([pos, h, head])
        elif head.startswith('UR'):
            items = head.split('.')
            bc, pos = int(items[1]), int(items[2])
            chain[bc].append([pos, 0, head])
        else:
            if len(seq) >= min_length:
                optional[head] = seq
    primary, alt = {}, {}
    pair = {}
    for bc, info in chain.items():
        info.sort(key = lambda x: (x[0], x[1]))
        st, en, prev = 0, 0, 0
        for item in info:
            if item[1] == 0:
                primary.setdefault(bc, []).append(item[2])
                st = en
                en = st + len(contigs[item[2]])
                prev = item[0]
            else:
                if len(contigs[item[2]]) >= min_length:
                    alt.setdefault(bc, []).append(item[2])
                    if item[0] == prev:
                        pair[item[2]] = [bc, st, en]

    # generate primary and alt contigs
    pi = 0
    mapping, name_mapping = {}, {}
    pri_contigs, alt_contigs = {}, {}
    with open(prefix + '.primary.fasta', 'w') as f:
        for bc, items in primary.items():
            name = '{:0>6d}F'.format(pi)
            mapping[bc] = name
            seq = []
            f.write('>{}\n'.format(name))
            for item in items:
                name_mapping[item] = name
                f.write(contigs[item])
                seq.append(contigs[item])
            f.write('\n')
            pri_contigs[name] = ''.join(seq)
            pi += 1
        pi += 10
        for head, seq in optional.items():
            name = '{:0>6d}F'.format(pi)
            name_mapping[head] = name
            f.write('>{}\n'.format(name))
            f.write(seq + '\n')
            pri_contigs[name] = seq
            pi += 1
    with open(prefix + '.alt.fasta', 'w') as f:
        for bc, items in alt.items():
            # if len(items) == 1: 
            #     pair.pop(items[0])
            #     continue
            ai = 1
            for item in items:
                name = '{}_{:0>3d}'.format(mapping[bc], ai)
                name_mapping[item] = name
                f.write('>{}\n'.format(name))
                f.write(contigs[item] + '\n')
                alt_contigs[name] = contigs[item]
                ai += 1
    with open(prefix + '.name_mapping', 'w') as f:
        for head, name in name_mapping.items():
            f.write('{}\t{}\n'.format(head, name))

    # generate pair file and homogolous file
    pair_new = []
    pri_cut, alt_cut, coll = {}, {}, {}
    with open(prefix + '.pair.cut.txt', 'w') as f:
        for alt, [pri, st, en] in pair.items():
            f.write('{}:{}-{}\t{}:{}-{}\n'.format(mapping[pri], st, en, name_mapping[alt], 0, len(contigs[alt])))
            pri_cut.setdefault(mapping[pri], []).append([st, en])
            alt_cut.setdefault(name_mapping[alt], []).append([0, len(contigs[alt])])
            pair_new.append([mapping[pri], st, en, name_mapping[alt], 0, len(contigs[alt])])
    
    for pri, items in pri_cut.items():
        items.sort(key = lambda x: x[0])
        prev = 0
        for st, en in items:
            if st > prev:
                coll.setdefault(pri, []).append([prev, st])
            prev = en
        if len(pri_contigs[pri]) > prev:
            coll.setdefault(pri, []).append([prev, len(pri_contigs[pri])])

    with open(prefix + '.hap1.cut.fasta', 'w') as f:
        for head, seq in alt_contigs.items():
            if head not in alt_cut: continue
            interval = alt_cut[head]
            interval.sort(key = lambda x: x[0])
            for st, en in interval:
                if st > en or st < 0 or en < 0:
                    logger.error('Invalid interval of {}:{}-{}'.format(head, st, en))
                    exit(1)
                f.write('>{}:{}-{}\n'.format(head, st, en))
                f.write(seq[st: en] + '\n')

    with open(prefix + '.hap2.cut.fasta', 'w') as fh, open(prefix + '.collapsed.cut.fasta', 'w') as fc:
        for head, seq in pri_contigs.items():
            if head not in pri_cut:
                length = len(seq)
                coll.setdefault(head, []).append([0, length])
                fc.write('>{}:0-{}\n'.format(head, length))
                fc.write(seq + '\n')
                continue
            interval = pri_cut[head]
            interval.sort(key = lambda x: x[0])
            for st, en in interval:
                if st > en or st < 0 or en < 0:
                    logger.error('Invalid interval of {}:{}-{}'.format(head, st, en))
                    exit(1)
                fh.write('>{}:{}-{}\n'.format(head, st, en))
                fh.write(seq[st: en] + '\n')
            if head not in coll: continue
            interval = coll[head]
            interval.sort(key = lambda x: x[0])
            for st, en in interval:
                if st >= en or st < 0 or en < 0:
                    logger.error('Invalid interval of {}:{}-{}'.format(head, st, en))
                    exit(1)
                fc.write('>{}:{}-{}\n'.format(head, st, en))
                fc.write(seq[st: en] + '\n')

    with open(prefix + '.collapsed.cut.bed', 'w') as f:
        for head, interval in sorted(coll.items(), key = lambda x: x[0]):
            interval.sort(key = lambda x: x[0])
            for st, en in interval:
                f.write('{}\t{}\t{}\n'.format(head, st, en))

    with open(prefix + '.hap1.cut.bed', 'w') as f:
        for head, interval in sorted(alt_cut.items(), key = lambda x: x[0]):
            interval.sort(key = lambda x: x[0])
            for st, en in interval:
                f.write('{}\t{}\t{}\n'.format(head, st, en))

    with open(prefix + '.hap2.cut.bed', 'w') as f:
        for head, interval in sorted(pri_cut.items(), key = lambda x: x[0]):
            interval.sort(key = lambda x: x[0])
            for st, en in interval:
                f.write('{}\t{}\t{}\n'.format(head, st, en))

    # mapping 
    for [pri, st, en, alt, ast, aen] in pair_new:
        file = os.path.join(dir, '{}.{}.{}.p.fasta'.format(pri, st, en))
        with open(file, 'w') as fw:
            fw.write('>{}:{}-{}\n'.format(pri, st, en))
            fw.write(pri_contigs[pri][st: en] + '\n')
        file = os.path.join(dir, '{}.a.fasta'.format(alt))
        with open(file, 'w') as fw:
            fw.write('>{}\n'.format(alt))
            fw.write(alt_contigs[alt][ast: aen] + '\n')

    pool = multiprocessing.Pool(threads)
    result = []
    for [pri, st, en, alt, ast, aen] in pair_new:
        result.append(pool.apply_async(map_separate, (os.path.join(dir, '{}.{}.{}.p.fasta'.format(pri, st, en)), os.path.join(dir, '{}.a.fasta'.format(alt)), 'paf', )))
    pool.close()
    pool.join()

    pattern = re.compile('([0-9]{6}F):([0-9]*)-([0-9]*)')
    with open(ofname, 'w') as fw:
        for rs in result:
            for r in rs.get():
                if len(r) == 0: continue
                items = r.strip().split()
                mat = pattern.match(items[5])
                if mat == None or len(mat.groups()) != 3:
                    logger.error('Invalid header of pri {}'.format(items[5]))
                    exit(1)
                items[5] = mat.group(1)
                items[6] = str(len(pri_contigs[items[5]]))
                items[7] = str(int(items[7]) + int(mat.group(2)))
                items[8] = str(int(items[8]) + int(mat.group(2)))
                fw.write('\t'.join(items) + '\n')

def shasta2unzip(ifname, prefix, min_length):
    contigs = {}
    ftype = _fa_or_fq(ifname)
    if ftype == None:
        logger.error('Unrecognized file type {}'.format(ifname))
        exit(1)
    for head, seq in eval('_load_' + ftype)(ifname):
        contigs[head] = seq
    logger.info('load {} segments in {}'.format(len(contigs), ifname))
    chain = defaultdict(list)
    optional = {}
    for head, seq in contigs.items():
        if head.startswith('PR'):
            items = head.split('.')
            bc, pos, _, h = int(items[1]), int(items[2]), items[3], int(items[4])
            chain[bc].append([pos, h, head])
        elif head.startswith('UR'):
            items = head.split('.')
            bc, pos = int(items[1]), int(items[2])
            chain[bc].append([pos, 0, head])
        else:
            if len(seq) >= min_length:
                optional[head] = seq
    primary, alt = {}, {}
    for bc, info in chain.items():
        info.sort(key = lambda x: (x[0], x[1]))
        for item in info:
            if item[1] == 0:
                primary.setdefault(bc, []).append(item[2])
            else:
                if len(contigs[item[2]]) >= min_length:
                    alt.setdefault(bc, []).append(item[2])
    
    pi = 0
    mapping, name_mapping = {}, {}
    with open(prefix + '.primary.fasta', 'w') as f:
        for bc, items in primary.items():
            name = '{:0>6d}F'.format(pi)
            mapping[bc] = name
            f.write('>{}\n'.format(name))
            for item in items:
                name_mapping[item] = name
                f.write(contigs[item])
            f.write('\n')
            pi += 1
        pi += 10
        for head, seq in optional.items():
            name = '{:0>6d}F'.format(pi)
            name_mapping[head] = name
            f.write('>{}\n'.format(name))
            f.write(seq + '\n')
            pi += 1
    with open(prefix + '.alt.fasta', 'w') as f:
        for bc, items in alt.items():
            # if len(items) == 1: continue
            ai = 1
            for item in items:
                name = '{}_{:0>3d}'.format(mapping[bc], ai)
                name_mapping[item] = name
                f.write('>{}\n'.format(name))
                f.write(contigs[item] + '\n')
                ai += 1
    with open(prefix + '.name_mapping', 'w') as f:
        for head, name in name_mapping.items():
            f.write('{}\t{}\n'.format(head, name))

def fix_switch2(sfname, afname, pfname, ofname):
    # sfname = args.sfname
    # afname = args.afname
    # pfname = args.pfname
    # ofname = args.ofname

    switch = {}
    name_mapping = {}
    if sfname != None:
        for line in open(sfname, 'r'):
            items = line.strip().split()
            if items[2] != 'S' and items[2] != 'F' and items[2] != 'C': continue
            switch[items[0]] = [int(items[1]), items[2]]
            name_mapping[items[0].split(':')[0]] = items[0]
    pair = {}
    for line in open(afname, 'r'):
        items = line.strip().split()
        if items[1] not in switch: continue
        pair[items[1]] = items[0]
        pair[items[0]] = items[1]
    pattern = re.compile(r'(\d+)([MIDNSHP=X])')
    fixed = {}
    for line in open(pfname, 'r'):
        items = line.strip().split()
        if len(items) < 11: continue
        if items[0] not in name_mapping: continue
        if name_mapping[items[0]] in fixed: continue
        alt = int(items[2])
        pri = int(items[7])
        if items[4] == '-':
            pri = int(items[6]) - int(items[8])
        cigar = ''
        for i in range(11, len(items)):
            if items[i].startswith('cg:Z:'):
                cigar = items[i][5:]
                break
        offset = int(pair[name_mapping[items[0]]].split(':')[1].split('-')[0])
        for mat in pattern.finditer(cigar):
            l, op = int(mat.group(1)), mat.group(2)
            if op == 'M' or op == 'X' or op == '=':
                if alt < switch[name_mapping[items[0]]][0] and alt + l > switch[name_mapping[items[0]]][0]:
                    fixed[name_mapping[items[0]]] = switch[name_mapping[items[0]]][0]
                    fixed[pair[name_mapping[items[0]]]] = pri + switch[name_mapping[items[0]]][0] - alt - offset
                    break
                alt += l
                pri += l
            elif op == 'I' or op == 'S':
                if alt < switch[name_mapping[items[0]]][0] and alt + l > switch[name_mapping[items[0]]][0]:
                    fixed[name_mapping[items[0]]] = alt + l
                    fixed[pair[name_mapping[items[0]]]] = pri - offset
                    break
                alt += l
            elif op == 'D' or op == 'N':
                pri += l
    def _offset(a, b):
        if a / b > 1.1 or a / b < 0.9:
            return True
    with open(ofname, 'w') as fw:
        for h, p in fixed.items():
            if _offset(p, switch[h][0]) or _offset(fixed[pair[h]], switch[pair[h]][0]):
                continue
            fw.write('{}\t{}\t{}\n'.format(h, p, switch[h][1]))

def fix_switch(args):
    fix_switch2(args.sfname, args.afname, args.pfname, args.ofname)

def group2(h1fname, h2fname, gfname, o1fname, o2fname):
    # h1fname = args.h1fname
    # h2fname = args.h2fname
    # gfname = args.gfname
    # o1fname = args.o1fname
    # o2fname = args.o2fname

    hap1, hap2 = {}, {}
    ftype = _fa_or_fq(h1fname)
    if ftype == None:
        logger.error('Unrecognized file type {}'.format(h1fname))
        exit(1)
    for head, seq in eval('_load_' + ftype)(h1fname):
        hap1[head] = seq
    ftype = _fa_or_fq(h2fname)
    if ftype == None:
        logger.error('Unrecognized file type {}'.format(h2fname))
        exit(1)
    for head, seq in eval('_load_' + ftype)(h2fname):
        hap2[head] = seq
    group1, group2 = [], []
    for line in open(gfname, 'r'):
        items = line.strip().split()
        group1.append(items[0])
        group2.append(items[1])
    with open(o1fname, 'w') as fw:
        for g in group1:
            if g in hap1:
                fw.write('>{}\n{}\n'.format(g, hap1[g]))
            elif g in hap2:
                fw.write('>{}\n{}\n'.format(g, hap2[g]))
    with open(o2fname, 'w') as fw:
        for g in group2:
            if g in hap1:
                fw.write('>{}\n{}\n'.format(g, hap1[g]))
            elif g in hap2:
                fw.write('>{}\n{}\n'.format(g, hap2[g]))

def group(args):
    group2(args.h1fname, args.h2fname, args.gfname, args.o1fname, args.o2fname)

def hapdup2unzip(ifname, prefix):
    ctg = {}
    ftype = _fa_or_fq(ifname)
    if ftype == None:
        logger.error('Unrecognized file type {}'.format(ifname))
        exit(1)
    for head, seq in eval('_load_' + ftype)(ifname):
        items = head.split('_')
        # if len(items) != 4:
        #     print(head, file = sys.stderr)
        ctg.setdefault(items[1], {})[items[-1]] = seq
    with open(prefix + '.fasta', 'w') as f:
        for head, seqs in ctg.items():
            f.write('>contig{}\n'.format(head))
            for id, seq in sorted(seqs.items(), key = lambda x: int(x[0])):
                f.write(seq)
            f.write('\n')

import uuid
def hap_map(pfname1, pfname2, threads, fmt, prefix):
    ctg1, ctg2 = {}, {}
    ftype = _fa_or_fq(pfname1)
    if ftype == None:
        logger.error('Unrecognized file type {}'.format(pfname1))
        exit(1)
    for head, seq in eval('_load_' + ftype)(pfname1):
        ctg1[head] = seq
    ftype = _fa_or_fq(pfname2)
    if ftype == None:
        logger.error('Unrecognized file type {}'.format(pfname2))
        exit(1)
    for head, seq in eval('_load_' + ftype)(pfname2):
        ctg2[head] = seq
    
    tmp_dir = str(uuid.uuid4())
    os.makedirs(tmp_dir, exist_ok = True)
    for head1, seq in ctg1.items():
        if head1[-1:] == '1':
            head2 = head1[:-1] + '2'
        elif head[-1:] == '2':
            head2 = head1[:-1] + '1'
        if head2 in ctg2:
            ref = os.path.join(tmp_dir, head1 + '.ref.fasta')
            qry = os.path.join(tmp_dir, head2 + '.qry.fasta')
            with open(ref, 'w') as fw:
                fw.write('>{}\n'.format(head1))
                fw.write(seq + '\n')
            with open(qry, 'w') as fw:
                fw.write('>{}\n'.format(head2))
                fw.write(ctg2[head] + '\n')

    results = []
    pool = multiprocessing.Pool(threads)
    for head1, seq in ctg1.items():
        if head1[-1:] == '1':
            head2 = head1[:-1] + '2'
        elif head[-1:] == '2':
            head2 = head1[:-1] + '1'
        if head2 in ctg2:
            ref = os.path.join(tmp_dir, head1 + '.ref.fasta')
            qry = os.path.join(tmp_dir, head2 + '.qry.fasta')
            results.append(pool.apply_async(map_separate, (ref, qry, fmt, )))
    pool.close()
    pool.join()

    sq, pg, record = [], [], []
    for rs in results:
        for r in rs.get():
            if r.startswith('@SQ'):
                sq.append(r)
            elif r.startswith('@PG'):
                pg.append(r)
            else:
                record.append(r)
    ofname = prefix + '.hapmap.' + fmt
    with open(ofname, 'w') as fw:
        for s in sq:
            fw.write(s + '\n')
        for p in pg:
            fw.write(p + '\n')
            break
        for r in record:
            fw.write(r + '\n')
    os.removedirs(tmp_dir)

def cut_contigs_dual(d1fname, d2fname, pair, prefix):
    l_dual1, l_dual2 = {}, {}
    # write pair information
    dual1, dual2, coll1, coll2 = {}, {}, {}, {}
    with open(prefix + '.pair.cut.txt', 'w') as fw:
        # [query, target]
        for [query, ql, qs, qe, strand, target, tl, ts, te] in pair:
            fw.write('{}:{}-{}\t{}:{}-{}\t{}\n'.format(query, qs, qe, target, ts, te, strand))
            dual1.setdefault(target, []).append([ts, te])
            dual2.setdefault(query, []).append([qs, qe])
            l_dual1[target] = tl
            l_dual2[query] = ql
    # get collasped regions of dual1 and dual2
    for ctg, items in dual1.items():
        items.sort(key = lambda x: x[0])
        prev = 0
        for st, en in items:
            if st > prev:
                coll1.setdefault(ctg, []).append([prev, st])
            prev = en
        if l_dual1[ctg] > prev:
            coll1.setdefault(ctg, []).append([prev, l_dual1[ctg]])
    for ctg, items in dual2.items():
        items.sort(key = lambda x: x[0])
        prev = 0
        for st, en in items:
            if st > prev:
                coll2.setdefault(ctg, []).append([prev, st])
            prev = en
        if l_dual2[ctg] > prev:
            coll2.setdefault(ctg, []).append([prev, l_dual2[ctg]])
    # write dual1 and dual2 contigs and collapsed contigs
    ftype = _fa_or_fq(d1fname)
    if ftype == None:
        logger.error('Unrecognized file type {}'.format(d1fname))
        exit(1)
    with open(prefix + '.hap1.cut.fasta', 'w') as fw, open(prefix + '.collapsed1.cut.fasta', 'w') as fw2:
        for head, seq in eval('_load_' + ftype)(d1fname):
            if head not in dual1: 
                fw2.write('>{}:{}-{}\n'.format(head, 0, len(seq)))
                fw2.write(seq + '\n')
                coll1.setdefault(head, []).append([0, len(seq)])
                continue
            interval = dual1[head]
            interval.sort(key = lambda x: x[0])
            for st, en in interval:
                if st >= en or st < 0 or en < 0:
                    logger.error('Invalid interval of {}:{}-{}'.format(head, st, en))
                    exit(1)
                fw.write('>{}:{}-{}\n'.format(head, st, en))
                fw.write(seq[st: en] + '\n')
            if head not in coll1: continue
            interval = coll1[head]
            interval.sort(key = lambda x: x[0])
            for st, en in interval:
                if st >= en or st < 0 or en < 0:
                    logger.error('Invalid interval of {}:{}-{}'.format(head, st, en))
                    exit(1)
                fw2.write('>{}:{}-{}\n'.format(head, st, en))
                fw2.write(seq[st: en] + '\n')
    ftype = _fa_or_fq(d2fname)
    if ftype == None:
        logger.error('Unrecognized file type {}'.format(d2fname))
        exit(1)
    with open(prefix + '.hap2.cut.fasta', 'w') as fw, open(prefix + '.collapsed2.cut.fasta', 'w') as fw2:
        for head, seq in eval('_load_' + ftype)(d2fname):
            if head not in dual2: 
                fw2.write('>{}:{}-{}\n'.format(head, 0, len(seq)))
                fw2.write(seq + '\n')
                coll2.setdefault(head, []).append([0, len(seq)])
                continue
            interval = dual2[head]
            interval.sort(key = lambda x: x[0])
            for st, en in interval:
                if st >= en or st < 0 or en < 0:
                    logger.error('Invalid interval of {}:{}-{}'.format(head, st, en))
                    exit(1)
                fw.write('>{}:{}-{}\n'.format(head, st, en))
                fw.write(seq[st: en] + '\n')
            if head not in coll2: continue
            interval = coll2[head]
            interval.sort(key = lambda x: x[0])
            for st, en in interval:
                if st >= en or st < 0 or en < 0:
                    logger.error('Invalid interval of {}:{}-{}'.format(head, st, en))
                    exit(1)
                fw2.write('>{}:{}-{}\n'.format(head, st, en))
                fw2.write(seq[st: en] + '\n')
    # write bed file
    with open(prefix + '.hap1.cut.bed', 'w') as fw:
        for head, interval in sorted(dual1.items(), key = lambda x: x[0]):
            interval.sort(key = lambda x: x[0])
            for st, en in interval:
                fw.write('{}\t{}\t{}\n'.format(head, st, en))
    with open(prefix + '.hap2.cut.bed', 'w') as fw:
        for head, interval in sorted(dual2.items(), key = lambda x: x[0]):
            interval.sort(key = lambda x: x[0])
            for st, en in interval:
                fw.write('{}\t{}\t{}\n'.format(head, st, en))
    with open(prefix + '.collapsed1.cut.bed', 'w') as fw:
        for head, interval in sorted(coll1.items(), key = lambda x: x[0]):
            interval.sort(key = lambda x: x[0])
            for st, en in interval:
                fw.write('{}\t{}\t{}\n'.format(head, st, en))
    with open(prefix + '.collapsed2.cut.bed', 'w') as fw:
        for head, interval in sorted(coll2.items(), key = lambda x: x[0]):
            interval.sort(key = lambda x: x[0])
            for st, en in interval:
                fw.write('{}\t{}\t{}\n'.format(head, st, en))

def process_one_ref(pafs, threshold, fw):
    # get most likely mapping contig
    pafs.sort(key = lambda x: (x.split()[0], int(x.split()[2])))
    covered = {}
    for paf in pafs:
        items = paf.strip().split()
        if items[5] not in covered:
            covered.setdefault(items[0] + items[4], 0)
        covered[items[0] + items[4]] += int(items[8]) - int(items[7])
    longest_covered = max(covered.values())
    cand, strand = '', ''
    for ctg, length in covered.items():
        if length == longest_covered:
            cand, strand = ctg[:-1], ctg[-1]
            break
    # get mapping interval
    pattern = re.compile(r'(\d+)([MIDNSHP=X])')
    interval1, interval2 = [], []
    query_len, target_len = 0, 0
    query, target = '', ''
    for paf in pafs:
        items = paf.strip().split()
        if items[0] != cand or items[4] != strand: continue
        fw.write(paf + '\n')
        query, ql, qs, qe = items[0], int(items[1]), int(items[2]), int(items[3])
        target, tl, ts, te = items[5], int(items[6]), int(items[7]), int(items[8])
        query_len, target_len = ql, tl
        for i in range(12, len(items)):
            if items[i].startswith('cg:Z:'):
                cg = items[i][5:]
                break
        curq, curt = qs, ts
        for m in pattern.finditer(cg):
            l, op = int(m.group(1)), m.group(2)
            if op == '=' or op == 'M':
                if l == ql:
                    continue
                if l >= threshold:
                    interval1.append([curq, curq + l])
                    interval2.append([curt, curt + l])
                curq, curt = curq + l, curt + l
            elif op == 'I':
                curq += l
            elif op == 'D':
                curt += l
            elif op == 'X':
                curq, curt = curq + l, curt + l
    # merge interval
    if len(interval1) == 0: return '', '', 0, 0, strand, [], []
    interval1_new, interval2_new = [], []
    index = 0
    if interval1[0][0] + interval2[0][0] != 0 and (interval1[0][0] < threshold or interval2[0][0] < threshold):
        index = 1
    if index == len(interval1):
        return '', '', 0, 0, strand, [], []
    interval1_new.append(interval1[index])
    interval2_new.append(interval2[index])
    index += 1
    if index != len(interval1):
        for i in range(index, len(interval1) - 1):
            if interval1[i][0] - interval1_new[-1][1] <= threshold:
                if interval1[i][1] - interval1[i][0] >= interval1_new[-1][1] - interval1_new[-1][0]:
                    interval1_new[-1] = interval1[i]
            elif interval2[i][0] - interval2_new[-1][1] <= threshold:
                if interval2[i][1] - interval2[i][0] >= interval2_new[-1][1] - interval2_new[-1][0]:
                    interval2_new[-1] = interval2[i]
            else:
                interval1_new.append(interval1[i])
                interval2_new.append(interval2[i])
        if (query_len - interval1[-1][1] >= threshold and target_len - interval2[-1][1] >= threshold) or (interval1[-1][1] == query_len and interval2[-1][1] == target_len):
            if interval1[-1][0] - interval1_new[-1][1] <= threshold or interval2[-1][0] - interval2_new[-1][1] <= threshold:
                if interval1[-1][1] - interval1[-1][0] >= interval1_new[-1][1] - interval1_new[-1][0]:
                    interval1_new[-1] = interval1[-1]
                if interval2[-1][1] - interval2[-1][0] >= interval2_new[-1][1] - interval2_new[-1][0]:
                    interval2_new[-1] = interval2[-1]
            else:
                interval1_new.append(interval1[-1])
                interval2_new.append(interval2[-1])
    if len(interval1_new) == 0: return '', '', 0, 0, strand, [], []
    interval1.clear()
    interval2.clear()
    if interval1_new[0][0] != 0 and interval2_new[0][0] != 0:
        interval1.append([0, interval1_new[0][0]])
        interval2.append([0, interval2_new[0][0]])
    for i in range(len(interval1_new) - 1):
        interval1.append([interval1_new[i][1], interval1_new[i + 1][0]])
        interval2.append([interval2_new[i][1], interval2_new[i + 1][0]])
    if interval1_new[-1][1] != query_len and interval2_new[-1][1] != target_len:
        interval1.append([interval1_new[-1][1], query_len])
        interval2.append([interval2_new[-1][1], target_len])
    return query, target, query_len, target_len, strand, interval1, interval2

def process_one_paf(paf, threshold):
    pattern = re.compile(r'(\d+)([MIDNSHP=X])')
    interval1, interval2 = [], []
    paf = paf.strip().split()
    ctg1, len1, st1, en1 = paf[0], int(paf[1]), int(paf[2]), int(paf[3])
    ctg2, len2, st2, en2 = paf[5], int(paf[6]), int(paf[7]), int(paf[8])
    for i in range(12, len(paf)):
        if paf[i].startswith('cg:Z:'):
            cg = paf[i][5:]
            break
    cur1, cur2 = st1, st2
    for m in pattern.finditer(cg):
        l, op = int(m.group(1)), m.group(2)
        if op == '=' or op == 'M':
            if l == len2:
                return ctg1, ctg2, len1, len2, [], []
            if l >= threshold:
                interval1.append([cur1, cur1 + l])
                interval2.append([cur2, cur2 + l])
                cur1, cur2 = cur1 + l, cur2 + l
            else:
                cur1, cur2 = cur1 + l, cur2 + l
        elif op == 'I':
            cur2 += l
        elif op == 'D':
            cur1 += l
        elif op == 'X':
            cur1, cur2 = cur1 + l, cur2 + l
    return ctg1, ctg2, len1, len2, interval1, interval2

def merge_intervals(interval1, interval2, len1, len2, threshold):
    interval1_new, interval2_new = [], []
    index = 0
    if interval1[0][0] < threshold or interval2[0][0] < threshold:
        index = 1
    if index == len(interval1):
        return [], []
    interval1_new.append(interval1[index])
    interval2_new.append(interval2[index])
    index += 1
    for i in range(index, len(interval1) - 1):
        if interval1[i][0] - interval1_new[-1][1] <= threshold:
            if interval1[i][1] - interval1[i][0] >= interval1_new[-1][1] - interval1_new[-1][0]:
                interval1_new[-1] = interval1[i]
        elif interval2[i][0] - interval2_new[-1][1] <= threshold:
            if interval2[i][1] - interval2[i][0] >= interval2_new[-1][1] - interval2_new[-1][0]:
                interval2_new[-1] = interval2[i]
        else:
            interval1_new.append(interval1[i])
            interval2_new.append(interval2[i])
    if len1 - interval1[-1][1] >= threshold and len2 - interval2[-1][1] >= threshold:
        interval1_new.append(interval1[-1])
        interval2_new.append(interval2[-1])
    return interval1_new, interval2_new

def dual_cut(d1fname, d2fname, pfname, threshold, threads, prefix):
    pafs = {}
    for line in open(pfname, 'r'):
        items = line.strip().split()
        pafs.setdefault(items[5], []).append(line.strip())
    pairs = []
    fw = open(prefix + '.dual.filtered.paf', 'w')
    for ctg, paf in pafs.items():
        query, target, query_len, target_len, strand, interval1, interval2 = process_one_ref(paf, threshold, fw)
        if len(interval1) == 0: continue
        for i in range(len(interval1)):
            pairs.append([query, query_len, interval1[i][0], interval1[i][1], strand, target, target_len, interval2[i][0], interval2[i][1]])
    fw.close()
    cut_contigs_dual(d1fname, d2fname, pairs, prefix)

def dual_cut1(pfname, threshold, threads, prefix):
    pafs = {}
    for line in open(pfname, 'r'):
        items = line.strip().split()
        pafs.setdefault(items[5], []).append(line.strip())
    
    fw = open(prefix + '.dual.filtered.paf', 'w')
    fw1 = open(prefix + '.homo.txt', 'w')
    for ctg, paf in pafs.items():
        query, target, query_len, target_len, strand, interval1, interval2 = process_one_ref(paf, threshold, fw)
        if len(interval1) == 0: continue
        for i in range(len(interval1)):
            fw1.write('{} {} {} {} {} {} {} {} {}\n'.format(query, query_len, interval1[i][0], interval1[i][1], strand, target, target_len, interval2[i][0], interval2[i][1]))
    fw.close()
    fw1.close()

def filter_short_reads(ifname, ofname, ffname, min_length):
    ftype = _fa_or_fq(ifname)
    if ftype == None:
        logger.error('Unrecognized file type {}'.format(ifname))
        exit(1)
    with open(ofname, 'w') as fw, open(ffname, 'w') as ff:
        for head, seq in eval('_load_' + ftype)(ifname):
            if len(seq) >= min_length:
                fw.write('>{}\n{}\n'.format(head, seq))
            else:
                ff.write('>{}\n{}\n'.format(head, seq))

def process_one_ref_advance(pafs, threshold):
    pattern = re.compile(r'(\d+)([MIDNSHP=X])')
    interval, homointer = [], []
    for paf in pafs:
        items = paf.strip().split()
        query, ql, qs, qe = items[0], int(items[1]), int(items[2]), int(items[3])
        target, tl, ts, te = items[5], int(items[6]), int(items[7]), int(items[8])
        strand = items[4]
        for i in range(12, len(items)):
            if items[i].startswith('cg:Z:'):
                cg = items[i][5:]
                break
        curq, curt = qs, ts
        if strand == '-':
            curq, qs = qe, qe
        homo = False
        for m in pattern.finditer(cg):
            l, op = int(m.group(1)), m.group(2)
            if op == '=' or op == 'M':
                if l == ql:
                    continue
                if strand == '+':
                    if l >= threshold and curq - qs >= threshold and curt - ts >= threshold:
                        interval.append([query, ql, qs, curq, strand, target, tl, ts, curt])
                        homointer.append([query, curq, curq + l, target, curt, curt + l, strand])
                        qs, ts = curq + l, curt + l
                        homo = True
                    curq, curt = curq + l, curt + l
                else:
                    if l >= threshold and qs - curq >= threshold and curt - ts >= threshold:
                        interval.append([query, ql, curq, qs, strand, target, tl, ts, curt])
                        homointer.append([query, curq - l, curq, target, curt, curt + l, strand])
                        qs, ts = curq - l, curt + l
                        homo = True
                    curq, curt = curq - l, curt + l
            elif op == 'I':
                if strand == '+':
                    curq += l
                else:
                    curq -= l
            elif op == 'D':
                curt += l
            elif op == 'X':
                if strand == '+':
                    curq, curt = curq + l, curt + l
                else:
                    curq, curt = curq - l, curt + l
        if strand == '+':
            if homo and curq - qs >= threshold and curt - ts >= threshold:
                interval.append([query, ql, qs, curq, strand, target, tl, ts, curt])
        else:
            if homo and qs - curq >= threshold and curt - ts >= threshold:
                interval.append([query, ql, curq, qs, strand, target, tl, ts, curt])
    return interval, homointer

def _contained_or_overlap(st1, en1, st2, en2):
    if st2 >= st1 and en2 <= en1:
        return True, False
    if st1 >= st2 and en1 <= en2:
        return True, False
    if st1 > st2 and st1 < en2 and en1 > en2:
        return False, True
    if st2 > st1 and st2 < en1 and en2 > en1:
        return False, True
    return False, False

def replace_mapping(pfname, threshold, prefix):
    mappings = {}
    with open(prefix + '.filtered.paf', 'w') as fw:
        for line in open(pfname, 'r'):
            items = line.strip().split()
            if int(items[3]) - int(items[2]) < threshold * 3 or int(items[8]) - int(items[7]) < threshold * 3: continue
            mappings.setdefault(items[5], []).append(line.strip())
            fw.write(line)
    
    interval = []
    for ref, pafs in mappings.items():
        pafs.sort(key = lambda x: int(x.split()[7]))
        inval, homo = process_one_ref_advance(pafs, threshold)
        interval.extend(inval)
    interval.sort(key = lambda x: (x[5], x[7]))

    contained, overlap = [0, 0], [0, 0]
    for i in range(len(interval)):
        for j in range(i + 1, len(interval)):
            if interval[j][5] == interval[i][5]:
                c, o = _contained_or_overlap(interval[j][7], interval[j][8], interval[i][7], interval[i][8])
                if c:
                    contained[0] += 1
                    print('target', interval[j], interval[i])
                if o:
                    overlap[0] += 1
            if interval[j][0] == interval[i][0]:
                c, o = _contained_or_overlap(interval[j][2], interval[j][3], interval[i][2], interval[i][3])
                if c:
                    contained[1] += 1
                    print('query', interval[j], interval[i])
                if o:
                    overlap[1] += 1
    print('contained: {} {} overlap: {} {}'.format(contained[0], contained[1], overlap[0], overlap[1]))

    with open(prefix + '.homo.txt', 'w') as fw:
        for [query, ql, qs, qe, strand, target, tl, ts, te] in interval:
            fw.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(query, ql, qs, qe, strand, target, tl, ts, te))

def filter_one_target(pafs, threshold):
    pattern = re.compile(r'(\d+)([MIDNSHP=X])')
    filtered = []
    for i in range(len(pafs)):
        contained = False
        itemsi = pafs[i].strip().split()
        for j in range(i):
            itemsj = pafs[j].strip().split()
            if int(itemsi[8]) <= int(itemsj[8]):
                contained = True
                break
        if not contained:
            for k in range(12, len(itemsi)):
                if itemsi[k].startswith('cg:Z:'):
                    cg = itemsi[k][5:]
            for m in pattern.finditer(cg):
                l, op = int(m.group(1)), m.group(2)
                if op == '=' or op == 'M':
                    if l == int(itemsi[2]):
                        continue
                    if l >= threshold:
                        filtered.append(pafs[i])
                        break
    return filtered

def filter_one_query(pafs):
    filtered = []
    for i in range(len(pafs)):
        contained = False
        itemsi = pafs[i].strip().split()
        for j in range(i):
            itemsj = pafs[j].strip().split()
            if int(itemsi[3]) <= int(itemsj[3]):
                contained = True
                break
        if not contained:
            filtered.append(pafs[i])
    return filtered

def chain_dual(d1fname, d2fname, pfname, percent, prefix):
    pattern = re.compile(r'(\d+)([MIDNSHP=X])')
    percentle = []
    for line in open(pfname):
        items = line.strip().split()
        for i in range(12, len(items)):
            if items[i].startswith('cg:Z:'):
                cg = items[i][5:]
                break
        for m in pattern.finditer(cg):
            l, op = int(m.group(1)), m.group(2)
            if op == '=' or op == 'M':
                if l == int(items[1]): continue
                if l >= 10000:
                    percentle.append(l)
    percentle.sort()
    # for i in range(101):
    #     print('percentile {} {}'.format(i, percentle[int(i * len(percentle) / 100) - 1]))
    # return
    if percent > 100:
        threshold = percent
    else:
        threshold = int((percent / 100) * len(percentle)) - 1 if percent > 1 else int(percent * len(percentle)) - 1
        threshold = percentle[threshold]
    logger.info('threshold {} for filtering paf'.format(threshold))

    target, query = {}, {}
    for line in open(pfname):
        items = line.strip().split()
        if int(items[3]) - int(items[2]) < threshold * 3 or int(items[8]) - int(items[7]) < threshold * 3: continue
        target.setdefault(items[5], []).append(line.strip())

    for ref, pafs in target.items():
        pafs.sort(key = lambda x: int(x.split()[7]))
        filtered = filter_one_target(pafs, threshold)
        for paf in filtered:
            items = paf.strip().split()
            query.setdefault(items[0], []).append(paf)
    target.clear()

    for _, pafs in query.items():
        pafs.sort(key = lambda x: int(x.split()[2]))
        filtered = filter_one_query(pafs)
        for paf in filtered:
            items = paf.strip().split()
            target.setdefault(items[5], []).append(paf)

    interval = []
    with open(prefix + '.filtered.paf', 'w') as fw, open(prefix + '.homo.txt', 'w') as fh:
        for _, pafs in target.items():
            pafs.sort(key = lambda x: int(x.split()[7]))
            inval, homo = process_one_ref_advance(pafs, threshold)
            interval.extend(inval)
            for paf in pafs:
                fw.write(paf + '\n')
            for [query, qs, qe, target, ts, te, strand] in homo:
                fh.write('{}:{}-{}\t{}:{}-{}\t{}\n'.format(query, qs, qe, target, ts, te, strand))
    cut_contigs_dual(d1fname, d2fname, interval, prefix)
    return
    contained, overlap = [0, 0], [0, 0]
    for i in range(len(interval)):
        for j in range(i + 1, len(interval)):
            if interval[j][5] == interval[i][5]:
                c, o = _contained_or_overlap(interval[j][7], interval[j][8], interval[i][7], interval[i][8])
                if c:
                    contained[0] += 1
                    print('target', interval[j], interval[i])
                if o:
                    overlap[0] += 1
            if interval[j][0] == interval[i][0]:
                c, o = _contained_or_overlap(interval[j][2], interval[j][3], interval[i][2], interval[i][3])
                if c:
                    contained[1] += 1
                    print('query', interval[j], interval[i])
                if o:
                    overlap[1] += 1
    print('contained: {} {} overlap: {} {}'.format(contained[0], contained[1], overlap[0], overlap[1]))

    with open(prefix + '.homo.txt', 'w') as fw:
        for [query, ql, qs, qe, strand, target, tl, ts, te] in interval:
            fw.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(query, ql, qs, qe, strand, target, tl, ts, te))

def merge_hap(pafs):
    if len(pafs) == 0: return
    result = [{'score': float('-inf'), 'start': -1, 'end': -1} for _ in range(len(pafs))]
    for i in range(len(pafs)):
        itemsi = pafs[i].strip().split()
        if result[i]['score'] < int(itemsi[3]) - int(itemsi[2]):
            result[i]['score'] = int(itemsi[3]) - int(itemsi[2])
            result[i]['start'] = -1
            result[i]['end'] = -1
        for j in range(i + 1, len(pafs)):
            itemsj = pafs[j].strip().split()
            if itemsi[0].split(':')[0] != itemsj[0].split(':')[0] and itemsi[4] != itemsj[4]: continue
            regioni = itemsi[0].split(':')[1].split('-')
            regionj = itemsj[0].split(':')[1].split('-')
            if regioni[0] != regionj[1] and regioni[1] != regionj[0]: continue
            result[j]['score'] = result[i]['score'] + int(itemsj[3]) - int(itemsj[2])
            result[j]['start'] = i
            result[j]['end'] = j
            result[i]['end'] = j
    max = float('-inf')
    best = -1
    for i in range(len(result)):
        if result[i]['score'] > max:
            max = result[i]['score']
            best = i
    index = [pafs[best].strip().split()[0]]
    prev = best
    while result[prev]['start'] != -1:
        index.append(pafs[result[prev]['start']].strip().split()[0])
        prev = result[prev]['start']
    index.reverse()
    regionf = index[0].split(':')[1].split('-')
    regionl = index[-1].split(':')[1].split('-')
    new_pair, indexed = [], {}
    if int(regionf[0]) < int(regionl[0]):
        new_pair.append('{}:{}-{}'.format(index[0].split(':')[0], regionf[0], regionl[1]))
    else:
        new_pair.append('{}:{}-{}'.format(index[0].split(':')[0], regionl[0], regionf[1]))
    for p in index:
        indexed[p] = new_pair[0]
    # if new_pair[0].startswith('h2tg000001l'):
    #     print('new pair {}, indexed {}'.format(new_pair, indexed))
    return new_pair, indexed
        
    merged = []
    query = {}
    for paf in pafs:
        items = paf.strip().split() #[0].split(':')[0]
        key = items[0].split(':')[0] + items[4]
        query[key] = query.get(key, 0) + 1
    ctg, maxc = '', 0
    for q, c in query.items():
        if c > maxc:
            ctg, maxc = q, c
    for paf in pafs:
        items = paf.strip().split()
        if items[0].split(':')[0] == ctg[:-1] and items[4] == ctg[-1]:
            merged.append(items[0])

    new_pair, indexd = [], {}
    region1 = merged[0].split(':')[1].split('-')
    region2 = merged[-1].split(':')[1].split('-')
    st = min(int(region1[0]), int(region2[0]))
    en = max(int(region1[1]), int(region2[1]))
    region_new = '{}:{}-{}'.format(merged[0].split(':')[0], st, en)
    new_pair.append(region_new)
    for p in merged:
        indexd[p] = region_new
    return new_pair, indexd

def chain_dual_advance(pafname, pfname, h1fname, h2fname, prefix):
    pair = {}
    for line in open(pafname, 'r'):
        items = line.strip().split()
        pair[items[0]] = [items[1], items[2]]
        pair[items[1]] = [items[0], items[2]]
        # pair[items[0]] = items[1]
        # pair[items[1]] = items[0]

    # target = {}
    new_pair, indexd = [], {}
    for p in pfname:
        target = {}
        for line in open(p, 'r'):
            items = line.strip().split()
            target.setdefault(items[5], []).append(line.strip())

        for _, pafs in target.items():
            pafs.sort(key = lambda x: int(x.split()[7]))
            new_p, idx = merge_hap(pafs)
            new_pair.extend(new_p)
            indexd.update(idx)

    pair_f = {}
    used = set()
    for p in pair.keys():
        if p not in indexd or pair[p][0] not in indexd: continue
        if indexd[p] in used or indexd[pair[p][0]] in used: continue
        pair_f[indexd[p]] = [indexd[pair[p][0]], pair[p][1]]
        used.add(indexd[p])
        used.add(indexd[pair[p][0]])
    
    pair_list = []
    for k, v in pair_f.items():
        itemsk = k.split(':')
        regionk = itemsk[1].split('-')
        itemsv = v[0].split(':')
        regionv = itemsv[1].split('-')
        pair_list.append([itemsk[0], int(regionk[0]), int(regionk[1]), itemsv[0], int(regionv[0]), int(regionv[1]), v[1]])
    
    def __contain(ast, aen, bst, ben):
        if ast <= bst and aen >= ben:
            return True
        return False
    filter = [0] * len(pair_list)
    for i in range(len(pair_list)):
        for j in range(i + 1, len(pair_list)):
            if pair_list[j][0] == pair_list[i][0]:
                if __contain(pair_list[i][1], pair_list[i][2], pair_list[j][1], pair_list[j][2]):
                    filter[j] = 1
                elif __contain(pair_list[j][1], pair_list[j][2], pair_list[i][1], pair_list[i][2]):
                    filter[i] = 1
            if pair_list[j][3] == pair_list[i][3]:
                if __contain(pair_list[i][4], pair_list[i][5], pair_list[j][4], pair_list[j][5]):
                    filter[j] = 1
                elif __contain(pair_list[j][4], pair_list[j][5], pair_list[i][4], pair_list[i][5]):
                    filter[i] = 1
    pair_f.clear()
    tmp_list = []
    for i in range(len(pair_list)):
        if filter[i] == 0:
            tmp_list.append(pair_list[i])
    
    pair_list.clear()
    P = [-1] * len(tmp_list)
    for i in range(len(tmp_list)):
        for j in range(i + 1, len(tmp_list)):
            if tmp_list[j][0] != tmp_list[i][0] or tmp_list[j][3] != tmp_list[i][3] or tmp_list[j][6] != tmp_list[i][6]: break
            if tmp_list[j][4] < tmp_list[i][5]:
                P[i] = j
    
    used = set()
    for i in range(len(tmp_list)):
        if i in used: continue
        if P[i] == -1:
            k = '{}:{}-{}'.format(tmp_list[i][0], tmp_list[i][1], tmp_list[i][2])
            v = '{}:{}-{}'.format(tmp_list[i][3], tmp_list[i][4], tmp_list[i][5])
            pair_f[k] = [v, tmp_list[i][6]]
            used.add(i)
        else:
            used.add(i)
            j = P[i]
            while P[j] != -1:
                used.add(j)
                j = P[j]
            used.add(j)
            if tmp_list[i][6] == '+':
                k = '{}:{}-{}'.format(tmp_list[i][0], tmp_list[i][1], tmp_list[j][2])
            else:
                k = '{}:{}-{}'.format(tmp_list[i][0], tmp_list[j][1], tmp_list[i][2])
            v = '{}:{}-{}'.format(tmp_list[i][3], tmp_list[i][4], tmp_list[j][5])
            pair_f[k] = [v, tmp_list[i][6]]
    tmp_list.clear()

    with open(prefix + '.new.pair.cut.txt', 'w') as fw:
        for p in pair_f.keys():
            fw.write('{}\t{}\t{}\n'.format(p, pair_f[p][0], pair_f[p][1]))

    # get length of dual1 and dual2
    l_dual1, l_dual2 = {}, {}
    ftype = _fa_or_fq(h1fname)
    if ftype == None:
        logger.error('Unrecognized file type {}'.format(h1fname))
        exit(1)
    for head, seq in eval('_load_' + ftype)(h1fname):
        l_dual1[head] = len(seq)
    ftype = _fa_or_fq(h2fname)
    if ftype == None:
        logger.error('Unrecognized file type {}'.format(h2fname))
        exit(1)
    for head, seq in eval('_load_' + ftype)(h2fname):
        l_dual2[head] = len(seq)
    
    # get bed information of hap1 and hab2
    dual1, dual2, coll1, coll2 = {}, {}, {}, {}
    for p in pair_f.keys():
        items = p.split(':')
        region = items[1].split('-')
        dual1.setdefault(items[0], []).append([int(region[0]), int(region[1])])
        items = pair_f[p][0].split(':')
        region = items[1].split('-')
        dual2.setdefault(items[0], []).append([int(region[0]), int(region[1])])
    # get collapsed regions of dual1 and dual2
    for ctg, items in dual1.items():
        items.sort(key = lambda x: x[0])
        prev = 0
        for st, en in items:
            if st > prev:
                coll1.setdefault(ctg, []).append([prev, st])
            prev = en
        if l_dual1[ctg] > prev:
            coll1.setdefault(ctg, []).append([prev, l_dual1[ctg]])
    for ctg, items in dual2.items():
        items.sort(key = lambda x: x[0])
        prev = 0
        for st, en in items:
            if st > prev:
                coll2.setdefault(ctg, []).append([prev, st])
            prev = en
        if l_dual2[ctg] > prev:
            coll2.setdefault(ctg, []).append([prev, l_dual2[ctg]])
    # write dual1 and dual2 contigs and collapsed regions
    with open(prefix + '.new.hap1.cut.fasta', 'w') as fw, open(prefix + '.new.collapsed1.cut.fasta', 'w') as fw2:
        for head, seq in eval('_load_' + ftype)(h1fname):
            if head not in dual1:
                fw2.write('>{}:{}-{}\n'.format(head, 0, len(seq)))
                fw2.write(seq + '\n')
                coll1.setdefault(head, []).append([0, len(seq)])
                continue
            interval = dual1[head]
            interval.sort(key = lambda x: x[0])
            for st, en in interval:
                if st >= en or st < 0 or en < 0:
                    logger.error('Invalid interval {}-{} for contig {}'.format(st, en, head))
                    exit(1)
                fw.write('>{}:{}-{}\n'.format(head, st, en))
                fw.write(seq[st: en] + '\n')
            if head not in coll1: continue
            interval = coll1[head]
            interval.sort(key = lambda x: x[0])
            for st, en in interval:
                if st >= en or st < 0 or en < 0:
                    logger.error('Invalid interval {}-{} for contig {}'.format(st, en, head))
                    exit(1)
                fw2.write('>{}:{}-{}\n'.format(head, st, en))
                fw2.write(seq[st: en] + '\n')
    with open(prefix + '.new.hap2.cut.fasta', 'w') as fw, open(prefix + '.new.collapsed2.cut.fasta', 'w') as fw2:
        for head, seq in eval('_load_' + ftype)(h2fname):
            if head not in dual2:
                fw2.write('>{}:{}-{}\n'.format(head, 0, len(seq)))
                fw2.write(seq + '\n')
                coll2.setdefault(head, []).append([0, len(seq)])
                continue
            interval = dual2[head]
            interval.sort(key = lambda x: x[0])
            for st, en in interval:
                if st >= en or st < 0 or en < 0:
                    logger.error('Invalid interval {}-{} for contig {}'.format(st, en, head))
                    exit(1)
                fw.write('>{}:{}-{}\n'.format(head, st, en))
                fw.write(seq[st: en] + '\n')
            if head not in coll2: continue
            interval = coll2[head]
            interval.sort(key = lambda x: x[0])
            for st, en in interval:
                if st >= en or st < 0 or en < 0:
                    logger.error('Invalid interval {}-{} for contig {}'.format(st, en, head))
                    exit(1)
                fw2.write('>{}:{}-{}\n'.format(head, st, en))
                fw2.write(seq[st: en] + '\n')
    # write bed file
    with open(prefix + '.new.hap1.cut.bed', 'w') as fw:
        for head, interval in sorted(dual1.items(), key = lambda x: x[0]):
            interval.sort(key = lambda x: x[0])
            for st, en in interval:
                fw.write('{}\t{}\t{}\n'.format(head, st, en))
    with open(prefix + '.new.hap2.cut.bed', 'w') as fw:
        for head, interval in sorted(dual2.items(), key = lambda x: x[0]):
            interval.sort(key = lambda x: x[0])
            for st, en in interval:
                fw.write('{}\t{}\t{}\n'.format(head, st, en))
    with open(prefix + '.new.collapsed1.cut.bed', 'w') as fw:
        for head, interval in sorted(coll1.items(), key = lambda x: x[0]):
            interval.sort(key = lambda x: x[0])
            for st, en in interval:
                fw.write('{}\t{}\t{}\n'.format(head, st, en))
    with open(prefix + '.new.collapsed2.cut.bed', 'w') as fw:
        for head, interval in sorted(coll2.items(), key = lambda x: x[0]):
            interval.sort(key = lambda x: x[0])
            for st, en in interval:
                fw.write('{}\t{}\t{}\n'.format(head, st, en))


def simulate_switch_error(args):
    pfname = args.pfname
    prifname = args.prifname
    altfname = args.altfname
    prefix = args.prefix

    pri_fasta, alt_fasta = {}, {}
    ftype = _fa_or_fq(prifname)
    if ftype == None:
        logger.error('Unrecognized file type {}'.format(prifname))
        exit(1)
    for head, seq in eval('_load_' + ftype)(prifname):
        pri_fasta[head] = seq
    ftype = _fa_or_fq(altfname)
    if ftype == None:
        logger.error('Unrecognized file type {}'.format(altfname))
        exit(1)
    for head, seq in eval('_load_' + ftype)(altfname):
        alt_fasta[head] = seq

    pafs = {}
    filtered = set()
    for line in open(pfname, 'r'):
        items = line.strip().split()
        if items[0] == '000007F_007' or items[0] == '000013F_12' or items[0] == '000023F_009':
            continue
        if items[0] in pafs:
            filtered.add(items[0])
        pafs[items[0]] = line.strip()
    switch = defaultdict()
    pattern = re.compile(r'(\d+)([MIDNSHP=X])')
    for h, paf in pafs.items():
        if len(switch) == 200:
            break
        if h in filtered:
            continue
        items = paf.split()
        if int(items[3]) - int(items[2]) < 5000000:
            continue
        alt, pri = int(items[2]), 0
        alt_seq = alt_fasta[items[0]][int(items[2]): int(items[3])]
        pri_seq = pri_fasta[items[5]][int(items[7]): int(items[8])]
        if items[4] == '-':
            pri_seq = _reverse(pri_seq)
        cigar = ''
        for i in range(12, len(items)):
            if items[i].startswith('cg:Z:'):
                cigar = items[i][5:]
                break
        for m in pattern.finditer(cigar):
            l, op = int(m.group(1)), m.group(2)
            if op == '=' or op == 'M' or op == 'X':
                alt += l
                pri += l
            elif op == 'I':
                alt += l
            elif op == 'D':
                pri += l
            if alt >= len(alt_seq) / 2:
                break
        alt_subseq1, alt_subseq2 = alt_seq[:alt], alt_seq[alt:]
        pri_subseq1, pri_subseq2 = pri_seq[:pri], pri_seq[pri:]
        if items[4] == '-':
            pri_subseq1, alt_subseq2 = _reverse(pri_subseq1), _reverse(alt_subseq2)
        alt_fasta[items[0]] = alt_fasta[items[0]][:int(items[2])] + alt_subseq1 + pri_subseq2 + alt_fasta[items[0]][int(items[3]):]
        pri_fasta[items[5]] = pri_fasta[items[5]][:int(items[7])] + pri_subseq1 + alt_subseq2 + pri_fasta[items[5]][int(items[8]):]
        switch[items[0] + ':' + str(alt)] = alt
        switch[items[5] + ':' + str(pri)] = pri
    with open(prefix + '.pri.fasta', 'w') as fw:
        for head, seq in pri_fasta.items():
            fw.write('>{}\n{}\n'.format(head, seq))
    with open(prefix + '.alt.fasta', 'w') as fw:
        for head, seq in alt_fasta.items():
            fw.write('>{}\n{}\n'.format(head, seq))
    with open(prefix + '.switch', 'w') as fw:
        for h, p in switch.items():
            fw.write('{}\t{}\n'.format(h, p))


def merge(args):
    chain_dual_advance(args.pafname, args.pfname, args.h1fname, args.h2fname, args.prefix)

def chain(args):
    chain_dual(args.d1fname, args.d2fname, args.pfname, args.percent, args.prefix)

def replace(args):
    replace_mapping(args.pfname, args.threshold, args.prefix)

def fsr(args):
    filter_short_reads(args.ifname, args.ofname, args.ffname, args.length)

def dualcut(args):
    dual_cut1(args.pfname, args.threshold, args.threads, args.prefix)

def shasta(args):
    shasta2unzip(args.ifname, args.prefix, args.min_length)

def hapdup(args):
    hapdup2unzip(args.ifname, args.prefix)

def hapmap(args):
    hap_map(args.pfname1, args.pfname2, args.threads, args.fmt, args.prefix)

def lachesis(args):
    create_scaffolding_lachesis(args.fname, args.dir)

def scaffolding(args):
    create_scaffolding(args.fname, args.dir, args.result, args.prefix)

def variants(args):
    merge_var(args.vf, args.bed1, args.bed2)

def place(args):
    place_alt(args.prefix, args.paf)

def pecat(args):
    pecat_to_unzip(args.pafname, args.aafname, args.ppfname, args.apfname, args.prefix)

def breakpoint(args):
    get_breakpoints(args.var, args.ri, args.ctg, args.prefix)

def view(args):
    view_contig(args.var, args.ri, args.ctg)

def filter(args):
    filter_alt(args.file, args.output)

def cut(args):
    cut_contigs(args.pri, args.alt, args.cfile, args.prefix)

def asm2unzip(args):
    asm_unzip(args.pri, args.alt, args.olp, args.prefix)

def phase(args):
    phase_contigs(args.res, args.pri, args.coll, args.mince, args.switch, args.prefix)

def map(args):
    map_unzip(args.pri, args.alt, args.map, args.out, args.threads, args.dir, args.fmt)

def scrub(args):
    if args.combine is not None:
        primary, alt = parse_combined(args.combine)
    elif args.pri is not None and args.alt is not None:
        primary = parse_separate(args.pri)
        alt = parse_separate(args.alt)
    else:
        logger.error('Please specify either a combined contig file (-c | --combine) or both p and a contig files (-p | --pri and -a | --alt)')
        exit(1)
    omitted = make_name_mapping(primary.keys(), alt.keys())
    with open(args.prefix + '.p_tigs.fasta', 'w') as fw:
        for head, seq in primary.items():
            fw.write('>{}\n'.format(head))
            fw.write(seq + '\n')
    with open(args.prefix + '.a_tigs.fasta', 'w') as fw:
        for head, seq in alt.items():
            if head not in omitted:
                fw.write('>{}\n'.format(head))
                fw.write(seq + '\n')

def clean_collesped(bfname1, bfname2, hfname1, hfname2, prefix):
    reads = {}
    for fname in [hfname1, hfname2]:
        ftype = _fa_or_fq(fname)
        if ftype == None:
            logger.error('Unrecognized file type {}'.format(fname))
            exit(1)
        for head, seq in eval('_load_' + ftype)(fname):
            reads[head] = seq
    bed1, bed2 = {}, {}
    prev = ''
    for line in open(bfname1):
        items = line.strip().split()
        name = items[0] + ':' + items[1] + '-' + items[2]
        if prev == '' or items[3] != prev:
            bed1.setdefault(items[3], [])
            prev = items[3]
        if name in reads:
            bed1[items[3]].append(name)
    prev = ''
    for line in open(bfname2):
        items = line.strip().split()
        name = items[0] + ':' + items[1] + '-' + items[2]
        if prev == '' or items[3] != prev:
            bed2.setdefault(items[3], [])
            prev = items[3]
        if name in reads:
            bed2[items[3]].append(name)
    with open(prefix + '.phased0.fasta', 'w') as fw:
        for name, ctgs in bed1.items():
            fw.write('>{}\n'.format(name))
            for c in ctgs:
                fw.write(reads[c])
            fw.write('\n')
    with open(prefix + '.phased1.fasta', 'w') as fw:
        for name, ctgs in bed2.items():
            fw.write('>{}\n'.format(name))
            for c in ctgs:
                fw.write(reads[c])
            fw.write('\n')
    with open(prefix + '.phased0.bed', 'w') as fw:
        for name, ctgs in bed1.items():
            for c in ctgs:
                fw.write('{}\t{}\n'.format(c, name))
    with open(prefix + '.phased1.bed', 'w') as fw:
        for name, ctgs in bed2.items():
            for c in ctgs:
                fw.write('{}\t{}\n'.format(c, name))

def clean(args):
    clean_collesped(args.bed1, args.bed2, args.hap1, args.hap2, args.prefix)

def cluster_ctgs(ffname, pfname, ofname):
    mapping = {}
    ftype = _fa_or_fq(ffname)
    if ftype == None:
        logger.error('Unrecognized file type {}'.format(ffname))
        exit(1)
    for head, seq in eval('_load_full_' + ftype)(ffname):
        r = re.search(r'.* chromosome ([0-9XY]+)', head)
        if r == None: continue
        items = head.split()
        mapping[items[0]] = r.group(1)

    cluster = {}
    prev = ''
    tmp = {}
    for line in open(pfname):
        items = line.strip().split()
        if items[0].endswith('1'): continue
        if prev == '' or items[0] != prev:
            if prev != '':
                ctg = sorted(tmp.items(), key = lambda x: x[1], reverse = True)[0][0]
                if ctg in mapping:
                    cluster.setdefault(mapping[ctg], []).append(prev)
                tmp.clear()
            prev = items[0]
        tmp[items[5]] = tmp.get(items[5], 0) + int(items[9])
    ctg = sorted(tmp.items(), key = lambda x: x[1], reverse = True)[0][0]
    if ctg in mapping:
        cluster.setdefault(mapping[ctg], []).append(prev)
    
    count = len(cluster)
    for c, _ in cluster.items():
        if not c.isdigit():
            count -= 1
    with open(ofname, 'w') as fw:
        for c, ctgs in sorted(cluster.items(), key = lambda x: x[0]):
            fw.write('\t'.join(ctgs) + '\n')
            index = c
            if not c.isdigit():
                index = count + 1
                count += 1
            fname = 'group{}.ordering'.format(int(index) - 1)
            out = open(fname, 'w')
            for ctg in ctgs:
                out.write('1\t{}\t1\t1\t.\n'.format(ctg))
            out.close()

def cluster(args):
    cluster_ctgs(args.ffname, args.pfname, args.ofname)

def main():
    parser = argparse.ArgumentParser('preprocessing contigs')
    subparser = parser.add_subparsers(title = 'subcommand', description = 'breakpoing\nview\n')

    parser_bp = subparser.add_parser('breakpoint', help = 'get breakpoint from variants')
    parser_bp.add_argument('-v', '--var', dest = 'var', default = None, metavar = 'path', help = 'variants file name')
    parser_bp.add_argument('-r', '--ri', dest = 'ri', default = None, metavar = 'path', help = 'readinfo file name')
    parser_bp.add_argument('-c', '--ctg', dest = 'ctg', default = None, metavar = 'str', help = 'contig file name')
    parser_bp.add_argument('-p', '--prefix', dest = 'prefix', default = 'phasing', metavar = 'str', help = 'prefix of output file')
    parser_bp.set_defaults(func = breakpoint)

    parser_view = subparser.add_parser('view', help = 'view variants in one region')
    parser_view.add_argument('-v', '--var', dest = 'var', default = None, metavar = 'path', help = 'variants file name')
    parser_view.add_argument('-r', '--ri', dest = 'ri', default = None, metavar = 'path', help = 'readinfo file name')
    parser_view.add_argument('-c', '--ctg', dest = 'ctg', default = None, metavar = 'str', help = 'region [chr:beg-end]')
    parser_view.set_defaults(func = view)

    parser_filter = subparser.add_parser('filter', help = 'filter nested and dumplicated alt contigs')
    parser_filter.add_argument('-f', '--file', dest = 'file', default = None, metavar = 'path', help = 'alt contigs placemnet file name')
    parser_filter.add_argument('-o', '--output', dest = 'output', default = None, metavar = 'path', help = 'output file name')
    parser_filter.set_defaults(func = filter)

    parser_cut = subparser.add_parser('cut', help = 'cut primary contigs')
    parser_cut.add_argument('-P', '--pri', dest = 'pri', default = None, metavar = 'path', help = 'primary contig file name')
    parser_cut.add_argument('-A', '--alt', dest = 'alt', default = None, metavar = 'path', help = 'alt contig file name')
    parser_cut.add_argument('-c', '--cfile', dest = 'cfile', default = None, metavar = 'str', help = 'cut file name')
    parser_cut.add_argument('-p', '--prefix', dest = 'prefix', default = 'phasing', metavar = 'str', help = 'prefix of output file')
    parser_cut.set_defaults(func = cut)

    parser_a2u = subparser.add_parser('asm2unzip', help = 'convert assembly to unzip format')
    parser_a2u.add_argument('-p', '--pri', dest = 'pri', default = None, metavar = 'path', help = 'primary assembly file name')
    parser_a2u.add_argument('-a', '--alt', dest = 'alt', default = None, metavar = 'path', help = 'alt assembly file name')
    parser_a2u.add_argument('-o', '--olp', dest = 'olp', default = None, metavar = 'path', help = 'overlap file name in paf format')
    parser_a2u.add_argument('-P', '--prefix', dest = 'prefix', default = 'phasing', metavar = 'str', help = 'prefix of output files')
    parser_a2u.set_defaults(func = asm2unzip)

    parser_phase = subparser.add_parser('phase', help = 'phase contigs according phased result')
    parser_phase.add_argument('-r', '--res', dest = 'res', default = None, metavar = 'path', help = 'phased result file name')
    parser_phase.add_argument('-P', '--pri', dest = 'pri', default = None, metavar = 'path', help = 'bed file name of haplotig on primary contigs')
    parser_phase.add_argument('-c', '--col', dest = 'coll', default = None, metavar = 'path', help = 'bed file name of collapsed contigs')
    parser_phase.add_argument('-m', '--mince', dest = 'mince', default = None, metavar = 'path', help = 'minced contigs file name')
    parser_phase.add_argument('-s', '--switch', dest = 'switch', default = None, metavar = 'path', help = 'switch error file name')
    parser_phase.add_argument('-p', '--prefix', dest = 'prefix', default = 'phasing', metavar = 'path', help = 'prefix of output files')
    parser_phase.set_defaults(func = phase)
    
    parser_map = subparser.add_parser('map', help = 'map falcon-unzip')
    parser_map.add_argument('-p', '--pri', dest = 'pri', default = None, metavar = 'path', help = 'primary contig file name')
    parser_map.add_argument('-a', '--alt', dest = 'alt', default = None, metavar = 'path', help = 'alt contig file name')
    parser_map.add_argument('-m', '--map', dest = 'map', default = None, metavar = 'path', help = 'name mapping file name')
    parser_map.add_argument('-o', '--out', dest = 'out', default = None, metavar = 'path', help = 'output file name')
    parser_map.add_argument('-d', '--dir', dest = 'dir', default = None, metavar = 'path', help = 'temporary directory')
    parser_map.add_argument('-t', '--threads', dest = 'threads', default = 1, type = int, metavar = 'int', help = 'number of threads')
    parser_map.add_argument('-f', '--fmt', dest = 'fmt', default = 'paf', choices = ['paf', 'sam'], type = str.lower, metavar = 'str', help = 'sam or paf format default paf')
    parser_map.set_defaults(func = map)

    parser_scrub = subparser.add_parser('scrub', help = 'scrub the contig names of Falcon-Unzip assembly')
    parser_scrub.add_argument('-p', '--pri', dest = 'pri', required = False, default = None, metavar = 'path', help = 'primary contig file name')
    parser_scrub.add_argument('-a', '--alt', dest = 'alt', required = False, default = None, metavar = 'path', help = 'alt contig file name')
    parser_scrub.add_argument('-c', '--combine', dest = 'combine', required = False, default = None, metavar = 'path', help = 'combined contig file name')
    parser_scrub.add_argument('-o', '--prefix', dest = 'prefix', default = 'phasing', metavar = 'str', help = 'prefix of output files')
    parser_scrub.set_defaults(func = scrub)

    parser_pecat = subparser.add_parser('pecat', help = 'convert the format of pecat to unzip')
    parser_pecat.add_argument('-pa', '--pafname', dest = 'pafname', default = None, metavar = 'path', help = 'primary assembly file name')
    parser_pecat.add_argument('-aa', '--aafname', dest = 'aafname', default = None, metavar = 'path', help = 'alt assembly file name')
    parser_pecat.add_argument('-pp', '--ppfname', dest = 'ppfname', default = None, metavar = 'path', help = 'primary polish file name')
    parser_pecat.add_argument('-ap', '--apfname', dest = 'apfname', default = None, metavar = 'path', help = 'alt polish file name')
    parser_pecat.add_argument('-p', '--prefix', dest = 'prefix', default = 'phasing', metavar = 'str', help = 'prefix of the output files')
    parser_pecat.set_defaults(func = pecat)

    parser_place = subparser.add_parser('place', help = 'place alt contigs along primary contigs according paf file')
    parser_place.add_argument('-f', '--paf', dest = 'paf', default = None, metavar = 'path', help = 'paf file name')
    parser_place.add_argument('-p', '--prefix', dest = 'prefix', default = 'phasing', metavar = 'str', help = 'prefix of output files')
    parser_place.set_defaults(func = place)

    parser_var = subparser.add_parser('variant', help = 'generate variants from first round phasing for second round phasing')
    parser_var.add_argument('-v', '--var', dest = 'vf', default = None, metavar = 'path', help = 'variants generated by first round phasing')
    parser_var.add_argument('-b1', '--bed1', dest = 'bed1', default = None, metavar = 'path', help = '*.phased0.bed')
    parser_var.add_argument('-b2', '--bed2', dest = 'bed2', default = None, metavar = 'path', help = '*.phased1.bed')
    parser_var.set_defaults(func = variants)

    parser_lachesis = subparser.add_parser('lachesis', help = 'generate scaffolds')
    parser_lachesis.add_argument('-f', '--file', dest = 'fname', default = None, metavar = 'path', help = 'fasta file')
    parser_lachesis.add_argument('-d', '--dir', dest = 'dir', default = None, metavar = 'path', help = 'output directory of lachesis')
    parser_lachesis.set_defaults(func = lachesis)

    parser_scaffold = subparser.add_parser('scaffold', help = 'generate scaffolds')
    parser_scaffold.add_argument('-f', '--file', dest = 'fname', default = None, metavar = 'path', help = 'fasta file')
    parser_scaffold.add_argument('-d', '--dir', dest = 'dir', default = None, metavar = 'path', help = 'output directory of lachesis')
    parser_scaffold.add_argument('-r', '--result', dest = 'result', default = None, metavar = 'path', help = 'phasing result file name')
    parser_scaffold.add_argument('-p', '--prefix', dest = 'prefix', default = 'phasing', metavar = 'str', help = 'prefix of output files')
    parser_scaffold.set_defaults(func = scaffolding)

    parser_clean = subparser.add_parser('clean')
    parser_clean.add_argument('-b1', '--bed1', dest = 'bed1', default = None, metavar = 'path')
    parser_clean.add_argument('-b2', '--bed2', dest = 'bed2', default = None, metavar = 'path')
    parser_clean.add_argument('-h1', '--hap1', dest = 'hap1', default = None, metavar = 'path')
    parser_clean.add_argument('-h2', '--hap2', dest = 'hap2', default = None, metavar = 'path')
    parser_clean.add_argument('-p', '--prefix', dest = 'prefix', default = 'phasing', metavar = 'str')
    parser_clean.set_defaults(func = clean)

    parser_cluster = subparser.add_parser('cluster')
    parser_cluster.add_argument('-f', '--ffname', dest = 'ffname', default = None, metavar = 'path', help = 'ctg file name')
    parser_cluster.add_argument('-p', '--pfname', dest = 'pfname', default = None, metavar = 'path', help = 'mapping file name, paf format')
    parser_cluster.add_argument('-o', '--ofname', dest = 'ofname', default = None, metavar = 'path', help = 'output file name')
    parser_cluster.set_defaults(func = cluster)

    parser_shasta = subparser.add_parser('shasta')
    parser_shasta.add_argument('-i', '--ifname', dest = 'ifname', default = None, metavar = 'path', help = 'shasta assembly file name')
    parser_shasta.add_argument('-p', '--prefix', dest = 'prefix', default = 'shasta', metavar = 'str', help = 'prefix of output files')
    parser_shasta.add_argument('-l', '--min_length', dest = 'min_length', default = 3000, type = int, metavar = 'int', help = 'minimum length of contigs')
    parser_shasta.set_defaults(func = shasta)

    parser_hapdup = subparser.add_parser('hapdup')
    parser_hapdup.add_argument('-i', '--ifname', dest = 'ifname', default = None, metavar = 'path', help = 'hapdup assembly file name')
    parser_hapdup.add_argument('-p', '--prefix', dest = 'prefix', default = 'hapdup', metavar = 'str', help = 'prefix of output files')
    parser_hapdup.set_defaults(func = hapdup)

    parser_hapmap = subparser.add_parser('hapmap')
    parser_hapmap.add_argument('-1', '--pfname1', dest = 'pfname1', default = None, metavar = 'path', help = 'hapdup assembly file name')
    parser_hapmap.add_argument('-2', '--pfname2', dest = 'pfname2', default = None, metavar = 'path', help = 'hapdup assembly file name')
    parser_hapmap.add_argument('-t', '--threads', dest = 'threads', default = 1, type = int, metavar = 'int', help = 'number of threads')
    parser_hapmap.add_argument('-f', '--fmt', dest = 'fmt', default = 'paf', choices = ['paf', 'sam'], type = str.lower, metavar = 'str', help = 'sam or paf format default paf')
    parser_hapmap.add_argument('-p', '--prefix', dest = 'prefix', default = 'hapmap', metavar = 'str', help = 'prefix of output files')
    parser_hapmap.set_defaults(func = hapmap)

    parser_dualcut = subparser.add_parser('dualcut')
    parser_dualcut.add_argument('-i', '--pfname', dest = 'pfname', default = None, metavar = 'path', help = 'hapdup assembly mapping file name')
    parser_dualcut.add_argument('-s', '--threshold', dest = 'threshold', default = 10000, type = int, metavar = 'int', help = 'threshold for cutting')
    parser_dualcut.add_argument('-p', '--prefix', dest = 'prefix', default = 'hapdup', metavar = 'str', help = 'prefix of output files')
    parser_dualcut.add_argument('-t', '--threads', dest = 'threads', default = 1, type = int, metavar = 'int', help = 'number of threads')
    parser_dualcut.set_defaults(func = dualcut)

    parser_fsr = subparser.add_parser('fsr')
    parser_fsr.add_argument('-i', '--ifname', dest = 'ifname', default = None, metavar = 'path', help = 'hapdup assembly mapping file name')
    parser_fsr.add_argument('-l', '--length', dest = 'length', default = 20000, type = int, metavar = 'int', help = 'threshold for cutting')
    parser_fsr.add_argument('-f', '--ffname', dest = 'ffname', default = None, metavar = 'path', help = 'hapdup assembly mapping file name')
    parser_fsr.add_argument('-o', '--ofname', dest = 'ofname', default = None, metavar = 'path', help = 'output file name')
    parser_fsr.set_defaults(func = fsr)

    parser_replace = subparser.add_parser('replace')
    parser_replace.add_argument('-i', '--pfname', dest = 'pfname', default = None, metavar = 'path', help = 'hapdup assembly mapping file name')
    parser_replace.add_argument('-s', '--threshold', dest = 'threshold', default = 10000, type = int, metavar = 'int', help = 'threshold for cutting')
    parser_replace.add_argument('-p', '--prefix', dest = 'prefix', default = 'phasing', metavar = 'str', help = 'prefix of output files')
    parser_replace.set_defaults(func = replace)

    parser_chain = subparser.add_parser('chain')
    parser_chain.add_argument('-1', '--d1fname', dest = 'd1fname', default = None, metavar = 'path', help = 'hapdup assembly mapping file name')
    parser_chain.add_argument('-2', '--d2fname', dest = 'd2fname', default = None, metavar = 'path', help = 'hapdup assembly mapping file name')
    parser_chain.add_argument('-i', '--pfname', dest = 'pfname', default = None, metavar = 'path', help = 'hapdup assembly mapping file name')
    parser_chain.add_argument('-s', '--percent', dest = 'percent', default = 80, type = float, metavar = 'float', help = 'percentile for filtering')
    parser_chain.add_argument('-p', '--prefix', dest = 'prefix', default = 'phasing', metavar = 'str', help = 'prefix of output files')
    parser_chain.set_defaults(func = chain)

    parser_merge = subparser.add_parser('merge')
    parser_merge.add_argument('-i', '--pfname', dest = 'pfname', nargs = '+', default = None, metavar = 'path', help = 'mapping file name')
    parser_merge.add_argument('-s', '--pafname', dest = 'pafname', default = None, metavar = 'path', help = 'pair file name')
    parser_merge.add_argument('-h1', '--hap1', dest = 'h1fname', default = None, metavar = 'path', help = 'dual1 file name')
    parser_merge.add_argument('-h2', '--hap2', dest = 'h2fname', default = None, metavar = 'path', help = 'dual2 file name')
    parser_merge.add_argument('-p', '--prefix', dest = 'prefix', default = 'phasing', metavar = 'str', help = 'prefix of output files')
    parser_merge.set_defaults(func = merge)

    parser_fix = subparser.add_parser('fix_phase')
    parser_fix.add_argument('-s', '--sfname', dest = 'sfname', default = None, metavar = 'path', help = 'switch file name')
    parser_fix.add_argument('-a', '--afname', dest = 'afname', default = None, metavar = 'path', help = 'pair file name')
    parser_fix.add_argument('-p', '--pfname', dest = 'pfname', default = None, metavar = 'path', help = 'paf file name')
    parser_fix.add_argument('-o', '--ofname', dest = 'ofname', default = None, metavar = 'path', help = 'output file name')
    parser_fix.set_defaults(func = fix_switch)

    parser_group = subparser.add_parser('group')
    parser_group.add_argument('--h1', dest = 'h1fname', default = None, metavar = 'path', help = 'hap1 file name')
    parser_group.add_argument('--h2', dest = 'h2fname', default = None, metavar = 'path', help = 'hap2 file name')
    parser_group.add_argument('--group', dest = 'gfname', default = None, metavar = 'path', help = 'group file name')
    parser_group.add_argument('--o1', dest = 'o1fname', default = None, metavar = 'path', help = 'output hap1 file name')
    parser_group.add_argument('--o2', dest = 'o2fname', default = None, metavar = 'path', help = 'output hap2 file name')
    parser_group.set_defaults(func = group)

    parser_sim = subparser.add_parser('sim')
    parser_sim.add_argument('--paf', dest = 'pfname', default = None, metavar = 'path', help = 'paf filename')
    parser_sim.add_argument('--pri', dest = 'prifname', default = None, metavar = 'path', help = 'primary filename')
    parser_sim.add_argument('--alt', dest = 'altfname', default = None, metavar = 'path', help = 'alternate filename')
    parser_sim.add_argument('--prefix', dest = 'prefix', default = None, metavar = 'str', help = 'prefix of output files')
    parser_sim.set_defaults(func = simulate_switch_error)

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()