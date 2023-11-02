#!/usr/bin/env python3

import re
import sys
import argparse
from collections import OrderedDict, defaultdict

def falconphase_ov_index(ffname, pfname, cfname):
    fai, index = {}, 0
    for line in open(ffname):
        items = line.strip().split()
        fai[items[0]] = index
        index += 1

    cluster, N_cluster = OrderedDict(), 0
    for line in open(cfname):
        if line.startswith('#'): continue
        items = line.strip().split()
        cluster['cluster' + str(N_cluster)] = items
        N_cluster += 1

    pair = {}
    for line in open(pfname):
        items = line.strip().split()
        pair[items[0]] = items[1]

    for c, items in cluster.items():
        line = c + '\t'
        for item in items:
            line += '{},'.format(fai[item])
        for item in items:
            line += '{},'.format(fai[pair[item]])
        line = line[:-1]
        line += '\t'
        for item in items:
            line += '{}:{},'.format(fai[item], fai[pair[item]])
        line = line[:-1]
        print(line, file = sys.stdout)

def falconphase(args):
    falconphase_ov_index(args.fai, args.pair, args.cluster)

def main():
    parser = argparse.ArgumentParser()
    subparser = parser.add_subparsers(title = 'subcommand', description = 'classify')
    
    parser_falcon = subparser.add_parser('falconphase', help = 'get ov_index.txt for falconphase')
    parser_falcon.add_argument('-f', '--fai', dest = 'fai', default = None, metavar = 'path', help = 'fai file')
    parser_falcon.add_argument('-p', '--pair', dest = 'pair', default = None, metavar = 'path', help = 'pair file')
    parser_falcon.add_argument('-c', '--cluster', dest = 'cluster', default = None, metavar = 'path', help = 'cluster file')
    parser_falcon.set_defaults(func = falconphase)
    
    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()