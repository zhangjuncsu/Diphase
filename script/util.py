import io
import gzip
import time

def _time():
    return time.asctime(time.localtime(time.time()))

def _fa_or_fq(fname):
    if fname.endswith('.gz'):
        fname = fname.replace('.gz', '')
    if fname.endswith('.fasta') or fname.endswith('.fa') or fname.endswith('.fna'):
        return 'fa'
    elif fname.endswith('fastq') or fname.endswith('fq'):
        return 'fq'
    return None

def _load_fq(fname):
    # yield one record from fastq file
    head, seq = '', ''
    if fname.endswith('.gz'):
        handle = io.TextIOWrapper(gzip.open(fname, 'r'))
    else:
        handle = open(fname, 'r')
    num = 0
    for line in handle:
        num += 1
        if num % 4 == 1 and line.startswith('@'):
            head = line[1:].split()[0]
        elif num % 4 == 2:
            seq = line.strip()
            yield head, seq

def _load_full_fq(fname):
    # yield one record from fastq file
    head, seq = '', ''
    if fname.endswith('.gz'):
        handle = io.TextIOWrapper(gzip.open(fname, 'r'))
    else:
        handle = open(fname, 'r')
    num = 0
    for line in handle:
        num += 1
        if num % 4 == 1 and line.startswith('@'):
            head = line[1:].strip()
        elif num % 4 == 2:
            seq = line.strip()
            yield head, seq

def _load_fa(fname):
    # yield one record from fasta file
    head, seq = '', ''
    if fname.endswith('.gz'):
        handle = io.TextIOWrapper(gzip.open(fname, 'r'))
    else:
        handle = open(fname, 'r')
    for line in handle:
        if line.startswith('>'):
            if head != '':
                yield head, seq
            head = line[1:].split()[0]
            seq = ''
        else:
            seq += line.strip()
    yield head, seq

def _load_full_fa(fname):
    # yield one record from fasta file
    head, seq = '', ''
    if fname.endswith('.gz'):
        handle = io.TextIOWrapper(gzip.open(fname, 'r'))
    else:
        handle = open(fname, 'r')
    for line in handle:
        if line.startswith('>'):
            if head != '':
                yield head, seq
            head = line[1:].strip()
            seq = ''
        else:
            seq += line.strip()
    yield head, seq

def _load_snp_matrix(mfname):
    header, matrix = [], []
    if mfname.endswith('.gz'):
        handle = io.TextIOWrapper(gzip.open(mfname, 'r'))
    else:
        handle = open(mfname, 'r')
    for line in handle:
        if len(line.strip()) == 0: continue
        if len(header) != 0:
            matrix.append(line.strip().split())
            if len(matrix[-1]) != len(header):
                raise ValueError('inconsistent matrix size')
            if len(matrix) == len(header):
                yield header, matrix
                header, matrix = [], []
        else:
            header = line.strip().split()
            matrix = []
    # yield header, matrix