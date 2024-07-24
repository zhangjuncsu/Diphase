#!/usr/bin/env python3

import os
import sys
import uuid
import time
import logging
import resource
import argparse
import subprocess

from preprocessing import fsa_rename, fsa_unzip, fsa_unzip1, parse_separate, make_name_mapping, map_unzip
from preprocessing import place_alt, filter_alt, cut_contigs, phase_contigs, setLogger, create_scaffolding
from preprocessing import shasta2unzip_and_mapping, dual_cut, phase_contigs_dual, fix_switch2, group2

logger = logging.getLogger('phasing pipeline')

def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program): return program
    else:
        for path in os.environ['PATH'].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file): return exe_file
    return None

def _check_program(*programs):
    for program in programs:
        prog = which(program)
        if prog is None:
            logger.error('Program {} not found!'.format(program))
            exit(1)
        else:
            logger.info('Program {} found in {}'.format(program, prog))

def _file_exist(*files):
    for f in files:
        if not os.access(f, os.R_OK):
            logger.error('Could not open {} for reading!'.format(f), exc_info = 1)
            exit(1)

def _test_file(dir, *files):
    for file in files:
        if not os.access(os.path.join(dir, file), os.R_OK):
            return False
    logger.info('File(s) {} already in {}'.format(' '.join(files), dir))
    return True

def _enable_logging(log_fname):
    # enable logging

    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('[%(asctime)s] [%(levelname)s] ==> %(message)s', '%Y-%m-%d %H:%M:%S')
    if log_fname is not None:
        fh = logging.FileHandler(log_fname)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

def pecat_to_unzip_format(pri_asm_fname, alt_asm_fname, pri_pol_fname, alt_pol_fname, prefix, dir):
    # pri_asm_fname     primary assembly filename
    # alt_asm_fname     alt assembly filename
    # pri_pol_fname     primary polish filename
    # alt_pol_fname     alt polish filename
    # prefix            prefix of output files
    # dir               working directory

    start_time = time.time()

    ## final files are ${prefix}.primary.clean.fasta and ${prefix}.alt.clean.fasta
    pri_clean_fname = prefix + '.primary.clean.fasta'
    alt_clean_fname = prefix + '.alt.clean.fasta'

    if _test_file(dir, pri_clean_fname, alt_clean_fname): return

    if pri_asm_fname is not None or alt_asm_fname is not None:
        _file_exist(pri_asm_fname, alt_asm_fname)
    _file_exist(pri_pol_fname, alt_pol_fname)

    _check_program('awk')

    try:
        os.chdir(dir)
    except FileNotFoundError:
        logger.error('Directory {} not exists!'.format(dir), exc_info = 1)
        exit(1)
    except PermissionError:
        logger.error('No permission to change to directory {}'.format(dir), exc_info = 1)
        exit(1)

    logger.info('{} Working directory {}'.format(sys._getframe().f_code.co_name, os.getcwd()))
    
    if pri_asm_fname is not None and alt_asm_fname is not None:
        # convert the name of polish sequences to the format of assembly
        pri_rename_fname = prefix + '.primary.rename.fasta'
        fsa_rename(pri_asm_fname, pri_pol_fname, pri_rename_fname)
        alt_rename_fname = prefix + '.alt.rename.fasta'
        fsa_rename(alt_asm_fname, alt_pol_fname, alt_rename_fname)

        # convert the name of PECAT to unzip
        pri_unzip_fname = prefix + '.primary.unzip.fasta'
        fsa_unzip(pri_rename_fname, pri_unzip_fname)
        alt_unzip_fname = prefix + '.alt.unzip.fasta'
        fsa_unzip(alt_rename_fname, alt_unzip_fname)
    else:
        pri_unzip_fname = prefix + '.primary.unzip.fasta'
        fsa_unzip1(pri_pol_fname, pri_unzip_fname)
        alt_unzip_fname = prefix + '.alt.unzip.fasta'
        fsa_unzip1(alt_pol_fname, alt_unzip_fname)

    # clean the name of the converted
    # pri_clean_fname = prefix + '.primary.clean.fasta'
    cmd_pri = 'awk \'{ if(NR%2==1) tmp=$1; if(NR%2==0 && !match(tmp,"_")) print tmp"\\n"$0 }\' '
    # cmd_pri += pri_unzip_fname + ' > ' + pri_clean_fname
    logger.info('Command for primary: {}'.format(cmd_pri))
    try:
        subprocess.call(cmd_pri, shell = True, 
                        stdin = open(pri_unzip_fname, 'r'), 
                        stdout = open(pri_clean_fname, 'w'), 
                        stderr = sys.stderr)
    except (subprocess.CalledProcessError, OSError) as e:
        logger.error('Error running awk for primary {}'.format(e.cmd), exc_info = 1)
        exit(1)
    # alt_clean_fname = prefix + '.alt.clean.fasta'
    cmd_alt = 'awk \'{ if(NR%2==1) tmp=$1; if(NR%2==0 && match(tmp,"_")) print tmp"\\n"$0 }\' '
    # cmd_alt += alt_unzip_fname + ' > ' + alt_clean_fname
    logger.info('Command for alt: {}'.format(cmd_alt))
    try:
        subprocess.call(cmd_alt, shell = True, 
                        stdin = open(alt_unzip_fname, 'r'), 
                        stdout = open(alt_clean_fname, 'w'), 
                        stderr = sys.stderr)
    except (subprocess.CalledProcessError, OSError) as e:
        logger.error('Error running awk for alt {}'.format(e.cmd), exc_info = 1)
        exit(1)
    
    end_time = time.time()
    logger.info('{} ends'.format(sys._getframe().f_code.co_name))
    logger.info('Wall clock time elapsed: {} s'.format(end_time - start_time))
    logger.info('CPU time elapsed: {} s'.format(time.process_time()))
    logger.info('Memory usage: {:>.3f} GB\n'.format(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024 / 1024))
    
def clean_unzip(pri_fname, alt_fname, prefix, dir):
    # pri_fname     primary filename
    # alt_fname     alt filename
    # prefix        prefix of output files
    # dir           working directory

    start_time = time.time()

    ## final files are ${prefix}.p_tigs.fasta and ${prefix}.a_tigs.fasta
    clean_pri_fname = prefix + '.p_tigs.fasta'
    clean_alt_fname = prefix + '.a_tigs.fasta'

    if _test_file(dir, clean_pri_fname, clean_alt_fname): return

    _file_exist(pri_fname, alt_fname)

    try:
        os.chdir(dir)
    except FileNotFoundError:
        logger.error('Directory {} not exists!'.format(dir), exc_info = 1)
        exit(1)
    except PermissionError:
        logger.error('No permission to change to directory {}'.format(dir), exc_info = 1)
        exit(1)
    
    logger.info('{} Working directory {}'.format(sys._getframe().f_code.co_name, os.getcwd()))

    primary = parse_separate(pri_fname)
    alt = parse_separate(alt_fname)
    omitted = make_name_mapping(primary.keys(), alt.keys())
    
    with open(clean_pri_fname, 'w') as fw:
        for head, seq in primary.items():
            fw.write('>{}\n'.format(head))
            fw.write(seq + '\n')
    with open(clean_alt_fname, 'w') as fw:
        for head, seq in alt.items():
            if head not in omitted:
                fw.write('>{}\n'.format(head))
                fw.write(seq + '\n')

    _check_program('samtools')
    try:
        subprocess.call('samtools faidx ' + clean_pri_fname, shell = True, stderr = sys.stderr)
    except (subprocess.CalledProcessError, OSError) as e:
        logger.error('Error when running samtools faidx for primary {}'.format(e.cmd), exc_info = 1)
        exit(1)
    try:
        subprocess.call('samtools faidx ' + clean_alt_fname, shell = True, stderr = sys.stderr)
    except (subprocess.CalledProcessError, OSError) as e:
        logger.error('Error when running samtools faidx for alt {}'.format(e.cmd), exc_info = 1)
        exit(1)

    end_time = time.time()
    logger.info('{} ends'.format(sys._getframe().f_code.co_name))
    logger.info('Wall clock time elapsed: {} s'.format(end_time - start_time))
    logger.info('CPU time elapsed: {} s'.format(time.process_time()))
    logger.info('Memory usage: {:>.3f} GB\n'.format(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024 / 1024))

def make_mapping_alt_to_primary(pri_fname, alt_fname, name_mapping, output, threads, fmt, dir):
    # pri_fname     primary filename
    # alt_fname     alt filename
    # name_mapping  name mapping filename
    # output        output file name
    # threads       number of threads
    # fmt           output format paf or sam
    # dir           working directory

    start_time = time.time()

    ## final file are ${output}.${fmt}
    check_fname = output + '.' + fmt
    if _test_file(dir, check_fname): return

    _file_exist(pri_fname, alt_fname, name_mapping)

    _check_program('minimap2')

    try:
        os.chdir(dir)
    except FileNotFoundError:
        logger.error('Directory {} not exists!'.format(dir), exc_info = 1)
        exit(1)
    except PermissionError:
        logger.error('No permission to change to directory {}'.format(dir), exc_info = 1)
        exit(1)

    logger.info('{} Working directory {}'.format(sys._getframe().f_code.co_name, os.getcwd()))

    if os.path.isfile(output):
        logger.warning('{} has exist, will be overwritten!'.format(output))
    if threads < 1:
        logger.info('threads {} must be >= 1, set to 1.'.format(threads))
        threads = 1
    if fmt.lower() != 'paf':
        logger.info('Format {} are not paf, set to paf.'.format(fmt))
        fmt = 'paf'

    tmp_dir = str(uuid.uuid4())
    while os.path.exists(tmp_dir):
        tmp_dir = str(uuid.uuid4())

    map_unzip(pri_fname, alt_fname, name_mapping, output, threads, tmp_dir, fmt)

    end_time = time.time()
    logger.info('{} ends'.format(sys._getframe().f_code.co_name))
    logger.info('Wall clock time elapsed: {} s'.format(end_time - start_time))
    logger.info('CPU time elapsed: {} s'.format(time.process_time()))
    logger.info('Memory usage: {:>.3f} GB\n'.format(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024 / 1024))

def place_alt_along_primary_and_filter_contained(paf_fname, prefix, dir):
    # paf_fname     mapping filename paf format
    # prefix        prefix of output files
    # dir           working directory

    start_time = time.time()

    ## final files are ${prefix}.unfiltered.alt2pri.txt ${prefix}.filtered.alt2pri.paf ${prefix}.filtered.alt2pri.txt
    unfiltered = prefix + '.unfiltered.alt2pri.txt'
    filtered_olp = prefix + '.filtered.alt2pri.paf'
    filtered = prefix + '.filtered.alt2pri.txt'
    if _test_file(dir, filtered_olp, filtered): return
    elif _test_file(dir, unfiltered, filtered_olp): 
        filter_alt(unfiltered, filtered)
        return

    _file_exist(paf_fname)
    
    try:
        os.chdir(dir)
    except FileNotFoundError:
        logger.error('Directory {} not exists!'.format(dir), exc_info = 1)
        exit(1)
    except PermissionError:
        logger.error('No permission to change to directory {}'.format(dir), exc_info = 1)
        exit(1)

    logger.info('{} Working directory {}'.format(sys._getframe().f_code.co_name, os.getcwd()))

    place_alt(paf_fname, prefix)
    
    unfiltered = prefix + '.unfiltered.alt2pri.txt'
    filtered = prefix + '.filtered.alt2pri.txt'
    _file_exist(unfiltered)
    filter_alt(unfiltered, filtered)

    end_time = time.time()
    logger.info('{} ends'.format(sys._getframe().f_code.co_name))
    logger.info('Wall clock time elapsed: {} s'.format(end_time - start_time))
    logger.info('CPU time elapsed: {} s'.format(time.process_time()))
    logger.info('Memory usage: {:>.3f} GB\n'.format(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024 / 1024))

def shasta_cut(shasta_fname, min_len, threads, prefix, dir):
    # shasta_fname  shasta filename
    # min_len       minimum length of contigs
    # prefix        prefix of output files
    # dir           working directory

    start_time = time.time()

    ## final files are ${prefix}.primary.fasta ${prefix}.alt.fasta name_mapping.txt
    ## ${prefix}.hap1.cut.fasta ${prefix}.hap1.cut.bed 
    ## ${prefix}.hap2.cut.fasta ${prefix}.hap2.cut.bed
    ## ${prefix}.collapsed.cut.fasta ${prefix}.collapsed.cut.bed ${prefix}.pair.cut.txt
    ## ${prefix}.filtered.alt2pri.txt ${prefix}.pair.cut.txt

    pri_fname, alt_fname = prefix + '.primary.fasta', prefix + '.alt.fasta'
    name_mapping = prefix + '.name_mapping'
    hap1_fname, hap1_bfname = prefix + '.hap1.cut.fasta', prefix + '.hap1.cut.bed'
    hap2_fname, hap2_bfname = prefix + '.hap2.cut.fasta', prefix + '.hap2.cut.bed'
    coll_fname, coll_bfname = prefix + '.collapsed.cut.fasta', prefix + '.collapsed.cut.bed'
    pair_fname = prefix + '.pair.cut.txt'
    alt2pri_fname = prefix + '.filtered.alt2pri.paf'

    if _test_file(dir, pri_fname, alt_fname, name_mapping, hap1_fname, hap1_bfname, hap2_fname, hap2_bfname, coll_fname, coll_bfname, pair_fname, alt2pri_fname): return
    
    _file_exist(shasta_fname)

    try:
        os.chdir(dir)
    except FileNotFoundError:
        logger.error('Directory {} not exists!'.format(dir), exc_info = 1)
        exit(1)
    except PermissionError:
        logger.error('No permission to change to directory {}'.format(dir), exc_info = 1)
        exit(1)

    logger.info('{} Working directory {}'.format(sys._getframe().f_code.co_name, os.getcwd()))

    # create a temporary directory
    tmp_dir = str(uuid.uuid4())
    while os.path.exists(tmp_dir):
        tmp_dir = str(uuid.uuid4())
    os.makedirs(tmp_dir, exist_ok = True)
    shasta2unzip_and_mapping(shasta_fname, prefix, alt2pri_fname, min_len, threads, tmp_dir)
    # remove the temporary directory
    os.removedirs(tmp_dir)

    end_time = time.time()
    logger.info('{} ends'.format(sys._getframe().f_code.co_name))
    logger.info('Wall clock time elapsed: {} s'.format(end_time - start_time))
    logger.info('CPU time elapsed: {} s'.format(time.process_time()))
    logger.info('Memory usage: {:>.3f} GB\n'.format(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024 / 1024))

def dual_contig_cut(d1fname, d2fname, threshold, threads, prefix, dir):
    # d1fname       dual1 filename
    # d2fname       dual2 filename
    # threshold     threshold of match
    # threads       number of threads
    # prefix        prefix of output files
    # dir           working directory

    start_time = time.time()

    ## final files are ${prefix}.dual.paf
    ## ${prefix}.hap1.cut.fasta ${prefix}.hap1.cut.bed 
    ## ${prefix}.hap2.cut.fasta ${prefix}.hap2.cut.bed
    ## ${prefix}.collapsed1.cut.fasta ${prefix}.collapsed1.cut.bed
    ## ${prefix}.collapsed2.cut.fasta ${prefix}.collapsed2.cut.bed ${prefix}.pair.cut.txt
    map_fname = prefix + '.dual.paf'
    hap1_fname, hap1_bfname = prefix + '.hap1.cut.fasta', prefix + '.hap1.cut.bed'
    hap2_fname, hap2_bfname = prefix + '.hap2.cut.fasta', prefix + '.hap2.cut.bed'
    coll1_fname, coll1_bfname = prefix + '.collapsed1.cut.fasta', prefix + '.collapsed1.cut.bed'
    coll2_fname, coll2_bfname = prefix + '.collapsed2.cut.fasta', prefix + '.collapsed2.cut.bed'
    pair_fname = prefix + '.pair.cut.txt'
    if _test_file(dir, hap1_fname, hap1_bfname, hap2_fname, hap2_bfname, coll1_fname, coll1_bfname, coll2_fname, coll2_bfname, pair_fname): return

    _file_exist(d1fname, d2fname)

    _check_program('minimap2')

    try:
        os.chdir(dir)
    except FileNotFoundError:
        logger.error('Directory {} not exists!'.format(dir), exc_info = 1)
        exit(1)
    except PermissionError:
        logger.error('No permission to change to directory {}'.format(dir), exc_info = 1)
        exit(1)

    logger.info('{} Working directory {}'.format(sys._getframe().f_code.co_name, os.getcwd()))

    if not _test_file(dir, map_fname):
        cmd_map = 'minimap2 -cx asm20 --eqx --secondary=no -t ' + str(threads) + ' ' + d1fname + ' ' + d2fname + ' -o ' + map_fname
        logger.info('Command for minimap2: {}'.format(cmd_map))
        try:
            subprocess.call(cmd_map, shell = True, stderr = sys.stderr)
        except (subprocess.CalledProcessError, OSError) as e:
            logger.error('Error running minimap2 {}'.format(e.cmd), exc_info = 1)
            exit(1)

    _file_exist(map_fname)

    dual_cut(d1fname, d2fname, map_fname, threshold, threads, prefix)

    end_time = time.time()
    logger.info('{} ends'.format(sys._getframe().f_code.co_name))
    logger.info('Wall clock time elapsed: {} s'.format(end_time - start_time))
    logger.info('CPU time elapsed: {} s'.format(time.process_time()))
    logger.info('Memory usage: {:>.3f} GB\n'.format(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024 / 1024))

def cut_primary(pri_fname, alt_fname, flt_fname, prefix, dir):
    # pri_fname     primary filename
    # alt_fname     alt filename
    # flt_fname     filtered filename
    # prefix        prefix of output files
    # dir           working directory

    start_time = time.time()

    ## final files are ${prefix}.hap1.cut.fasta ${prefix}.hap1.cut.bed 
    ## ${prefix}.hap2.cut.fasta ${prefix}.hap2.cut.bed
    ## ${prefix}.collapsed.cut.fasta ${prefix}.collapsed.cut.bed ${prefix}.pair.cut.txt
    hap1_fname, hap1_bfname = prefix + '.hap1.cut.fasta', prefix + '.hap1.cut.bed'
    hap2_fname, hap2_bfname = prefix + '.hap2.cut.fasta', prefix + '.hap2.cut.bed'
    coll_fname, coll_bfname = prefix + '.collapsed.cut.fasta', prefix + '.collapsed.cut.bed'
    pair_fname = prefix + '.pair.cut.txt'
    if _test_file(dir, hap1_fname, hap1_bfname, hap2_fname, hap2_bfname, coll_fname, coll_bfname, pair_fname): return

    _file_exist(pri_fname, alt_fname, flt_fname)

    try:
        os.chdir(dir)
    except FileNotFoundError:
        logger.error('Directory {} not exists!'.format(dir), exc_info = 1)
        exit(1)
    except PermissionError:
        logger.error('No permission to change to directory {}'.format(dir), exc_info = 1)
        exit(1)

    logger.info('{} Working directory {}'.format(sys._getframe().f_code.co_name, os.getcwd()))

    cut_contigs(pri_fname, alt_fname, flt_fname, prefix)

    end_time = time.time()
    logger.info('{} ends'.format(sys._getframe().f_code.co_name))
    logger.info('Wall clock time elapsed: {} s'.format(end_time - start_time))
    logger.info('CPU time elapsed: {} s'.format(time.process_time()))
    logger.info('Memory usage: {:>.3f} GB\n'.format(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024 / 1024))

def mince_contigs(h1fname, h2fname, cfname, output, dir, c2fname = None):
    # h1fname   hap1 filename
    # h2fname   hap2 filename
    # cfname    collapsed filename
    # output    output filename
    # dir       working directory
    # c2fname   collapsed2 filename can be None

    start_time = time.time()

    ## final files are ${output} and ${output}.fai
    fai_fname = output + '.fai'
    if _test_file(dir, output, fai_fname): return

    _check_program('samtools')

    _file_exist(h1fname, h2fname, cfname)

    try:
        os.chdir(dir)
    except FileNotFoundError:
        logger.error('Directory {} not exists!'.format(dir), exc_info = 1)
        exit(1)
    except PermissionError:
        logger.error('No permission to change to directory {}'.format(dir), exc_info = 1)
        exit(1)

    logger.info('{} Working directory {}'.format(sys._getframe().f_code.co_name, os.getcwd()))

    if os.path.isfile(output):
        logger.warning('{} has exist, will be overwritten!'.format(output))
    
    cmd = 'cat ' + h1fname + ' ' + h2fname + ' ' + cfname + ' > ' + output
    if c2fname is not None:
        cmd = 'cat ' + h1fname + ' ' + h2fname + ' ' + cfname + ' ' + c2fname + ' > ' + output
    logger.info('Command for cat: {}'.format(cmd))
    try:
        subprocess.call(cmd, shell = True, stderr = sys.stderr)
    except (subprocess.CalledProcessError, OSError) as e:
        logger.error('Error running cat {}'.format(e.cmd), exc_info = 1)
        exit(1)

    cmd_faidx = 'samtools faidx ' + output
    logger.info('Command for samtools faidx: {}'.format(cmd_faidx))
    try:
        subprocess.call(cmd_faidx, shell = True, stderr = sys.stderr)
    except (subprocess.CalledProcessError, OSError) as e:
        logger.error('Error running samtools faidx {}'.format(e.cmd), exc_info = 1)
        exit(1)

    end_time = time.time()
    logger.info('{} ends'.format(sys._getframe().f_code.co_name))
    logger.info('Wall clock time elapsed: {} s'.format(end_time - start_time))
    logger.info('CPU time elapsed: {} s'.format(time.process_time()))
    logger.info('Memory usage: {:>.3f} GB\n'.format(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024 / 1024))

def run_bwa(rfname, p1fname, p2fname, output, threads, dir):
    # rfname    reference filename
    # p1fname   pair 1 filename
    # p2fname   pair 2 filename
    # ofname    output file name
    # threads   number of threads
    # dir       working directory

    start_time = time.time()

    ## final files are ${output}
    if _test_file(dir, output): return

    # _file_exist(rfname, p1fname, p2fname)

    _check_program('bwa', 'samtools')

    try:
        os.chdir(dir)
    except FileNotFoundError:
        logger.error('Directory {} not exists!'.format(dir), exc_info = 1)
        exit(1)
    except PermissionError:
        logger.error('No permission to change to directory {}'.format(dir), exc_info = 1)
        exit(1)

    logger.info('{} Working directory {}'.format(sys._getframe().f_code.co_name, os.getcwd()))

    if os.path.isfile(output):
        logger.warning('{} has exist, will be overwritten!'.format(output))
    if threads < 1:
        logger.info('threads {} must be >= 1, set to 1.'.format(threads))
        threads = 1

    cmd_index = 'bwa index ' + rfname
    logger.info('Command for bwa index: {}'.format(cmd_index))
    try:
        subprocess.call(cmd_index, shell = True, stderr = open('bwa.log', 'w'))
    except (subprocess.CalledProcessError, OSError) as e:
        logger.error('Error running bwa index {}'.format(e.cmd), exc_info = 1)
        exit(1)

    cmd_map = 'bwa mem -5 -t ' + str(threads) + ' ' + rfname + ' ' + p1fname + ' ' + p2fname + ' | '
    cmd_map += 'samtools view -bS -h -F 2316 -@ 47 -o ' + output
    logger.info('Command for bwa mem: {}'.format(cmd_map))
    try:
        subprocess.call(cmd_map, shell = True, stderr = open('bwa.log', 'a'))
    except (subprocess.CalledProcessError, OSError) as e:
        logger.error('Error running bwa mem {}'.format(e.cmd), exc_info = 1)
        exit(1)

    end_time = time.time()
    logger.info('{} ends'.format(sys._getframe().f_code.co_name))
    logger.info('Wall clock time elapsed: {} s'.format(end_time - start_time))
    logger.info('CPU time elapsed: {} s'.format(time.process_time()))
    logger.info('Memory usage: {:>.3f} GB\n'.format(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024 / 1024))

def run_clair3(pri_fname, alt_fname, rdfname, model, stype, threads, prefix, dir):
    # pri_fname     primary filename
    # alt_fname     alt filename
    # rdfname       reads filename
    # model         model name used by clair3
    # stype         sequence platform [clr | hifi | ont]
    # threads       number of threads
    # prefix        prefix of output files
    # dir           working directory

    ## final files are clair3pri/merge_output.vcf.gz and clair3alt/merge_output.vcf.gz
    vcf_pri = 'clair3pri/merge_output.vcf.gz'
    vcf_alt = 'clair3alt/merge_output.vcf.gz'
    if _test_file(dir, vcf_pri, vcf_alt): return

    _file_exist(pri_fname, alt_fname, rdfname, model)

    _check_program('minimap2', 'samtools', 'run_clair3.sh')

    try:
        os.chdir(dir)
    except FileNotFoundError:
        logger.error('Directory {} not exists!'.format(dir), exc_info = 1)
        exit(1)
    except PermissionError:
        logger.error('No permission to change to directory {}'.format(dir), exc_info = 1)
        exit(1)

    logger.info('{} Working directory {}'.format(sys._getframe().f_code.co_name, os.getcwd()))

    if threads < 1:
        logger.info('threads {} must be >= 1, set to 1.'.format(threads))
        threads = 1

    threads_samtools = threads - 1
    cmd_map = 'minimap2 -ax map-'
    if stype == 'clr':
        cmd_map += 'pb '
    elif stype == 'hifi':
        cmd_map += 'hifi'
    elif stype == 'ont':
        cmd_map += 'ont'
    else:
        logger.error('Unregonized sequence type {}, [clr | hifi | ont]'.format(stype), exc_info = 1)
        exit(1)
    cmd_map += ' --secondary=no -t ' + str(threads)

    # mapping reads to primary 
    mapped_pri = prefix + '.rd2pri.sorted.bam'
    pri_err_fname = 'minimap2.pri.log'
    # cmd_map_pri = cmd_map + ' ' + pri_fname + ' ' + rdfname + ' | samtools view -bSh -F 2316 -@ ' + str(threads_samtools)
    # cmd_map_pri += ' | samtools sort -@ ' + str(threads_samtools) + ' - -o ' + mapped_pri
    # print('\t{}'.format(cmd_map_pri), file = sys.stderr)
    # try:
    #     subprocess.call(cmd_map_pri, shell = True, stderr = sys.stderr)
    # except (subprocess.CalledProcessError, OSError) as e:
    #     print('Error mapping reads to primary {}'.format(e.cmd), file = sys.stderr)
    #     exit(1)
    # try:
    #     subprocess.call('samtools index ' + mapped_pri + ' -@ ' + str(threads_samtools), shell = True, 
    #                     stderr = sys.stderr)
    # except (subprocess.CalledProcessError, OSError) as e:
    #     print('Error running samtools index {}'.format(e.cmd), file = sys.stderr)
    #     exit(1)
    
    # mapping reads to alt
    mapped_alt = prefix + '.rd2alt.sorted.bam'
    alt_err_fname = 'minimap2.alt.log'
    # cmd_map_alt = cmd_map + ' ' + alt_fname + ' ' + rdfname + ' | samtools view -bSh -F 2316 -@ ' + str(threads_samtools)
    # cmd_map_alt += ' | samtools sort -@ ' + str(threads_samtools) + ' - -o ' + mapped_alt
    # print('\t{}'.format(cmd_map_alt), file = sys.stderr)
    # try:
    #     subprocess.call(cmd_map_alt, shell = True, stderr = sys.stderr)
    # except (subprocess.CalledProcessError, OSError) as e:
    #     print('Error mapping reads to alt {}'.format(e.cmd), file = sys.stderr)
    #     exit(1)
    # try:
    #     subprocess.call('samtools index ' + mapped_alt + ' -@ ' + str(threads_samtools), shell = True, 
    #                     stderr = sys.stderr)
    # except (subprocess.CalledProcessError, OSError) as e:
    #     print('Error running samtools index {}'.format(e.cmd), file = sys.stderr)
    #     exit(1)
    
    def _map_reads_(ref_fname, read_fname, out_fname, err_fname):

        start_time = time.time()

        cmd_map_x = cmd_map + ' ' + ref_fname + ' ' + read_fname
        cmd_map_x += ' | samtools view -bSh -F 2316 -@ ' + str(threads_samtools)
        cmd_map_x += ' | samtools sort -@ ' + str(threads_samtools) + ' - -o ' + out_fname
        logger.info('Command for mapping reads: {}'.format(cmd_map_x))
        try:
            subprocess.call(cmd_map_x, shell = True, stderr = open(err_fname, 'w'))
        except (subprocess.CalledProcessError, OSError) as e:
            logger.error('Error mapping reads {}'.format(e.cmd), exc_info = 1)
            exit(1)
        
        end_time = time.time()
        logger.info('{} ends'.format(sys._getframe().f_code.co_name))
        logger.info('Wall clock time elapsed: {} s'.format(end_time - start_time))
        logger.info('CPU time elapsed: {} s'.format(time.process_time()))
        logger.info('Memory usage: {:>.3f} GB\n'.format(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024 / 1024))

    cmd_clair3 = 'run_clair3.sh --threads=' + str(threads) + ' --platform='
    if stype == 'clr':
        cmd_clair3 += 'ont'
    elif stype == 'hifi':
        cmd_clair3 += 'hifi'
    elif stype == 'ont':
        cmd_clair3 += 'ont'
    else:
        logger.error('Unregonized sequence type {}, [clr | hifi | ont]'.format(stype), exc_info = 1)
        exit(1)
    cmd_clair3 += ' --model_path=' + model + ' --call_snp_only --include_all_ctgs --output=clair3'

    pri_clair3_err_fname = 'clair3.pri.log'
    # cmd_clair3_pri = cmd_clair3 + 'pri --bam_fn=' + mapped_pri + ' --ref_fn=' + pri_fname
    # print('\t{}'.format(cmd_clair3_pri), file = sys.stderr)
    # try:
    #     subprocess.call(cmd_clair3_pri, shell = True, stderr = sys.stderr)
    # except (subprocess.CalledProcessError, OSError) as e:
    #     print('Error running clair3 for primary {}'.format(e.cmd), file = sys.stderr)
    #     exit(1)
    
    alt_clair3_err_fname = 'clair3.alt.log'
    # cmd_clair3_alt = cmd_clair3 + 'alt --bam_fn=' + mapped_alt + ' --ref_fn=' + alt_fname
    # print('\t{}'.format(cmd_clair3_alt), file = sys.stderr)
    # try:
    #     subprocess.call(cmd_clair3_alt, shell = True, stderr = sys.stderr)
    # except (subprocess.CalledProcessError, OSError) as e:
    #     print('Error running clair3 for alt {}'.format(e.cmd), file = sys.stderr)
    #     exit(1)
    
    def _run_clair3_x_(ref_fname, bam_fname, out_pre, err_fname):

        start_time = time.time()

        cmd_clair3_x = cmd_clair3 + out_pre + ' --bam_fn=' + bam_fname + ' --ref_fn=' + ref_fname
        logger.info('Command for clair3: {}'.format(cmd_clair3_x))
        if not os.path.exists(ref_fname + '.fai'):
            try:
                subprocess.call('samtools faidx ' + ref_fname, shell = True, stderr = sys.stderr)
            except (subprocess.CalledProcessError, OSError) as e:
                logger.error('Error running samtools faidx {}'.format(e.cmd), exc_info = 1)
                exit(1)
            
        if not os.path.exists(bam_fname + '.bai'):
            try:
                subprocess.call('samtools index ' + bam_fname + ' -@ ' + str(threads_samtools), 
                                shell = True, stderr = sys.stderr)
            except (subprocess.CalledProcessError, OSError) as e:
                logger.error('Error running samtools index {}'.format(e.cmd), exc_info = 1)
                exit(1)
            
        try:
            subprocess.call(cmd_clair3_x, shell = True, stderr = open(err_fname, 'w'))
        except (subprocess.CalledProcessError, OSError) as e:
            logger.error('Error running clair3 {}'.format(e.cmd), exc_info = 1)
            exit(1)

        end_time = time.time()
        logger.info('{} ends'.format(sys._getframe().f_code.co_name))
        logger.info('Wall clock time elapsed: {} s'.format(end_time - start_time))
        logger.info('CPU time elapsed: {} s'.format(time.process_time()))
        logger.info('Memory usage: {:>.3f} GB\n'.format(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024 / 1024))

    if not _test_file(dir, vcf_pri):
        if not _test_file(dir, mapped_pri):
            _map_reads_(pri_fname, rdfname, mapped_pri, pri_err_fname)
        _run_clair3_x_(pri_fname, mapped_pri, 'pri', pri_clair3_err_fname)
    if not _test_file(dir, vcf_alt):
        if not _test_file(dir, mapped_alt):
            _map_reads_(alt_fname, rdfname, mapped_alt, alt_err_fname)
        _run_clair3_x_(alt_fname, mapped_alt, 'alt', alt_clair3_err_fname)

def run_phasing(pair_fname, fai_fname, clair_pri_fname, clair_alt_fname, alt2pri_fname, hic_fname, mat_fname, var_fname, 
              threads, iter, seed, mapq, prefix, filtered, porec, 
              dump_var, dump_filtered, dump_mat, print_mat, dir):
    # pair_fname         pair filename
    # fai_fname         fai filename
    # clair_pri_fname   clair3 primary filename
    # clair_alt_fname   clair3 alt filename
    # alt2pri_fname     mapping filename of alt to pri
    # hic_fname         hic mapping filename
    # mat_fname         matrix filenme
    # var_fname         snp filename
    # threads           number of threads to be used
    # iter              iteration to phase
    # seed              seed for random function
    # mapq              mapping quality threshold
    # prefix            prefix of output files
    # filtered          whether hic mapping are filtered
    # porec             whether the reads is porec
    # dump_var          dump snp into ${prefix}.var.txt
    # dump_filtered     dump filtered hic mapping into ${prefix}.hic.filtered.bam
    # dump_mat          dump the matrix into ${prefix}.binmat
    # print_mat         print the matrix to stdout
    # dir               working directory

    start_time = time.time()

    ## final file is ${prefix}.result.txt
    final_file = prefix + '.result.txt'
    if _test_file(dir, final_file): return

    _file_exist(pair_fname, fai_fname)

    _check_program('phasing')
    
    try:
        os.chdir(dir)
    except FileNotFoundError:
        logger.error('Directory {} not exists!'.format(dir), exc_info = 1)
        exit(1)
    except PermissionError:
        logger.error('No permission to change to directory {}'.format(dir), exc_info = 1)
        exit(1)

    logger.info('{} Working directory {}'.format(sys._getframe().f_code.co_name, os.getcwd()))
    
    if mat_fname is not None:
        _file_exist(mat_fname)
    else:
        if hic_fname is not None:
            _file_exist(hic_fname)
        else:
            logger.error('Please specify a hic mapping file or a matrix file!', exc_info = 1)
            exit(1)
    
    if var_fname is not None:
        _file_exist(var_fname)
    else: 
        if alt2pri_fname is not None:
            _file_exist(alt2pri_fname)
        else:
            logger.error('Please specify a snp file or a mapping file of alt to pri!', exc_info = 1)
            exit(1)
        if clair_pri_fname is not None and clair_alt_fname is not None:
            _file_exist(clair_pri_fname)
            _file_exist(clair_alt_fname)
        elif clair_pri_fname is None and clair_alt_fname is None:
            logger.info('No clair3 called snp file specified, use the mapping file {} to get snp infomation.'.format(alt2pri_fname))
        else:
            logger.error('Please specify both clair3 called snp file!', exc_info = 1)
            exit(1)

    if threads < 1:
        logger.info('threads {} must be >= 1, set to 1.'.format(threads))
        threads = 1
    if iter < 100:
        logger.info('iter {} should be >= 100, set to 100'.format(iter))
        iter = 100

    cmd_phase = '/usr/bin/time -v phasing -a ' + pair_fname + ' -c ' + fai_fname + ' -t ' + str(threads) + ' -n ' 
    cmd_phase += str(iter) + ' -s ' + str(seed) + ' -q ' + str(mapq) + ' -p ' + prefix
    if mat_fname is not None:
        cmd_phase += ' -m ' + mat_fname
    else:
        cmd_phase += ' -i ' + hic_fname
    if var_fname is not None:
        cmd_phase += ' -v ' + var_fname
    else:
        cmd_phase += ' -b ' + alt2pri_fname
        if clair_pri_fname is not None and clair_alt_fname is not None:
            cmd_phase += ' --priclair ' + clair_pri_fname + ' --altclair ' + clair_alt_fname
    if filtered:
        cmd_phase += ' --filtered'
    if porec:
        cmd_phase += ' --porec'
    if dump_var:
        cmd_phase += ' --dump-var'
    if dump_filtered:
        cmd_phase += ' --dump-filt'
    if dump_mat:
        cmd_phase += ' --dump-mat'
    if print_mat:
        ofname = prefix + '.snp.count.txt'
        cmd_phase += ' --print-mat > ' + ofname
    
    logger.info(cmd_phase)
    try:
        subprocess.call(cmd_phase, shell = True, stderr = open('phasing.log', 'w'))
    except (subprocess.CalledProcessError, OSError) as e:
        logger.error('Error running phasing {}'.format(e.cmd), exc_info = 1)
        exit(1)

    end_time = time.time()
    logger.info('{} ends'.format(sys._getframe().f_code.co_name))
    logger.info('Wall clock time elapsed: {} s'.format(end_time - start_time))
    logger.info('CPU time elapsed: {} s'.format(time.process_time()))
    logger.info('Memory usage: {:>.3f} GB\n'.format(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024 / 1024))

def finalize(rfname, h2fname, cfname, mfname, sfname, prefix, dir, afname, pfname):
    # rfname    phasing result filename
    # h2fname   hap2 filename
    # cfname    collapsed filename
    # mfname    minced filename
    # sfname    switch filename
    # prefix    prefix of output files
    # dir       working directory
    # afname    pair filename
    # pfname    paf filename

    # hap2 and collapsed are both from primary

    start_time = time.time()

    ## final files are ${prefix}.phased0.fasta ${prefix}.phased1.fasta ${prefix}.phased0.bed ${prefix}.phased1.bed
    hap1_fname, hap1_bfname = prefix + '.phased0.fasta', prefix + '.phased0.bed'
    hap2_fname, hap2_bfname = prefix + '.phased1.fasta', prefix + '.phased1.bed'
    fix_sfname = prefix + '.fixed.switch'
    if _test_file(dir, hap1_fname, hap2_fname, hap1_bfname, hap2_bfname): return

    _file_exist(rfname, h2fname, cfname, mfname)

    try:
        os.chdir(dir)
    except FileNotFoundError:
        logger.error('Directory {} not exists!'.format(dir), exc_info = 1)
        exit(1)
    except PermissionError:
        logger.error('No permission to change to directory {}'.format(dir), exc_info = 1)
        exit(1)

    logger.info('{} Working directory {}'.format(sys._getframe().f_code.co_name, os.getcwd()))

    fix_switch2(sfname, afname, pfname, fix_sfname)
    phase_contigs(rfname, h2fname, cfname, mfname, fix_sfname, prefix)
    # phase_contigs(rfname, h2fname, cfname, mfname, sfname, prefix)

    end_time = time.time()
    logger.info('{} ends'.format(sys._getframe().f_code.co_name))
    logger.info('Wall clock time elapsed: {} s'.format(end_time - start_time))
    logger.info('CPU time elapsed: {} s'.format(time.process_time()))
    logger.info('Memory usage: {:>.3f} GB\n'.format(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024 / 1024))

def group(bfname, b1fname, b2fname, h1fname, h2fname, dir):
    # bfname    bam filename
    # b1fname   bed1 filename
    # b2fname   bed2 filename
    # h1fname   hap1 filename
    # h2fname   hap2 filename
    # dir       working directory

    start_time = time.time()
    o1fname, o2fname = h1fname, h2fname
    try:
        subprocess.call('mv ' + h1fname + ' ' + h1fname + '.bak', shell = True)
        subprocess.call('mv ' + h2fname + ' ' + h2fname + '.bak', shell = True)
    except (subprocess.CalledProcessError, OSError) as e:
        logger.error('Error running mv {}'.format(e.cmd), exc_info = 1)
        exit(1)
    gfname = 'group.txt'
    if _test_file(dir, gfname): return
    _file_exist(bfname, b1fname, b2fname, h1fname, h2fname)
    
    try:
        os.chdir(dir)
    except FileNotFoundError:
        logger.error('Directory {} not exists!'.format(dir), exc_info = 1)
        exit(1)
    except PermissionError:
        logger.error('No permission to change to directory {}'.format(dir), exc_info = 1)
        exit(1)

    logger.info('{} Working directory {}'.format(sys._getframe().f_code.co_name, os.getcwd()))

    cmd_group = 'group -b ' + bfname + ' --b1 ' + b1fname + ' --b2 ' + b2fname + ' -o ' + gfname + ' -i 100'
    logger.info(cmd_group)

    try:
        subprocess.call(cmd_group, shell = True, stderr = open('group.log', 'w'))
    except (subprocess.CalledProcessError, OSError) as e:
        logger.error('Error running phasing {}'.format(e.cmd), exc_info = 1)
        exit(1)
    group2(h1fname, h2fname, gfname, o1fname, o2fname)

    end_time = time.time()
    logger.info('{} ends'.format(sys._getframe().f_code.co_name))
    logger.info('Wall clock time elapsed: {} s'.format(end_time - start_time))
    logger.info('CPU time elapsed: {} s'.format(time.process_time()))
    logger.info('Memory usage: {:>.3f} GB\n'.format(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024 / 1024))

def finalize_dual(rfname, h1fname, c1fname, h2fname, c2fname, mfname, sfname, prefix, dir, afname, pfname):
    # rfname    phasing result filename
    # h1fname   hap1 bed filename
    # c1fname   collapsed1 bed filename
    # h2fname   hap2 bed filename
    # c2fname   collapsed2 bed filename
    # mfname    minced filename
    # sfname    switch filename
    # prefix    prefix of output files
    # dir       working directory
    # afname    pair filename
    # pfname    paf filename

    # hap2 and collapsed are both from primary

    start_time = time.time()

    ## final files are ${prefix}.phased0.fasta ${prefix}.phased1.fasta ${prefix}.phased0.bed ${prefix}.phased1.bed
    hap1_fname, hap1_bfname = prefix + '.phased0.fasta', prefix + '.phased0.bed'
    hap2_fname, hap2_bfname = prefix + '.phased1.fasta', prefix + '.phased1.bed'
    fix_sfname = prefix + '.fixed.switch'
    if _test_file(dir, hap1_fname, hap2_fname, hap1_bfname, hap2_bfname): return

    _file_exist(rfname, h1fname, c1fname, h2fname, c2fname, mfname)

    try:
        os.chdir(dir)
    except FileNotFoundError:
        logger.error('Directory {} not exists!'.format(dir), exc_info = 1)
        exit(1)
    except PermissionError:
        logger.error('No permission to change to directory {}'.format(dir), exc_info = 1)
        exit(1)

    logger.info('{} Working directory {}'.format(sys._getframe().f_code.co_name, os.getcwd()))

    fix_switch2(sfname, afname, pfname, fix_sfname)
    phase_contigs_dual(rfname, h1fname, h2fname, c1fname, c2fname, mfname, fix_sfname, prefix)
    # phase_contigs_dual(rfname, h1fname, h2fname, c1fname, c2fname, mfname, sfname, prefix)

    end_time = time.time()
    logger.info('{} ends'.format(sys._getframe().f_code.co_name))
    logger.info('Wall clock time elapsed: {} s'.format(end_time - start_time))
    logger.info('CPU time elapsed: {} s'.format(time.process_time()))
    logger.info('Memory usage: {:>.3f} GB\n'.format(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024 / 1024))

# logging all the arguments
def log_arguments(args):
    args_vars = vars(args)
    logger.info('==================== Arguments ====================')
    for arg in args_vars.keys():
        if arg != 'func' and args_vars[arg] is not None:
            logger.info('{:<20s} {}'.format(arg, args_vars[arg]))
    logger.info('==================== Arguments ====================\n')

def run_phase(args):
    # phasing for unzip format

    args.pri_fname = os.path.abspath(args.pri_fname)
    args.alt_fname = os.path.abspath(args.alt_fname)
    # if args.porec is True, then args.porec_fname must be not None
    if args.porec:
        if args.hic1fname is not None or args.hic2fname is not None:
            logger.warning('Hi-C mapping files will not be used in porec process!')
        if args.porec_fname is None:
            logger.error('Please specify the porec filename!')
            exit(1)
        args.porec_fname = os.path.abspath(args.porec_fname)
    elif args.hic1fname is not None and args.hic2fname is not None:
        args.hic1fname = os.path.abspath(args.hic1fname)
        args.hic2fname = os.path.abspath(args.hic2fname)
    else:
        logger.error('Please specify the Hi-C mapping files!')
        exit(1)
    # if args.model is not None, then args.rdfname and args.type must be not None
    if args.model is not None:
        if args.rdfname is None:
            logger.error('Please specify the reference filename!')
            exit(1)
        if args.type is None:
            logger.error('Please specify the reference type!')
            exit(1)
        args.rdfname = os.path.abspath(args.rdfname)
    
    if args.threads < 1:
        logger.warning('Thread number must be greater than 0! Set to 1.')
        args.thread = 1
    if args.iter < 10000:
        logger.warning('Iteration number must be greater than 10000! Set to 10000.')
        args.iter = 10000
    if args.mapq < 1:
        logger.info('mapq {} should be >= 1, set to 1'.format(args.mapq))
        args.mapq = 1
    
    log_arguments(args)

    logger.info('{} Working directory {}'.format(sys._getframe().f_code.co_name, os.getcwd()))

    try:
        os.chdir(args.dir)
    except FileNotFoundError:
        logger.error('Directory {} not exist!'.format(args.dir), exc_info = 1)
        exit(1)
    except PermissionError:
        logger.error('No permissions to change to {}!'.format(args.dir), exc_info = 1)
        exit(1)
    
    args.matfname = os.path.abspath(args.matfname) if args.matfname is not None else None
    args.varfname = os.path.abspath(args.varfname) if args.varfname is not None else None

    
    # sub directory to store cut and bwa result
    dir_cut = os.path.join(args.dir, 'minced')
    # sub directory to store clair3 result
    dir_clair3 = os.path.join(args.dir, 'clair3')
    # sub directory to store phasing result
    dir_phase = os.path.join(args.dir, 'phase')

    # clean primary and alt, results store in ${args.dir}
    clean_unzip(args.pri_fname, args.alt_fname, args.prefix, args.dir)

    # mapping alt to primary, results store in ${args.dir}
    clean_pri_fname = os.path.join(args.dir, args.prefix + '.p_tigs.fasta')
    clean_alt_fname = os.path.join(args.dir, args.prefix + '.a_tigs.fasta')
    name_mapping_fname = os.path.join(args.dir, 'name_mapping.txt')
    alt2pri_fname = os.path.join(args.dir, args.prefix + '.alt2pri')
    make_mapping_alt_to_primary(clean_pri_fname, clean_alt_fname, name_mapping_fname, alt2pri_fname, args.threads, 'paf', args.dir)

    # place alt along with primary, results store in ${args.dir}
    place_alt_along_primary_and_filter_contained(alt2pri_fname + '.paf', args.prefix, args.dir)

    # cut primary, results store in ${dir_cut}
    if os.path.exists(dir_cut):
        logger.warning('The sub directory minced exist in {}, some files will be overwritten!'.format(args.dir))
    else:
        os.makedirs(dir_cut, exist_ok = True)
    filtered_placed_file = os.path.join(args.dir, args.prefix + '.filtered.alt2pri.txt')
    cut_primary(clean_pri_fname, clean_alt_fname, filtered_placed_file, args.prefix, dir_cut)

    # mince, results store in ${dir_cut}
    hap1_fname = os.path.join(dir_cut, args.prefix + '.hap1.cut.fasta')
    hap1_bfname = os.path.join(dir_cut, args.prefix + '.hap1.cut.bed')
    hap2_fname = os.path.join(dir_cut, args.prefix + '.hap2.cut.fasta')
    hap2_bfname = os.path.join(dir_cut, args.prefix + '.hap2.cut.bed')
    collapsed_fname = os.path.join(dir_cut, args.prefix + '.collapsed.cut.fasta')
    collapsed_bfname = os.path.join(dir_cut, args.prefix + '.collapsed.cut.bed')
    minced_fname = os.path.join(dir_cut, args.prefix + '.minced.fasta')
    mince_contigs(hap1_fname, hap2_fname, collapsed_fname, minced_fname, dir_cut)

    # run bwa, results store in ${dir_cut}
    if not args.porec:
        hic_mapping_fname = os.path.join(dir_cut, args.prefix + '.bwa.filtered2316.bam')
        run_bwa(minced_fname, args.hic1fname, args.hic2fname, hic_mapping_fname, args.threads, dir_cut)
    else:
        # porec process
        hic_mapping_fname = ''

    # run clair3 if model provided, results store in ${dir_clair3}
    if args.model is not None:
        if os.path.exists(dir_clair3):
            logger.warning('The sub directory clair3 exist in {}, some files may be overwritten!'.format(args.dir))
        else:
            os.makedirs(dir_clair3, exist_ok = True)

        run_clair3(clean_pri_fname, clean_alt_fname, args.rdfname, args.model, args.type, args.threads, args.prefix, dir_clair3)
    else:
        logger.warning('No model provided, clair3 will not run. the accuracy of phasing will be affected!')
    
    # run phasing, results store in ${dir_phase}
    if os.path.exists(dir_phase):
        logger.warning('The sub directory clair3 exist in {}, some files will be overwritten!'.format(args.dir))
    else:
        os.makedirs(dir_phase, exist_ok = True)
    
    pair_fname = os.path.join(dir_cut, args.prefix + '.pair.cut.txt')
    fai_fname = os.path.join(dir_cut, args.prefix + '.minced.fasta.fai')
    if args.model is not None:
        pri_clair3_fname = os.path.join(dir_clair3, 'clair3pri/merge_output.vcf.gz')
        alt_clair3_fname = os.path.join(dir_clair3, 'clair3alt/merge_output.vcf.gz')
    else:
        pri_clair3_fname, alt_clair3_fname = None, None
    
    filtered_paf_file = os.path.join(args.dir, args.prefix + '.filtered.alt2pri.paf')

    run_phasing(pair_fname, fai_fname, pri_clair3_fname, alt_clair3_fname, filtered_paf_file, hic_mapping_fname, 
            args.matfname, args.varfname, args.threads, args.iter, args.seed, args.mapq, args.prefix, 
            args.filtered, args.porec, args.dump_var, args.dump_filtered, args.dump_mat, args.print_mat, dir_phase)

    result_fname = os.path.join(dir_phase, args.prefix + '.result.txt')
    switch_fname = os.path.join(dir_phase, args.prefix + '.switch')
    finalize(result_fname, hap2_bfname, collapsed_bfname, minced_fname, switch_fname, args.prefix, args.dir, pair_fname, filtered_paf_file)

def run_phase_pecat(args):
    # phasing for pecat

    if args.pri_asm_fname is not None or args.alt_asm_fname is not None:
        args.pri_asm_fname = os.path.abspath(args.pri_asm_fname)
        args.alt_asm_fname = os.path.abspath(args.alt_asm_fname)
    args.pri_pol_fname = os.path.abspath(args.pri_pol_fname)
    args.alt_pol_fname = os.path.abspath(args.alt_pol_fname)
    # if args.porec is True, then args.porec_fname must be not None
    if args.porec:
        if args.hic1fname is not None or args.hic2fname is not None:
            logger.warning('Hi-C mapping files will not be used in porec process!')
        if args.porec_fname is None:
            logger.error('Please specify the porec filename!')
            exit(1)
        args.porec_fname = os.path.abspath(args.porec_fname)
    elif args.hic1fname is None or args.hic2fname is None:
        logger.error('Please specify the hic mapping files!')
        exit(1)
    else:
        args.hic1fname = os.path.abspath(args.hic1fname)
        args.hic2fname = os.path.abspath(args.hic2fname)
    # if args.model is not None, then args.rdfname and args.type must be not None
    if args.model:
        if args.rdfname is None:
            logger.error('Please specify the raw reads filename!')
            exit(1)
        if args.type is None:
            logger.error('Please specify the type of the raw sequence!')
            exit(1)
        args.rdfname = os.path.abspath(args.rdfname)
        args.model = os.path.abspath(args.model)
    # args.thread must be >= 1
    if args.threads < 1:
        logger.info('threads {} must be >= 1, set to 1.'.format(args.threads))
        args.threads = 1
    # args.iter must be >= 100
    if args.iter < 10000:
        logger.info('iter {} should be >= 100, set to 100'.format(args.iter))
        args.iter = 10000
    if args.mapq < 1:
        logger.info('mapq {} should be >= 1, set to 1'.format(args.mapq))
        args.mapq = 1
    
    log_arguments(args)
    
    logger.info('{} Working directory {}'.format(sys._getframe().f_code.co_name, os.getcwd()))

    # if args.dir is None:
    #     args.dir = os.getcwd()
    # else:
    #     args.dir = os.path.abspath(args.dir)
    
    # # create the working directory
    # if not os.path.exists(args.dir):
    #     os.makedirs(args.dir, exist_ok = True)

    try:
        os.chdir(args.dir)
    except FileNotFoundError:
        logger.error('Directory {} not exist!'.format(args.dir), exc_info = 1)
        exit(1)
    except PermissionError:
        logger.error('No permissions to change to {}!'.format(args.dir), exc_info = 1)
        exit(1)
    
    args.matfname = os.path.abspath(args.matfname) if args.matfname is not None else None
    args.varfname = os.path.abspath(args.varfname) if args.varfname is not None else None

    # sub directory to store preprocessing result
    dir_pre = os.path.join(args.dir, 'preprocessing')
    # sub directory to store cut and bwa result
    dir_cut = os.path.join(args.dir, 'minced')
    # sub directory to store clair3 result
    dir_clair3 = os.path.join(args.dir, 'clair3')
    # sub directory to store phasing result
    dir_phase = os.path.join(args.dir, 'phase')

    # process the format of pecat to unzip, results store in ${dir_pre}
    if os.path.exists(dir_pre):
        logger.warning('The sub directory preprocessing exists in {}, some files may be overwritten!'.format(args.dir))
    else:
        os.makedirs(dir_pre, exist_ok = True)

    pecat_to_unzip_format(args.pri_asm_fname, args.alt_asm_fname, args.pri_pol_fname, args.alt_pol_fname, args.prefix, dir_pre)

    # clean primary and alt, results store in ${args.dir}
    pre_pri_fname = os.path.join(dir_pre, args.prefix + '.primary.clean.fasta')
    pre_alt_fname = os.path.join(dir_pre, args.prefix + '.alt.clean.fasta')
    clean_unzip(pre_pri_fname, pre_alt_fname, args.prefix, args.dir)

    # mapping alt to primary, results store in ${args.dir}
    clean_pri_fname = os.path.join(args.dir, args.prefix + '.p_tigs.fasta')
    clean_alt_fname = os.path.join(args.dir, args.prefix + '.a_tigs.fasta')
    name_mapping_fname = os.path.join(args.dir, 'name_mapping.txt')
    alt2pri_fname = os.path.join(args.dir, args.prefix + '.alt2pri')
    make_mapping_alt_to_primary(clean_pri_fname, clean_alt_fname, name_mapping_fname, alt2pri_fname, args.threads, 'paf', args.dir)

    # place alt along with primary, results store in ${args.dir}
    place_alt_along_primary_and_filter_contained(alt2pri_fname + '.paf', args.prefix, args.dir)

    # cut primary, results store in ${dir_cut}
    if os.path.exists(dir_cut):
        logger.warning('The sub directory minced exist in {}, some files will be overwritten!'.format(args.dir))
    else:
        os.makedirs(dir_cut, exist_ok = True)
    filtered_placed_file = os.path.join(args.dir, args.prefix + '.filtered.alt2pri.txt')
    cut_primary(clean_pri_fname, clean_alt_fname, filtered_placed_file, args.prefix, dir_cut)

    # mince, results store in ${dir_cut}
    hap1_fname = os.path.join(dir_cut, args.prefix + '.hap1.cut.fasta')
    hap1_bfname = os.path.join(dir_cut, args.prefix + '.hap1.cut.bed')
    hap2_fname = os.path.join(dir_cut, args.prefix + '.hap2.cut.fasta')
    hap2_bfname = os.path.join(dir_cut, args.prefix + '.hap2.cut.bed')
    collapsed_fname = os.path.join(dir_cut, args.prefix + '.collapsed.cut.fasta')
    collapsed_bfname = os.path.join(dir_cut, args.prefix + '.collapsed.cut.bed')
    minced_fname = os.path.join(dir_cut, args.prefix + '.minced.fasta')
    mince_contigs(hap1_fname, hap2_fname, collapsed_fname, minced_fname, dir_cut)

    # run bwa, results store in ${dir_cut}
    if not args.porec:
        hic_mapping_fname = os.path.join(dir_cut, args.prefix + '.bwa.filtered2316.bam')
        run_bwa(minced_fname, args.hic1fname, args.hic2fname, hic_mapping_fname, args.threads, dir_cut)
    else:
        # porec process
        hic_mapping_fname = ''

    # run clair3 if model provided, results store in ${dir_clair3}
    if args.model is not None:
        if os.path.exists(dir_clair3):
            logger.warning('The sub directory clair3 exist in {}, some files may be overwritten!'.format(args.dir))
        else:
            os.makedirs(dir_clair3, exist_ok = True)

        run_clair3(clean_pri_fname, clean_alt_fname, args.rdfname, args.model, args.type, args.threads, args.prefix, dir_clair3)
    else:
        logger.warning('No model provided, clair3 will not run. the accuracy of phasing will be affected!')
    
    # run phasing, results store in ${dir_phase}
    if os.path.exists(dir_phase):
        logger.warning('The sub directory clair3 exist in {}, some files will be overwritten!'.format(args.dir))
    else:
        os.makedirs(dir_phase, exist_ok = True)
    
    pair_fname = os.path.join(dir_cut, args.prefix + '.pair.cut.txt')
    fai_fname = os.path.join(dir_cut, args.prefix + '.minced.fasta.fai')
    if args.model is not None:
        pri_clair3_fname = os.path.join(dir_clair3, 'clair3pri/merge_output.vcf.gz')
        alt_clair3_fname = os.path.join(dir_clair3, 'clair3alt/merge_output.vcf.gz')
    else:
        pri_clair3_fname, alt_clair3_fname = None, None
    
    filtered_paf_file = os.path.join(args.dir, args.prefix + '.filtered.alt2pri.paf')

    run_phasing(pair_fname, fai_fname, pri_clair3_fname, alt_clair3_fname, filtered_paf_file, hic_mapping_fname, 
            args.matfname, args.varfname, args.threads, args.iter, args.seed, args.mapq, args.prefix, 
            args.filtered, args.porec, args.dump_var, args.dump_filtered, args.dump_mat, args.print_mat, dir_phase)

    result_fname = os.path.join(dir_phase, args.prefix + '.result.txt')
    switch_fname = os.path.join(dir_phase, args.prefix + '.switch')
    finalize(result_fname, hap2_bfname, collapsed_bfname, minced_fname, switch_fname, args.prefix, args.dir, pair_fname, filtered_paf_file)

def run_phase_shasta(args):
    # phasing for shasta

    args.shasta_fname = os.path.abspath(args.shasta_fname)
    # if args.porec is True, then args.porec_fname must be not None
    if args.porec:
        if args.hic1fname is not None or args.hic2fname is not None:
            logger.warning('Hi-C mapping files will not be used in porec process!')
        if args.porec_fname is None:
            logger.error('Please specify the porec filename!')
            exit(1)
        args.porec_fname = os.path.abspath(args.porec_fname)
    elif args.hic1fname is not None and args.hic2fname is not None:
        args.hic1fname = os.path.abspath(args.hic1fname)
        args.hic2fname = os.path.abspath(args.hic2fname)
    else:
        logger.error('Please specify the Hi-C mapping files!')
        exit(1)
    # if args.model is not None, then args.rdfname and args.type must be not None
    if args.model is not None:
        if args.rdfname is None:
            logger.error('Please specify the reference filename!')
            exit(1)
        if args.type is None:
            logger.error('Please specify the reference type!')
            exit(1)
        args.rdfname = os.path.abspath(args.rdfname)
    
    if args.threads < 1:
        logger.warning('Thread number must be greater than 0! Set to 1.')
        args.thread = 1
    if args.iter < 10000:
        logger.warning('Iteration number must be greater than 10000! Set to 10000.')
        args.iter = 10000
    if args.mapq < 1:
        logger.info('mapq {} should be >= 1, set to 1'.format(args.mapq))
        args.mapq = 1
    
    log_arguments(args)

    logger.info('{} Working directory {}'.format(sys._getframe().f_code.co_name, os.getcwd()))

    try:
        os.chdir(args.dir)
    except FileNotFoundError:
        logger.error('Directory {} not exist!'.format(args.dir), exc_info = 1)
        exit(1)
    except PermissionError:
        logger.error('No permissions to change to {}!'.format(args.dir), exc_info = 1)
        exit(1)
    
    args.matfname = os.path.abspath(args.matfname) if args.matfname is not None else None
    args.varfname = os.path.abspath(args.varfname) if args.varfname is not None else None

    
    # sub directory to store cut and bwa result
    dir_cut = os.path.join(args.dir, 'minced')
    # sub directory to store clair3 result
    dir_clair3 = os.path.join(args.dir, 'clair3')
    # sub directory to store phasing result
    dir_phase = os.path.join(args.dir, 'phase')

    # clean primary and alt, results store in ${args.dir}
    # clean_unzip(args.pri_fname, args.alt_fname, args.prefix, args.dir)

    # mapping alt to primary, results store in ${args.dir}
    # clean_pri_fname = os.path.join(args.dir, args.prefix + '.p_tigs.fasta')
    # clean_alt_fname = os.path.join(args.dir, args.prefix + '.a_tigs.fasta')
    # name_mapping_fname = os.path.join(args.dir, 'name_mapping.txt')
    # alt2pri_fname = os.path.join(args.dir, args.prefix + '.alt2pri')
    # make_mapping_alt_to_primary(clean_pri_fname, clean_alt_fname, name_mapping_fname, alt2pri_fname, args.threads, 'paf', args.dir)

    # place alt along with primary, results store in ${args.dir}
    # place_alt_along_primary_and_filter_contained(alt2pri_fname + '.paf', args.prefix, args.dir)

    # cut primary, results store in ${dir_cut}
    if os.path.exists(dir_cut):
        logger.warning('The sub directory minced exist in {}, some files will be overwritten!'.format(args.dir))
    else:
        os.makedirs(dir_cut, exist_ok = True)
    # filtered_placed_file = os.path.join(args.dir, args.prefix + '.filtered.alt2pri.txt')
    # cut_primary(clean_pri_fname, clean_alt_fname, filtered_placed_file, args.prefix, dir_cut)
    shasta_cut(args.shasta_fname, args.min_len, args.threads, args.prefix, dir_cut)

    # mince, results store in ${dir_cut}
    hap1_fname = os.path.join(dir_cut, args.prefix + '.hap1.cut.fasta')
    hap1_bfname = os.path.join(dir_cut, args.prefix + '.hap1.cut.bed')
    hap2_fname = os.path.join(dir_cut, args.prefix + '.hap2.cut.fasta')
    hap2_bfname = os.path.join(dir_cut, args.prefix + '.hap2.cut.bed')
    collapsed_fname = os.path.join(dir_cut, args.prefix + '.collapsed.cut.fasta')
    collapsed_bfname = os.path.join(dir_cut, args.prefix + '.collapsed.cut.bed')
    minced_fname = os.path.join(dir_cut, args.prefix + '.minced.fasta')
    clean_pri_fname = os.path.join(dir_cut, args.prefix + '.primary.fasta')
    clean_alt_fname = os.path.join(dir_cut, args.prefix + '.alt.fasta')
    mince_contigs(hap1_fname, hap2_fname, collapsed_fname, minced_fname, dir_cut)

    # run bwa, results store in ${dir_cut}
    if not args.porec:
        hic_mapping_fname = os.path.join(dir_cut, args.prefix + '.bwa.filtered2316.bam')
        run_bwa(minced_fname, args.hic1fname, args.hic2fname, hic_mapping_fname, args.threads, dir_cut)
    else:
        # porec process
        hic_mapping_fname = ''

    # run clair3 if model provided, results store in ${dir_clair3}
    if args.model is not None:
        if os.path.exists(dir_clair3):
            logger.warning('The sub directory clair3 exist in {}, some files may be overwritten!'.format(args.dir))
        else:
            os.makedirs(dir_clair3, exist_ok = True)

        run_clair3(clean_pri_fname, clean_alt_fname, args.rdfname, args.model, args.type, args.threads, args.prefix, dir_clair3)
    else:
        logger.warning('No model provided, clair3 will not run. the accuracy of phasing will be affected!')
    
    # run phasing, results store in ${dir_phase}
    if os.path.exists(dir_phase):
        logger.warning('The sub directory clair3 exist in {}, some files will be overwritten!'.format(args.dir))
    else:
        os.makedirs(dir_phase, exist_ok = True)
    
    pair_fname = os.path.join(dir_cut, args.prefix + '.pair.cut.txt')
    fai_fname = os.path.join(dir_cut, args.prefix + '.minced.fasta.fai')
    if args.model is not None:
        pri_clair3_fname = os.path.join(dir_clair3, 'clair3pri/merge_output.vcf.gz')
        alt_clair3_fname = os.path.join(dir_clair3, 'clair3alt/merge_output.vcf.gz')
    else:
        pri_clair3_fname, alt_clair3_fname = None, None
    
    # filtered_paf_file = os.path.join(args.dir, args.prefix + '.filtered.alt2pri.paf')
    filtered_paf_file = os.path.join(dir_cut, args.prefix + '.filtered.alt2pri.paf')

    run_phasing(pair_fname, fai_fname, pri_clair3_fname, alt_clair3_fname, filtered_paf_file, hic_mapping_fname, 
            args.matfname, args.varfname, args.threads, args.iter, args.seed, args.mapq, args.prefix, 
            args.filtered, args.porec, args.dump_var, args.dump_filtered, args.dump_mat, args.print_mat, dir_phase)

    result_fname = os.path.join(dir_phase, args.prefix + '.result.txt')
    switch_fname = os.path.join(dir_phase, args.prefix + '.switch')
    finalize(result_fname, hap2_bfname, collapsed_bfname, minced_fname, switch_fname, args.prefix, args.dir, pair_fname, filtered_paf_file)

def run_phase_dual(args):
    # phasing for dual

    args.dual1_fname = os.path.abspath(args.dual1_fname)
    args.dual2_fname = os.path.abspath(args.dual2_fname)
    # if args.porec is True, then args.porec_fname must be not None
    if args.porec:
        if args.hic1fname is not None or args.hic2fname is not None:
            logger.warning('Hi-C mapping files will not be used in porec process!')
        if args.porec_fname is None:
            logger.error('Please specify the porec filename!')
            exit(1)
        args.porec_fname = os.path.abspath(args.porec_fname)
    elif args.hic1fname is not None and args.hic2fname is not None:
        args.hic1fname = os.path.abspath(args.hic1fname)
        args.hic2fname = os.path.abspath(args.hic2fname)
    else:
        logger.error('Please specify the Hi-C mapping files!')
        exit(1)
    # if args.model is not None, then args.rdfname and args.type must be not None
    if args.model is not None:
        if args.rdfname is None:
            logger.error('Please specify the reference filename!')
            exit(1)
        if args.type is None:
            logger.error('Please specify the reference type!')
            exit(1)
        args.rdfname = os.path.abspath(args.rdfname)
    
    if args.threads < 1:
        logger.warning('Thread number must be greater than 0! Set to 1.')
        args.thread = 1
    if args.iter < 10000:
        logger.warning('Iteration number must be greater than 10000! Set to 10000.')
        args.iter = 10000
    if args.threshold < 10000:
        logger.info('mapq {} should be >= 10000, set to 10000'.format(args.mapq))
        args.mapq = 10000
    
    log_arguments(args)

    logger.info('{} Working directory {}'.format(sys._getframe().f_code.co_name, os.getcwd()))

    try:
        os.chdir(args.dir)
    except FileNotFoundError:
        logger.error('Directory {} not exist!'.format(args.dir), exc_info = 1)
        exit(1)
    except PermissionError:
        logger.error('No permissions to change to {}!'.format(args.dir), exc_info = 1)
        exit(1)
    
    args.matfname = os.path.abspath(args.matfname) if args.matfname is not None else None
    args.varfname = os.path.abspath(args.varfname) if args.varfname is not None else None

    
    # sub directory to store cut and bwa result
    dir_cut = os.path.join(args.dir, 'minced')
    # sub directory to store clair3 result
    dir_clair3 = os.path.join(args.dir, 'clair3')
    # sub directory to store phasing result
    dir_phase = os.path.join(args.dir, 'phase')

    # cut dual, results store in ${dir_cut}
    if os.path.exists(dir_cut):
        logger.warning('The sub directory minced exist in {}, some files will be overwritten!'.format(args.dir))
    else:
        os.makedirs(dir_cut, exist_ok = True)
    
    dual_contig_cut(args.dual1_fname, args.dual2_fname, args.threshold, args.threads, args.prefix, dir_cut)

    # mince, results store in ${dir_cut}
    hap1_fname = os.path.join(dir_cut, args.prefix + '.hap1.cut.fasta')
    hap1_bfname = os.path.join(dir_cut, args.prefix + '.hap1.cut.bed')
    hap2_fname = os.path.join(dir_cut, args.prefix + '.hap2.cut.fasta')
    hap2_bfname = os.path.join(dir_cut, args.prefix + '.hap2.cut.bed')
    collapsed1_fname = os.path.join(dir_cut, args.prefix + '.collapsed1.cut.fasta')
    collapsed1_bfname = os.path.join(dir_cut, args.prefix + '.collapsed1.cut.bed')
    collapsed2_fname = os.path.join(dir_cut, args.prefix + '.collapsed2.cut.fasta')
    collapsed2_bfname = os.path.join(dir_cut, args.prefix + '.collapsed2.cut.bed')
    minced_fname = os.path.join(dir_cut, args.prefix + '.minced.fasta')
    
    mince_contigs(hap1_fname, hap2_fname, collapsed1_fname, minced_fname, dir_cut, collapsed2_fname)

    # run bwa, results store in ${dir_cut}
    if not args.porec:
        hic_mapping_fname = os.path.join(dir_cut, args.prefix + '.bwa.filtered2316.bam')
        run_bwa(minced_fname, args.hic1fname, args.hic2fname, hic_mapping_fname, args.threads, dir_cut)
    else:
        # porec process
        hic_mapping_fname = ''

    # run clair3 if model provided, results store in ${dir_clair3}
    if args.model is not None:
        if os.path.exists(dir_clair3):
            logger.warning('The sub directory clair3 exist in {}, some files may be overwritten!'.format(args.dir))
        else:
            os.makedirs(dir_clair3, exist_ok = True)

        run_clair3(args.dual1_fname, args.dual2_fname, args.rdfname, args.model, args.type, args.threads, args.prefix, dir_clair3)
    else:
        logger.warning('No model provided, clair3 will not run. the accuracy of phasing will be affected!')
    
    # run phasing, results store in ${dir_phase}
    if os.path.exists(dir_phase):
        logger.warning('The sub directory clair3 exist in {}, some files will be overwritten!'.format(args.dir))
    else:
        os.makedirs(dir_phase, exist_ok = True)
    
    pair_fname = os.path.join(dir_cut, args.prefix + '.pair.cut.txt')
    fai_fname = os.path.join(dir_cut, args.prefix + '.minced.fasta.fai')
    if args.model is not None:
        pri_clair3_fname = os.path.join(dir_clair3, 'clair3pri/merge_output.vcf.gz')
        alt_clair3_fname = os.path.join(dir_clair3, 'clair3alt/merge_output.vcf.gz')
    else:
        pri_clair3_fname, alt_clair3_fname = None, None
    
    filtered_paf_file = os.path.join(dir_cut, args.prefix + '.dual.filtered.paf')

    run_phasing(pair_fname, fai_fname, pri_clair3_fname, alt_clair3_fname, filtered_paf_file, hic_mapping_fname, 
            args.matfname, args.varfname, args.threads, args.iter, args.seed, args.mapq, args.prefix, 
            args.filtered, args.porec, args.dump_var, args.dump_filtered, args.dump_mat, args.print_mat, dir_phase)

    result_fname = os.path.join(dir_phase, args.prefix + '.result.txt')
    switch_fname = os.path.join(dir_phase, args.prefix + '.switch')
    finalize_dual(result_fname, hap1_bfname, collapsed1_bfname, hap2_bfname, collapsed2_bfname, minced_fname, switch_fname, args.prefix, args.dir, pair_fname, filtered_paf_file)

def phasing_lachesis(pair_fname, fai_fname, dir, bed1_fname, bed2_fname, 
                     cluster_fname, hic_fname, iter, seed, prefix, dump_mat, print_mat):
    # pair_fname     pair filename
    # fai_fname      fai filename
    # dir            working directory
    # bed1_fname     bed filename of hap1
    # bed2_fname     bed filename of hap2
    # cluster_fname  cluster filename
    # hic_fname      hic mapping filename
    # iter           iteration to phase
    # seed           seed for random function
    # prefix         prefix of output files
    # dump_mat       dump the matrix into ${prefix}.binmat
    # print_mat      print the matrix to stdout

    start = time.time()

    # final file is ${prefix}.result.txt
    final_file = prefix + '.result.txt'
    if _test_file(dir, final_file): return

    _file_exist(pair_fname, fai_fname, bed1_fname, bed2_fname, cluster_fname, hic_fname)

    _check_program('phasing')

    try:
        os.chdir(dir)
    except FileNotFoundError:
        logger.error('Directory {} not exists!'.format(dir), exc_info = 1)
        exit(1)
    except PermissionError:
        logger.error('No permission to change to directory {}'.format(dir), exc_info = 1)
        exit(1)

    logger.info('{} Working directory {}'.format(sys._getframe().f_code.co_name, os.getcwd()))
    
    cmd = 'phasing -a ' + pair_fname + ' -c ' + fai_fname + ' --bed1 ' + bed1_fname 
    cmd += ' --bed2 ' + bed2_fname + ' -i ' + hic_fname + ' -p ' + prefix + ' -n ' + str(iter)
    cmd += ' -s ' + str(seed) + ' --cluster ' + cluster_fname + ' --scaffold'
    if dump_mat:
        cmd += ' --dump-mat'
    if print_mat:
        cmd += ' --print-mat > ' + prefix + '.snp.count.txt'

    try:
        subprocess.call(cmd, shell = True, stderr = open('phasing.log', 'w'))
    except (subprocess.CalledProcessError, OSError) as e:
        logger.error('Error running phasing {}'.format(e.cmd), exc_info = 1)
        exit(1)

    end = time.time()
    logger.info('{} ends'.format(sys._getframe().f_code.co_name))
    logger.info('Wall clock time elapsed: {} s'.format(end - start))
    logger.info('CPU time elapsed: {} s'.format(time.process_time()))
    logger.info('Memory usage: {:>.3f} GB\n'.format(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024 / 1024))

def finalize_lachesis(mfname, ldir, dir, rfname, prefix):
    # mfname    minced filename
    # ldir      lachesis directory
    # dir       working directory
    # rfname    phasing result filename
    # prefix    prefix of output files

    start = time.time()

    # final files are ${prefix}.phased0.fasta ${prefix}.phased1.fasta ${prefix}.phased0.bed ${prefix}.phased1.bed
    hap1_fname, hap1_bfname = prefix + 'phased0.fasta', prefix + 'phased0.bed'
    hap2_fname, hap2_bfname = prefix + 'phased1.fasta', prefix + 'phased1.bed'
    if _test_file(dir, hap1_fname, hap2_fname, hap1_bfname, hap2_bfname): return

    _file_exist(mfname, rfname)

    try:
        os.chdir(dir)
    except FileNotFoundError:
        logger.error('Directory {} not exists!'.format(dir), exc_info = 1)
        exit(1)
    except PermissionError:
        logger.error('No permission to change to directory {}'.format(dir), exc_info = 1)
        exit(1)

    logger.info('{} Working directory {}'.format(sys._getframe().f_code.co_name, os.getcwd()))

    create_scaffolding(mfname, ldir, rfname, prefix)

    end = time.time()
    logger.info('{} ends'.format(sys._getframe().f_code.co_name))
    logger.info('Wall clock time elapsed: {} s'.format(end - start))
    logger.info('CPU time elapsed: {} s'.format(time.process_time()))
    logger.info('Memory usage: {:>.3f} GB\n'.format(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024 / 1024))

def run_phase_lachesis(args):

    if args.iter < 100:
        logger.info('iter {} should be >= 100, set to 100'.format(args.iter))
        args.iter = 100
    if args.ldir is None or not os.path.isdir(args.ldir):
        logger.error('Please specify the lachesis directory!')
        exit(1)

    log_arguments(args)
    try:
        os.chdir(args.dir)
    except FileNotFoundError:
        logger.error('Directory {} not exists!'.format(args.dir), exc_info = 1)
        exit(1)
    except PermissionError:
        logger.error('No permission to change to directory {}'.format(args.dir), exc_info = 1)
        exit(1)

    logger.info('{} Working directory {}'.format(sys._getframe().f_code.co_name, os.getcwd()))

    dir_scaffold = os.path.join(args.dir, 'scaffold')
    if os.path.exists(dir_scaffold):
        logger.warning('The sub directory scaffold exists in {}, some files may be overwritten!'.format(args.dir))
    else:
        os.makedirs(dir_scaffold, exist_ok = True)

    # phasing scaffolding
    phasing_lachesis(args.pair_fname, args.fai_fname, dir_scaffold, args.bed1_fname, args.bed2_fname, 
                    args.cluster_fname, args.hic_fname, args.iter, args.seed, args.prefix, 
                    args.dump_mat, args.print_mat)
    # finalize scaffolding
    result_fname = os.path.join(dir_scaffold, args.prefix + '.result.txt')
    finalize_lachesis(args.mfname, args.ldir, args.dir, result_fname, args.prefix)

def main():
    parser = argparse.ArgumentParser('pipeline for phasing')
    subparser = parser.add_subparsers(title = 'subcommand')

    parser_pecat = subparser.add_parser('pecat', help = 'phasing PECAT assembly')
    parser_pecat.add_argument('--pri_asm', dest = 'pri_asm_fname', type = str, 
                              required = True, default = None, metavar = 'file', 
                              help = 'primary assembly filename')
    parser_pecat.add_argument('--pri_pol', dest = 'pri_pol_fname', type = str, 
                              required = True, default = None, metavar = 'file', 
                              help = 'primary polish filename')
    parser_pecat.add_argument('--alt_asm', dest = 'alt_asm_fname', type = str, 
                              required = True, default = None, metavar = 'file', 
                              help = 'alt assembly filename')
    parser_pecat.add_argument('--alt_pol', dest = 'alt_pol_fname', type = str, 
                              required = True, default = None, metavar = 'file', 
                              help = 'alt polish filename')
    parser_pecat.add_argument('--rdfname', dest = 'rdfname', type = str, default = None, 
                              metavar = 'file', help = 'raw reads filename')
    parser_pecat.add_argument('--model', dest = 'model', type = str, default = None, 
                              metavar = 'file', help = 'model filename used by clair3')
    parser_pecat.add_argument('--hic1', dest = 'hic1fname', type = str, default = None, 
                              metavar = 'file', help = 'pair 1 filename of Hi-C')
    parser_pecat.add_argument('--hic2', dest = 'hic2fname', type = str, default = None, 
                              metavar = 'file', help = 'pair 2 filename of Hi-C')
    parser_pecat.add_argument('--porec_fname', dest = 'porec_fname', type = str, default = None, 
                              metavar = 'file', help = 'Pore-C filename')
    parser_pecat.add_argument('--mat_fname', dest = 'matfname', type = str, default = None, 
                              metavar = 'file', help = 'matrix filename')
    parser_pecat.add_argument('--var_fname', dest = 'varfname', type = str, default = None, 
                              metavar = 'file', help = 'snp filename')
    
    parser_pecat.add_argument('-d', '--dir', dest = 'dir', type = str, required = True, 
                              default = '.', metavar = 'file', help = 'working directory')
    parser_pecat.add_argument('-p', '--prefix', dest = 'prefix', type = str, default = 'phasing', 
                              metavar = 'str', help = 'prefix of output files')
    parser_pecat.add_argument('-t', '--threads', dest = 'threads', type = int, default = 1, 
                              metavar = 'int', help = 'number of threads to use')
    parser_pecat.add_argument('-q', '--mapq', dest = 'mapq', type = int, default = 1,
                              metavar = 'int', help = 'mapq cutoff for hic mapping')
    parser_pecat.add_argument('--type', dest = 'type', choices = ['clr', 'hifi', 'ont'], 
                              metavar = 'str', help = 'type of the raw sequence')
    parser_pecat.add_argument('--iter', dest = 'iter', type = int, default = 10000, 
                              metavar = 'int', help = 'number of iteration to phase')
    parser_pecat.add_argument('--seed', dest = 'seed', type = int, default = 10000, 
                              metavar = 'int', help = 'seed for random function')
    parser_pecat.add_argument('--filtered', dest = 'filtered', action = 'store_true', 
                              help = 'whether the hic mapping are filtered')
    parser_pecat.add_argument('--porec', dest = 'porec', action = 'store_true', 
                              help = 'whether the reads is porec')
    parser_pecat.add_argument('--dump_var', dest = 'dump_var', action = 'store_true', 
                              help = 'dump snp into ${prefix}.var.txt')
    parser_pecat.add_argument('--dump_filtered', dest = 'dump_filtered', action = 'store_true', 
                              help = 'dump filtered hic mapping into ${prefix}.hic.filtered.bam')
    parser_pecat.add_argument('--dump_mat', dest = 'dump_mat', action = 'store_true', 
                              help = 'dump the matrix into ${prefix}.binmat')
    parser_pecat.add_argument('--print_mat', dest = 'print_mat', action = 'store_true', 
                              help = 'print the matrix to stdout')
    
    parser_pecat.set_defaults(func = run_phase_pecat)

    parser_phase = subparser.add_parser('phase', help = 'phasing unzip format assembly')
    parser_phase.add_argument('--pri', dest = 'pri_fname', type = str, 
                              required = True, default = None, metavar = 'file', 
                              help = 'primary assembly filename')
    parser_phase.add_argument('--alt', dest = 'alt_fname', type = str, 
                              required = True, default = None, metavar = 'file', 
                              help = 'alt assembly filename')
    parser_phase.add_argument('--rdfname', dest = 'rdfname', type = str, default = None, 
                              metavar = 'file', help = 'raw reads filename')
    parser_phase.add_argument('--model', dest = 'model', type = str, default = None, 
                              metavar = 'file', help = 'model filename used by clair3')
    parser_phase.add_argument('--hic1', dest = 'hic1fname', type = str, default = None, 
                              metavar = 'file', help = 'pair 1 filename of Hi-C')
    parser_phase.add_argument('--hic2', dest = 'hic2fname', type = str, default = None, 
                              metavar = 'file', help = 'pair 2 filename of Hi-C')
    parser_phase.add_argument('--porec_fname', dest = 'porec_fname', type = str, default = None, 
                              metavar = 'file', help = 'Pore-C filename')
    parser_phase.add_argument('--mat_fname', dest = 'matfname', type = str, default = None, 
                              metavar = 'file', help = 'matrix filename')
    parser_phase.add_argument('--var_fname', dest = 'varfname', type = str, default = None, 
                              metavar = 'file', help = 'snp filename')
    
    parser_phase.add_argument('-d', '--dir', dest = 'dir', type = str, required = True, 
                              default = '.', metavar = 'file', help = 'working directory')
    parser_phase.add_argument('-p', '--prefix', dest = 'prefix', type = str, default = 'phasing', 
                              metavar = 'str', help = 'prefix of output files')
    parser_phase.add_argument('-t', '--threads', dest = 'threads', type = int, default = 1, 
                              metavar = 'int', help = 'number of threads to use')
    parser_phase.add_argument('-q', '--mapq', dest = 'mapq', type = int, default = 1,
                              metavar = 'int', help = 'mapq cutoff for hic mapping')
    parser_phase.add_argument('--type', dest = 'type', choices = ['clr', 'hifi', 'ont'], 
                              metavar = 'str', help = 'type of the raw sequence')
    parser_phase.add_argument('--iter', dest = 'iter', type = int, default = 10000, 
                              metavar = 'int', help = 'number of iteration to phase')
    parser_phase.add_argument('--seed', dest = 'seed', type = int, default = 10000, 
                              metavar = 'int', help = 'seed for random function')
    parser_phase.add_argument('--filtered', dest = 'filtered', action = 'store_true', 
                              help = 'whether the hic mapping are filtered')
    parser_phase.add_argument('--porec', dest = 'porec', action = 'store_true', 
                              help = 'whether the reads is porec')
    parser_phase.add_argument('--dump_var', dest = 'dump_var', action = 'store_true', 
                              help = 'dump snp into ${prefix}.var.txt')
    parser_phase.add_argument('--dump_filtered', dest = 'dump_filtered', action = 'store_true', 
                              help = 'dump filtered hic mapping into ${prefix}.hic.filtered.bam')
    parser_phase.add_argument('--dump_mat', dest = 'dump_mat', action = 'store_true', 
                              help = 'dump the matrix into ${prefix}.binmat')
    parser_phase.add_argument('--print_mat', dest = 'print_mat', action = 'store_true', 
                              help = 'print the matrix to stdout')
    
    parser_phase.set_defaults(func = run_phase)

    parser_shasta = subparser.add_parser('shasta', help = 'phasing shasta assembly')
    parser_shasta.add_argument('--assembly', dest = 'shasta_fname', type = str,
                                required = True, default = None, metavar = 'file',
                                help = 'shasta assembly filename')
    parser_shasta.add_argument('--rdfname', dest = 'rdfname', type = str, default = None,
                                metavar = 'file', help = 'raw reads filename')
    parser_shasta.add_argument('--model', dest = 'model', type = str, default = None,
                                metavar = 'file', help = 'model filename used by clair3')
    parser_shasta.add_argument('--hic1', dest = 'hic1fname', type = str, default = None,
                                metavar = 'file', help = 'pair 1 filename of Hi-C')
    parser_shasta.add_argument('--hic2', dest = 'hic2fname', type = str, default = None,
                                metavar = 'file', help = 'pair 2 filename of Hi-C')
    parser_shasta.add_argument('--porec_fname', dest = 'porec_fname', type = str, default = None,
                                metavar = 'file', help = 'Pore-C filename')
    parser_shasta.add_argument('--mat_fname', dest = 'matfname', type = str, default = None,
                                metavar = 'file', help = 'matrix filename')
    parser_shasta.add_argument('--var_fname', dest = 'varfname', type = str, default = None,
                                metavar = 'file', help = 'snp filename')
    
    parser_shasta.add_argument('--min_len', dest = 'min_len', type = int, default = 3000,
                                metavar = 'int', help = 'minimum length of contigs')
    parser_shasta.add_argument('-d', '--dir', dest = 'dir', type = str, required = True,
                                default = '.', metavar = 'file', help = 'working directory')
    parser_shasta.add_argument('-p', '--prefix', dest = 'prefix', type = str, default = 'phasing',
                                metavar = 'str', help = 'prefix of output files')
    parser_shasta.add_argument('-t', '--threads', dest = 'threads', type = int, default = 1,
                                metavar = 'int', help = 'number of threads to use')
    parser_shasta.add_argument('-q', '--mapq', dest = 'mapq', type = int, default = 1,
                                metavar = 'int', help = 'mapq cutoff for hic mapping')
    parser_shasta.add_argument('--type', dest = 'type', choices = ['clr', 'hifi', 'ont'],
                                metavar = 'str', help = 'type of the raw sequence')
    parser_shasta.add_argument('--iter', dest = 'iter', type = int, default = 10000,
                                metavar = 'int', help = 'number of iteration to phase')
    parser_shasta.add_argument('--seed', dest = 'seed', type = int, default = 10000,
                                metavar = 'int', help = 'seed for random function')
    parser_shasta.add_argument('--filtered', dest = 'filtered', action = 'store_true',
                                help = 'whether the hic mapping are filtered')
    parser_shasta.add_argument('--porec', dest = 'porec', action = 'store_true',
                                help = 'whether the reads is porec')
    parser_shasta.add_argument('--dump_var', dest = 'dump_var', action = 'store_true',
                                help = 'dump snp into ${prefix}.var.txt')
    parser_shasta.add_argument('--dump_filtered', dest = 'dump_filtered', action = 'store_true',
                                help = 'dump filtered hic mapping into ${prefix}.hic.filtered.bam')
    parser_shasta.add_argument('--dump_mat', dest = 'dump_mat', action = 'store_true',
                                help = 'dump the matrix into ${prefix}.binmat')
    parser_shasta.add_argument('--print_mat', dest = 'print_mat', action = 'store_true',
                                help = 'print the matrix to stdout')
    
    parser_shasta.set_defaults(func = run_phase_shasta)

    parser_dual = subparser.add_parser('dual', help = 'phasing dual assembly')
    parser_dual.add_argument('--dual1', dest = 'dual1_fname', type = str,
                                required = True, default = None, metavar = 'file',
                                help = 'dual assembly 1 filename')
    parser_dual.add_argument('--dual2', dest = 'dual2_fname', type = str,
                                required = True, default = None, metavar = 'file',
                                help = 'dual assembly 2 filename')
    parser_dual.add_argument('--rdfname', dest = 'rdfname', type = str, default = None,
                                metavar = 'file', help = 'raw reads filename')
    parser_dual.add_argument('--model', dest = 'model', type = str, default = None,
                                metavar = 'file', help = 'model filename used by clair3')
    parser_dual.add_argument('--hic1', dest = 'hic1fname', type = str, default = None,
                                metavar = 'file', help = 'pair 1 filename of Hi-C')
    parser_dual.add_argument('--hic2', dest = 'hic2fname', type = str, default = None,
                                metavar = 'file', help = 'pair 2 filename of Hi-C')
    parser_dual.add_argument('--porec_fname', dest = 'porec_fname', type = str, default = None,
                                metavar = 'file', help = 'Pore-C filename')
    parser_dual.add_argument('--mat_fname', dest = 'matfname', type = str, default = None,
                                metavar = 'file', help = 'matrix filename')
    parser_dual.add_argument('--var_fname', dest = 'varfname', type = str, default = None,
                                metavar = 'file', help = 'snp filename')
    
    parser_dual.add_argument('-d', '--dir', dest = 'dir', type = str, required = True,
                                default = '.', metavar = 'file', help = 'working directory')
    parser_dual.add_argument('-p', '--prefix', dest = 'prefix', type = str, default = 'phasing',
                                metavar = 'str', help = 'prefix of output files')
    parser_dual.add_argument('-t', '--threads', dest = 'threads', type = int, default = 1,
                                metavar = 'int', help = 'number of threads to use')
    parser_dual.add_argument('-q', '--mapq', dest = 'mapq', type = int, default = 1,
                                metavar = 'int', help = 'mapq cutoff for hic mapping')
    parser_dual.add_argument('-s', '--threshold', dest = 'threshold', type = int, default = 15000,
                                metavar = 'int', help = 'mapq cutoff for hic mapping')
    parser_dual.add_argument('--type', dest = 'type', choices = ['clr', 'hifi', 'ont'],
                                metavar = 'str', help = 'type of the raw sequence')
    parser_dual.add_argument('--iter', dest = 'iter', type = int, default = 10000,
                                metavar = 'int', help = 'number of iteration to phase')
    parser_dual.add_argument('--seed', dest = 'seed', type = int, default = 10000,
                                metavar = 'int', help = 'seed for random function')
    parser_dual.add_argument('--filtered', dest = 'filtered', action = 'store_true',
                                help = 'whether the hic mapping are filtered')
    parser_dual.add_argument('--porec', dest = 'porec', action = 'store_true',
                                help = 'whether the reads is porec')
    parser_dual.add_argument('--dump_var', dest = 'dump_var', action = 'store_true',
                                help = 'dump snp into ${prefix}.var.txt')
    parser_dual.add_argument('--dump_filtered', dest = 'dump_filtered', action = 'store_true',
                                help = 'dump filtered hic mapping into ${prefix}.hic.filtered.bam')
    parser_dual.add_argument('--dump_mat', dest = 'dump_mat', action = 'store_true',
                                help = 'dump the matrix into ${prefix}.binmat')
    parser_dual.add_argument('--print_mat', dest = 'print_mat', action = 'store_true',
                                help = 'print the matrix to stdout')
    
    parser_dual.set_defaults(func = run_phase_dual)

    parser_lachesis = subparser.add_parser('lachesis', help = 'phasing lachesis scaffold')
    parser_lachesis.add_argument('--pair', dest = 'pair_fname', type = str, 
                                required = True, default = None, metavar = 'file', 
                                help = 'pair filename')
    parser_lachesis.add_argument('--fai', dest = 'fai_fname', type = str, 
                                required = True, default = None, metavar = 'file', 
                                help = 'fai filename')
    parser_lachesis.add_argument('--bed1', dest = 'bed1_fname', type = str, 
                                required = True, default = None, metavar = 'file', 
                                help = 'bed filename of hap1')
    parser_lachesis.add_argument('--bed2', dest = 'bed2_fname', type = str,
                                required = True, default = None, metavar = 'file', 
                                help = 'bed filename of hap2')
    parser_lachesis.add_argument('--cluster', dest = 'cluster_fname', type = str,
                                required = True, default = None, metavar = 'file', 
                                help = 'cluster filename')
    parser_lachesis.add_argument('--hic', dest = 'hic_fname', type = str,
                                required = True, default = None, metavar = 'file', 
                                help = 'hic mapping filename')
    parser_lachesis.add_argument('--minced', dest = 'mfname', type = str,
                                required = True, default = None, metavar = 'file', 
                                help = 'minced filename')
    parser_lachesis.add_argument('--ldir', dest = 'ldir', type = str,
                                required = True, default = None, metavar = 'file', 
                                help = 'lachesis directory')
    parser_lachesis.add_argument('-d', '--dir', dest = 'dir', type = str, required = True, 
                                default = '.', metavar = 'file', help = 'working directory')
    parser_lachesis.add_argument('-p', '--prefix', dest = 'prefix', type = str, default = 'phasing',
                                metavar = 'str', help = 'prefix of output files')
    parser_lachesis.add_argument('--iter', dest = 'iter', type = int, default = 10000,
                                metavar = 'int', help = 'number of iteration to phase')
    parser_lachesis.add_argument('--seed', dest = 'seed', type = int, default = 10000,
                                metavar = 'int', help = 'seed for random function')
    parser_lachesis.add_argument('--dump_mat', dest = 'dump_mat', action = 'store_true',
                                help = 'dump the matrix into ${prefix}.binmat')
    parser_lachesis.add_argument('--print_mat', dest = 'print_mat', action = 'store_true',
                                help = 'print the matrix to stdout')
    parser_lachesis.set_defaults(func = run_phase_lachesis)

    args = parser.parse_args()

    args.dir = os.getcwd() if args.dir is None else os.path.abspath(args.dir)
    if not os.path.isdir(args.dir):
        try:
            os.makedirs(args.dir, exist_ok = True)
        except PermissionError:
            logger.error('No permissions to create directory {}'.format(args.dir), exc_info = 1)
            exit(1)

    log_fname = os.path.join(args.dir, 'pipeline.log')
    # _enable_logging(log_fname)
    setLogger(log_fname)

    args.func(args)

if __name__ == '__main__':
    main()

