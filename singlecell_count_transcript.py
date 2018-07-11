#!/usr/bin/env python
'''
Copyright 2018 Matt Settles
Created July 11, 2018

Read in mapped reads with BWA to transcriptome and produces a counts table.
'''

import os.path
import sys
import time
import argparse
import traceback
import errno
from subprocess import Popen, PIPE, STDOUT
from collections import Counter
import numpy

def median(lst):
    return numpy.median(numpy.array(lst))

def sp_bam_read(file, bufsize=-1):
    if os.path.isfile(file):
        p = Popen('samtools view'.split() + [file], stdout=PIPE, stderr=STDOUT, bufsize=bufsize)
        return p.stdout
    else:
        sys.stderr.write("Cannot open file " + file + '\n')
## NEED TO FIX, does not work
def sp_bam_write(file, bufsize=-1):
    filep = open(file, 'wb')
    p = Popen('samtools view', stdin=PIPE, stdout=filep, shell=True, bufsize=bufsize)
    return p.stdin

def main(infile, output_dir, verbose):
    # Set up the global variables
    global read_count
    global stime

    read_assigned_fwd = 0
    read_assigned_rvs = 0
    read_umi_fwd = 0
    read_umi_rvs = 0
    read_ambiguous = 0
    read_unmapped = 0

    barcode_list = Counter()
    umi_Counter_fwd = {}
    umi_Counter_rvs = {}

    if infile == 'stdin':
        # reading from stdin
        insam = sys.stdin
    elif infile.split(".")[-1] == "bam":
        insam = sp_bam_read(infile)
    else:
        insam = open(infile, 'r')

    try:
        while 1:
            line = insam.next()
            # Comment/header lines start with @
            if line[0] != "@" and len(line.strip().split()) > 2:
                read_count += 1

                line2 = line.strip().split()
                # Handle TAG:
                # get the final concatenatd tag
                tag = line2[0]
                rid = tag.split(':')[0:2] # [cell_bc, umi]
                barcode_list[rid[0]] += 1

                flag = line2[1]
                if (flag & 0x100):  # secondary alignment, ignore
                    continue
                if (flag & 0x4): # segment unmapped
                    read_unmapped += 1
                    continue
                transcript_id = line2[2]
                gene_id = '_'.join(transcript_id.split("_")[:-1])

                # Look for multi-mappers
                # Looking for tags XA (From bwa)
                # XA : Alternative hits; format: (chr,pos,CIGAR,NM;)*
                xa = [x for x in line2 if x.startswith('XA')]
                if len(xa) != 0:
                    xa = xa[0][5:]
                    xa = xa.split(";")[:-1]
                    xa_transcript = [x.split(",")[0] for x in xa]
                    xa_gene = ['_'.join(x.split("_")[:-1]) for x in xa_transcript]
                    xa_gene.append(gene_id) # add the mapped Gene_ID
                    num_hits = len(set(xa_gene))
                    if num_hits > 1:
                    read_ambiguous +=1
                    continue
                if (flag & 0x10):
                    if gene_id in umi_Counter_rvs:
                        if rid[0] in umi_Counter_rvs[gene_id]:
                            if rid[1] in umi_Counter_rvs[gene_id][rid[0]]:
                                read_umi_fwd += 1
                            else:
                                read_assigned_fwd += 1
                            umi_Counter_rvs[gene_id][rid[0]][rid[1]] += 1
                        else:
                            umi_Counter_rvs[gene_id][rid[0]] = Counter([rid[1]])
                            read_assigned_fwd += 1
                    else:
                        umi_Counter_rvs[gene_id] = {}
                        umi_Counter_rvs[gene_id][rid[0]] = Counter([rid[1]])
                        read_assigned_fwd += 1
                else:
                    if gene_id in umi_Counter_rvs:
                        if rid[0] in umi_Counter_rvs[gene_id]:
                            if rid[1] in umi_Counter_rvs[gene_id][rid[0]]:
                                read_umi_rvs += 1
                            else:
                                read_assigned_rvs += 1
                            umi_Counter_rvs[gene_id][rid[0]][rid[1]] += 1
                        else:
                            umi_Counter_rvs[gene_id][rid[0]] = Counter([rid[1]])
                            read_assigned_rvs += 1
                    else:
                        umi_Counter_rvs[gene_id] = {}
                        umi_Counter_rvs[gene_id][rid[0]] = Counter([rid[1]])
                        read_assigned_rvs += 1

    except StopIteration:
        base = output_dir
        with open(output_dir + 'fwd_counts.txt', 'w') as f:
            bc_keys = barcode_list.keys()
            f.write('Gene_ID\t' + '\t'.join(bc_keys) + '\n')

            for gene in umi_Counter_fwd:
                txt = gene
                for bc in bc_keys:
                    if bc in umi_Counter_fwd[gene]:
                        txt = "\t".join([txt,str(len(umi_Counter_fwd[gene][bc]))])
                    else:
                        txt = "\t".join([txt,str(0)])
                f.write(txt + '\n')
        with open(output_dir + 'rvs_counts.txt', 'w') as f:
            bc_keys = barcode_list.keys()
            f.write('Gene_ID\t' + '\t'.join(bc_keys) + '\n')

            for gene in umi_Counter_rvs:
                txt = gene
                for bc in bc_keys:
                    if bc in umi_Counter_rvs[gene]:
                        txt = "\t".join([txt,str(len(umi_Counter_rvs[gene][bc]))])
                    else:
                        txt = "\t".join([txt,str(0)])
                f.write(txt + '\n')

        if verbose:
            sys.stderr.write("PROCESS\tREADS\treads analyzed:%i|reads/sec:%i|barcodes:%i|reads/barcode:%f\n" % (read_count, round(read_count / (time.time() - stime), 0), len(barcode_list), median(barcode_list.values())))
            sys.stderr.write("PROCESS\tREADS\tASSIGNED FWD: %i (%.2f%%)\n" % (read_assigned_fwd, (float(read_assigned_fwd) / read_count) * 100))
            sys.stderr.write("PROCESS\tREADS\tUMI FWD: %i (%.2f%%)\n" % (read_umi_fwd, (float(read_umi_fwd) / read_count) * 100))
            sys.stderr.write("PROCESS\tREADS\tASSIGNED RVS: %i (%.2f%%)\n" % (read_assigned_rvs, (float(read_assigned_rvs) / read_count) * 100))
            sys.stderr.write("PROCESS\tREADS\tUMI RVS: %i (%.2f%%)\n" % (read_umi_rvs, (float(read_umi_rvs) / read_count) * 100))
            sys.stderr.write("PROCESS\tREADS\tAMBIGUOUS: %i (%.2f%%)\n" % (read_ambiguous, (float(read_ambiguous) / read_count) * 100))
            sys.stderr.write("PROCESS\tREADS\tUNMAPPED: %i (%.2f%%)\n" % (read_unmapped, (float(read_unmapped / read_count) * 100))
        pass
    except (KeyboardInterrupt, SystemExit):
        sys.exit("PROCESS\tERROR\t%s unexpectedly terminated\n" % (__name__))
    except Exception:
        sys.stderr.write("".join(traceback.format_exception(*sys.exc_info())))
        sys.exit("PROCESS\tERROR\tAn unknown fatal error was encountered.\n")


#####################################
# Parse options and setup #
version_num = "0.0.1"
parser = argparse.ArgumentParser(description='singlecell_count.py, to process mapped and assigned reads, outputing a counts table',
                                 epilog='For questions or comments, please contact Matt Settles <settles@ucdavis.edu>\n%(prog)s version: ' + version_num, add_help=True)
parser.add_argument('--version', action='version', version="%(prog)s version: " + version_num)

parser.add_argument('--quiet', help="turn off verbose output",
                    action="store_false", dest="verbose", default=True)

parser.add_argument('-o', '--output', help="Directory + prefix to output reads, [default: %(default)s]",
                    action="store", type=str, dest="output_dir", default="stdout")

group = parser.add_argument_group("Inputs", "sam/bam files to input.")

group.add_argument('-s', '--sam', metavar="sam", dest='sam', help='sam file (can also be bam), can also use stdin',
                   action='store', type=str)

options = parser.parse_args()

output_dir = options.output_dir

infile = options.sam
if infile is None:
    sys.stderr.write("PROCESS\tERROR\tRead sam file is missing\n")
    sys.exit(1)

verbose = options.verbose

# need to check, can write to output folder

# global variables
read_count = 0

stime = time.time()

main(infile, output_dir, verbose)

sys.exit(0)
