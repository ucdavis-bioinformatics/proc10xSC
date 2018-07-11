#!/usr/bin/env python
'''
Copyright 2018 Matt Settles
Created June 28, 2018

Read in mapped reads with STAR assigned to gene with FeatureCounts and produce counts table.
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

    read_assigned = 0
    read_umi = 0
    read_ambiguous = 0
    read_unknown = 0

    barcode_list = Counter()
    umi_Counter = {}

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
                # Looking for tags XS, XN and XT
                # XS: Assignment status
                # XN: Number of targets
                # XT: Comma separated target list
                xs = [x for x in line2 if x.startswith('XS')]

                if xs[0][5:] != 'Assigned':
                    read_unknown += 1
                    continue
                xn = [x for x in line2 if x.startswith('XN')]
                if int(xn[0][5:]) != 1:
                    read_ambiguous +=1
                    continue
                xt = [x for x in line2 if x.startswith('XT')]
                xt = xt[0][5:].split(',')[0]

                if xt in umi_Counter:
                    if rid[0] in umi_Counter[xt]:
                        if rid[1] in umi_Counter[xt][rid[0]]:
                            read_umi += 1
                        else:
                            read_assigned += 1
                        umi_Counter[xt][rid[0]][rid[1]] += 1
                    else:
                        umi_Counter[xt][rid[0]] = Counter([rid[1]])
                        read_assigned += 1
                else:
                    umi_Counter[xt] = {}
                    umi_Counter[xt][rid[0]] = Counter([rid[1]])
                    read_assigned += 1

    except StopIteration:
        with open(output_dir + 'counts.txt', 'w') as f:
            bc_keys = barcode_list.keys()
            f.write('Gene_ID\t' + '\t'.join(bc_keys) + '\n')

            for gene in umi_Counter:
                txt = gene
                for bc in bc_keys:
                    if bc in umi_Counter[gene]:
                        txt = "\t".join([txt,str(len(umi_Counter[gene][bc]))])
                    else:
                        txt = "\t".join([txt,str(0)])
                f.write(txt + '\n')

        if verbose:
            sys.stderr.write("PROCESS\tREADS\treads analyzed:%i|reads/sec:%i|barcodes:%i|reads/barcode:%f\n" % (read_count, round(read_count / (time.time() - stime), 0), len(barcode_list), median(barcode_list.values())))
            sys.stderr.write("PROCESS\tREADS\tASSIGNED: %i (%.2f%%)\n" % (read_assigned, (float(read_assigned) / read_count) * 100))
            sys.stderr.write("PROCESS\tREADS\tUMI: %i (%.2f%%)\n" % (read_umi, (float(read_umi) / read_count) * 100))
            sys.stderr.write("PROCESS\tREADS\tAMBIGUOUS: %i (%.2f%%)\n" % (read_ambiguous, (float(read_ambiguous) / read_count) * 100))
            sys.stderr.write("PROCESS\tREADS\tUNKNOWN: %i (%.2f%%)\n" % (read_unknown, (float(read_unknown) / read_count) * 100))
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
