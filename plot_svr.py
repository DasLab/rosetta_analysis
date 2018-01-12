#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none'
import sys

import argparse

def floatable(st):
    try:
        float(st)
        return True
    except:
        return False

def detect_named_col(lines, name_list, default):
    for line in lines:
        if 'description' in line:
            # This is the header line. Good! We can split it and find
            # our desired col.
            for name in name_list:
                if name in line.split():
                    return line.split().index(name)
    return default

def detect_score_col(filename):
    scoreidx = 1
    scoreidx = detect_named_col(filename, ["score", "total_score"], scoreidx)
    return scoreidx

def detect_rms_col(filename):
    rmsidx = 24
    rmsidx = detect_named_col(filename, ["rms", "rms_fill", "rms_stem"], rmsidx)
    return rmsidx


def plot(filename, scoreidx, rmsidx, dump, outfile_format): 
    scores = []
    rmsds = []
    lines = []

    with open(filename) as f:
        lines = [ line for line in f if line[0:5] == "SCORE"]
    
    if scoreidx == -1: scoreidx = detect_score_col(lines)
    if rmsidx == -1: rmsidx = detect_rms_col(lines)

    scores = [float(line.split()[scoreidx]) for line in lines if floatable(line.split()[scoreidx])]
    rmsds = [float(line.split()[rmsidx]) for line in lines if floatable(line.split()[rmsidx])]
    print len(scores), len(rmsds)
    plt.scatter(rmsds, scores)

    # Arbitrarily set ymax to be 0 to exclude FARFAR outliers.
    # Later develop a smart strategy having to do with actual outlier detection
    plt.ylim([min(scores)-10, max([score for score in scores if score < 0])+10])

    if dump:
        plt.savefig("%s.%s" % (filename, outfile_format))
    else:
        plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Make a score-vs.-RMSD plot from a Rosetta score or silent file.')
    parser.add_argument('filename', type=str, #nargs='1',
            help='the silent or score file to plot')
    parser.add_argument('--score_col', type=int, nargs='?',
            help='the column (zero-indexed) where the file contains score data, default detect (1 if failed)', default=-1)
    parser.add_argument('--rmsd_col', type=int, nargs='?',
            help='the column (zero-indexed) where the file contains rms data, default detect (24 if failed)', default=-1)
    parser.add_argument('--dump', type=bool, nargs='?',
            help='dump to file? if false, show the plot in-line. default true', default=True)
    parser.add_argument('--format', type=str, nargs='?',
            help='what format output file? default svg, can be png or pdf', default='svg')

    args = parser.parse_args()

    #plot(filename, scoreidx, rmsidx, dump, outfile_format)
    plot(args.filename, args.score_col, args.rmsd_col, args.dump, args.format)
