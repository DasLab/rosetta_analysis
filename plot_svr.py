#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none'
import sys
import mypy
from typing import List
import seaborn as sns

import argparse

def floatable(st: str) -> bool:
    try:
        float(st)
        return True
    except:
        return False

def detect_named_col(lines: List[str], name_list: List[str], default: int) -> int:
    for line in lines:
        if 'description' in line:
            # This is the header line. Good! We can split it and find
            # our desired col.
            for name in name_list:
                if name in line.split():
                    return line.split().index(name)
    return default

def detect_score_col(filename: str) -> int:
    scoreidx = 1
    scoreidx = detect_named_col(filename, ["score", "total_score", "silent_score"], scoreidx)
    return scoreidx

def detect_rms_col(filename: str) -> int:
    rmsidx = 24
    rmsidx = detect_named_col(filename, ["new_rms", "rms_fill", "rms",  "rms_stem"], rmsidx)
    return rmsidx


def plot(filename: List[str], scoreidx: int, rmsidx: int, dump: bool, outfile_format: str, xmax: float) -> None: 
    scores = []
    rmsds = []
    lines_arr = []

    for fn in filename:
        with open(fn) as f:
            lines_arr.append([ line for line in f if line[0:5] == "SCORE"])
   
    mins, maxes = [], []
    for lines in lines_arr:
        if scoreidx == -1: scoreidx = detect_score_col(lines)
        if rmsidx == -1: rmsidx = detect_rms_col(lines)

        scores = [float(line.split()[scoreidx]) for line in lines if floatable(line.split()[scoreidx])]
        rmsds = [float(line.split()[rmsidx]) for line in lines if floatable(line.split()[rmsidx])]
        print(len(scores), len(rmsds))

        # How large should markersize be? I'd say, start at 25 and reduce every order of magnitude.
        markersize = (6 - np.log10(len(rmsds))) ** 2
        #print((scores), (rmsds))
        plt.scatter(rmsds, scores, s=markersize, alpha=0.5)
        #sns.kdeplot(rmsds, scores, cmap="Blues", shade=True)
        
        pctilescore = sorted(scores)[int(.95*len(scores))]
        print(pctilescore)
        mins.append(min(scores))
        maxes.append(min(max([score for score in scores if score <= pctilescore]), min(scores)+100))


    # Arbitrarily set ymax to be 0 to exclude FARFAR outliers.
    # Later develop a smart strategy having to do with actual outlier detection
    # Strategy: get range of min and max, expand by 5% either way.
    ymin, ymax = min(mins), max(maxes)
    ymin2, ymax2 = ymin - (ymax - ymin)*0.05, ymax + (ymax - ymin)*0.05
    #plt.ylim([min(mins)-10, max(maxes)+10])
    plt.ylim([ymin2, ymax2])
    #plt.show()
    #plt.ylim([-15,-3])
    plt.xlim([0, xmax])

    if dump:
        plt.savefig("%s.%s" % ("_".join(filename), outfile_format))
    else:
        plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Make a score-vs.-RMSD plot from a Rosetta score or silent file.')
    parser.add_argument('filename', type=str, nargs='+', # any nargs ok -- many series
            help='the silent or score file to plot')
    parser.add_argument('--score_col', type=int, nargs='?',
            help='the column (zero-indexed) where the file contains score data, default detect (1 if failed)', default=-1)
    parser.add_argument('--rmsd_col', type=int, nargs='?',
            help='the column (zero-indexed) where the file contains rms data, default detect (24 if failed)', default=-1)
    parser.add_argument('--dump', type=bool, nargs='?',
            help='dump to file? if false, show the plot in-line. default true', default=True)
    parser.add_argument('--format', type=str, nargs='?',
            help='what format output file? default svg, can be png or pdf', default='svg')
    parser.add_argument('--xmax', type=float, nargs='?',
            help='xmax? default 16.', default=16)

    args = parser.parse_args()

    #plot(filename, scoreidx, rmsidx, dump, outfile_format)
    plot(args.filename, args.score_col, args.rmsd_col, args.dump, args.format, args.xmax)
