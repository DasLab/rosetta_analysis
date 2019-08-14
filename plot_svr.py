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
    """
    Allows filtering column contents by numeric-ness.
    """
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

def detect_score_col(lines: List[str]) -> int:
    """
    From the given lines, find the score col.
    """
    return detect_named_col(lines, ["score", "total_score", "silent_score"], 1)

def detect_rms_col(lines: List[str]) -> int:
    """
    From the given lines, find the rms col.
    """
    return detect_named_col(lines, ["new_rms", "rms_fill", "rms",  "rms_stem"], 24)

def output_name(filenames: List[str]):
	"""
	Returns _-joined by default, but if there's a common dir prefix, removes most of it.
	"""
	if '/' not in filenames[0]:
		return '_'.join(filenames)
	else:
		print("\"{}\"".format(filenames[0][:filenames[0].index('/')]))
		print(list(map(lambda f: f[f.index('/')+1:], filenames)))

		return filenames[0][:filenames[0].index('/')]+'/'+'_'.join(list(map(lambda f: f[f.index('/')+1:], filenames)))

def plot_simrna(filename: List[str], dump: bool, outfile_format: str, xmax: float, ymax_: float) -> None: 

    print("In simrna mode")
    scores = []
    rmsds = []
    lines_arr = []
    scoreidx = 1
    rmsidx = 0

    for fn in filename:
        with open(fn) as f:
            lines_arr.append([ line for line in f])
   
    mins, maxes = [], []
    for lines in lines_arr:

        scores = [float(line.split()[scoreidx]) for line in lines if floatable(line.split()[scoreidx])]
        rmsds = [float(line.split()[rmsidx]) for line in lines if floatable(line.split()[rmsidx])]
        print(len(scores), len(rmsds))

        # How large should markersize be? I'd say, start at 25 and reduce every order of magnitude.
        markersize = (6 - np.log10(len(rmsds))) ** 2
        #print((scores), (rmsds))
        plt.scatter(rmsds, scores, s=markersize, alpha=min(0.1,1/(np.log(len(rmsds))-2.5)))
        #sns.kdeplot(rmsds, scores, cmap="Blues", shade=True)
        
        pctilescore = sorted(scores)[int(.95*len(scores))]
        print(pctilescore)
        mins.append(min(scores))
        #maxes.append(min(max([score for score in scores if score <= pctilescore]), min(scores)+100))
        maxes.append(max([score for score in scores if score <= pctilescore]))

    # Arbitrarily set ymax to be 0 to exclude FARFAR outliers.
    # Later develop a smart strategy having to do with actual outlier detection
    # Strategy: get range of min and max, expand by 5% either way.
    ymin, ymax = min(mins), max(maxes)
    #ymin2, ymax2 = ymin - (ymax - ymin)*0.05, -2100 #ymax + (ymax - ymin)*0.05
    ymin2, ymax2 = ymin - (ymax - ymin)*0.05, ymax + (ymax - ymin)*0.05
    #plt.ylim([min(mins)-10, max(maxes)+10])
    if ymax_ != -1:
        plt.ylim([ymin2, ymax_])
    else:
        plt.ylim([ymin2, ymax2])
    #plt.show()
    #plt.ylim([-90,0])
    plt.xlim([0, xmax])

    if dump:
        plt.savefig("%s.%s" % ("_".join(filename), outfile_format))
    else:
        plt.show()

def plot_svs(filename: List[str], dump: bool, outfile_format: str) -> None:
    scores = []
    rmsds = []
    lines_arr = []

    for fn in filename:
        with open(fn) as f:
            lines_arr.append([ line for line in f if line[0:5] == "SCORE"])
   
    mins, maxes = [], []
    for lines in lines_arr:
        if scoreidx == -1 or len(lines_arr) > 1: scoreidx = detect_score_col(lines)
        if rmsidx == -1 or len(lines_arr) > 1: rmsidx = detect_rms_col(lines)

        scores = [float(line.split()[scoreidx]) for line in lines if floatable(line.split()[scoreidx])]
        rmsds = [float(line.split()[rmsidx]) for line in lines if floatable(line.split()[rmsidx])]
        print(len(scores), len(rmsds))

        # How large should markersize be? I'd say, start at 25 and reduce every order of magnitude.
        markersize = (6 - np.log10(len(rmsds))) ** 2
        #print((scores), (rmsds))
        plt.scatter(rmsds, scores, s=markersize, alpha=min(0.1,1/(np.log(len(rmsds))-2.5)))
        plt.xlabel("silent_score")
        plt.ylabel("score")
        #sns.kdeplot(rmsds, scores, cmap="Blues", shade=True)
        
        pctilescore = sorted(scores)[int(.95*len(scores))]
        print(pctilescore)
        mins.append(min(scores))
        #maxes.append(min(max([score for score in scores if score <= pctilescore]), min(scores)+100))
        maxes.append(max([score for score in scores if score <= pctilescore]))


    scoremax = sorted(scores)[int(0.95*len(scores))]
    # Arbitrarily set ymax to be 0 to exclude FARFAR outliers.
    # Later develop a smart strategy having to do with actual outlier detection
    # Strategy: get range of min and max, expand by 5% either way.
    ymin, ymax = min(mins), max(maxes)
    ymin2, ymax2 = ymin - (ymax - ymin)*0.05, ymax + (ymax - ymin)*0.05
    #plt.ylim([min(mins)-10, max(maxes)+10])
    plt.ylim([ymin2, ymax2])
    #plt.show()
    #plt.ylim([-15,-3])
    #plt.xlim([0, xmax])

    if dump:
        plt.savefig("%s.%s" % ("_".join(filename), outfile_format))
    else:
        plt.show()

def plot(filename: List[str], scoreidx: int, rmsidx: int, dump: bool, outfile_format: str, xmax: float, ymax_: float) -> None: 
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

        scores = [float(line.split()[scoreidx]) for line in lines if floatable(line.split()[scoreidx])]# and float(line.split()[scoreidx]) < 0]
        rmsds = [float(line.split()[rmsidx]) for line in lines if floatable(line.split()[rmsidx])]# and float(line.split()[scoreidx]) < 0 ]
        print(len(scores), len(rmsds))

        # How large should markersize be? I'd say, start at 25 and reduce every order of magnitude.
        markersize = (6 - np.log10(len(rmsds))) ** 2
        #print((scores), (rmsds))
        plt.scatter(rmsds, scores, s=markersize)#, alpha=min(0.1,1/(np.log(len(rmsds))-2.5)))
        #sns.kdeplot(rmsds, scores, cmap="Blues", shade=True)
        
        pctilescore = sorted(scores)[int(.95*len(scores))]
        print(pctilescore)
        mins.append(min(scores))
        #maxes.append(min(max([score for score in scores if score <= pctilescore]), min(scores)+100))
        maxes.append(max([score for score in scores if score <= pctilescore]))


    # Arbitrarily set ymax to be 0 to exclude FARFAR outliers.
    # Later develop a smart strategy having to do with actual outlier detection
    # Strategy: get range of min and max, expand by 5% either way.
    ymin, ymax = min(mins), max(maxes)
    ymin2, ymax2 = ymin - (ymax - ymin)*0.05, ymax + (ymax - ymin)*0.05
    #ymin2, ymax2 = -90, 0
    #plt.ylim([min(mins)-10, max(maxes)+10])
    if ymax_ != -1:
        plt.ylim([ymin2, ymax_])
    else:
        plt.ylim([ymin2, ymax2])
    #plt.show()
    #plt.ylim([-15,-3])
    plt.xlim([0, xmax])

    if dump:
        plt.savefig("%s.%s" % (output_name(filename), outfile_format))
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
    parser.add_argument('--ymax', type=float, nargs='?',
            help='ymax? default -1 (set from data).', default=-1)
    parser.add_argument('--simrna', action='store_true', 
            help='just do simrna mode')
    parser.add_argument('--svs', action='store_true',
            help='score_vs_score mode')

    args = parser.parse_args()

    if args.simrna:
        plot_simrna(args.filename, args.dump, args.format, args.xmax, args.ymax)
    elif args.svs:
        plot_svs(args.filename, args.dump, args.format)
    else:
        #plot(filename, scoreidx, rmsidx, dump, outfile_format)
        plot(args.filename, args.score_col, args.rmsd_col, args.dump, args.format, args.xmax, args.ymax)
