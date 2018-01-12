#!/usr/bin/env python2.7

import sys
import pandas
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
#import statistics

def no_bracketed(st):
    """
    Remove bracketed substrings from a string -- used to convert to a 
    one letter sequence:
    acgu[URA:VirtualPhosphate]acgu ==> acguacgu
    """

    ret = ""
    incl = True
    for c in st:
        if incl:
            if c == '[': incl = False
            else: ret += c
        else:
            if c == ']': incl = True
            else: pass
    return ret

def median(arr):
    if len(arr) % 2:
        return (sorted(arr)[len(arr)/2] + sorted(arr)[len(arr)/2-1])/2
    else:
        return sorted(arr)[len(arr)/2]

def most_common(seq_to_score_dict, num_limit, best_median_scores=None):
    """
    Obtain the `num_limit` most common sequences -- which are the keys
    of the first argument (a dict from sequence to a list of scores).
    If `best_median_scores` is provided, it is a list to which any
    sequences selected must also belong.
    """

    seq_to_freq = { seq: len(seq_to_score_dict[seq]) for seq in seq_to_score_dict.keys() }
    i = 1
    most_populated_sequences = []
    print seq_to_freq
    for seq, pop in sorted(seq_to_freq.iteritems(), key=lambda (k,v): (v,k), reverse=True):
        if best_median_scores is not None and not seq in best_median_scores: continue
        most_populated_sequences.append(seq)
        if i >= num_limit: break
        i += 1
    print "most_populated_sequences",most_populated_sequences
    return most_populated_sequences

def best_medians(seq_to_score_dict, num_limit):
    """
    Obtain the `num_limit` sequences with the lowest median scores.
    """
    seq_to_freq = { seq: median(sorted(seq_to_score_dict[seq])) for seq in seq_to_score_dict.keys() }
    i = 1
    best_scoring = []
    for seq, pop in sorted(seq_to_freq.iteritems(), key=lambda (k,v): (v,k)):
        if len(seq_to_score_dict[seq]) == 1: continue
        best_scoring.append(seq)
        if i >= num_limit: break
        i += 1
    return best_scoring

def combine_dicts(dict1, dict2, which_values_are_key=1):
    """
    Suppose you have two dicts with a common key that you want to 'factor out'
    to form a new dict where one of the parent's values is the new key and the
    values are lists of the other dict's values. There is probably a concise
    name for this operation, but to be clear I'll give the concrete example
    from this application:

    Suppose you have been reading through a flat text file. Reading separate 
    lines at separate times, you discover a bunch of mappings from tag => score
    and tag => sequence. The tag is just an incidental common key though -- we
    care literally not at all about S_000001 for our purposes. What we really
    want and need is to be able to get 'all scores for a certain sequence.'
    """

    if which_values_are_key == 2:
        return combine_dicts(dict2, dict1, which_values_are_key=1)

    combined = {}
    for tag in dict2.keys():
        if dict1[tag] in combined:
            combined[dict1[tag]].append(dict2[tag])
        else:
            combined[dict1[tag]] = [dict2[tag]]
    return combined

def make_lambda(seq):
    """
    Make a function that will slice an input string in the same way that you would
    slice `seq` to obtain all of its 'n's.
    """

    indices = []
    for i,c in enumerate(seq):
        if c == 'n': indices.append(i)
    return lambda(st): "".join([c for i,c in enumerate(st) if i in indices])

def data_frame_input_format(tag_score, tag_sequence, include_sequences):
    """
    Take two dicts -- from tag to score and tag to sequence -- as well as a list
    to filter by sequences -- and turn them into the 'list-of-identically-keyed-dicts'
    input format for a pandas DataFrame constructor.
    """
    protodf = []
    for tag in tag_score.keys():
        if not tag in tag_sequence: continue
        if not tag_sequence[tag] in include_sequences: continue
        protodf.append({'tag': tag, 'sequence': tag_sequence[tag], 'score': tag_score[tag]})
    return protodf

def process_silent(outfile):
    """
    Take a silent file and extract the mappings of
    interest: from tag to score and from tag to sequence.
    """

    lines = []
    with open(outfile) as f: lines = f.readlines()
    tag_to_score = {}
    tag_to_sequence = {}
    description_id = -1
    score_id = -1
    missing_id = -1
    seq_slice = None
    for line in lines:
        objs = line.strip().split()
        if objs[0] == "SEQUENCE:":
            seq_slice = make_lambda(objs[1])

        if description_id == -1 and 'description' in objs:
            description_id = objs.index('description')

        if score_id == -1 and 'score' in objs:
            score_id = objs.index('score')
        
        if missing_id == -1 and 'missing' in objs:
            missing_id = objs.index('missing')

        if 'SCORE:' in objs and not 'description' in objs and objs[missing_id] == '0':
            tag = objs[description_id]
            tag_to_score[tag] = float(objs[score_id])

        if 'ANNOTATED_SEQUENCE:' in objs:
            tag_to_sequence[objs[-1]] = seq_slice(no_bracketed(objs[1]))
    return tag_to_score, tag_to_sequence

def sorted_by_median(most_populated_sequences, seq_to_score_dict):
    seq_to_freq = { seq: median(sorted(seq_to_score_dict[seq])) for seq in seq_to_score_dict.keys() }
    sorted_list = []
    for seq, pop in sorted(seq_to_freq.iteritems(), key=lambda (k,v): (v,k)):
        if seq not in most_populated_sequences: continue
        sorted_list.append(seq)
    return sorted_list

def sequence_recovery(outfile, image_format, num_limit, consider_best_medians):

    tag_to_score, tag_to_sequence = process_silent(outfile)
    sequence_to_scores = combine_dicts(tag_to_score, tag_to_sequence, which_values_are_key=2)

    # Filter for the top -- say -- 16 most often recovered sequences.
    if consider_best_medians:
        best_median_scores = best_medians(sequence_to_scores, num_limit*2)
        most_populated_sequences = most_common(sequence_to_scores, num_limit, best_median_scores)
    else:
        most_populated_sequences = most_common(sequence_to_scores, num_limit)

    # convert to dataframe
    protodf = data_frame_input_format(tag_to_score, tag_to_sequence, most_populated_sequences)
    df = pandas.DataFrame.from_dict(protodf)
    g = sns.violinplot(x='sequence', y='score', data=df, 
            #order=sorted(sequence_to_scores.keys()),
            #order=sorted(most_populated_sequences),
            # sort by median score
            order=sorted_by_median(most_populated_sequences, sequence_to_scores),
            dodge=False)
    for item in g.get_xticklabels(): item.set_rotation(45)

    # if image format is 'complex' do each 
    image_formats = image_format.split(',')
    for fmt in image_formats: plt.savefig("%s.%s" % (outfile, fmt))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process a design simulation silent file for sequence recovery, assuming FULL SEQUENCE/ ANNOTATED_SEQUENCE are present.')
    parser.add_argument('silent', metavar='silent', type=str, help='Silent file to be parsed')
    parser.add_argument('--image_format', dest='fmt', metavar='format', default='png', type=str, nargs='?', help='File format for output. svg, png are perhaps preferred; separate multiple formats with commas.') 
    parser.add_argument('--num_sequences', dest='num_limit', metavar='limit', default=0, type=int, nargs='?', help='Number of sequences to retain for plotting (default: all).') 
    parser.add_argument('--best_medians', dest='median_filter', metavar='limit', default=True, type=bool, nargs='?', help='Consider only the nlimit*2 sequences with the best median scores.') 
    args = parser.parse_args()

    sequence_recovery(args.silent, args.fmt, args.num_limit, args.median_filter)
