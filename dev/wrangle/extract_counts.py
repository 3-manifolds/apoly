import os, sys, re

outfile = open('rep_counts', 'w')
for line in open('knot_reps_out'):
    parts = line.split('\t')
    outfile.write('%s\t%s\n' % (parts[0], parts[3]))
