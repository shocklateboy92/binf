#!/usr/bin/env python3

import os, sys
import fileinput
import re
import statistics

means = []

for f in sys.argv[1:]:
    fnum = int(re.search('results(\d+)\.txt', f).group(1))
    vals = []
    run = 0
    for line in open(f):
        res = re.search('user\D+(\d+)m([\d.]+)s', line)
        if res:
            time = int(res.group(1)) * 60 + float(res.group(2))
            vals.append(time)
            run += 1

    means.append((fnum, statistics.mean(vals)))

print(r"\addplot +[] coordinates {{{}}};".format(" ".join([str(v) for v in means])))

#assert means[0][0] == 1
#
#st = means[0][1]
#speedups = [st / t for (x, t) in means]
#
#print (speedups)
