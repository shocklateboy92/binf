#!/usr/bin/env python3

import os, sys
import fileinput
import re

for f in sys.argv[1:]:
    fnum = int(re.match('results(\d)\.txt', f).group(1))
    vals = []
    run = 0
    for line in open(f):
        res = re.search('real\D+(\d+)m([\d.]+)s', line)
        if res:
            time = int(res.group(1)) * 60 + float(res.group(2))
            vals.append("({}, {})".format(run, time))
            run += 1

    print(r"\addplot +[mark=none] coordinates {{{}}};".format(" ".join(vals)))
    print(r"\addlegendentry{{n = {}}}".format(fnum))


