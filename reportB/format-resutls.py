#!/usr/bin/env python3

import os, sys
import fileinput
import re
import statistics

for f in ["i1-results.txt", "i2-results.txt"]:
    vals = []

    for line in open(f):
        res = re.search('real\D+(\d+)m([\d.]+)s', line)
        if res:
            time = int(res.group(1)) * 60 + float(res.group(2))
            vals.append(time)

    vals.sort()

    mid = round(len(vals) / 2)
    begin = vals[0:mid]
    end = vals[mid:len(vals)]

    output = r"""
        \addplot+[
            boxplot prepared={{
                lower whisker={},
                lower quartile={},
                median={},
                upper quartile={},
                upper whisker={},
            }},
        ]
        coordinates{{}};
    """.format(
            min(vals),
            statistics.median(vals),
            statistics.median(begin),
            statistics.median(end),
            max(vals)
        )

    print (output)
