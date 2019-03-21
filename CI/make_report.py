#!/usr/bin/env python3

"""
This produces a structured normalized report json file based on warnings generated by another tool.
Currently implemented is clang-tidy warnings.
"""

import argparse
import re
from collections import namedtuple
from itertools import groupby
import os
import html
from fnmatch import fnmatch
import json
import sys


from codereport import CodeReport, ReportItem

def parse_clang_tidy_item(itemstr):

    try:
        m = re.match(r"(?P<file>[/.\-+\w]+):(?P<line>\d+):(?P<col>\d+): (?P<sev>.*?):(?P<msg>[\s\S]*?)\[(?P<code>.*)\]\n(?P<info>[\s\S]*)", itemstr)

        item = ReportItem(
            path=m.group("file"),
            line=int(m.group("line")),
            col=int(m.group("col")),
            message=m.group("msg").strip(),
            code=m.group("code"),
            severity=m.group("sev")
        )
        
        print(repr(item))

        return item
    except:
        print("Failed parsing clang-tidy item:")
        print("-"*20)
        print(itemstr)
        print("-"*20)
        raise

def parse_clang_tidy_output(output):

    # cleanup
    itemstr = output
    itemstr = re.sub(r"Enabled checks:\n[\S\s]+?\n\n", "", itemstr)
    itemstr = re.sub(r"clang-tidy-\d\.\d.*\n?", "", itemstr)
    itemstr = re.sub(r"clang-apply-.*", "", itemstr)
    itemstr = re.sub(r".*-header-filter.*", "", itemstr)

    items = []
    prevstart = 0

    matches = list(re.finditer(r"([\w/.\-+]+):(\d+):(\d+): (?:(?:warning)|(?:error)):", itemstr))
    for idx, m in enumerate(matches):
        # print(m)
        start, end = m.span()
        if idx > 0:
            item = itemstr[prevstart:start]
            items.append(item)
        if idx+1 == len(matches):
            item = itemstr[start:]
            items.append(item)
        prevstart = start

    items = set(map(parse_clang_tidy_item, sorted(items)))

    return items


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode", choices=("clang-tidy",),
                   help="Type of input warnings")
    p.add_argument("inputfile",
                   help="The input file containing the warnings")
    p.add_argument("output", default="codereport_clang_tidy.json",
                   help="The resulting JSON file")
    p.add_argument("--exclude", "-e", action="append", default=[],
                   help="Exclude files that match any of these patterns")
    p.add_argument("--filter", action="append", default=[],
                   help="Only include files that match any of these patterns")

    args = p.parse_args()

    if args.mode == "clang-tidy":
        with open(args.inputfile, "r", encoding="utf-8") as f:
            inputstr = f.read()
        items = parse_clang_tidy_output(inputstr)

        def select(item):
            accept = True
            if len(args.filter) > 0:
                accept = accept and all(fnmatch(item.path, e) for e in args.filter)

            accept = accept and not any(fnmatch(item.path, e) for e in args.exclude)
            return accept

        items = filter(select, items)


        data = [i.dict() for i in items]
        print("Write to", args.output)
        with open(args.output, "w+") as jf:
            json.dump(data, jf, indent=2)


if "__main__" == __name__:
    main()


