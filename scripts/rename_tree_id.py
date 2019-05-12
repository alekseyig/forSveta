#!/usr/bin/env python

from __future__ import print_function

import re
import sys

re_desc = re.compile("(?P<id>[0-9.]+)\s+(?P<name>.*)")


def main(args):
    desc_file_name, tree_file_name = args[:2]

    desc = {}
    with open(desc_file_name) as f:
        for ln in f:
            ln = ln.strip('\n')
            match = re.match(re_desc, ln)
            if match is not None:
                desc[match.group('id')] = match.group('name').strip()

    with open(tree_file_name) as f:
        for ln in f:
            ln = ln.strip('\n')
            fields = ln.split(' ')
            new_fields = []
            for i in fields:
                if ':' in i:
                    k = i.rstrip(':')
                    if k in desc:
                        new_fields.append('{}:'.format(desc[k]))
                    else:
                        new_fields.append(i)
                else:
                    new_fields.append(i)

            print(' '.join(new_fields))


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
