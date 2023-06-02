import argparse

import re
from typing import List, TextIO


def filter_chromosomes(interval_list: TextIO, chromosomes: List[str]) -> None:
    res = [f'\\b{c}\\b|{c}_' for c in chromosomes]
    print(res)
    chromosome_patterns = [re.compile(r) for r in res]
    for line in interval_list:
        if re.search(r"^@HD", line):
            print(line, end='')
        else:
            if any([re.search(pattern, line) for pattern in chromosome_patterns]):
                print(line, end='')


def main(full_interval_list: str, chromosomes: List[str]):
    with open(full_interval_list, 'r') as interval_list:
        filter_chromosomes(interval_list, chromosomes)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False,
                                     description='Filters interval list for the specified chromosomes')

    parser.add_argument('--full-interval-list', help='Full interval list', required=True)
    parser.add_argument('--chromosome', nargs='+', required=True)

    args = parser.parse_args()
    main(args.full_interval_list, args.chromosome)
