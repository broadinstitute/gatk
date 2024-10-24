import argparse

import re
from typing import List, TextIO


def filter_chromosomes(output_path: str, full_interval_list_path: str, *chromosomes: str) -> None:
    # Chromosome names surrounded by word boundaries or trailed by an underscore.
    with open(output_path, 'w') as output, open(full_interval_list_path, 'r') as full_interval_list:
        res = [f'\\b{c}\\b|{c}_' for c in chromosomes]
        chromosome_patterns = [re.compile(r) for r in res]
        for line in full_interval_list:
            if re.search(r"^@HD|^@PG", line):
                output.write(line)
            else:
                if any([re.search(pattern, line) for pattern in chromosome_patterns]):
                    output.write(line)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False,
                                     description='Filters interval list for the specified chromosomes')

    parser.add_argument('--input-interval-list', help='Input full interval list', required=True)
    parser.add_argument('--output-interval-list', help='Output filtered interval list', required=True)
    parser.add_argument('--chromosome', action='append', required=True)

    args = parser.parse_args()
    filter_chromosomes(args.output_interval_list, args.input_interval_list, *args.chromosome)
