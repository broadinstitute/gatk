from datetime import datetime
from subprocess import run
import argparse
import subprocess


def build_tag(args):
    current_date = datetime.now()
    date_string = current_date.strftime('%Y-%m-%d')

    return f"{date_string}-alpine-{args.image_id}"


def build_argument_parser():
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Build a tag for Variants Docker image')
    parser.add_argument('-i', '--image-id', required=True, help='Docker image ID')
    parser.add_argument('-t', '--image-type', required=True, help='Docker image type, should be "alpine" or "slim".')
    return parser


if __name__ == '__main__':
    arg_parser = build_argument_parser()
    parsed_args = arg_parser.parse_args()
    tag = build_tag(parsed_args)
    print(tag, end='')
